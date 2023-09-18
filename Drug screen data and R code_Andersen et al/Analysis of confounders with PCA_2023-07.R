library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(drc)
library(data.table)
library(fastDummies)
library(survival)
library(survminer)
library(readxl)

# Path ----
path <- "/mnt/CommonStorageRAI/Aram/Drug screen - clinical forecasting/Clinical forecasting manuscript - revision 2023-07"

# Load ----
load(file = paste0(path,"/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(file = paste0(path,"/Drug sensitivity QC results and batch data_2023-07.RData"))
load(file = paste0(path,"/Survival, clinical features and ELN2022 classifications_2023-07.RData"))

# Analysis of confounding components ----

X_biological <- cbind(df_genetics, df_prognostics_dummy[,grep("Primary|Age|Sex|FAB|ELN",colnames(df_prognostics_dummy))])

mat_rAUC <- reshape2::dcast(df_scores %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = "rAUC")
rownames(mat_rAUC) <- as.character(mat_rAUC$Patient.num); mat_rAUC$Patient.num <- NULL
mat_rAUC <- mat_rAUC[match(rownames(Y), rownames(mat_rAUC)),]
SDs <- matrixStats::colSds(data.matrix(mat_rAUC))

mat_DSS3 <- reshape2::dcast(df_scores %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = "DSS3")
rownames(mat_DSS3) <- as.character(mat_DSS3$Patient.num); mat_DSS3$Patient.num <- NULL
mat_DSS3 <- mat_DSS3[match(rownames(Y), rownames(mat_DSS3)),]
SDs3 <- matrixStats::colSds(data.matrix(mat_DSS3))

df_Z.prime_avg1 <- df_Z.prime_avg %>% subset(!grepl("re", Patient.ID)) %>% as.data.frame()
df_Z.prime_avg1$Mean_DMSO <- log(df_Z.prime_avg1$Mean_DMSO)
df_Z.prime_avg1$Mean_BZCL <- log(df_Z.prime_avg1$Mean_BZCL)
df_Z.prime_avg1$CV_DMSO <- (df_Z.prime_avg1$SD_DMSO/df_Z.prime_avg1$Mean_DMSO)
df_Z.prime_avg1$CV_BZCL <- (df_Z.prime_avg1$SD_BZCL/df_Z.prime_avg1$Mean_BZCL)
rownames(df_Z.prime_avg1) <- df_Z.prime_avg1$Patient.num
df_Z.prime_avg1 <-  df_Z.prime_avg1[match(rownames(Y), df_Z.prime_avg1$Patient.num),]

df_control_densities <- df_Z.prime_avg1[,c("Mean_DMSO","Mean_BZCL")] %>% as.data.frame()

df_control_noise <- df_Z.prime_avg1[,c("CV_DMSO","CV_BZCL", "Z_score")] %>% as.data.frame()

df_batch <- df_Z.prime_avg1[,c("Period","Instrument","Source", "Seeding_number")] %>% as.data.frame()

df_batch <- dummy_cols(df_batch)
df_batch <- df_batch[,-c(1:4)]
df_batch$Period_1 <- NULL
df_batch$Instrument_EnVision <- NULL
df_batch$Source_Blood <- NULL
df_batch$`Seeding_number_<10000` <- NULL
rownames(df_batch) <- rownames(df_Z.prime_avg1) 

df_fit_uncertainty1 <- df_fit_uncertainty %>% subset(!grepl("re", Patient.ID)) %>% as.data.frame()
df_fit_uncertainty1 <-  df_fit_uncertainty1[match(rownames(Y), df_fit_uncertainty1$Patient.num),]

df_non_responders <- df_fit_uncertainty1[,c("DSS3_0", "EC50_max","TEC50_max")] %>% as.data.frame()
rownames(df_non_responders) <- df_fit_uncertainty1$Patient.num

df_fit <- df_fit_uncertainty1[,c("believe_DSS",#"SE_EC50_40",  "MAE_20","Max_residual_30",
                                "SE_EC50_mean", "MAE_mean","Max_residual_mean")] %>% as.data.frame()
rownames(df_fit) <- df_fit_uncertainty1$Patient.num

df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$rAUC <- 1-df_scores$rAUC

df_r2 <- c()
df_PC <- c()
for(t in c(quantile(SDs, 1-c(100)/ncol(mat_rAUC)),0)){
  for(pt.st in c(T, F)){
  for(metric in colnames(df_scores)[6:12]){
    
    X <- reshape2::dcast(df_scores %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = metric)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
    X <- X[match(rownames(Y), rownames(X)),]
    X <- data.matrix(X)

    X <- X[,SDs>=t]
    
    if(pt.st){
      X <- (X - rowMeans(X))/matrixStats::rowSds(X)
    }
    
    x_sd <- matrixStats::colSds(X)
    x_mean <- colMeans(X)
    X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
    if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
    svdz <- svd(t(X))
    df_PCA <- data.frame(svdz$v[,1:50]); colnames(df_PCA) <- paste0("PC", 1:50); rownames(df_PCA) <- rownames(mat_rAUC)
    
    m <- data.frame(Percent_of_variance=((svdz$d^2)/sum(svdz$d^2))[1:50], PC = colnames(df_PCA))
    m$Metric <- metric
    m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
    m$Pre_selection_N <- ncol(X)
    df_PC <- rbind(df_PC, m)
    
    for(i in colnames(df_PCA)[1:50]){
      
      
      dat <- cbind(df_PCA[,i],df_control_densities[match(rownames(df_PCA), rownames(df_control_densities)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Control densities")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      dat <- cbind(df_PCA[,i],df_control_noise[match(rownames(df_PCA), rownames(df_control_noise)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Control noise")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      dat <- cbind(df_PCA[,i],df_non_responders[match(rownames(df_PCA), rownames(df_non_responders)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Curvefit non-responders")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      
      dat <- cbind(df_PCA[,i],
                   df_fit[match(rownames(df_PCA), rownames(df_fit)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Curvefit error")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      dat <- cbind(df_PCA[,i],df_batch[match(rownames(df_PCA), rownames(df_batch)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Batch covariates")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      
      dat <- cbind(df_PCA[,i],
                   X_biological[match(rownames(df_PCA), rownames(X_biological)),])
      
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Biological/clinical covariates")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
    }
  }
}
}
df_PC$PC <- factor(df_PC$PC, levels = unique(df_r2$PC))
df_r2$PC <- factor(df_r2$PC, levels = unique(df_r2$PC))

df_r2$Percent_of_variance <- df_PC$Percent_of_variance[match(paste(df_r2$Metric,df_r2$Pre_selection_N, df_r2$Patient_stdz,as.character(df_r2$PC)),paste(df_PC$Metric,df_PC$Pre_selection_N, df_PC$Patient_stdz,as.character(df_PC$PC)))]
#df_r2$adj.r.squared[df_r2$adj.r.squared<0] <- 0
df_r2 <- df_r2 %>%  mutate(Variance_explained=adj.r.squared*Percent_of_variance)
df_r2 <- df_r2 %>% group_by(Variable,Pre_selection_N,Metric,Patient_stdz) %>% mutate(Cumulative_variance_explained=cumsum(Variance_explained))
df_r2$Variable <- factor(df_r2$Variable, levels = unique(df_r2$Variable))

df_PC$Metric <- factor(df_PC$Metric, levels=unique(df_PC$Metric))
df_PC$Metric2 <- paste0(gsub("_"," ",df_PC$Metric),"\n",ifelse(df_PC$Patient_stdz=="", "", "z-score"))
df_PC <- df_PC[order(df_PC$Metric, df_PC$Patient_stdz),]
df_PC$Metric2 <- factor(df_PC$Metric2, levels=unique(df_PC$Metric2))

df_r2$Metric <- factor(df_r2$Metric, levels=unique(df_r2$Metric))
df_r2$Metric2 <- paste0(gsub("_"," ",df_r2$Metric),"\n",ifelse(df_r2$Patient_stdz=="", "", "z-score"))
df_r2 <- df_r2[order(df_r2$Metric, df_r2$Patient_stdz),]
df_r2$Metric2 <- factor(df_r2$Metric2, levels=unique(df_r2$Metric2))

ggplot(df_r2 %>% subset(Pre_selection_N == 349))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Variable), alpha=1)+
  theme_bw()+
  scale_color_manual(values=pal_d3()(9)[c(1:3,6:7,9)])+
  facet_wrap(~Metric2 )+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=6.5, height=5.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_cumulative variance explained_349_1.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N != 349))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Variable), alpha=1)+
  theme_bw()+
  scale_color_manual(values=pal_d3()(9)[c(1:3,6:7,9)])+
  facet_wrap(~Metric2 )+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=6.5, height=5.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_cumulative variance explained_100_1.pdf"))


ggplot(df_PC %>% subset(Pre_selection_N == 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_jco()+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_349_1.pdf"))


ggplot(df_PC %>% subset(Pre_selection_N == 349) %>%
         subset(as.numeric(PC) %in% 1:10) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_349_2.pdf"))

ggplot(df_PC %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_jco()+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_100_1.pdf"))


ggplot(df_PC %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_100_2.pdf"))




ggplot(df_r2 %>% subset(Pre_selection_N == 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(9))+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_349_1.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N == 349) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_349_2.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(9))+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_100_1.pdf"))

ggplot(df_r2 %>% subset(Pre_selection_N != 349) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_100_2.pdf"))

df_PC.1 <- c()
df_r2.1 <- c()
for(t in c(quantile(SDs3, 1-c(100)/ncol(mat_rAUC)))){
  for(pt.st in c(T, F)){
    for(metric in colnames(df_scores)[6:12]){
      
      X <- reshape2::dcast(df_scores %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = metric)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
      X <- X[match(rownames(Y), rownames(X)),]
      X <- data.matrix(X)
      
      X <- X[,SDs3>=t]
      
      if(pt.st){
        X <- (X - rowMeans(X))/matrixStats::rowSds(X)
      }
      
      x_sd <- matrixStats::colSds(X)
      x_mean <- colMeans(X)
      X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
      if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
      svdz <- svd(t(X))
      df_PCA <- data.frame(svdz$v[,1:50]); colnames(df_PCA) <- paste0("PC", 1:50); rownames(df_PCA) <- rownames(mat_rAUC)
      
      m <- data.frame(Percent_of_variance=((svdz$d^2)/sum(svdz$d^2))[1:50], PC = colnames(df_PCA))
      m$Metric <- metric
      m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
      m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
      df_PC.1 <- rbind(df_PC.1, m)
      
      for(i in colnames(df_PCA)[1:50]){
        
        
        dat <- cbind(df_PCA[,i],df_control_densities[match(rownames(df_PCA), rownames(df_control_densities)),])
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Control densities")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
        dat <- cbind(df_PCA[,i],df_control_noise[match(rownames(df_PCA), rownames(df_control_noise)),])
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Control noise")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
        dat <- cbind(df_PCA[,i],df_non_responders[match(rownames(df_PCA), rownames(df_non_responders)),])
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Curvefit non-responders")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
        
        dat <- cbind(df_PCA[,i],
                     df_fit[match(rownames(df_PCA), rownames(df_fit)),])
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Curvefit error")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
        dat <- cbind(df_PCA[,i],df_batch[match(rownames(df_PCA), rownames(df_batch)),])
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Batch covariates")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
        
        dat <- cbind(df_PCA[,i],
                     X_biological[match(rownames(df_PCA), rownames(X_biological)),])
        
        colnames(dat)[1] <- "PC"
        lm.fit <- summary(lm(PC~.,dat))
        m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Biological/clinical covariates")
        m$Metric <- metric
        m$Patient_stdz <- ifelse(pt.st,"Patient standardization", "")
        m$Pre_selection_N <- paste0(" ","DSS3",ncol(X))
        df_r2.1 <- rbind(df_r2.1, m); rm(lm.fit, m, dat)
        
      }
    }
  }
}
df_PC.1$PC <- factor(df_PC.1$PC, levels = unique(df_r2.1$PC))
df_r2.1$PC <- factor(df_r2.1$PC, levels = unique(df_r2.1$PC))

df_PC.1$PC <- factor(df_PC.1$PC, levels = unique(df_r2.1$PC))
df_r2.1$PC <- factor(df_r2.1$PC, levels = unique(df_r2.1$PC))

df_r2.1$Percent_of_variance <- df_PC.1$Percent_of_variance[match(paste(df_r2.1$Metric,df_r2.1$Pre_selection_N, df_r2.1$Patient_stdz,as.character(df_r2.1$PC)),paste(df_PC.1$Metric,df_PC.1$Pre_selection_N, df_PC.1$Patient_stdz,as.character(df_PC.1$PC)))]
#df_r2.1$adj.r.squared[df_r2.1$adj.r.squared<0] <- 0
df_r2.1 <- df_r2.1 %>%  mutate(Variance_explained=adj.r.squared*Percent_of_variance)
df_r2.1 <- df_r2.1 %>% group_by(Variable,Pre_selection_N,Metric,Patient_stdz) %>% mutate(Cumulative_variance_explained=cumsum(Variance_explained))
df_r2.1$Variable <- factor(df_r2.1$Variable, levels = unique(df_r2.1$Variable))

df_PC.1$Metric <- factor(df_PC.1$Metric, levels=unique(df_PC.1$Metric))
df_PC.1$Metric2 <- paste0(gsub("_"," ",df_PC.1$Metric),"\n",ifelse(df_PC.1$Patient_stdz=="", "", "z-score"))
df_PC.1 <- df_PC.1[order(df_PC.1$Metric, df_PC.1$Patient_stdz),]
df_PC.1$Metric2 <- factor(df_PC.1$Metric2, levels=unique(df_PC.1$Metric2))

df_r2.1$Metric <- factor(df_r2.1$Metric, levels=unique(df_r2.1$Metric))
df_r2.1$Metric2 <- paste0(gsub("_"," ",df_r2.1$Metric),"\n",ifelse(df_r2.1$Patient_stdz=="", "", "z-score"))
df_r2.1 <- df_r2.1[order(df_r2.1$Metric, df_r2.1$Patient_stdz),]
df_r2.1$Metric2 <- factor(df_r2.1$Metric2, levels=unique(df_r2.1$Metric2))



ggplot(df_r2.1 %>% subset(Pre_selection_N != 349))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Variable), alpha=1)+
  theme_bw()+
  scale_color_manual(values=pal_d3()(9)[c(1:3,6:7,9)])+
  facet_wrap(~Metric2 )+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=6.5, height=5.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_cumulative variance explained_100 DSS3_1.pdf"))

ggplot(df_PC.1 %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_jco()+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_100 DSS3_1.pdf"))


ggplot(df_PC.1 %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(~Pre_selection_N+Patient_stdz)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="", y = "Percent of variance")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=2.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance_100 DSS3_2.pdf"))



ggplot(df_r2.1 %>% subset(Pre_selection_N != 349) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(9))+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_100 DSS3_1.pdf"))

ggplot(df_r2.1 %>% subset(Pre_selection_N != 349) %>% subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")) %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(Variable~Patient_stdz)+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=6.5, height=5.75, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_PC variance explained_100 DSS3_2.pdf"))


ggplot(bind_rows(df_r2 %>% mutate(Pre_selection_N = as.character(Pre_selection_N)),df_r2.1 %>% mutate(Pre_selection_N=gsub("DSS3","DSS3 ",Pre_selection_N))))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric), alpha=1)+
  theme_bw()+
  scale_color_jco()+
  facet_grid(Variable ~Pre_selection_N +Patient_stdz, scales="free" )+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=7, height=5.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_cumulative variance explained_1.pdf"))

ggplot(bind_rows(df_r2 %>% mutate(Pre_selection_N = as.character(Pre_selection_N)),df_r2.1 %>% mutate(Pre_selection_N=gsub("DSS3","DSS3 ",Pre_selection_N))) %>% 
         subset(Metric %in% c("DSS2", "DSS3", "rAUC", "rAUC_log2")) )+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric), alpha=1)+
  theme_bw()+
  scale_color_jco()+
  scale_color_manual(values=pal_jco()(10)[c(1:2,6:7)])+
  facet_grid(Variable ~Pre_selection_N +Patient_stdz, scales="free" )+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=7, height=5.5, 
       filename = paste0(path, "/Results_PCA/PCA of confounding factors_cumulative variance explained_2.pdf"))


