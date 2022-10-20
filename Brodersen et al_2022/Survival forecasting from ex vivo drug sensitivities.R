# Testing and analysis of survival predictions from ex vivo drug sensitivities for https://doi.org/10.1101/2022.10.11.509866
# Author: Aram N. Andersen
# Date: 20220711
#----Libraries----
rm(list=ls())
library(CoxTools)
library(glmnet)
library(glmnetUtils)
library(reshape)
library(dplyr)
library(doSNOW)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(ggsci)
library(matrixStats)

#----Comparing AUC/DSS/EC50 survival prediction----
path <- "D:/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Drug sensitivity data and surivival data - NOR_20220711.RData"))

df_scores$AUC <- df_scores$AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100

hist(df_scores$AUC, breaks=100)
hist(df_scores$DSS1, breaks=100)
hist(df_scores$DSS2, breaks=100)
hist(df_scores$DSS3, breaks=100)
hist(log10(df_scores$EC50), breaks=100)
hist(log10(df_scores$TEC50), breaks=100)



list_res_DS_metrics <- list()
for(m in colnames(df_scores)[6:11]){
  X <- cast(data = df_scores, Patient.num ~ drug , value=m, fun.aggregate = NULL)
  rownames(X) <-X$Patient.num
  X$Patient.num<-NULL
  if(grepl("EC50",m)){
    X <- log10(X)
    X <- X - mean(unlist(X))
  }
  list_res_DS_metrics[[m]] <-  Cox_forecasting_glmnet(X_data=X,
                                                          y_data=data.matrix(Y),
                                                          alpha=0,
                                                          lambda=c(exp(seq(-4,6, 0.1))),
                                                          free_cores = 20,
                                                          test.n= c(6,4),
                                                          nfolds = nrow(Y),
                                                          iter=200,
                                                          log_AUC=2,
                                                          Patient.Z=1:2,
                                                          Drug.Z =1:2,
                                                          RCPC=0:10)
}

list_res_AUC.raw <- Cox_forecasting_glmnet(X_data=X.AUC.raw,
                                               y_data=data.matrix(Y),
                                               alpha=0,
                                               lambda=c(exp(seq(-4,6, 0.1))),
                                               free_cores = 20,
                                               test.n= c(6,4),
                                               nfolds = nrow(Y),
                                               iter=200,
                                               log_AUC=1:2,
                                               Patient.Z=1:2,
                                               Drug.Z =1:2,
                                               RCPC=0:10)

save(list_res_DS_metrics, list_res_AUC.raw, file=paste0(path,"/Results 20220711/Drug sensitivity score comparison - NOR_20220711.RData"))


load(file=paste0(path,"/Results 20220711/Drug sensitivity score comparison - NOR_20220711.RData"))


names(list_res_DS_metrics)[1] <- "DSS_AUC"
list_res_DS_metrics[["AUC_raw"]] <- list_res_AUC.raw
df_C.index.alldata <- list_res_DS_metrics
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
  df_C.index.alldata[[m]]$ID <- gsub("AUC", m, df_C.index.alldata[[m]]$ID)
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}

df_C.index.alldata <- bind_rows(df_C.index.alldata)
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient stdz", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug stdz", "")
df_C.index.alldata$Metric <- gsub("\\/.*","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub("AUC_raw"," AUC_raw",df_C.index.alldata$Metric)
df_C.index.alldata$Metric <- gsub("log2"," log2",df_C.index.alldata$Metric, fixed=TRUE)
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, RCPC) %>% summarise(Median=median(C_index_test),
                                                                           Mean=mean(C_index_test))

ggexport(
  ggarrange(ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_test))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              #geom_jitter(alpha=0.2, size=0.1, width=0.2)+
              geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Testing: n=",10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank()),

            ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_train))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              #geom_jitter(alpha=0.2, size=0.1, width=0.2)+
              geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Training: n=",74-10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank())),
  width=6.5, height=4,
  filename = paste0(path,"/Results 20220711/Drug sensitivity score comparison_NOR_1.pdf")
)

ggexport(
  ggplot(df_C.index.alldata %>% subset(), aes(x=RCPC, y=C_index_test))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    #geom_jitter(alpha=0.2, size=0.1, width=0.2)+
    geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
    facet_grid(Patient_stdz+Drug_stdz~Metric)+
    ylim(0,1)+
    theme_bw()+
    labs(title=paste0("Testing: n=",10), x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  width=9, height=6,
  filename = paste0(path,"/Results 20220711/Drug sensitivity score comparison_NOR_2.pdf")
)

ggexport(
  ggplot(df_C.index.alldata %>% subset(), aes(x=RCPC, y=C_index_train))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    #geom_jitter(alpha=0.2, size=0.1, width=0.2)+
    geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
    facet_grid(Patient_stdz+Drug_stdz~Metric)+
    ylim(0,1)+
    theme_bw()+
    labs(title=paste0("Training: n=",74-10), x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  width=9, height=6,
  filename = paste0(path,"/Results 20220711/Drug sensitivity score comparison_NOR_3.pdf")
)



#----Comparing Ridge, lasso, elastic net and variable pre-selection----
path <- "D:/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Drug sensitivity data and surivival data - NOR_20220711.RData"))

list_res_AUC.raw.CVA <- Cox_forecasting_glmnet_CVA(X_data=X.AUC.raw,
                                                   y_data=data.matrix(Y),
                                                   alpha="CVA",
                                                   lambda=c(exp(seq(-4,6, 0.1))),
                                                   free_cores = 2,
                                                   test.n= c(6,4),
                                                   nfolds = nrow(Y),
                                                   iter=50,
                                                   log_AUC=1,
                                                   Patient.Z=1,
                                                   Drug.Z =2,
                                                   RCPC=0)

list_res_AUC.raw <- Cox_forecasting_glmnet(X_data=X.AUC.raw,
                                           y_data=data.matrix(Y),
                                           alpha=list_res_AUC.raw.CVA$CVA_results$alpha,
                                           lambda=c(exp(seq(-4,6, 0.1))),
                                           free_cores = 2,
                                           test.n= c(6,4),
                                           nfolds = nrow(Y),
                                           iter=50,
                                           log_AUC=1,
                                           Patient.Z=1,
                                           Drug.Z =2,
                                           RCPC=0)

list_res_AUC.raw_list0 <- list()

for(t in c(0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.170,0.175,0.18,0.185,0.190,0.195,0.2)){
  list_res_AUC.raw_list0[[paste("sd cutoff",t)]] <- Cox_forecasting_glmnet(X_data=X.AUC.raw[,colSds(data.matrix(X.AUC.raw), na.rm = T)>t],
                                                                           y_data=data.matrix(Y),
                                                                           alpha=list_res_AUC.raw.CVA$CVA_results$alpha,
                                                                           lambda=c(exp(seq(-4,6, 0.1))),
                                                                           free_cores = 2,
                                                                           test.n= c(6,4),
                                                                           nfolds = nrow(Y),
                                                                           iter=50,
                                                                           log_AUC=1,
                                                                           Patient.Z=1,
                                                                           Drug.Z =2,
                                                                           RCPC=0)
}
list_res_AUC.raw_list0 <- list()

save(list_res_AUC.raw,list_res_AUC.raw.CVA, list_res_AUC.raw_list0,
     file=paste0(path,"/Results 20220711/Elnet comp_with variable pre-selection from AUC_Nor.RData"))


df_C.index.alldata <- list_res_AUC.raw_list0
for(i in 1:length(df_C.index.alldata)){
  df_C.index.alldata[[i]] <- df_C.index.alldata[[i]]$C_index_results
}


df_C.index.alldata <- bind_rows(list_res_AUC.raw$C_index_results %>% mutate(Data = "full"),
                                bind_rows(df_C.index.alldata, .id="Data"))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient stdz", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug stdz", "")
df_C.index.alldata$Penalty <- gsub(".*_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub("\\/.*","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub("AUC_raw"," AUC_raw",df_C.index.alldata$Metric)
df_C.index.alldata$Metric <- gsub("log2"," log2",df_C.index.alldata$Metric, fixed=TRUE)
df_C.index.alldata$Pre_selection <- gsub(".* ", "", df_C.index.alldata$Data)
df_C.index.alldata$Pre_selection[df_C.index.alldata$Pre_selection=="full"] <- "0"
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Pre_selection, Penalty) %>% summarise(Median=median(C_index_test),
                                                                                             Mean=mean(C_index_test),
                                                                                             MAD=mad(C_index_test),
                                                                                             SD=sd(C_index_test))

ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = Pre_selection, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="Pre-selection (SD threshold)", fill="Mean C-index")+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=5, filename = paste0(path,"/Results 20220711/Pre-selection and penalty comparison_1.pdf"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = Pre_selection, fill = Median)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="Pre-selection (SD threshold)", fill="Median C-index")+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=5, filename = paste0(path,"/Results 20220711/Pre-selection and penalty comparison_2.pdf"))



#----Comparing clinical data with drug sensitivity data----
load(file=paste0(path,"/Results 20220711/Clinical data_curated.RDAta"))
Y_C75 <- Y

Z <- data.frame(clinical_74[,6:ncol(clinical_74)])
Z <- Z[match(rownames(X.AUC.raw), rownames(Z)),]

Y_extra <- data.frame(clinical_220[,2:3])
colnames(Y_extra) <- colnames(Y_C75)
pt_extra <- as.character(clinical_220$PatientCode_JElist)
Z_extra <- data.frame(clinical_220[,6:ncol(clinical_220)])

Z_data_extra <- data.matrix(Z_data_extra)
y_data_extra <- data.matrix(Y_extra)
y_data <- data.matrix(Y_C75)


df_C.index.C220 <- c()
for(a in c(0,0.4,1)){
  i=paste0("Clinical/genetic", "/Penalty_", a)
  X <- Z_data_extra
  Y <- y_data_extra

  cores <- parallel::detectCores()-free_cores
  cluster.cores<-makeCluster(cores)
  registerDoSNOW(cluster.cores)
  pb <- txtProgressBar(max=iter, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
    set.seed(b)
    setTxtProgressBar(pb,b)
    X.0 <- X[Y[,2]==0,]
    X.1 <- X[Y[,2]==1,]
    Y.0 <- Y[Y[,2]==0,]
    Y.1 <- Y[Y[,2]==1,]

    ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
    ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)

    train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
    #train <- rbind(train, X0)
    test <- rbind(X.0[ind.0,], X.1[ind.1,])

    Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
    #Y.train.loop <- rbind(Y.train.loop, Y0)
    Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])

    model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
    nfolds=nrow(train)
    model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)

    c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
    ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))

    rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)

    return(c(c,ct))
  }
  stopCluster(cluster.cores)
  closeAllConnections()
  rm(cluster.cores)
  rm(cores)
  gc()

  df_C_index <- do.call(rbind, list_C_index)
  df_C_index <- data.frame(df_C_index)
  colnames(df_C_index) <- c("C_index_test", "C_index_train")
  df_C_index$Iteration <- 1:iter
  df_C_index$ID <- i
  df_C_index$Data <- paste0("Clinical/genetic_C", nrow(Z_data_extra))

  df_C.index.C220 <- rbind(df_C.index.C220, df_C_index)

  cat("\nAnalysis completed for: ", i,"\n")
  rm(X,Y,  df_C_index, list_C_index, opts, pb, p)


}
boxplot(df_C.index.C220$C_index_test)
boxplot(df_C.index.C220$C_index_train)

list_res_Cox_combination <- Cox_forecasting_glmnet_combination(X_data=X.AUC.raw,
                                                               y_data=data.matrix(Y_C75),
                                                               Z_data=data.matrix(Z),
                                                               alpha=c(0, 0.4, 1),
                                                               lambda=exp(seq(-4,6, 0.1)),
                                                               free_cores = 2,
                                                               test.n= c(6,4),
                                                               iter=200,
                                                               log_AUC=1,
                                                               Patient.Z=1,
                                                               Drug.Z =2,
                                                               RCPC=0)

list_res_Cox_combination_sd0.15 <- Cox_forecasting_glmnet_combination(X_data=X.AUC.raw[,colSds(data.matrix(X.AUC.raw), na.rm = T)>0.15],
                                                                      y_data=data.matrix(Y_C75),
                                                                      Z_data=data.matrix(Z),
                                                                      alpha=c(0, 0.4, 1),
                                                                      lambda=exp(seq(-4,6, 0.1)),
                                                                      free_cores = 2,
                                                                      test.n= c(6,4),
                                                                      iter=200,
                                                                      log_AUC=1,
                                                                      Patient.Z=1,
                                                                      Drug.Z =2,
                                                                      RCPC=0)


save(df_C.index.C220,list_res_Cox_combination, list_res_Cox_combination_sd0.15,
     file = paste0(path,"/Results 20220711/Clinical data and Elnet comp_with variable pre-selection from AUC_Nor.RData"))



df_C.index.combination <- bind_rows(df_C.index.C220 %>%
                                      mutate(Data=paste("",Data),Pre_selection = "") %>%
                                      mutate(Data=gsub("_C220", " (C220)", Data),Pre_selection = ""),
                                    list_res_Cox_combination$C_index_results %>%
                                      mutate(Pre_selection = ""),
                                    list_res_Cox_combination_sd0.15$C_index_results %>%
                                      subset(grepl("Drug|Combined",Data)) %>%
                                      mutate(Pre_selection = "Pre-selection | "))


df_C.index.combination$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient stdz", "")
df_C.index.combination$Drug_stdz <- ifelse(grepl("Drug", df_C.index.combination$ID), "Drug stdz", "")
df_C.index.combination$Penalty <- gsub(".*_","", df_C.index.combination$ID)
df_C.index.combination$Metric <- gsub("\\/.*","", df_C.index.combination$ID)
df_C.index.combination$Metric <- gsub("AUC_raw"," AUC_raw",df_C.index.combination$Metric)
df_C.index.combination$Metric <- gsub("log2"," log2",df_C.index.combination$Metric, fixed=TRUE)
df_C.index.stat.combination <- df_C.index.combination %>% group_by(ID, Pre_selection,Data, Penalty) %>% summarise(Median=median(C_index_test),
                                                                                                                  Mean=mean(C_index_test),
                                                                                                                  MAD=mad(C_index_test),
                                                                                                                  SD=sd(C_index_test))
df_C.index.combination$Penalty[df_C.index.combination$Penalty==0] <- " Ridge"
df_C.index.combination$Penalty[df_C.index.combination$Penalty==0.4] <- "Elastic net"
df_C.index.combination$Penalty[df_C.index.combination$Penalty==1] <- "Lasso"

df_C.index.combination_m <- reshape2::melt(data.frame(df_C.index.combination), id.vars=c("Pre_selection", "Data", "Penalty"), measure.vars = c("C_index_test", "C_index_train"))
df_C.index.combination_m$variable <- ifelse(grepl("test", df_C.index.combination_m$variable), "Test", "Training")

ggexport(
  ggplot(df_C.index.combination_m, aes(x=paste(Pre_selection, Data), y=value, fill=Penalty))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_boxplot(alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
    ylim(0,1)+
    facet_grid(variable~.)+
    theme_bw()+
    labs(x="", y="C-index")+
    scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  width = 4.5, height=4, filename = paste0(path,"/Results 20220711/Clinical data combination and models comparison.pdf")
)

#----

#----Iterative model reduction----
path <- "D:/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Drug sensitivity data and surivival data - NOR_20220711.RData"))

X <- X.AUC.raw[,]
list_res_iterations <- Cox_forecasting_drug_withdrawal_2(X[,],
                                                         data.matrix(Y),
                                                         Reduce=c("Iterative"),
                                                         alpha=0,
                                                         lambda=c(exp(seq(-4,6, 0.1))),
                                                         free_cores = 5,
                                                         test.n= c(6,4),
                                                         iter=50,
                                                         log_AUC=1,
                                                         Patient.Z=1,
                                                         Drug.Z =2,
                                                         RCPC=0,
                                                         path=paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722"))
save(list_res_iterations, file = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722.RData"))

list_red <- list()
for(i in 1:length(list_res_iterations)){
  list_red[[i]] <- list_res_iterations[[i]]$WD_optimization$Naive_reduction %>%
    subset(Data %in% list_res_iterations[[i]]$WD_optimization$Naive_reduction_summary$Data)
  list_red[[i]]$Iter <- i-1
}


df_reduction <- bind_rows(list_red)
df_reduction$Removal <- paste(df_reduction$Iter,df_reduction$Data)
df_reduction$Removal <- factor(df_reduction$Removal, levels=unique(df_reduction$Removal))


df_reduction_sum <- df_reduction %>%
  subset(!grepl("full", Data)) %>%
  group_by(Iter, Data, Removal) %>%
  summarise_if(is.numeric, list(Mean=mean,
                                Median=median,
                                SD=sd))
df_reduction_sum$Drug <- gsub(".*_","",df_reduction_sum$Data)
df_reduction_sum$Drug <- factor(df_reduction_sum$Drug, levels = unique(df_reduction_sum$Drug))
df_reduction_sum$Iters <- ifelse(1:180 <= 60, "1-60", "61-120")
df_reduction_sum$Iters[1:180 > 120] <- "121-180"

ggexport(ggplot(df_reduction_sum, aes(x=Drug, y=C_index_test_Mean))+
           geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
           geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
           geom_pointrange(aes(ymin=C_index_test_Mean-C_index_test_SD/sqrt(50)*1.96,ymax=C_index_test_Mean+C_index_test_SD/sqrt(50)*1.96),size=0.1)+
           ylim(0,1)+
           theme_bw()+
           facet_wrap(~Iters, nrow=3,ncol=1, scales = "free_x")+
           labs(title=paste0("Testing: n=",10), x="Removal", y="C-index")+
           theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 strip.background = element_blank()),
  width = 7, height = 7.5,
  filename = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722_drug removals_1.pdf"))

ggexport(ggplot(df_reduction_sum, aes(x=Drug, y=C_index_train_Mean))+
           geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
           geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
           geom_pointrange(aes(ymin=C_index_train_Mean-C_index_train_SD/sqrt(50)*1.96,ymax=C_index_train_Mean+C_index_train_SD/sqrt(50)*1.96),size=0.1)+
           ylim(0,1)+
           theme_bw()+
           facet_wrap(~Iters, nrow=3,ncol=1, scales = "free_x")+
           labs(title=paste0("Training: n=",74-10), x="Removal", y="C-index")+
           theme(axis.text.x = element_text(angle = 60, hjust = 1, size=6),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 strip.background = element_blank()),
         width = 7, height = 7.5,
         filename = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722_drug removals_2.pdf"))

ggexport(ggarrange(ggplot(df_reduction %>% subset(grepl("full", Data)), aes(x=factor(Iter), y=C_index_test, group=Iter))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(fill = "darkblue", alpha=0.3, outlier.size = NULL)+
                     geom_point(data=df_reduction %>%
                                  subset(grepl("full", Data)) %>%
                                  group_by(Removal,Iter,Data) %>%
                                  summarise_if(is.numeric, mean),
                                size=1)+
                     ylim(0,1)+
                     theme_bw()+
                     scale_x_discrete(breaks=c("0","60","120"))+
                     labs(title=paste0("Testing: n=",10), x="Iteration", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank()),
                   ggplot(df_reduction %>% subset(grepl("full", Data)), aes(x=factor(Iter), y=C_index_train, group=Iter))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(fill = "darkblue", alpha=0.3, outlier.size = NULL)+
                     geom_point(data=df_reduction  %>%
                                  subset(grepl("full", Data)) %>%
                                  group_by(Removal,Iter,Data) %>%
                                  summarise_if(is.numeric, mean),
                                size=1)+
                     ylim(0,1)+
                     theme_bw()+
                     scale_x_discrete(breaks=c("0","60","120"))+
                     labs(title=paste0("Training: n=",74-10), x="Iteration", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank()), nrow=2, ncol=1)
         ,
         width = 8, height = 5.5,
         filename = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722_model iterations_1.pdf"))


df_reduction_sum=df_reduction %>%
  subset(grepl("full", Data) & Iter %in% c(0,20,40,60,80,100,120,134)) %>%
  group_by(Removal,Iter,Data) %>%
  summarise_if(is.numeric, median)

ggexport(ggarrange(ggplot(df_reduction %>% subset(grepl("full", Data) & Iter %in% c(0,20,40,60,80,100,120,134)), aes(x=factor(Iter), y=C_index_test, group=Iter))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(fill = "darkblue", alpha=0.3, outlier.size = NULL)+
                     ylim(0,1)+
                     theme_bw()+
                     scale_x_discrete(breaks=c("0","60","120"))+
                     labs(title=paste0("Testing: n=",10), x="Iteration", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank()),
                   ggplot(df_reduction %>% subset(grepl("full", Data) & Iter %in% c(0,20,40,60,80,100,120,134)), aes(x=factor(Iter), y=C_index_train, group=Iter))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(fill = "darkblue", alpha=0.3, outlier.size = NULL)+
                     ylim(0,1)+
                     theme_bw()+
                     scale_x_discrete(breaks=c("0","60","120"))+
                     labs(title=paste0("Training: n=",74-10), x="Iteration", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, hjust = 1),
                           panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank()), nrow=1, ncol=2)
         ,
         width = 4, height = 3.5,
         filename = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722_model iterations_2.pdf"))


#----Coefficient bootstrapping----
path <- "/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Drug sensitivity data and surivival data - NOR_20220711.RData"))
laod(file = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722.RData"))

Results_Cox_model_reduction <-list_res_iterations


X_initial <- Results_Cox_model_reduction[[1]]$X_initial
X_iter_final <- Results_Cox_model_reduction[[length(Results_Cox_model_reduction)]]$X_final

list_X <- list(X_initial=X_initial,
               X_iter_final=X_iter_final)


list_res_X_boot_CV <- list()
for(m in names(list_X)){
  X <- list_X[[m]]
  list_res_X_boot_CV[[m]] <-  Cox_bootstrapping(X_data=X,
                                                y_data=data.matrix(Y),
                                                alpha=0,
                                                lambda=c(exp(seq(-4,6, 0.1))),
                                                pre.CV=FALSE,
                                                lambda_opt = 0,
                                                free_cores = 5,
                                                iter=200,
                                                log_AUC=2,
                                                Patient.Z=2,
                                                Drug.Z =2,
                                                RCPC=0)
}

save(list_X, list_res_X_boot_CV,
     file = paste0(path,"/Results 20220711/Drug sensitivity model reduction X data bootstrapping and testing_20220722.RData"))

#----Withdrawal and bootstrapping results----
path <- "/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Drug sensitivity data and surivival data - NOR_20220711.RData"))
laod(file = paste0(path,"/Results 20220711/Drug sensitivity model reduction_20220722.RData"))
load(file = paste0(path,"/Results 20220711/Drug sensitivity model reduction X data bootstrapping and testing_20220722.RData"))

df_WD_intial <- Results_Cox_model_reduction[[1]]$WD_initial$df_loss
df_WD_intial$SE  <- matrixStats::rowSds(t(Results_Cox_model_reduction[[1]]$WD_initial$mat.C.index.test-Results_Cox_model_reduction[[1]]$df_parent$C_index_test))/sqrt(50)
df_WD_intial$Lower_95CI <- df_WD_intial$C_loss_test - df_WD_intial$SE*1.96; df_WD_intial$Upper_95CI <- df_WD_intial$C_loss_test + df_WD_intial$SE*1.96
df_WD_final <- Results_Cox_model_reduction[[length(Results_Cox_model_reduction)]]$WD_initial$df_loss
df_WD_final$SE  <- matrixStats::rowSds(t(Results_Cox_model_reduction[[length(Results_Cox_model_reduction)]]$WD_initial$mat.C.index.test-Results_Cox_model_reduction[[length(Results_Cox_model_reduction)]]$df_parent$C_index_test))/sqrt(50)
df_WD_final$Lower_95CI <- df_WD_final$C_loss_test - df_WD_final$SE*1.96; df_WD_final$Upper_95CI <- df_WD_final$C_loss_test + df_WD_final$SE*1.96


df_betas_initial <- list_res_X_boot_CV[[1]]$Betas_statistics
df_betas_final <- list_res_X_boot_CV[[2]]$Betas_statistics



df_drugs <- data.frame(Drug=colnames(list_X$X_initial), Label= NA)
df_drugs$Label[df_drugs$Drug %in% c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")] <- "Standard treatment"
df_drugs$Label[df_drugs$Drug %in% c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
                                    "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")] <- "Alternate treament"
df_drugs <- na.omit(df_drugs)

df_loss <- df_WD_intial
df_betas <- df_betas_initial
df_loss$Rank <- rank(df_loss$C_loss_test,ties.method="random")


df_betas_lab <- subset(df_betas, abs(t.stat) > 1.64 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(abs(df_betas_lab$t.stat) > 1.96, "*", ".")
df_betas_lab$Sign[abs(df_betas_lab$t.stat)<1.64] <- ""
df_betas_lab$drug <- gsub("\\(.*|HCl|2HCl|sodium salt", "", df_betas_lab$drug)
df_loss_lab <- subset(df_loss, (C_loss_test/SE) < -1.64 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse((df_loss_lab$C_loss_test/df_loss_lab$SE) < -1.96, "*", ".")
df_loss_lab$Sign[(df_loss_lab$C_loss_test/df_loss_lab$SE)> -1.64] <- ""
df_loss_lab$Var <- gsub("\\(.*|HCl|2HCl|sodium salt", "", df_loss_lab$Var)



ggexport(ggarrange(ggplot(df_betas, aes(x=Rank, y=Mean))+
                     geom_errorbar(aes(ymin=Lower_95CI , ymax=Upper_95CI),size=0.1,alpha=1,col="#7eb593")+
                     geom_point(size=0.9,col="#7eb593")+
                     geom_point(data=df_betas_lab, aes(x=Rank, y=Mean, col=Label),size=1.2)+
                     geom_text_repel(data=subset(df_betas_lab, Mean > 0 ) ,

                                     aes(x=Rank, y=Mean, label=paste0(drug, Sign), col=Label),
                                     nudge_y      = 0.02,
                                     nudge_x=-60,
                                     direction    = "y",
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     angle        = 0,
                                     vjust        = 0.2,
                                     segment.size = 0.2,
                                     size=2.5)+
                     geom_text_repel(data=subset(df_betas_lab, Mean < 0 ) ,
                                     aes(x=Rank, y=Mean, label=paste0(drug, Sign), col=Label),
                                     nudge_y      = -0.06,
                                     nudge_x = 70,
                                     direction    = "y",
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     angle        = 0,
                                     vjust        = 0.7,
                                     segment.size = 0.2,
                                     size=2.5)+
                     scale_color_manual(values=c("black", "darkblue", "darkred"))+
                     labs(y= "Coefficient (95% CI)", x = "t-statistic rank", title="log(HR) association")+
                     theme_bw(),
                   ggplot(df_loss, aes(Rank,C_loss_test))+
                     geom_errorbar(aes(ymin=Lower_95CI , ymax=Upper_95CI),size=0.1,alpha=1,col="#7eb593")+
                     geom_point(alpha=0.7, col="#7eb593",size=0.9)+
                     geom_point(data=df_loss_lab,
                                aes(Rank,C_loss_test, col=Label),size=1.2,alpha=1)+
                     geom_text_repel(data=df_loss_lab,
                                     aes(Rank,C_loss_test,label=paste0(Var, Sign), col=Label),
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     size=2.5,
                                     segment.size = 0.2,
                                     nudge_y      = -0.01,
                                     nudge_x = 60,
                                     direction    = "y")+
                     scale_color_manual(values=c("black", "darkblue", "darkred"))+
                     labs(x="Rank",y="C-index difference (95% CI)", title=paste0("Drug withdrawal"))+
                     #ylim(min(df_loss$C_loss_test),max(df_loss$C_loss_test))+
                     theme_bw(), common.legend = T),
         width = 8, height = 4,
         filename = paste0(path,"/Results 20220711/Drug survival association and prediction importance_full matrix_labeled.pdf"))




df_loss <- df_WD_final
df_betas <- df_betas_final
df_loss$Rank <- rank(df_loss$C_loss_test,ties.method="random")


df_betas_lab <- subset(df_betas, abs(t.stat) > 1.64 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(abs(df_betas_lab$t.stat) > 1.96, "*", ".")
df_betas_lab$Sign[abs(df_betas_lab$t.stat)<1.64] <- ""
df_betas_lab$drug <- gsub("\\(.*|HCl|2HCl|sodium salt", "", df_betas_lab$drug)

df_loss_lab <- subset(df_loss, (C_loss_test/SE) < -2.58 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse((df_loss_lab$C_loss_test/df_loss_lab$SE) < -1.96, "**", "*")
df_loss_lab$Var <- gsub("\\(.*|HCl|2HCl|sodium salt", "", df_loss_lab$Var)



ggexport(ggarrange(ggplot(df_betas, aes(x=Rank, y=Mean))+
                     geom_errorbar(aes(ymin=Lower_95CI , ymax=Upper_95CI),size=0.1,alpha=1,col="#7eb593")+
                     geom_point(size=0.9,col="#7eb593")+
                     geom_point(data=df_betas_lab, aes(x=Rank, y=Mean, col=Label),size=1.2)+
                     geom_text_repel(data=subset(df_betas_lab, Mean > 0 ) ,

                                     aes(x=Rank, y=Mean, label=paste0(drug, Sign), col=Label),
                                     nudge_y      = 0.02,
                                     nudge_x=-10,
                                     direction    = "y",
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     angle        = 0,
                                     vjust        = 0.1,
                                     segment.size = 0.2,
                                     size=2.5)+
                     geom_text_repel(data=subset(df_betas_lab, Mean < 0 ) ,
                                     aes(x=Rank, y=Mean, label=paste0(drug, Sign), col=Label),
                                     nudge_y      = -0.04,
                                     nudge_x = 12,
                                     direction    = "y",
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     angle        = 0,
                                     vjust        = 0.2,
                                     segment.size = 0.2,
                                     size=2.5)+
                     scale_color_manual(values=c("black", "darkblue", "darkred"))+
                     labs(y= "Coefficient (95% CI)", x = "t-statistic rank", title="log(HR) association")+
                     theme_bw(),
                   ggplot(df_loss, aes(Rank,C_loss_test))+
                     geom_errorbar(aes(ymin=Lower_95CI , ymax=Upper_95CI),size=0.1,alpha=1,col="#7eb593")+
                     geom_point(alpha=0.7, col="#7eb593",size=0.9)+
                     geom_point(data=df_loss_lab,
                                aes(Rank,C_loss_test, col=Label),size=1.2,alpha=1)+
                     geom_text_repel(data=df_loss_lab,
                                     aes(Rank,C_loss_test,label=paste0(Var, Sign), col=Label),
                                     max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                                     size=2.5,
                                     segment.size = 0.2,
                                     nudge_y      = -0.01,
                                     nudge_x = 10,
                                     direction    = "y")+
                     scale_color_manual(values=c("black", "darkblue", "darkred"))+
                     labs(x="Rank",y="C-index difference (95% CI)", title=paste0("Drug withdrawal"))+
                     #ylim(min(df_loss$C_loss_test),max(df_loss$C_loss_test))+
                     theme_bw(), common.legend = T),
         width = 8, height = 4,
         filename = paste0(path,"/Results 20220711/Drug survival association and prediction importance_iterative reduction_labeled.pdf"))

#----

#----Prepare association matrix----
path <- "/Aram/Drug screen - clinical forecasting"
load(file = paste0(path,"/Results 20220711/Drug sensitivity data and surivival data - NOR_20220711.RData"))
load(file = paste0(path,"/Results 20220711/Drug sensitivity model reduction X data bootstrapping and testing_20220722.RData"))

#Naively reduced matrix
X <- list_X$X_initial
X_initial <- X

set.seed(1)
model <- glmnet(X, data.matrix(Y), family = "cox", alpha = 0, standardize = FALSE,lambda=c(exp(seq(-4,6, 0.1))),  type.measure = "deviance")
nfolds=nrow(X)
model.cv <- cv.glmnet(X, data.matrix(Y), family = "cox", alpha = 0, standardize = FALSE,lambda=c(exp(seq(-4,6, 0.1))),  type.measure = "deviance", nfolds=nfolds)
model.coefs <- coefficients(model, s=model.cv$lambda.min)

X_association_initial <- X %*% diag(as.numeric(model.coefs))
colnames(X_association_initial) <- colnames(X)
colnames(X_initial) <- colnames(X)


#Iterative reduced matrix
X_iter <- X_initial[,colnames(X_initial) %in% colnames(list_X$X_iter_final)]

set.seed(1)
model <- glmnet(X_iter, data.matrix(Y), family = "cox", alpha = 0, standardize = FALSE,lambda=c(exp(seq(-4,6, 0.1))),  type.measure = "deviance")
nfolds=nrow(X_iter)
model.cv <- cv.glmnet(X_iter, data.matrix(Y), family = "cox", alpha = 0, standardize = FALSE,lambda=c(exp(seq(-4,6, 0.1))),  type.measure = "deviance", nfolds=nfolds)
model.coefs <- coefficients(model, s=model.cv$lambda.min)

X_association_iter <- X_iter %*% diag(as.numeric(model.coefs))
colnames(X_association_iter) <- colnames(X_iter)

save(X_association_iter, X_association_initial,
     X_iter, X_initial,
     Y, file=paste0(path,"/Results 20220711/Drug sensitivity survival association matrices.RData"))


library("devtools")
library("heatmap3")
library("gplots")
library("RColorBrewer")
library("pheatmap")
library("tagcloud")
library("grid")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
load(file=paste0(path,"/Results 20220711/Clinical data_curated.RDAta"))

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels,
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)),
                    length.out = sum(d.select)),
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4,
                                     l = 4
  )

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # plot result
  grid.newpage()
  grid.draw(heatmap)

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

#----Heatmap color labels----
#predicted risk, death and gender

pred <- predict(model, newx=data.matrix(X_iter), s=model.cv$lambda.min, type="response")
#HR <- rank(pred)
#HR <- viridis::viridis_pal(option = "B")(length(HR))[HR]
HR <- smoothPalette(log(pred), pal= viridis::viridis_pal(option = "B")(length(pred)) )
names(HR) <- rownames(X_iter)
HR <- HR[match(clinical_74$PatientCode_JElist, names(HR))]

survival.legend <- c( "#946877","#C0D0DF")
names(survival.legend) <- c("Alive","Dead")
survival <- clinical_74$death
names(survival) <- paste0(clinical_74$PatientCode_JElist)
survival[survival=="0" & clinical_74$time>365] <- survival.legend[1]
survival[survival=="1"] <- survival.legend[2]
survival[which(!(survival %in% survival.legend))] <- "white"

gender.legend <- c( "#e0c8a2","#70aae6")
names(gender.legend) <- c("Female","Male")
gender <- clinical_74$Gender
names(gender) <- paste0(clinical_74$PatientCode_JElist)
gender[gender=="0"] <- gender.legend[1]
gender[gender=="1"] <- gender.legend[2]
gender[which(!(gender %in% gender.legend))] <- "white"

# Prognosis and classifications

WHO_class.legend <- c(brewer.pal(3, "Set3"))
names(WHO_class.legend) <- c("WHO class 0","WHO class 1", "WHO class 2")
WHO_class<-clinical_74$WHO_class_0
names(WHO_class)<-paste0(clinical_74$PatientCode_JElist)
WHO_class[clinical_74$WHO_class_0=="1"] <- WHO_class.legend[1]
WHO_class[clinical_74$WHO_class_1=="1"] <- WHO_class.legend[2]
WHO_class[clinical_74$WHO_class_2=="1"] <- WHO_class.legend[3]
WHO_class[which(!(WHO_class %in% WHO_class.legend))] <- "white"

AML_fab.legend <- c(brewer.pal(5, "Pastel2"))
names(AML_fab.legend) <- c("AML fab class 0","AML fab class 1 to 2", "AML fab class 3", "AML fab class 4 to 5", "AML fab class 6")
AML_fab<-clinical_74$AMLFAB_class_AML_zero
names(AML_fab)<-paste0(clinical_74$PatientCode_JElist)
AML_fab[clinical_74$AMLFAB_class_AML_zero=="1"] <- AML_fab.legend[1]
AML_fab[clinical_74$AMLFAB_class_AML_one_to_two=="1"] <- AML_fab.legend[2]
AML_fab[clinical_74$AMLFAB_class_AML_three=="1"] <- AML_fab.legend[3]
AML_fab[clinical_74$AMLFAB_class_AML_four_to_five=="1"] <- AML_fab.legend[4]
AML_fab[clinical_74$AMLFAB_class_AML_six=="1"] <- AML_fab.legend[5]
AML_fab[which(!(AML_fab %in% AML_fab.legend))] <- "white"

Prognosis.legend <- c(brewer.pal(3, "Pastel1"))
names(Prognosis.legend) <- c("Average","Bad", "Good")
Prognosis<-clinical_74$WHO_class_0
names(Prognosis)<-paste0(clinical_74$PatientCode_JElist)
Prognosis[clinical_74$Cytogenetics_Average_Prognosis=="1"] <- Prognosis.legend[1]
Prognosis[clinical_74$Cytogenetics_Bad_Prognosis=="1"] <- Prognosis.legend[2]
Prognosis[clinical_74$Cytogenetics_Good_Prognosis=="1"] <- Prognosis.legend[3]
Prognosis[which(!(Prognosis %in% Prognosis.legend))] <- "white"

# Mutations
mut.legend <- c( "#ff8c00","white")
names(mut.legend) <- c("Positive","Negative")

NPM1 <- clinical_74$NPM1
names(NPM1) <- paste0(clinical_74$PatientCode_JElist)
NPM1[NPM1=="1"] <- mut.legend[1]
NPM1[NPM1=="0"] <- mut.legend[2]
NPM1[which(!(NPM1 %in% mut.legend))] <- "white"

FLT3ITD <- clinical_74$FLT3ITD
names(FLT3ITD) <- paste0(clinical_74$PatientCode_JElist)
FLT3ITD[FLT3ITD=="1"] <- mut.legend[1]
FLT3ITD[FLT3ITD=="0"] <- mut.legend[2]
FLT3ITD[which(!(FLT3ITD %in% mut.legend))] <- "white"

CEBPA <- clinical_74$CEBPA
names(CEBPA) <- paste0(clinical_74$PatientCode_JElist)
CEBPA[CEBPA=="1"] <- mut.legend[1]
CEBPA[CEBPA=="0"] <- mut.legend[2]
CEBPA[which(!(CEBPA %in% mut.legend))] <- "white"

EVI1 <- clinical_74$EVI1
names(EVI1) <- paste0(clinical_74$PatientCode_JElist)
EVI1[EVI1=="1"] <- mut.legend[1]
EVI1[EVI1=="0"] <- mut.legend[2]
EVI1[which(!(EVI1 %in% mut.legend))] <- "white"

RUNX1 <- clinical_74$RUNX1
names(RUNX1) <- paste0(clinical_74$PatientCode_JElist)
RUNX1[RUNX1=="1"] <- mut.legend[1]
RUNX1[RUNX1=="0"] <- mut.legend[2]
RUNX1[which(!(RUNX1 %in% mut.legend))] <- "white"

ASXL1 <- clinical_74$ASXL1
names(ASXL1) <- paste0(clinical_74$PatientCode_JElist)
ASXL1[ASXL1=="1"] <- mut.legend[1]
ASXL1[ASXL1=="0"] <- mut.legend[2]
ASXL1[which(!(ASXL1 %in% mut.legend))] <- "white"

PMLRARA <- clinical_74$PMLRARA
names(PMLRARA) <- paste0(clinical_74$PatientCode_JElist)
PMLRARA[PMLRARA=="1"] <- mut.legend[1]
PMLRARA[PMLRARA=="0"] <- mut.legend[2]
PMLRARA[which(!(PMLRARA %in% mut.legend))] <- "white"

FLT3TKD <- clinical_74$FLT3TKD
names(FLT3TKD) <- paste0(clinical_74$PatientCode_JElist)
FLT3TKD[FLT3TKD=="1"] <- mut.legend[1]
FLT3TKD[FLT3TKD=="0"] <- mut.legend[2]
FLT3TKD[which(!(FLT3TKD %in% mut.legend))] <- "white"

CBFMYH11 <- clinical_74$CBFMYH11
names(CBFMYH11) <- paste0(clinical_74$PatientCode_JElist)
CBFMYH11[CBFMYH11=="1"] <- mut.legend[1]
CBFMYH11[CBFMYH11=="0"] <- mut.legend[2]
CBFMYH11[which(!(CBFMYH11 %in% mut.legend))] <- "white"

inv16 <- clinical_74$inv16
names(inv16) <- paste0(clinical_74$PatientCode_JElist)
inv16[inv16=="1"] <- mut.legend[1]
inv16[inv16=="0"] <- mut.legend[2]
inv16[which(!(inv16 %in% mut.legend))] <- "white"

bcrabl <- clinical_74$bcrabl
names(bcrabl) <- paste0(clinical_74$PatientCode_JElist)
bcrabl[bcrabl=="1"] <- mut.legend[1]
bcrabl[bcrabl=="0"] <- mut.legend[2]
bcrabl[which(!(bcrabl %in% mut.legend))] <- "white"


#----Heatmap----
clab=cbind(NPM1, FLT3ITD, CEBPA, EVI1, RUNX1, ASXL1, PMLRARA, FLT3TKD, CBFMYH11, inv16, BCR_ABL=bcrabl, WHO_class, AML_fab,
           Prognosis, Gender=gender, Survival=survival, log_HR=HR)


df_column_col <- clab[match(rownames(X_association_iter),rownames(clab)),]
my_colour=list(NPM1=mut.legend, FLT3ITD=mut.legend, CEBPA=mut.legend, EVI1=mut.legend,
               RUNX1=mut.legend, ASXL1=mut.legend, PMLRARA=mut.legend, FLT3TKD=mut.legend,
               CBFMYH11=mut.legend, inv16=mut.legend, BCR_ABL=mut.legend,
               WHO_class=WHO_class.legend, AML_fab=AML_fab.legend,Prognosis=Prognosis.legend, Gender=gender.legend,
               Survival=survival.legend, log_HR=viridis::viridis_pal(option = "B")(length(HR)))
colnames(df_column_col)
df_column_col[,1:11] <- apply(df_column_col[,1:11],2,function(x) names(mut.legend)[match(x,mut.legend)])
df_column_col <- data.frame(df_column_col)
df_column_col$WHO_class <- names(WHO_class.legend)[match(df_column_col$WHO_class,WHO_class.legend)]
df_column_col$AML_fab <- names(AML_fab.legend)[match(df_column_col$AML_fab,AML_fab.legend)]
df_column_col$Prognosis <- names(Prognosis.legend)[match(df_column_col$Prognosis,Prognosis.legend)]
df_column_col$Gender <- names(gender.legend)[match(df_column_col$Gender,gender.legend)]
df_column_col$Survival <- names(survival.legend)[match(df_column_col$Survival,survival.legend)]
df_column_col$log_HR <- unlist(log(as.numeric(pred)))
df_column_col$log_HR <- (df_column_col$log_HR- mean(df_column_col$log_HR))/sd(df_column_col$log_HR)
rownames(df_column_col) <- paste0("Patient", rownames(df_column_col))
X_association_iter.t <- t(X_association_iter)
colnames(X_association_iter.t) <- paste0("Patient", colnames(X_association_iter.t))


HC_r <- hclust(dist(t(X_association_iter), method="euclidean"), method="ward.D")
HC_c <- hclust(dist((X_association_iter), method="maximum"), method="ward.D")


heat <- pheatmap(X_association_iter.t[],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)),#bluered(100),
                 breaks=seq(-2.5,2.5,length.out= 100),
                 cutree_cols = 7,

                 labels_col = rep("", ncol(X_association_iter)),
                 annotation_col = df_column_col[,],

                 annotation_colors = my_colour,

                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)


add.flag(heat,
         kept.labels = colnames(X_association_iter[,matrixStats::colSds(data.matrix(X_association_iter))>0.2]),
         repel.degree = 1)



pdf(width=8.5,height=7,
    file = file = paste0(path,"/Results 20220711/Patient risk group clustering_heatmap_cut_1.pdf"))
add.flag(heat,
         kept.labels = colnames(X_association_iter[,matrixStats::colSds(data.matrix(X_association_iter))>0.2]),
         repel.degree = 1)
dev.off()


X_iter.t <- t(X_iter)
colnames(X_iter.t) <- paste0("Patient", colnames(X_iter.t))
heat <- pheatmap(X_iter.t,
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)),#bluered(100),
                 breaks=seq(-2.5,2.5,length.out= 100),
                 cutree_cols = 7,

                 labels_col = rep("", ncol(X_association_iter)),
                 annotation_col = df_column_col[,],

                 annotation_colors = my_colour,

                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)


add.flag(heat,
         kept.labels = colnames(X_association_iter[,matrixStats::colSds(data.matrix(X_association_iter))>0.2]),
         repel.degree = 1)

#----

#----Kaplan-Meyer plots of drug sensitivity groups----
library(survival)
library(ggfortify)
library(survminer)


hits <- colnames(X_initial)[grepl("Daunorubicin|Idarubicin|Cytarabine|Navitoclax|Geldanamycin|YM155|Mitoxantrone|Triptolide|Cyt387|Estrone|Ganetespib|CEP33779|Ubenimex|Mesna|ABT???888|MLN2238",colnames(X_initial))]
for(Q in c(0.25,0.5,0.75)){
  #Q = 0.75
  X.hits <- X_initial[ ,hits]
  X.hits <- melt(X.hits)
  X.hits <- X.hits %>%
    group_by(Var2) %>%
    #mutate(k = kmeans(value, centers=2, iter.max = 100, algorithm = "MacQueen")$cluster)
    #mutate(k = kmod(data.frame(x=value, y=1), k = 2, l = 6, i_max = 100, conv_method = "delta_C")$XC_dist_sqr_assign[,2])
    mutate(rank = rank(value)/74) %>%
    mutate(k = ifelse(rank>Q, 1, 0))
  X.hits <- X.hits %>%
    group_by(Var2, k) %>%
    mutate(k.avg = mean(value))
  X.hits <- X.hits %>%
    group_by(Var2) %>%
    mutate(Sensitive = ifelse(k.avg - min(k.avg) > 0, 1,0 ))
  X.hits$Sensitive <- as.character(X.hits$Sensitive)

  ggexport(
    ggplot(X.hits %>% mutate(Drug = Var2) %>% mutate(Drug = factor(Drug)))+
      geom_histogram(aes(value), bins =20, fill="white", col="black")+
      geom_point(aes(x=value, y=-2.5, col=Sensitive), size=0.5)+
      facet_wrap(~Drug, scales="free_x")+
      labs(x="-log2(AUC) ")+
      scale_color_manual(values=c("#999999", "#E69F00"))+
      theme_minimal()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks = element_line(),
            axis.text.x = element_text(angle = 45, hjust = 1)),
    width = 8.5, height = 8.5,
    filename = paste0(path,"/Results 20220711/Patient risk group stratification_AUC distribution_q",Q,".pdf")
  )


  plot.list<-list()
  for(i in 1:length(hits)){
    fit.drug <- survfit(Surv(time, death) ~ drug.response, data = Y %>%
                          mutate(death=status,drug.response = X.hits[which(X.hits$Var2 == hits[i]),]$Sensitive=="1") )

    plot.list[[i]] <-  autoplot(fit.drug,  conf.int = F)+
      theme_minimal()+
      labs(title= paste(hits[i],"\np-value =",round(surv_pvalue(fit.drug)[2],3)), x="Time (days)", y="Survival", colour = paste0("Top ",(1-Q)*100,"%"," responders"))+
      ylim(0,1)+
      scale_y_continuous(labels = scales::percent, limits=c(0,1))+
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks = element_line(),
            plot.title = element_text(size=8),
            axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8))
  }
  ggexport(
    ggarrange(plotlist=plot.list, common.legend = TRUE),
  width = 8.5, height = 8.5, filename = paste0(path,"/Results 20220711/Patient risk group stratification_Kaplan-Meyer plot_q",Q,".pdf")
  )

}

#----


