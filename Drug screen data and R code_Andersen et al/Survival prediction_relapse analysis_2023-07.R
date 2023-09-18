library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(survival)
library(survminer)
library(CoxTools)
library(glmnet)
library(glmnetUtils)
library(doSNOW)
path <- "C:/Users/EnserinkLab2019/Desktop/Clinical forecasting - revision 2023-07"

# Load ----
load(paste0(path, "/Drug sensitivity QC results and batch data_2023-07.RData"))
load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(paste0(path, "/Survival, clinical features and ELN2022 classifications_2023-07.RData"))
df_scores <- df_scores %>% subset(!grepl("remission", Patient.ID))

X_rAUC <- dcast(data = df_scores %>% subset(!grepl("re", Patient.ID)), Patient.num ~ drug , value.var="rAUC")
rownames(X_rAUC) <- as.character(X_rAUC$Patient.num); X_rAUC$Patient.num <- NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]
SDs <- matrixStats::colSds(data.matrix(X_rAUC))

df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$rAUC <- 1-df_scores$rAUC

df_scores_re <- df_scores %>% subset(Patient.num %in% df_scores$Patient.num[grepl("re", df_scores$Patient.ID)])
df_scores1 <- df_scores %>% subset(!(Patient.num %in% df_scores$Patient.num[grepl("re", df_scores$Patient.ID)]))

df_relapse <- df_batch %>% subset(Sample_ID %in% df_scores_re$Patient.ID & grepl("re",Sample_ID))
df_relapse <- cbind(df_relapse[,1:5], df_survival[match(df_relapse$ID,df_survival$ID),12:17])
df_relapse$Relapse_time <- difftime(df_relapse$Sample_date, df_relapse$`Date of diagnosis`, units = "days")
df_relapse$Relapse_time <- round(as.numeric(df_relapse$Relapse_time))

ggplot(df_relapse %>% mutate(Rank = rank(-as.numeric(Time))),aes(Rank, Time))+
  geom_linerange(aes(ymin=0, ymax=Time))+
  geom_point(aes(y=Relapse_time),pch=17, size=3, col=pal_jco()(2)[2])+
  geom_point(aes(y=0),pch=16, size=3, col=pal_jco()(2)[1])+
  geom_point(pch=18, size=3)+
  labs(x="Patients")+
  theme_classic()+
  coord_flip()
ggsave(filename = paste0(path,"/Results_relapse analysis/Relapse times.pdf"), 
         width=3, height=2.75)

# Relapse prediction----
df_C.index_relapse <- c()
df_CoxPH_relapse <- c()
for(m in c("DSS3","rAUC_log2")){
  for(pt_stz in c(TRUE, FALSE)){
    X <- dcast(data = df_scores, Patient.ID ~ drug , value.var=m)
    rownames(X) <-X$Patient.ID; X$Patient.ID<-NULL

    if(pt_stz){X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))}
    
    X_re <- X[which(rownames(X) %in% df_scores_re$Patient.ID),]
    X <- X[-which(rownames(X) %in% df_scores_re$Patient.ID),]
    
    X_re0 <- X_re[-grep("re", rownames(X_re)),]
    X_re <- X_re[grep("re", rownames(X_re)),]


    Y_re <- Y[match(rownames(X_re0), paste("Patient",rownames(Y))),]
    Y1 <- Y[match(gsub(".* ","",rownames(X)), rownames(Y)),]
    
    for(tr in c(quantile(SDs, 1-c(40, 100)/ncol(X_rAUC)),0)){
      X1 <- X[,which(SDs>=tr)]
      for(a in c(0,0.4,1)){
        set.seed(1)
        model <- glmnet(data.matrix(X1), data.matrix(Y1), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
        nfolds=nrow(X1)
        model.cv <- cv.glmnet(data.matrix(X1), data.matrix(Y1),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
        
        
        reference <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X1),type = "response")
        #y <- cbind(predict(model, s = model.cv$lambda.min, newx= data.matrix(X1[,which(SDs>=tr)]),type = "response"), Y1)
        pr0 <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re0[,which(SDs>=tr)]),type = "response")
        pr <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re[,which(SDs>=tr)]),type = "response")
        c0=Cindex(pr0, y=data.matrix(Y_re))
        c=Cindex(pr, y=data.matrix(Y_re))
        df_C.index_relapse <- rbind(df_C.index_relapse,
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", C_index=c0),
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", C_index=c))
        df_CoxPH_relapse <- rbind(df_CoxPH_relapse,
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", 
                                             HR=pr0[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)),
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", 
                                             HR=pr[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)))
      }
    }
  }
}



df_CoxPH_relapse <- df_CoxPH_relapse %>% 
  group_by(Model, Metric, Standardize, Pre_selection_N) %>% 
  mutate(pval=pairwise.wilcox.test(HR,Sample,paired = TRUE)[[3]])
df_CoxPH_relapse$pval <- ifelse(df_CoxPH_relapse$pval<0.05,
                                       scales::pvalue(df_CoxPH_relapse$pval ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")

ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_HR of median_rAUC_1.pdf"), 
  width=12, height=6)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_HR of mean_rAUC_1.pdf"), 
  width=12, height=6)



ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_HR of median_DSS3_1.pdf"), 
  width=12, height=6)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_HR of mean_DSS3_1.pdf"), 
  width=12, height=6)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Pre_selection_N), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_bar(alpha=1, outlier.size = NULL, outlier.alpha = 0, stat="identity", position=position_dodge())+
    facet_grid(ifelse(Standardize, "Patient standardized", "")~Model, scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_C-index_1.pdf"), 
  width=7, height=4)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Sample), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1, outlier.alpha = 0, position=position_dodge())+
    geom_line(aes(group=paste(Model, Pre_selection_N, Metric)), col="gray")+
    geom_point(aes(pch=factor(Model), size=factor(Pre_selection_N), group=Sample),alpha=0.5)+
    scale_size_manual(values=rev(c(0.5,0.75,1)*3))+
    facet_grid(~ifelse(Standardize, "z-score", ""), scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_C-index_2.pdf"), 
  width=5.5, height=3.25)


# Relapse prediction - sts----
Y2 <- Y
Y2$status[Y2$time>365*2] <- 0
Y2$time[Y2$time>365*2] <- 365*2

df_C.index_relapse <- c()
df_CoxPH_relapse <- c()
for(m in c("DSS3","rAUC_log2")){
  for(pt_stz in c(TRUE, FALSE)){
    X <- dcast(data = df_scores, Patient.ID ~ drug , value.var=m)
    rownames(X) <-X$Patient.ID; X$Patient.ID<-NULL
    
    if(pt_stz){X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))}
    
    X_re <- X[which(rownames(X) %in% df_scores_re$Patient.ID),]
    X <- X[-which(rownames(X) %in% df_scores_re$Patient.ID),]
    
    X_re0 <- X_re[-grep("re", rownames(X_re)),]
    X_re <- X_re[grep("re", rownames(X_re)),]
    
    
    Y_re <- Y2[match(rownames(X_re0), paste("Patient",rownames(Y2))),]
    Y1 <- Y2[match(gsub(".* ","",rownames(X)), rownames(Y2)),]
    
    for(tr in c(quantile(SDs, 1-c(40, 100)/ncol(X_rAUC)),0)){
      X1 <- X[,which(SDs>=tr)]
      for(a in c(0,0.4,1)){
        set.seed(1)
        model <- glmnet(data.matrix(X1), data.matrix(Y1), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
        nfolds=nrow(X1)
        model.cv <- cv.glmnet(data.matrix(X1), data.matrix(Y1),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
        
        
        reference <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X1),type = "response")
        #y <- cbind(predict(model, s = model.cv$lambda.min, newx= data.matrix(X1[,which(SDs>=tr)]),type = "response"), Y1)
        pr0 <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re0[,which(SDs>=tr)]),type = "response")
        pr <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re[,which(SDs>=tr)]),type = "response")
        c0=Cindex(pr0, y=data.matrix(Y_re))
        c=Cindex(pr, y=data.matrix(Y_re))
        df_C.index_relapse <- rbind(df_C.index_relapse,
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", C_index=c0),
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", C_index=c))
        df_CoxPH_relapse <- rbind(df_CoxPH_relapse,
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", 
                                             HR=pr0[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)),
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", 
                                             HR=pr[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)))
      }
    }
  }
}

df_CoxPH_relapse$time0 <- Y$time[match(df_CoxPH_relapse$ID,rownames(Y))]

df_CoxPH_relapse <- df_CoxPH_relapse %>% 
  group_by(Model, Metric, Standardize, Pre_selection_N) %>% 
  mutate(pval=pairwise.wilcox.test(HR,Sample,paired = TRUE)[[3]])
df_CoxPH_relapse$pval <- ifelse(df_CoxPH_relapse$pval<0.05,
                                scales::pvalue(df_CoxPH_relapse$pval ,
                                               accuracy = 0.01,
                                               decimal.mark = ".", 
                                               add_p = TRUE),"ns")
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_HR of median_rAUC_1.pdf"), 
  width=12, height=6)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_HR of mean_rAUC_1.pdf"), 
  width=12, height=6)



ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_HR of median_DSS3_1.pdf"), 
  width=12, height=6)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    scale_y_continuous(trans="log2")+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N+pval, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_HR of mean_DSS3_1.pdf"), 
  width=12, height=6)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Pre_selection_N), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_bar(alpha=1, outlier.size = NULL, outlier.alpha = 0, stat="identity", position=position_dodge())+
    facet_grid(ifelse(Standardize, "Patient standardized", "")~Model, scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_C-index_1.pdf"), 
  width=7, height=4)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Sample), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1, outlier.alpha = 0, position=position_dodge())+
    geom_line(aes(group=paste(Model, Pre_selection_N, Metric)), col="gray")+
    geom_point(aes(pch=factor(Model), size=factor(Pre_selection_N), group=Sample),alpha=0.5)+
    scale_size_manual(values=rev(c(0.5,0.75,1)*3))+
    facet_grid(~ifelse(Standardize, "z-score", ""), scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_sts_C-index_2.pdf"), 
  width=5.5, height=3.25)

# Relapse prediction - treatment naive only----
df_C.index_relapse <- c()
df_CoxPH_relapse <- c()
for(m in c("DSS3","rAUC_log2")){
  for(pt_stz in c(TRUE, FALSE)){
    X <- dcast(data = df_scores %>% subset(Patient.num %in% df_survival$ID[which(df_survival$`Treatment Naive` == "Yes")]), Patient.ID ~ drug , value.var=m)
    rownames(X) <-X$Patient.ID; X$Patient.ID<-NULL
    
    if(pt_stz){X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))}
    
    X_re <- X[which(rownames(X) %in% df_scores_re$Patient.ID),]
    X <- X[-which(rownames(X) %in% df_scores_re$Patient.ID),]
    
    X_re0 <- X_re[-grep("re", rownames(X_re)),]
    X_re <- X_re[grep("re", rownames(X_re)),]
    
    
    Y_re <- Y[match(rownames(X_re0), paste("Patient",rownames(Y))),]
    Y1 <- Y[match(gsub(".* ","",rownames(X)), rownames(Y)),]
    
    for(tr in c(quantile(SDs, 1-c(40, 100)/ncol(X_rAUC)),0)){
      X1 <- X[,which(SDs>=tr)]
      for(a in c(0,0.4,1)){
        set.seed(1)
        model <- glmnet(data.matrix(X1), data.matrix(Y1), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
        nfolds=nrow(X1)
        model.cv <- cv.glmnet(data.matrix(X1), data.matrix(Y1),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
        
        
        reference <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X1),type = "response")
        #y <- cbind(predict(model, s = model.cv$lambda.min, newx= data.matrix(X1[,which(SDs>=tr)]),type = "response"), Y1)
        pr0 <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re0[,which(SDs>=tr)]),type = "response")
        pr <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re[,which(SDs>=tr)]),type = "response")
        c0=Cindex(pr0, y=data.matrix(Y_re))
        c=Cindex(pr, y=data.matrix(Y_re))
        df_C.index_relapse <- rbind(df_C.index_relapse,
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", C_index=c0),
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", C_index=c))
        df_CoxPH_relapse <- rbind(df_CoxPH_relapse,
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", 
                                             HR=pr0[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)),
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", 
                                             HR=pr[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)))
      }
    }
  }
}


ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_HR of median_rAUC_1.pdf"), 
  width=10, height=5)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_HR of mean_rAUC_1.pdf"), 
  width=10, height=5)



ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_HR of median_DSS3_1.pdf"), 
  width=10, height=5)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_HR of mean_DSS3_1.pdf"), 
  width=10, height=5)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Pre_selection_N), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_bar(alpha=1, outlier.size = NULL, outlier.alpha = 0, stat="identity", position=position_dodge())+
    facet_grid(ifelse(Standardize, "Patient standardized", "")~Model, scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_C-index_1.pdf"), 
  width=7, height=4)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Sample), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1, outlier.alpha = 0, position=position_dodge())+
    geom_line(aes(group=paste(Model, Pre_selection_N, Metric)), col="gray")+
    geom_point(aes(pch=factor(Model), size=factor(Pre_selection_N), group=Sample),alpha=0.5, pch=16)+
    scale_size_manual(values=c(0.5,0.75,1)*3)+
    facet_grid(~ifelse(Standardize, "z-score", ""), scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only_C-index_2.pdf"), 
  width=5.5, height=3.25)

# Relapse prediction - treatment naive only + sts----
Y2 <- Y
Y2$status[Y2$time>365*2] <- 0
Y2$time[Y2$time>365*2] <- 365*2

df_C.index_relapse <- c()
df_CoxPH_relapse <- c()
for(m in c("DSS3","rAUC_log2")){
  for(pt_stz in c(TRUE, FALSE)){
    X <- dcast(data = df_scores %>% subset(Patient.num %in% df_survival$ID[which(df_survival$`Treatment Naive` == "Yes")]), Patient.ID ~ drug , value.var=m)
    rownames(X) <-X$Patient.ID; X$Patient.ID<-NULL
    
    if(pt_stz){X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))}
    
    X_re <- X[which(rownames(X) %in% df_scores_re$Patient.ID),]
    X <- X[-which(rownames(X) %in% df_scores_re$Patient.ID),]
    
    X_re0 <- X_re[-grep("re", rownames(X_re)),]
    X_re <- X_re[grep("re", rownames(X_re)),]
    
    
    Y_re <- Y2[match(rownames(X_re0), paste("Patient",rownames(Y2))),]
    Y1 <- Y2[match(gsub(".* ","",rownames(X)), rownames(Y2)),]
    
    for(tr in c(quantile(SDs, 1-c(40, 100)/ncol(X_rAUC)),0)){
      X1 <- X[,which(SDs>=tr)]
      for(a in c(0,0.4,1)){
        set.seed(1)
        model <- glmnet(data.matrix(X1), data.matrix(Y1), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
        nfolds=nrow(X1)
        model.cv <- cv.glmnet(data.matrix(X1), data.matrix(Y1),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
        
        
        reference <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X1),type = "response")
        #y <- cbind(predict(model, s = model.cv$lambda.min, newx= data.matrix(X1[,which(SDs>=tr)]),type = "response"), Y1)
        pr0 <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re0[,which(SDs>=tr)]),type = "response")
        pr <- predict(model, s = model.cv$lambda.min, newx= data.matrix(X_re[,which(SDs>=tr)]),type = "response")
        c0=Cindex(pr0, y=data.matrix(Y_re))
        c=Cindex(pr, y=data.matrix(Y_re))
        df_C.index_relapse <- rbind(df_C.index_relapse,
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", C_index=c0),
                                    data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", C_index=c))
        df_CoxPH_relapse <- rbind(df_CoxPH_relapse,
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Treatment_naive", 
                                             HR=pr0[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)),
                                  data.frame(Metric = m,Model=a, Standardize=pt_stz, Pre_selection_N=ncol(X1), Sample="Relapsed", 
                                             HR=pr[,1],Ref_mean=mean(reference),Ref_mean_log=exp(mean(log(reference))), Ref_median=median(reference), 
                                             time=Y_re$time, status=Y_re$status, ID=rownames(Y_re)))
      }
    }
  }
}

df_CoxPH_relapse$time0 <- Y$time[match(df_CoxPH_relapse$ID,rownames(Y))]

ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_HR of median_rAUC_1.pdf"), 
  width=10, height=5)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="rAUC_log2") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_HR of mean_rAUC_1.pdf"), 
  width=10, height=5)



ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_median, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_HR of median_DSS3_1.pdf"), 
  width=10, height=5)
ggexport(
  ggplot(df_CoxPH_relapse %>% subset(Metric=="DSS3") %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample)), 
         aes(x=Sample, y=HR/Ref_mean_log, fill=Sample))+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_hline(yintercept = 1, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1,outlier.alpha = 0, width=0.5)+
    geom_line(aes(group=ID), col="gray")+
    geom_point(aes(size=as.numeric(time0), group=Sample), alpha=0.5, pch=16)+
    facet_wrap(ifelse(Standardize, "Patient standardized", "")~Model+Pre_selection_N, scales = "free_y", ncol=9, nrow=2)+
    scale_fill_jco()+
    scale_shape_manual(values=c(1,16))+
    theme_bw()+
    labs(x="", y="Hazard ratio")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_HR of mean_DSS3_1.pdf"), 
  width=10, height=5)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Pre_selection_N), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_bar(alpha=1, outlier.size = NULL, outlier.alpha = 0, stat="identity", position=position_dodge())+
    facet_grid(ifelse(Standardize, "Patient standardized", "")~Model, scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_C-index_1.pdf"), 
  width=7, height=4)

ggexport(
  ggplot(df_C.index_relapse %>% 
           mutate(Sample=gsub("Tr", " Tr", Sample),
                  Metric=gsub("rAUC"," rAUC",Metric),
                  Pre_selection_N = gsub("349"," 349",Pre_selection_N)), 
         aes(x=paste(Metric, Sample), y=C_index, fill=Sample))+
    geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
    geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
    geom_boxplot(alpha=1, outlier.alpha = 0, position=position_dodge())+
    geom_line(aes(group=paste(Model, Pre_selection_N, Metric)), col="gray")+
    geom_point(aes(pch=factor(Model), size=factor(Pre_selection_N), group=Sample),alpha=0.5, pch=16)+
    scale_size_manual(values=c(0.5,0.75,1)*3)+
    facet_grid(~ifelse(Standardize, "z-score", ""), scales = "free")+
    scale_fill_jco()+
    ylim(0,1)+
    theme_bw()+
    labs(x="", y="C-index")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank()),
  filename = paste0(path,"/Results_relapse analysis/Relapse sample predictions_naive only + sts_C-index_2.pdf"), 
  width=5.5, height=3.25)