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

# Path ----
path <- "/mnt/CommonStorageRAI/Aram/Drug screen - clinical forecasting/Clinical forecasting manuscript - revision 2023-07"

# Load ----
load(file = paste0(path,"/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(file = paste0(path,"/Survival, clinical features and ELN2022 classifications_2023-07.RData"))

df_scores <- df_scores %>% subset(!grepl("re", Patient.ID))

X_rAUC <- dcast(data = df_scores, Patient.num ~ drug , value.var="rAUC")
rownames(X_rAUC) <- as.character(X_rAUC$Patient.num); X_rAUC$Patient.num <- NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]


df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$rAUC <- 1-df_scores$rAUC
#df_scores$DSS_AUC_log2 <- NULL

table(Y$status)/nrow(Y)*10

# Comparing AUC/DSS/EC50 survival prediction ----
list_res_DS_metrics <- list()
for(m in colnames(df_scores)[9]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  if(grepl("EC50",m)){
    X <- log10(X)
    X[,] <- apply(X,2, function(x) x - mean(x))
  }
  list_res_DS_metrics[[m]] <-  Cox_forecasting_glmnet(X_data=X,
                                                      y_data=data.matrix(Y),
                                                      alpha=c(0),
                                                      lambda=c(exp(seq(-8,6, 0.1))),
                                                      free_cores = 5,
                                                      test.n= c(5,5),
                                                      nfolds = nrow(Y),
                                                      iter=200,
                                                      log_AUC=2,
                                                      Patient.Z=1:2,
                                                      Drug.Z =1:2,
                                                      RCPC=0)
}

save(list_res_DS_metrics, file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Ridge.RData"))



df_C.index.alldata <- list_res_DS_metrics
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}

df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient standardization", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug standardization", "")
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))

ggexport(
  ggarrange(ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_test))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Testing: n=",10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank()),
            
            ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_train))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkorange", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank())),
  width=7, height=4,
  filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Ridge.pdf")
)



list_res_DS_metrics_L1 <- list()
for(m in colnames(df_scores)[6:14]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  if(grepl("EC50",m)){
    X <- log10(X)
    X[,] <- apply(X,2, function(x) x - mean(x))
  }
  list_res_DS_metrics_L1[[m]] <-  Cox_forecasting_glmnet(X_data=X,
                                                         y_data=data.matrix(Y),
                                                         alpha=c(1,0.4),
                                                         lambda=c(exp(seq(-8,6, 0.1))),
                                                         free_cores = 5,
                                                         test.n= c(5,5),
                                                         nfolds = nrow(Y),
                                                         iter=200,
                                                         log_AUC=2,
                                                         Patient.Z=1:2,
                                                         Drug.Z =1:2,
                                                         RCPC=0)
}

save(list_res_DS_metrics_L1, file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_elnet and lasso.RData"))



df_C.index.alldata <- list_res_DS_metrics_L1
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- data.frame(df_C.index.alldata[[m]]$C_index_results)
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Penalty <- as.factor(as.numeric(gsub(".*Penalty_","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}

df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient standardization", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug standardization", "")
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))


ggexport(
  ggarrange(ggplot(subset(df_C.index.alldata, Penalty == 1), aes(x=Metric, y=C_index_test))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkblue", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Testing: n=",10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank()),
            
            ggplot(subset(df_C.index.alldata, Penalty == 1), aes(x=Metric, y=C_index_train))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkblue", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank())),
  width=7, height=4,
  filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Lasso.pdf")
)

ggexport(
  ggarrange(ggplot(subset(df_C.index.alldata, Penalty != 1), aes(x=Metric, y=C_index_test))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkgreen", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Testing: n=",10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank()),
            
            ggplot(subset(df_C.index.alldata, Penalty != 1), aes(x=Metric, y=C_index_train))+
              geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
              geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
              geom_boxplot(fill = "darkgreen", alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
              facet_grid(Drug_stdz~Patient_stdz)+
              ylim(0,1)+
              theme_bw()+
              labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.background = element_blank())),
  width=7, height=4,
  filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Elastic net.pdf")
)


# Pre-selection and alpha testing ----


list_res_rAUC_cva <- Cox_forecasting_glmnet_CVA(X_data=X_rAUC,
                                                y_data=data.matrix(Y),
                                                alpha="CVA",
                                                lambda=c(exp(seq(-4,6, 0.1))),
                                                free_cores = 5,
                                                test.n= c(6,4),
                                                nfolds = nrow(Y),
                                                iter=50,
                                                log_AUC=1,
                                                Patient.Z=1,
                                                Drug.Z =2,
                                                RCPC=0)


list_res_model_testing_quant <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X))
  #X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(t in c(quantile(SDs, 1-c(25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,125,150,200,250, 300)/ncol(X)),0)){
    list_res_model_testing_quant[[paste(m,"sd cutoff",t,ncol(X[,SDs>t]))]] <- Cox_forecasting_glmnet(X_data=X[,SDs>t],
                                                                                                     y_data=data.matrix(Y),
                                                                                                     alpha=list_res_rAUC_cva$CVA_results$alpha,
                                                                                                     lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                     free_cores = 2,
                                                                                                     test.n= c(5,5),
                                                                                                     nfolds = nrow(Y),
                                                                                                     iter=50,
                                                                                                     log_AUC=2,
                                                                                                     Patient.Z=1,
                                                                                                     Drug.Z =2,
                                                                                                     RCPC=0)
    cat("\n ", paste(m,"sd cutoff",t,ncol(X[,SDs>t])), "completed \n ")
  }
}

save(list_res_rAUC_cva, list_res_model_testing_quant,
     file=paste0(path,"/Results_survival predictions/Alpha comp_with variable pre-selection from SD quantiles.RData"))

list_res_model_testing_quant_rAUC <- list()
for(m in colnames(df_scores)[c(6,7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  #X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(t in c(quantile(SDs, 1-c(25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,125,150,200,250,300)/ncol(X)),0)){
    list_res_model_testing_quant_rAUC[[paste(m,"sd cutoff",t,ncol(X[,SDs>t]))]] <- Cox_forecasting_glmnet(X_data=X[,SDs>t],
                                                                                                          y_data=data.matrix(Y),
                                                                                                          alpha=list_res_rAUC_cva$CVA_results$alpha,
                                                                                                          lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                          free_cores = 5,
                                                                                                          test.n= c(5,5),
                                                                                                          nfolds = nrow(Y),
                                                                                                          iter=50,
                                                                                                          log_AUC=2,
                                                                                                          Patient.Z=1,
                                                                                                          Drug.Z =2,
                                                                                                          RCPC=0)
    cat("\n ", paste(m,"sd cutoff",t,ncol(X[,SDs>t])), "completed \n ")
  }
}

save(list_res_rAUC_cva, list_res_model_testing_quant_rAUC,
     file=paste0(path,"/Results_survival predictions/Alpha comp_with variable pre-selection from SD quantiles_rAUC.RData"))


df_C.index.alldata <- list_res_model_testing_quant
for(i in 1:length(df_C.index.alldata)){
  df_C.index.alldata[[i]] <- df_C.index.alldata[[i]]$C_index_results %>% mutate(Selection = "rAUC_log2/DSS3")
}


df_C.index.alldata1 <- list_res_model_testing_quant_rAUC
for(i in 1:length(df_C.index.alldata1)){
  df_C.index.alldata1[[i]] <- df_C.index.alldata1[[i]]$C_index_results %>% mutate(Selection = "rAUC")
}


df_C.index.alldata <- rbind(bind_rows(df_C.index.alldata, .id="Data"),
                            bind_rows(df_C.index.alldata1, .id="Data"))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient stdz", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug stdz", "")
df_C.index.alldata$Penalty <- gsub(".*_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub(" .*","", df_C.index.alldata$Data)
df_C.index.alldata$Pre_selection <- gsub(".* sd cutoff ", "", df_C.index.alldata$Data)

df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Pre_selection, Penalty, Metric, Selection) %>% summarise(Median=median(C_index_test),
                                                                                                                Mean=mean(C_index_test),
                                                                                                                MAD=mad(C_index_test),
                                                                                                                SD=sd(C_index_test),
                                                                                                                q975=quantile(C_index_test, 0.975),
                                                                                                                q95=quantile(C_index_test, 0.95),
                                                                                                                q05=quantile(C_index_test, 0.05),
                                                                                                                q025=quantile(C_index_test, 0.025))
df_C.index.stat$Pre_selection_SD <- round(as.numeric(gsub(" .*","",df_C.index.stat$Pre_selection)),3)
df_C.index.stat$Pre_selection_N <- round(as.numeric(gsub(".* ","",df_C.index.stat$Pre_selection)))
df_C.index.stat$Metric <- factor(df_C.index.stat$Metric, levels=c("rAUC", "rAUC_log2", "DSS3"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_1.pdf"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Median)) +
           geom_tile(color = "black") +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_2.pdf"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = (Mean-0.5)/SD)) +
           geom_tile(color = "black") +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="(Mean-0.5)/SD C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_3.pdf"))


ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           geom_text(aes(label=ifelse((Mean-0.5)/SD>2, "*", "")),color = "black", size=2.5) +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_4.pdf"))


ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Median)) +
           geom_tile(color = "black") +
           geom_text(aes(label=ifelse(q025>0.5, "*", "")),color = "black", size=2.5) +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_5.pdf"))


ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Median)) +
           geom_tile(color = "black") +
           geom_text(aes(label=ifelse(q05>0.5, "*", "")),color = "black", size=2.5) +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_6.pdf"))



ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           geom_text(aes(label=ifelse(q025>0.5, "*", "")),color = "black", size=2.5) +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_5b.pdf"))


ggexport(ggplot(df_C.index.stat, aes(x = factor(round(as.numeric(Penalty),4)), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           geom_text(aes(label=ifelse(q05>0.5, "*", "")),color = "black", size=2.5) +
           #scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn")))+
           #scale_fill_viridis_c(option="A") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x=expression(paste("Penalty mixture factor (", alpha,")")), y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Selection~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=7, filename = paste0(path,"/Results_survival predictions/Pre-selection and penalty comparison_quantile selection_6b.pdf"))


# Comparing Hill AUC survival prediction ----

load(file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))
df.hillAUC <- df.hillAUC %>% subset(!grepl("re", Patient.ID))

df.hillAUC$rAUC <- 1-df.hillAUC$rAUC
df.hillAUC$hill_AUC <- 1-df.hillAUC$hill_AUC
df.hillAUC$hill0_AUC <- 1-df.hillAUC$hill0_AUC
df.hillAUC$hillB_AUC <- 1-df.hillAUC$hillB_AUC

list_res_hill_metrics <- list()
for(m in colnames(df.hillAUC)[c(6:8,10:12)]){
  X <- dcast(data = df.hillAUC, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  
  list_res_hill_metrics[[m]] <-  Cox_forecasting_glmnet(X_data=X,
                                                         y_data=data.matrix(Y),
                                                         alpha=c(1,0.4,0),
                                                         lambda=c(exp(seq(-8,6, 0.1))),
                                                         free_cores = 5,
                                                         test.n= c(5,5),
                                                         nfolds = nrow(Y),
                                                         iter=200,
                                                         log_AUC=2,
                                                         Patient.Z=1:2,
                                                         Drug.Z =2,
                                                         RCPC=0)
}
save(list_res_hill_metrics, file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Hill AUC.RData"))




load(file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Ridge.RData"))
load(file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_elnet and lasso.RData"))

df_C.index.alldata_ref <- c(list_res_DS_metrics[1:2], list_res_DS_metrics_L1[1:2])
for(m in 1:4){
  df_C.index.alldata_ref[[m]] <- df_C.index.alldata_ref[[m]]$C_index_results
  df_C.index.alldata_ref[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata_ref[[m]]$ID)))
  df_C.index.alldata_ref[[m]]$Penalty <- as.factor(as.numeric(gsub(".*Penalty_","",df_C.index.alldata_ref[[m]]$ID)))
  df_C.index.alldata_ref[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata_ref[[m]]$ID))
}

df_C.index.alldata_ref <- bind_rows(df_C.index.alldata_ref, .id="Metric")
df_C.index.alldata_ref$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata_ref$ID), "Patient standardization", "")
df_C.index.alldata_ref$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata_ref$ID), "Drug standardization", "")
df_C.index.alldata_ref$Metric <- factor(df_C.index.alldata_ref$Metric, levels=unique(df_C.index.alldata_ref$Metric))
df_C.index.stat_ref <- df_C.index.alldata_ref %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))


df_C.index.alldata <- list_res_hill_metrics
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- data.frame(df_C.index.alldata[[m]]$C_index_results)
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Penalty <- as.factor(as.numeric(gsub(".*Penalty_","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}

df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient standardization", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug standardization", "")
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))

df_C.index.alldata$log2 <- ifelse(grepl("log2", df_C.index.alldata$Metric), "log2", "")
df_C.index.alldata_ref$log2 <- ifelse(grepl("log2", df_C.index.alldata_ref$Metric), "log2", "")
df_C.index.alldata$Penalty_type <-ifelse(df_C.index.alldata$Penalty==0, "Ridge", "Lasso")
df_C.index.alldata$Penalty_type[df_C.index.alldata$Penalty==0.4] <- "Elastic net"
df_C.index.alldata_ref$Penalty_type <-ifelse(df_C.index.alldata_ref$Penalty==0, "Ridge", "Lasso")
df_C.index.alldata_ref$Penalty_type[df_C.index.alldata_ref$Penalty==0.4] <- "Elastic net"

df_C.index.alldata$C_index_test_ref <- df_C.index.alldata_ref$C_index_test[match(paste(df_C.index.alldata$ID,df_C.index.alldata$Iteration, df_C.index.alldata$log2),
                                                                             paste(df_C.index.alldata_ref$ID,df_C.index.alldata_ref$Iteration, df_C.index.alldata_ref$log2))]
df_C.index.alldata$C_index_train_ref <- df_C.index.alldata_ref$C_index_train[match(paste(df_C.index.alldata$ID,df_C.index.alldata$Iteration, df_C.index.alldata$log2),
                                                                               paste(df_C.index.alldata_ref$ID,df_C.index.alldata_ref$Iteration, df_C.index.alldata_ref$log2))]


df_C.index.alldata$Penalty_type <- factor(df_C.index.alldata$Penalty_type, levels=c("Ridge", "Elastic net", "Lasso"))
ggarrange(ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_test))+
            geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
            geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
            geom_boxplot(aes(fill=Penalty_type), alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
            facet_grid(~Patient_stdz)+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Testing: n=",10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),
          
          ggplot(subset(df_C.index.alldata, RCPC == 0), aes(x=Metric, y=C_index_train))+
            geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
            geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
            geom_boxplot(aes(fill=Penalty_type), alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
            facet_grid(~Patient_stdz)+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),common.legend = T)
ggsave(  width=8, height=3,
         filename = paste0(path,"/Results_survival predictions/Drug sensitivity score_Hill fitting comparison.pdf")
)
ggarrange(ggplot(df_C.index.alldata, aes(x=Metric, y=C_index_test - C_index_test_ref))+
            geom_hline(yintercept = 0, lty=2, alpha=0.5)+
            geom_boxplot(aes(fill=Penalty_type), alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
            facet_grid(~Patient_stdz)+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            #ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Testing: n=",10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),
          
          ggplot(df_C.index.alldata, aes(x=Metric, y=C_index_train - C_index_train_ref))+
            geom_hline(yintercept = 0, lty=2, alpha=0.5)+
            geom_boxplot(aes(fill=Penalty_type), alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
            facet_grid(~Patient_stdz)+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            #ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),common.legend = T)
ggsave(width=8, height=3,
       filename = paste0(path,"/Results_survival predictions/Drug sensitivity score_Hill fitting effect_1.pdf"))

ggarrange(ggplot(df_C.index.alldata %>% group_by(Patient_stdz, Metric, Penalty_type) %>% 
                   summarise(Mean=mean(C_index_test - C_index_test_ref),
                             SD=sd(C_index_test - C_index_test_ref)), 
                 aes(x=Metric, y=Mean))+
            geom_hline(yintercept = 0, lty=2, alpha=0.5)+
            geom_linerange(aes(col=Penalty_type, ymin=Mean-SD, ymax=Mean+SD),  position = position_dodge2(width = 0.75), col="black")+
            geom_point(aes(fill=Penalty_type), pch=21, position = position_dodge2(width = 0.75), size=2)+
            facet_grid(~Patient_stdz)+
            scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            #ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Testing: n=",10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),
          
          ggplot(df_C.index.alldata %>% group_by(Patient_stdz, Metric, Penalty_type) %>% 
                   summarise(Mean=mean(C_index_train - C_index_train_ref),
                             SD=sd(C_index_train - C_index_train_ref)), 
                 aes(x=Metric, y=Mean))+
            geom_hline(yintercept = 0, lty=2, alpha=0.5)+
            geom_linerange(aes(col=Penalty_type, ymin=Mean-SD, ymax=Mean+SD),  position = position_dodge2(width = 0.75), col="black")+
            geom_point(aes(fill=Penalty_type), pch=21, position = position_dodge2(width = 0.75), size=2)+
            facet_grid(~Patient_stdz)+
            scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
            #ylim(0,1)+
            theme_bw()+
            labs(title=paste0("Training: n=",69-10), x="", y="C-index")+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank()),common.legend = T)
ggsave(width=8, height=3,
         filename = paste0(path,"/Results_survival predictions/Drug sensitivity score_Hill fitting effect_2.pdf"))









# De-confounding ----

list_res_RCPC_testing <- list()
for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  #X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(t in c(quantile(SDs, 1-c(40,60,80,100)/ncol(X_rAUC)),0)){
    list_res_RCPC_testing[[paste(m,"sd cutoff",t,ncol(X[,SDs>t]))]] <- Cox_forecasting_glmnet(X_data=X[,SDs>=t],
                                                                                              y_data=data.matrix(Y),
                                                                                              alpha=c(0,0.4,1),
                                                                                              lambda=c(exp(seq(-8,6, 0.1))),
                                                                                              free_cores = 5,
                                                                                              test.n= c(5,5),
                                                                                              nfolds = nrow(Y),
                                                                                              iter=50,
                                                                                              log_AUC=2,
                                                                                              Patient.Z=1:2,
                                                                                              Drug.Z =2,
                                                                                              RCPC=0:8)
    cat("\n ", paste(m,"sd cutoff",t,ncol(X[,SDs>t])), "completed \n ")
  }
}

list_res_RCPC_testing_rd <- list()
for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  
  for(t in 1:50){
    set.seed(t)
    ids <- sample(1:ncol(X), 100, replace = FALSE, prob = SDs)
    list_res_RCPC_testing_rd[[paste(m,t)]] <- Cox_forecasting_glmnet(X_data=X[,ids],
                                                                     y_data=data.matrix(Y),
                                                                     alpha=1,
                                                                     lambda=c(exp(seq(-8,6, 0.1))),
                                                                     free_cores = 5,
                                                                     test.n= c(5,5),
                                                                     nfolds = nrow(Y),
                                                                     iter=50,
                                                                     log_AUC=2,
                                                                     Patient.Z=1:2,
                                                                     Drug.Z =2,
                                                                     RCPC=0:8)
    cat("\n ", paste(m,t), "completed \n ")
  }
}

save(list_res_RCPC_testing,list_res_RCPC_testing_rd,
     file=paste0(path,"/Results_survival predictions/RCPC model testing with variable pre-selection_with weighted sampling.RData"))

load(file=paste0(path,"/Results_survival predictions/RCPC model testing with variable pre-selection_with weighted sampling.RData"))



df_C.index.alldata <- list_res_RCPC_testing
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Data")
df_C.index.alldata$RCPC <- as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata$ID))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "z-score", "")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub(" .*","",df_C.index.alldata$Data)
df_C.index.alldata$Pre_selection_SD <- as.numeric(gsub(" .*","",gsub(".*cutoff ","",df_C.index.alldata$Data)))
df_C.index.alldata$Pre_selection_N <- as.numeric(gsub(".* ","",df_C.index.alldata$Data))

df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))
df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))


df_C.index.stat <- df_C.index.alldata %>% 
  group_by(RCPC, Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  summarise(Median=median(C_index_test),
            Mean=mean(C_index_test),
            MAD=mad(C_index_test),
            SD=sd(C_index_test))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) 

lm_test <- lm(RCPC~Metric*Patient_stdz+Pre_selection_N+Model,data=df_C.index_optimPC)
summary(lm_test)

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]




ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_optimal RCPC_1.pdf"))

kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""


ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Pre_selection_N)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_optimal RCPC_2.pdf"))


pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "Mean") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval1 <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


df_C.index_optimPC$pval1_pair <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval1_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval1_pair <- ifelse(df_C.index_optimPC$pval1_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval1_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval1_pair[is.na(df_C.index_optimPC$pval1_pair)] <- ""

ggexport(ggarrange(ggplot(df_C.index_optimPC, 
                          aes(x=Patient_stdz, y=Mean))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     facet_grid(~Metric+paste(scales::pvalue(pval1,
                                                             accuracy = 0.01,
                                                             decimal.mark = ".", 
                                                             add_p = TRUE)))+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(df_C.index_optimPC, 
                          aes(x=Metric, y=Mean))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=df_C.index_optimPC %>% 
                                 group_by(Patient_stdz,Metric) %>% top_n(1,Mean),
                               aes(label=pval1_pair, y=Mean+0.025),size=3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_optimal RCPC_3.pdf"))






df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) 


lm_test <- lm(RCPC~Metric*Patient_stdz+Pre_selection_N+Model,data=df_C.index_optimPC)
summary(lm_test)

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]




ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_median_optimal RCPC_1.pdf"))


kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Pre_selection_N)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_median_optimal RCPC_2.pdf"))

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "Mean") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval1 <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


df_C.index_optimPC$pval1_pair <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz==""], 
                                                      df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                      p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval1_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                           df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                           p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval1_pair <- ifelse(df_C.index_optimPC$pval1_pair<0.05,
                                        scales::pvalue(df_C.index_optimPC$pval1_pair ,
                                                       accuracy = 0.01,
                                                       decimal.mark = ".", 
                                                       add_p = TRUE),"ns")
df_C.index_optimPC$pval1_pair[is.na(df_C.index_optimPC$pval1_pair)] <- ""
ggexport(ggarrange(ggplot(df_C.index_optimPC, 
                          aes(x=Patient_stdz, y=Median))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     facet_grid(~Metric+paste(scales::pvalue(pval1,
                                                             accuracy = 0.01,
                                                             decimal.mark = ".", 
                                                             add_p = TRUE)))+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(df_C.index_optimPC, 
                          aes(x=Metric, y=Median))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=df_C.index_optimPC %>% 
                                 group_by(Patient_stdz,Metric) %>% top_n(1,Mean),
                               aes(label=pval1_pair, y=Mean+0.05),size=3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6, filename = paste0(path,"/Results_RCPC/RCPC testing_median_optimal RCPC_3.pdf"))




ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Metric+Patient_stdz~Model)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_1.pdf"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Median)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Metric+Patient_stdz~Model)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_median_1.pdf"))



df_C.index.stat <- df_C.index.stat[order(df_C.index.stat$Model, df_C.index.stat$Penalty, df_C.index.stat$Metric, df_C.index.stat$Pre_selection_N,df_C.index.stat$Patient_stdz),]
df_C.index.stat$Median_minus <- c(NA, df_C.index.stat$Median[-nrow(df_C.index.stat)])
df_C.index.stat$Median_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Median_delta <- df_C.index.stat$Median - df_C.index.stat$Median_minus
df_C.index.stat$Mean_minus <- c(NA,df_C.index.stat$Mean[-nrow(df_C.index.stat)])
df_C.index.stat$Mean_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Mean_delta <- df_C.index.stat$Mean - df_C.index.stat$Mean_minus

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Mean_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]),mid="white",high=scales::muted(pal_jco()(2)[2]),midpoint=0) +
  labs(x="RCPC", y="# features\nSD pre-selection", fill="Mean C-index change")+
  facet_grid(Metric+Patient_stdz~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_mean change_1.pdf"))
ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Median_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]),mid="white",high=scales::muted(pal_jco()(2)[2]),midpoint=0) +
  labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index change")+
  facet_grid(Metric+Patient_stdz~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_median change_1.pdf"))




df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Mean, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Mean=max(Mean))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_optimal RCPC_4.pdf"))


ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_mean_optimal RCPC_4b.pdf"))



df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Median, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Median=max(Median))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_median_optimal RCPC_4.pdf"))




ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_median_optimal RCPC_4b.pdf"))





df_C.index.alldata <- list_res_RCPC_testing_rd
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Data")
df_C.index.alldata$RCPC <- as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata$ID))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "z-score", "")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub(" .*","",df_C.index.alldata$Data)

df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))

df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))
df_C.index.alldata$Data <- gsub(".* ","",df_C.index.alldata$Data)

df_C.index.stat <- df_C.index.alldata %>% 
  group_by(RCPC, Model, Penalty, Metric, Data, Patient_stdz) %>% 
  summarise(Median=median(C_index_test),
            Mean=mean(C_index_test),
            MAD=mad(C_index_test),
            SD=sd(C_index_test))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC)

wilcox_data <- dcast(df_C.index_optimPC, Data~Metric+Patient_stdz, value.var = "RCPC")
df_wilcox_res <- wilcox.test(wilcox_data[,2],wilcox_data[,3], paired = TRUE, alternative = "two.sided")
df_wilcox_res$p.value
wilcox_data <- wilcox_data[,-1]
wilcox_optim_RCPC_mean <- outer(wilcox_data[,], 
                                wilcox_data[,], 
                                FUN = Vectorize(function(x,y) wilcox.test(x,y, paired = TRUE, alternative = "two.sided")$p.value))

write.csv(wilcox_optim_RCPC_mean, 
          file = paste0(path,"/Results_RCPC/RCPC testing_random 100_mean_optimal RCPC_wilcoxon test.csv"))

pvals = dcast(df_C.index_optimPC, Metric+Model+Data ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_mean_optimal RCPC_1.pdf"))

kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric , y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Data)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_mean_optimal RCPC_2.pdf"))

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Mean))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
           facet_grid(Patient_stdz ~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=10, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_mean_1.pdf"))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Mean, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Mean=max(Mean))


ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_mean_optimal RCPC_3.pdf"))


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC)

wilcox_data <- dcast(df_C.index_optimPC, Data~Metric+Patient_stdz, value.var = "RCPC")
df_wilcox_res <- wilcox.test(wilcox_data[,2],wilcox_data[,3], paired = TRUE, alternative = "two.sided")
df_wilcox_res$p.value
wilcox_data <- wilcox_data[,-1]
wilcox_optim_RCPC_median <- outer(wilcox_data[,], 
                                  wilcox_data[,], 
                                  FUN = Vectorize(function(x,y) wilcox.test(x,y, paired = TRUE, alternative = "two.sided")$p.value))


write.csv(wilcox_optim_RCPC_median, 
          file = paste0(path,"/Results_RCPC/RCPC testing_random 100_median_optimal RCPC_wilcoxon test.csv"))


pvals = dcast(df_C.index_optimPC, Metric+Model+Data ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]



ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_median_optimal RCPC_1.pdf"))


kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Data)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+facet_grid(~Patient_stdz)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_median_optimal RCPC_2.pdf"))

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Median))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Median)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Patient_stdz ~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=10, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_median_1.pdf"))


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Median, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Median=max(Median))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_random 100_median_optimal RCPC_3.pdf"))



df_C.index.stat <- df_C.index.stat[order(df_C.index.stat$Model, df_C.index.stat$Penalty, df_C.index.stat$Metric, df_C.index.stat$Data,df_C.index.stat$Patient_stdz),]
df_C.index.stat$Median_minus <- c(NA, df_C.index.stat$Median[-nrow(df_C.index.stat)])
df_C.index.stat$Median_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Median_delta <- df_C.index.stat$Median - df_C.index.stat$Median_minus
df_C.index.stat$Mean_minus <- c(NA,df_C.index.stat$Mean[-nrow(df_C.index.stat)])
df_C.index.stat$Mean_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Mean_delta <- df_C.index.stat$Mean - df_C.index.stat$Mean_minus


m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Mean))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Mean_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
  labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
  facet_grid(Patient_stdz ~Metric)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Median))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Median_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
  labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
  facet_grid(Patient_stdz ~Metric)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  summarise(RCPC_min.dMedian = RCPC[which.min(Median_delta)],
            RCPC_max.dMedian = RCPC[which.max(Median_delta)],
            RCPC_min.dMean = RCPC[which.min(Mean_delta)],
            RCPC_max.dMean = RCPC[which.max(Mean_delta)])
df_C.index_optimPC <- melt(df_C.index_optimPC, measure.vars=colnames(df_C.index_optimPC)[6:9])

ggplot(df_C.index_optimPC, 
       aes(x=Patient_stdz, y=value))+
  geom_violin(fill="black", alpha=0.3, size=0.5)+
  geom_line(aes(group=paste(Model,Data)),fill="black", alpha=0.3, size=0.5)+
  geom_point(fill="black", pch=21, alpha=0.5)+
  facet_grid(variable~Metric)+
  theme_bw()+
  labs(x="", y="Optimal RCPC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")


# De-confounding - short term survival----
Y1 <- Y
Y1$status[Y1$time>365*2] <- 0
Y1$time[Y1$time>365*2] <- 365*2

res.cox <- coxph(Surv(time, status) ~ ., data =  Y1)
fit <- survfit(res.cox)
survplot <- ggsurvplot(fit, Y1,conf.int = TRUE, palette = "jco", 
                       censor = TRUE, surv.median.line = "hv")
survplot$plot+
  scale_y_continuous(labels = scales::percent, limits=c(0,1))+
  labs(y="Survival", x="Time (days)")+
  theme(legend.position = "none")


table(Y1$status)/nrow(Y1)

list_res_RCPC_testing_sts <- list()
for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  #X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(t in c(quantile(SDs, 1-c(40,60,80,100)/ncol(X_rAUC)),0)){
    list_res_RCPC_testing_sts[[paste(m,"sd cutoff",t,ncol(X[,SDs>t]))]] <- Cox_forecasting_glmnet(X_data=X[,SDs>=t],
                                                                                              y_data=data.matrix(Y1),
                                                                                              alpha=c(0,0.4,1),
                                                                                              lambda=c(exp(seq(-8,6, 0.1))),
                                                                                              free_cores = 5,
                                                                                              test.n= c(6,4),
                                                                                              nfolds = nrow(Y),
                                                                                              iter=50,
                                                                                              log_AUC=2,
                                                                                              Patient.Z=1:2,
                                                                                              Drug.Z =2,
                                                                                              RCPC=0:8)
    cat("\n ", paste(m,"sd cutoff",t,ncol(X[,SDs>t])), "completed \n ")
  }
}

list_res_RCPC_testing_rd_sts <- list()
for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  
  for(t in 1:50){
    set.seed(t)
    ids <- sample(1:ncol(X), 100, replace = FALSE, prob = SDs)
    list_res_RCPC_testing_rd_sts[[paste(m,t)]] <- Cox_forecasting_glmnet(X_data=X[,ids],
                                                                     y_data=data.matrix(Y1),
                                                                     alpha=1,
                                                                     lambda=c(exp(seq(-8,6, 0.1))),
                                                                     free_cores = 5,
                                                                     test.n= c(6,4),
                                                                     nfolds = nrow(Y),
                                                                     iter=50,
                                                                     log_AUC=2,
                                                                     Patient.Z=1:2,
                                                                     Drug.Z =2,
                                                                     RCPC=0:8)
    cat("\n ", paste(m,t), "completed \n ")
  }
}

save(list_res_RCPC_testing_sts,list_res_RCPC_testing_rd_sts,
     file=paste0(path,"/Results_survival predictions/RCPC model testing with variable pre-selection_with weighted sampling_STS.RData"))

load(file=paste0(path,"/Results_survival predictions/RCPC model testing with variable pre-selection_with weighted sampling_STS.RData"))



df_C.index.alldata <- list_res_RCPC_testing_sts
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Data")
df_C.index.alldata$RCPC <- as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata$ID))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "z-score", "")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub(" .*","",df_C.index.alldata$Data)
df_C.index.alldata$Pre_selection_SD <- as.numeric(gsub(" .*","",gsub(".*cutoff ","",df_C.index.alldata$Data)))
df_C.index.alldata$Pre_selection_N <- as.numeric(gsub(".* ","",df_C.index.alldata$Data))

df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% group_by(ID, Metric) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))
df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))


df_C.index.stat <- df_C.index.alldata %>% 
  group_by(RCPC, Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  summarise(Median=median(C_index_test),
            Mean=mean(C_index_test),
            MAD=mad(C_index_test),
            SD=sd(C_index_test))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) 

lm_test <- lm(RCPC~Metric*Patient_stdz+Pre_selection_N+Model,data=df_C.index_optimPC)
summary(lm_test)

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]




ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_optimal RCPC_1.pdf"))

kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""


ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Pre_selection_N)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_optimal RCPC_2.pdf"))


pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "Mean") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval1 <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


df_C.index_optimPC$pval1_pair <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz==""], 
                                                      df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                      p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval1_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                           df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                           p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval1_pair <- ifelse(df_C.index_optimPC$pval1_pair<0.05,
                                        scales::pvalue(df_C.index_optimPC$pval1_pair ,
                                                       accuracy = 0.01,
                                                       decimal.mark = ".", 
                                                       add_p = TRUE),"ns")
df_C.index_optimPC$pval1_pair[is.na(df_C.index_optimPC$pval1_pair)] <- ""

ggexport(ggarrange(ggplot(df_C.index_optimPC, 
                          aes(x=Patient_stdz, y=Mean))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     facet_grid(~Metric+paste(scales::pvalue(pval1,
                                                             accuracy = 0.01,
                                                             decimal.mark = ".", 
                                                             add_p = TRUE)))+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(df_C.index_optimPC, 
                          aes(x=Metric, y=Mean))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=df_C.index_optimPC %>% 
                                 group_by(Patient_stdz,Metric) %>% top_n(1,Mean),
                               aes(label=pval1_pair, y=Mean+0.025),size=3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_optimal RCPC_3.pdf"))






df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) 


lm_test <- lm(RCPC~Metric*Patient_stdz+Pre_selection_N+Model,data=df_C.index_optimPC)
summary(lm_test)

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]




ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_optimal RCPC_1.pdf"))


kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="black", alpha=0.2)+
           geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
           geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Pre_selection_N)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 8.5, height=3.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_optimal RCPC_2.pdf"))

pvals = dcast(df_C.index_optimPC, Metric+ Pre_selection_N+Model ~Patient_stdz, value.var = "Mean") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval1 <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


df_C.index_optimPC$pval1_pair <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz==""], 
                                                      df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                      p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval1_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$Mean[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                           df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                           p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval1_pair <- ifelse(df_C.index_optimPC$pval1_pair<0.05,
                                        scales::pvalue(df_C.index_optimPC$pval1_pair ,
                                                       accuracy = 0.01,
                                                       decimal.mark = ".", 
                                                       add_p = TRUE),"ns")
df_C.index_optimPC$pval1_pair[is.na(df_C.index_optimPC$pval1_pair)] <- ""
ggexport(ggarrange(ggplot(df_C.index_optimPC, 
                          aes(x=Patient_stdz, y=Median))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     facet_grid(~Metric+paste(scales::pvalue(pval1,
                                                             accuracy = 0.01,
                                                             decimal.mark = ".", 
                                                             add_p = TRUE)))+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(df_C.index_optimPC, 
                          aes(x=Metric, y=Median))+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=df_C.index_optimPC %>% 
                                 group_by(Patient_stdz,Metric) %>% top_n(1,Mean),
                               aes(label=pval1_pair, y=Mean+0.05),size=3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_optimal RCPC_3.pdf"))




ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Mean C-index")+
           facet_grid(Metric+Patient_stdz~Model)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_1.pdf"))

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Median)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Metric+Patient_stdz~Model)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_1.pdf"))



df_C.index.stat <- df_C.index.stat[order(df_C.index.stat$Model, df_C.index.stat$Penalty, df_C.index.stat$Metric, df_C.index.stat$Pre_selection_N,df_C.index.stat$Patient_stdz),]
df_C.index.stat$Median_minus <- c(NA, df_C.index.stat$Median[-nrow(df_C.index.stat)])
df_C.index.stat$Median_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Median_delta <- df_C.index.stat$Median - df_C.index.stat$Median_minus
df_C.index.stat$Mean_minus <- c(NA,df_C.index.stat$Mean[-nrow(df_C.index.stat)])
df_C.index.stat$Mean_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Mean_delta <- df_C.index.stat$Mean - df_C.index.stat$Mean_minus

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Mean_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]),mid="white",high=scales::muted(pal_jco()(2)[2]),midpoint=0) +
  labs(x="RCPC", y="# features\nSD pre-selection", fill="Mean C-index change")+
  facet_grid(Metric+Patient_stdz~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_mean change_1.pdf"))
ggplot(df_C.index.stat, aes(x = factor(RCPC), y = factor(Pre_selection_N), fill = Median_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]),mid="white",high=scales::muted(pal_jco()(2)[2]),midpoint=0) +
  labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index change")+
  facet_grid(Metric+Patient_stdz~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 7, height=8, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median change_1.pdf"))




df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Mean, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Mean=max(Mean))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_optimal RCPC_4.pdf"))


ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_mean_optimal RCPC_4b.pdf"))



df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Pre_selection_N, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Median, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Median=max(Median))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_optimal RCPC_4.pdf"))




ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N),fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     #geom_hline(yintercept=0.5,lty=2)+
                     geom_violin(fill="black", col="black", alpha=0.2)+
                     geom_line(aes(group=paste(Model,Pre_selection_N), col=Model), alpha=0.3, size=0.5)+
                     geom_point(aes(size=factor(Pre_selection_N), fill=Model), pch=21, alpha=0.5)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     scale_size_manual(values = c(0.2,0.4,0.6,0.8,1)*3)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
                     theme_bw()+
                     labs(x="", y="Best C-index median")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_median_optimal RCPC_4b.pdf"))





df_C.index.alldata <- list_res_RCPC_testing_rd_sts
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Data")
df_C.index.alldata$RCPC <- as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata$ID))
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "z-score", "")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Metric <- gsub(" .*","",df_C.index.alldata$Data)

df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))

df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))
df_C.index.alldata$Data <- gsub(".* ","",df_C.index.alldata$Data)

df_C.index.stat <- df_C.index.alldata %>% 
  group_by(RCPC, Model, Penalty, Metric, Data, Patient_stdz) %>% 
  summarise(Median=median(C_index_test),
            Mean=mean(C_index_test),
            MAD=mad(C_index_test),
            SD=sd(C_index_test))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC)

wilcox_data <- dcast(df_C.index_optimPC, Data~Metric+Patient_stdz, value.var = "RCPC")
df_wilcox_res <- wilcox.test(wilcox_data[,2],wilcox_data[,3], paired = TRUE, alternative = "two.sided")
df_wilcox_res$p.value
wilcox_data <- wilcox_data[,-1]
wilcox_optim_RCPC_mean <- outer(wilcox_data[,], 
                                wilcox_data[,], 
                                FUN = Vectorize(function(x,y) wilcox.test(x,y, paired = TRUE, alternative = "two.sided")$p.value))

write.csv(wilcox_optim_RCPC_mean, 
          file = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_mean_optimal RCPC_wilcoxon test.csv"))

pvals = dcast(df_C.index_optimPC, Metric+Model+Data ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]


ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_mean_optimal RCPC_1.pdf"))

kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric , y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Data)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_mean_optimal RCPC_2.pdf"))

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Mean))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
           facet_grid(Patient_stdz ~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=10, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_mean_1.pdf"))

df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Mean) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Mean, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Mean=max(Mean))


ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Mean+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Mean))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Mean+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_mean_optimal RCPC_3.pdf"))


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC)

wilcox_data <- dcast(df_C.index_optimPC, Data~Metric+Patient_stdz, value.var = "RCPC")
df_wilcox_res <- wilcox.test(wilcox_data[,2],wilcox_data[,3], paired = TRUE, alternative = "two.sided")
df_wilcox_res$p.value
wilcox_data <- wilcox_data[,-1]
wilcox_optim_RCPC_median <- outer(wilcox_data[,], 
                                  wilcox_data[,], 
                                  FUN = Vectorize(function(x,y) wilcox.test(x,y, paired = TRUE, alternative = "two.sided")$p.value))


write.csv(wilcox_optim_RCPC_median, 
          file = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_median_optimal RCPC_wilcoxon test.csv"))


pvals = dcast(df_C.index_optimPC, Metric+Model+Data ~Patient_stdz, value.var = "RCPC") %>% 
  group_by(Metric) %>% 
  summarise(pval_t=t.test(Var.4, `z-score`, paired=T)$p.value,
            pval_wilcox=wilcox.test(Var.4, `z-score`, paired=T)$p.value)
df_C.index_optimPC$pval <- pvals$pval_wilcox[match(df_C.index_optimPC$Metric, pvals$Metric)]



ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Patient_stdz, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           facet_grid(~Metric+paste(scales::pvalue(pval,
                                                   accuracy = 0.01,
                                                   decimal.mark = ".", 
                                                   add_p = TRUE)))+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_median_optimal RCPC_1.pdf"))


kruskal.test(RCPC ~ Metric , data = df_C.index_optimPC)


df_C.index_optimPC$pval_pair <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz==""], 
                                                     df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz==""],
                                                     p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric, pvals$Metric)]
df_C.index_optimPC$pval_pair[df_C.index_optimPC$Patient_stdz!=""] <- pairwise.wilcox.test(df_C.index_optimPC$RCPC[df_C.index_optimPC$Patient_stdz!=""], 
                                                                                          df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""],
                                                                                          p.adjust.method = "BH", paired=T)[[3]][c(1,7:12)][match(df_C.index_optimPC$Metric[df_C.index_optimPC$Patient_stdz!=""], pvals$Metric)]
df_C.index_optimPC$pval_pair <- ifelse(df_C.index_optimPC$pval_pair<0.05,
                                       scales::pvalue(df_C.index_optimPC$pval_pair ,
                                                      accuracy = 0.01,
                                                      decimal.mark = ".", 
                                                      add_p = TRUE),"ns")
df_C.index_optimPC$pval_pair[is.na(df_C.index_optimPC$pval_pair)] <- ""

ggexport(ggplot(df_C.index_optimPC, 
                aes(x=Metric, y=RCPC))+
           geom_violin(fill="black", col="white", alpha=0.2)+
           geom_jitter(fill=pal_jco()(2)[1],
                       position = position_jitter(w = 0.25, h = 0.15),
                       pch=21, alpha=0.9)+
           geom_text(data=df_C.index_optimPC %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,RCPC) %>% 
                       group_by(Patient_stdz,Metric) %>% top_n(1,paste(Model,Data)),
                     aes(label=pval_pair, y=RCPC+1),size=3)+facet_grid(~Patient_stdz)+
           facet_grid(~Patient_stdz)+
           scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
           theme_bw()+
           labs(x="", y="Optimal RCPC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),
         width = 10, height=3, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_median_optimal RCPC_2.pdf"))

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Median))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggexport(ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Median)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]),"white",scales::muted(pal_jco()(2)[2]))) +
           labs(x="RCPC", y="# features\nSD pre-selection", fill="Median C-index")+
           facet_grid(Patient_stdz ~Metric)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),
         width = 5, height=10, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_median_1.pdf"))


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, Median) %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  top_n(1, -RCPC) %>% 
  mutate(RCPC = "Optimal")

df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


df_C.index_1st <- df_C.index.stat %>% subset(RCPC == 0) %>% mutate(RCPC="0")


pvals <- rbind(df_C.index_optimPC, df_C.index_1st) %>% group_by(Patient_stdz,Metric) %>% 
  summarise(pval=pairwise.wilcox.test(Median, 
                                      RCPC,
                                      p.adjust.method = "BH", paired=T)[[3]][1,1],
            RCPC=max(RCPC),
            Median=max(Median))

ggexport(ggarrange(ggplot(rbind(df_C.index_optimPC, df_C.index_1st), 
                          aes(x=paste(Patient_stdz, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals,aes(y=Median+0.025,
                                              label=scales::pvalue(pval,
                                                                   accuracy = 0.01,
                                                                   decimal.mark = ".", 
                                                                   add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Metric)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"),
                   ggplot(rbind(df_C.index_optimPC, df_C.index_1st) %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)), 
                          aes(x=paste(Metric, RCPC), y=Median))+
                     geom_hline(yintercept=0.5,lty=2)+
                     geom_boxplot(aes(fill=Metric), alpha=1, width=0.5, outlier.shape = NA)+
                     geom_text(data=pvals %>% mutate(Metric = gsub("rAUC"," rAUC", Metric)),aes(y=Median+0.025,
                                                                                                label=scales::pvalue(pval,
                                                                                                                     accuracy = 0.01,
                                                                                                                     decimal.mark = ".", 
                                                                                                                     add_p = TRUE)), size=3, nudge_x =-0.25)+
                     facet_grid(~Patient_stdz)+
                     scale_fill_manual(values=pal_jco()(9))+
                     theme_bw()+
                     labs(x="", y="Best C-index mean")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank(),
                           legend.position = "top"), nrow=2, common.legend = T),
         width = 8.5, height=6.5, filename = paste0(path,"/Results_RCPC/RCPC testing_sts_random 100_median_optimal RCPC_3.pdf"))



df_C.index.stat <- df_C.index.stat[order(df_C.index.stat$Model, df_C.index.stat$Penalty, df_C.index.stat$Metric, df_C.index.stat$Data,df_C.index.stat$Patient_stdz),]
df_C.index.stat$Median_minus <- c(NA, df_C.index.stat$Median[-nrow(df_C.index.stat)])
df_C.index.stat$Median_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Median_delta <- df_C.index.stat$Median - df_C.index.stat$Median_minus
df_C.index.stat$Mean_minus <- c(NA,df_C.index.stat$Mean[-nrow(df_C.index.stat)])
df_C.index.stat$Mean_minus[df_C.index.stat$RCPC==0] <- NA
df_C.index.stat$Mean_delta <- df_C.index.stat$Mean - df_C.index.stat$Mean_minus


m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Mean))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Mean_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
  labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
  facet_grid(Patient_stdz ~Metric)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")

m <- df_C.index.stat %>% group_by(Data) %>% summarise(Mean=mean(Median))
df_C.index.stat$Data <- factor(df_C.index.stat$Data, levels = m$Data[order(m$Mean)])

ggplot(df_C.index.stat, aes(x = factor(RCPC), y = Data, fill = Median_delta)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = rev(hcl.colors(20, "RdYlGn"))) +
  labs(x="RCPC", y="Random sub-sample", fill="Mean C-index")+
  facet_grid(Patient_stdz ~Metric)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")


df_C.index_optimPC <- df_C.index.stat %>% 
  group_by(Model, Penalty, Metric, Data, Patient_stdz) %>% 
  summarise(RCPC_min.dMedian = RCPC[which.min(Median_delta)],
            RCPC_max.dMedian = RCPC[which.max(Median_delta)],
            RCPC_min.dMean = RCPC[which.min(Mean_delta)],
            RCPC_max.dMean = RCPC[which.max(Mean_delta)])
df_C.index_optimPC <- melt(df_C.index_optimPC, measure.vars=colnames(df_C.index_optimPC)[6:9])

ggplot(df_C.index_optimPC, 
       aes(x=Patient_stdz, y=value))+
  geom_violin(fill="black", alpha=0.3, size=0.5)+
  geom_line(aes(group=paste(Model,Data)),fill="black", alpha=0.3, size=0.5)+
  geom_point(fill="black", pch=21, alpha=0.5)+
  facet_grid(variable~Metric)+
  theme_bw()+
  labs(x="", y="Optimal RCPC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")

