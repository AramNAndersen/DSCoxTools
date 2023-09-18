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
path <- "C:/Users/EnserinkLab2019/Desktop/Clinical forecasting - revision 2023-07"

# Load ----
load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(paste0(path, "/Survival, clinical features and ELN2022 classifications_2023-07.RData"))

# Combination testing with clinical variables----
df_scores <- df_scores %>% subset(!grepl("re", Patient.ID))
X_rAUC <- dcast(data = df_scores, Patient.num ~ drug , value.var="rAUC")
X <- X_rAUC[match(rownames(Y),X_rAUC$Patient.num),]; X $Patient.num <- NULL
rownames(X_rAUC) <-X_rAUC$Patient.num;X_rAUC$Patient.num<-NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]

plot(unlist(X[,1:10]), unlist(X_rAUC[,1:10]))

df_scores$DSS3 <- df_scores$DSS3/100

df_clinical1 <- df_genetics[match(rownames(Y), rownames(df_genetics)),]
df_clinical1.1 <- df_prognostics_dummy[match(rownames(Y), rownames(df_prognostics_dummy)),grep("ELN|Age|Sex",colnames(df_prognostics_dummy))]
df_clinical1.2 <- df_prognostics_dummy[match(rownames(Y), rownames(df_prognostics_dummy)),grep("ELN|Age|Sex|Primary",colnames(df_prognostics_dummy))]
df_clinical2 <- cbind(df_clinical1, df_clinical1.1)
df_clinical2.1 <- cbind(df_clinical1, df_clinical1.2)

list_res_combination <- list()
list_res_clinical <- list()
list_res_drug <- list()
for(clin in c(1,1.1,1.2,2,2.1)){
  if(clin == 1){
    X2 <- df_clinical1
  }else if(clin == 1.1){
    X2 <- df_clinical1.1
  }else if(clin == 1.2){
    X2 <- df_clinical1.2
  }else if(clin == 2){
    X2 <- df_clinical2
  }else if(clin == 2.1){
    X2 <- df_clinical2.1
  }
  
  list_res_clinical[[paste0("clinical_",clin)]] <- Cox_forecasting_glmnet(X_data=X2,
                                                                          y_data=data.matrix(Y),
                                                                          alpha=c(0,0.4,1),
                                                                          lambda=c(exp(seq(-8,6, 0.1))),
                                                                          free_cores = 1,
                                                                          test.n= c(5,5),
                                                                          nfolds = nrow(Y),
                                                                          iter=200,
                                                                          log_AUC=2,
                                                                          Patient.Z=2,
                                                                          Drug.Z =2,
                                                                          RCPC=0)
  for(m in colnames(df_scores)[c(7,12)]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y), rownames(X)),]
    
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    
    for(tr in c(quantile(SDs, 1-c(40,60,80,100)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      list_res_combination[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=cbind(X2,X1),
                                                                                                                           y_data=data.matrix(Y),
                                                                                                                           alpha=c(0,0.4,1),
                                                                                                                           lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                           free_cores = 1,
                                                                                                                           test.n= c(5,5),
                                                                                                                           nfolds = nrow(Y),
                                                                                                                           iter=200,
                                                                                                                           log_AUC=2,
                                                                                                                           Patient.Z=2,
                                                                                                                           Drug.Z =2,
                                                                                                                           RCPC=0)
      
      cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
}

for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
  X <- X[match(rownames(Y), rownames(X)),]
  
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40,60,80,100)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_res_drug[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=X1,
                                                                                  y_data=data.matrix(Y),
                                                                                  alpha=c(0,0.4,1),
                                                                                  lambda=c(exp(seq(-8,6, 0.1))),
                                                                                  free_cores = 1,
                                                                                  test.n= c(5,5),
                                                                                  nfolds = nrow(Y),
                                                                                  iter=200,
                                                                                  log_AUC=2,
                                                                                  Patient.Z=2,
                                                                                  Drug.Z =2,
                                                                                  RCPC=0)
    
    
    cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
  }
}

save(list_res_combination, list_res_clinical, list_res_drug,
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction.RData"))
load( file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction.RData"))

df_C.index.alldata <- c(list_res_clinical,list_res_drug,list_res_combination)
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_C.index.alldata$Metric)
df_C.index.alldata$Pre_selection_N[grepl("clinical",df_C.index.alldata$Pre_selection_N)] <- ""
df_C.index.alldata$Data <- gsub(" cutoff.*|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric[grepl("clinical",df_C.index.alldata$Metric)] <- ""
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))
df_C.index.alldata$Pre_selection_N <- factor(df_C.index.alldata$Pre_selection_N, levels=unique(df_C.index.alldata$Pre_selection_N)[c(1,6:2)])

df_C.index.alldata$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1.1", "Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1", "Genetics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1_drug_", "Genetics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_C.index.alldata$Data)




df_C.index.alldata$Data <- factor(df_C.index.alldata$Data, 
                                  levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                           "Genetics/Diagnostics",
                                           "Genetics/Diagnostics+Clinical",
                                           
                                           "rAUC_log2","Genetics/rAUC_log2",
                                           "Diagnostics/rAUC_log2","Diagnostics+Clinical/rAUC_log2",
                                           "Genetics/Diagnostics/rAUC_log2","Genetics/Diagnostics+Clinical/rAUC_log2",
                                           
                                           "DSS3","Genetics/DSS3",
                                           "Diagnostics/DSS3","Diagnostics+Clinical/DSS3",
                                           "Genetics/Diagnostics/DSS3","Genetics/Diagnostics+Clinical/DSS3"))
df_C.index.alldata <- df_C.index.alldata[!grepl("Clinical", df_C.index.alldata$Data),]
df_C.index.stat <- df_C.index.alldata %>% group_by(Model, Metric, Data, Pre_selection_N) %>% summarise(Median=median(C_index_test),
                                                                                                       Mean=mean(C_index_test))
# ggplot(df_C.index.alldata %>% subset(!(Pre_selection_N %in% c(100,80,60))), 
#        aes(x=Data, y=C_index_test, fill=Model))+
#   geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#   geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#   geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
#   facet_grid(~Metric+Pre_selection_N, scales = "free")+
#   scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#   ylim(0,1)+
#   theme_bw()+
#   labs(x="", y="C-index")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank())
# 
# ggexport(ggarrange(ggplot(df_C.index.alldata %>% subset(!(Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_test, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5, outlier.alpha = 0)+
#                      facet_grid(~Metric+Pre_selection_N, scales = "free")+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()),
#                    ggplot(df_C.index.alldata %>% subset(!(Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_train, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5, outlier.alpha = 0)+
#                      facet_grid(~Metric+Pre_selection_N, scales = "free")+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
#          width = 10, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_1.pdf"))
# 
# 
# ggexport(ggarrange(ggplot(df_C.index.alldata %>% subset((Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_test, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
#                      facet_grid(~Metric+Pre_selection_N, scales = "free")+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()),
#                    ggplot(df_C.index.alldata %>% subset((Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_train, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
#                      facet_grid(~Metric+Pre_selection_N, scales = "free")+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
#          width = 11.5, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_2.pdf"))
# 
# 
# 
# 
# 
# ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                             subset(!(Pre_selection_N %in% c(100,80,60))) %>%
#                             subset((Data %in% c("Treatment","Genetics","Diagnostics",
#                                                 "Genetics/Diagnostics",
#                                                 "Genetics/Diagnostics/Treatment",
#                                                 "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
#                                                 "Genetics/Diagnostics/rAUC_log2",
#                                                 "Treatment/Genetics/Diagnostics/rAUC_log2",
#                                                 "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
#                                                 "Genetics/Diagnostics/DSS3",
#                                                 "Treatment/Genetics/Diagnostics/DSS3"))), 
#                           aes(x=Data, y=C_index_test, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
#                      facet_grid(Pre_selection_N~.)+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()),
#                    ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                             subset(!(Pre_selection_N %in% c(100,80,60))) %>%
#                             subset((Data %in% c("Treatment","Genetics","Diagnostics",
#                                                 "Genetics/Diagnostics",
#                                                 "Genetics/Diagnostics/Treatment",
#                                                 "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
#                                                 "Genetics/Diagnostics/rAUC_log2",
#                                                 "Treatment/Genetics/Diagnostics/rAUC_log2",
#                                                 "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
#                                                 "Genetics/Diagnostics/DSS3",
#                                                 "Treatment/Genetics/Diagnostics/DSS3"))), 
#                           aes(x=Data, y=C_index_train, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
#                      facet_grid(Pre_selection_N~.)+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
#          width = 6, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_3.pdf"))
# 
# 
# ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                             subset(!(Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_test, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5, outlier.alpha = 0)+
#                      facet_grid(Pre_selection_N~.)+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()),
#                    ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                             subset(!(Pre_selection_N %in% c(100,80,60))), 
#                           aes(x=Data, y=C_index_train, fill=Model))+
#                      geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#                      geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#                      geom_boxplot(alpha=0.5, outlier.alpha = 0)+
#                      facet_grid(Pre_selection_N~.)+
#                      scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
#                      ylim(0,1)+
#                      theme_bw()+
#                      labs(x="", y="C-index")+
#                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                            panel.grid.major = element_blank(),
#                            panel.grid.minor = element_blank(),
#                            strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
#          width = 7, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_4.pdf"))
# 
# 
# 
# ggexport(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                   subset(!(Pre_selection_N %in% c(100,80,60))) %>% group_by(Model, Data, Pre_selection_N) %>%
#                   summarise(Mean=mean(C_index_test),
#                             SD=sd(C_index_test)), 
#                 aes(x=Data, y=Mean, col=Pre_selection_N))+
#            ylim(0,1)+
#            geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#            geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#            geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
#            facet_grid(~Model)+
#            scale_color_manual(values = pal_jama()(7)[c(1,4)])+
#            theme_bw()+
#            labs(x="", y="C-index")+
#            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  strip.background = element_blank()),
#          width = 12.5, height=4, filename = paste0(path,"/Results_survival predictions/Combination data testing_5.pdf"))
# 
# 
# ggexport(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
#                   subset(!(Pre_selection_N %in% c(100,80,60))) %>% group_by(Model, Data, Pre_selection_N) %>%
#                   summarise(Mean=mean(C_index_train),
#                             SD=sd(C_index_train)), 
#                 aes(x=Data, y=Mean, col=Pre_selection_N))+
#            ylim(0,1)+
#            geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
#            geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
#            geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
#            facet_grid(~Model)+
#            scale_color_manual(values = pal_jama()(7)[c(1,4)])+
#            theme_bw()+
#            labs(x="", y="C-index")+
#            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  strip.background = element_blank()),
#          width = 12.5, height=4, filename = paste0(path,"/Results_survival predictions/Combination data testing_6.pdf"))


# Testing for other outcomes ----
rownames(df_survival) <- df_survival$ID
df_survival <- df_survival[match(rownames(Y), rownames(df_survival)),]
df_outcomes <- df_outcomes[match(rownames(Y), rownames(df_outcomes)),]


df_Y <- cbind(df_survival[,c("Persist._leukemia_post_first_ind._treatment", "Relapse")],
              df_outcomes[,c("Status", "one_year_survival","two_year_survival","three_year_survival","four_year_survival","five_year_survival")])

df_Y$Persist._leukemia_post_first_ind._treatment[df_Y$Persist._leukemia_post_first_ind._treatment=="Yes"] <- 1
df_Y$Persist._leukemia_post_first_ind._treatment[df_Y$Persist._leukemia_post_first_ind._treatment=="No"] <- 0
df_Y$Persist._leukemia_post_first_ind._treatment <- as.numeric(df_Y$Persist._leukemia_post_first_ind._treatment)
df_Y$Relapse[df_Y$Relapse=="Yes"] <- 1
df_Y$Relapse[df_Y$Relapse=="No"] <- 0
df_Y$Relapse <- as.numeric(df_Y$Relapse)
# df_Y$Allogenic_BM_transplant[df_Y$Allogenic_BM_transplant=="Yes"] <- 1
# df_Y$Allogenic_BM_transplant[df_Y$Allogenic_BM_transplant=="No"] <- 0
# df_Y$Allogenic_BM_transplant <- as.numeric(df_Y$Allogenic_BM_transplant)
# df_Y$Was_CR_achieved[df_Y$Was_CR_achieved=="Yes"] <- 1
# df_Y$Was_CR_achieved[df_Y$Was_CR_achieved=="No"] <- 0
# df_Y$Was_CR_achieved <- as.numeric(df_Y$Was_CR_achieved)


Binomial_forecasting_glmnet_kfold <- function(X_data,
                                              y_data,
                                              alpha=0,
                                              lambda=c(exp(seq(-4,6, 0.1))),
                                              kfold=7,
                                              log_AUC=c(1:2),
                                              Patient.Z=c(1:2),
                                              Drug.Z =c(1:2),
                                              RCPC=c(0,1,2,3,4)){
  
  df_accuracy.alldata <- c()
  
  for(a in alpha){
    for(Transform in c("log2(AUC)", "AUC")[log_AUC]){
      for(pt.st in c(TRUE, FALSE)[Patient.Z] ){
        for(cpc in RCPC){
          for(drug.st in c(TRUE, FALSE)[Drug.Z]){
            i=paste0(Transform, c("/Patient_stdz")[pt.st], c("/Drug_stdz")[drug.st],"/RCPC_", cpc, "/Penalty_", a)
            set.seed(1)
            X <- X_data
            l=min(data.frame(replace(X, X == 0, 1)))
            if(Transform=="log2(AUC)"){
              if(length(which(X <= 0))){
                X <- data.matrix(-log2(X + l))
              }else{
                X <- data.matrix(-log2(X))
              }
            }else{
              X <- data.matrix(X)
            }
            if(pt.st){
              X <- (X - rowMeans(X))/matrixStats::rowSds(X)
            }
            if(cpc>0){
              x_sd <- matrixStats::colSds(X)
              x_mean <- colMeans(X)
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
              svdz <- svd(t(X))
              X <- svdz$u[,-c(1:cpc)] %*% diag(svdz$d[-c(1:cpc)]) %*% t(svdz$v[,-c(1:cpc)])
              X <- X*x_sd + x_mean
              X <- t(X)
              rm(svdz, x_mean, x_sd)
            }
            
            if(drug.st){
              X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
              if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
            }
            
            Y <- y_data
            
            set.seed(1)
            ids <- sample(nrow(X),nrow(X))
            Y <- Y[ids]
            X <- X[ids,]
            ListX <- split(ids,1:kfold) 
            
            list_C_index <- list()
            for(b in 1:kfold){
              set.seed(b)
              
              train <- X[-ListX[[b]],]
              test <- X[ListX[[b]],]
              
              Y.train.loop <- Y[-ListX[[b]]]
              Y.test.loop <- Y[ListX[[b]]]
              
              
              model.loop <- glmnet(data.matrix(train), Y.train.loop, family = "binomial", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
              nfolds=nrow(train)
              model.loop.cv <- cv.glmnet(data.matrix(train), Y.train.loop,family = "binomial", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
              
              accuracy <- assess.glmnet(
                object=model.loop,
                newx = train,
                newy=Y.train.loop,
                s = model.loop.cv$lambda.min,
                family = c("binomial"))
              accuracy <- data.frame(bind_cols(accuracy))
              confusion <- confusion.glmnet(
                object=model.loop,
                newx = train,
                newy=Y.train.loop,
                s = model.loop.cv$lambda.min,
                family = c("binomial"))
              confusion <- data.frame(confusion);confusion$True <- as.character(confusion$True);confusion$Predicted <- as.character(confusion$Predicted)
              confusion <- reshape2::dcast(confusion, Predicted~True, value.var = "Freq")
              if(nrow(confusion)<2){
                if(confusion$Predicted==0){
                  confusion <- rbind(confusion[,2:3], 
                                     c(0,0))
                }else{
                  confusion <- rbind(c(0,0),
                                     confusion[,2:3])
                }
              }
              confusion <- data.matrix(confusion)
              accuracy$accuracy <- sum(diag(confusion))/sum(confusion)
              accuracy$TPR <- confusion[2,2]/sum(confusion[,2])
              accuracy$TNR <- confusion[1,1]/sum(confusion[,1])
              accuracy$FPR <- confusion[1,2]/sum(confusion[,2])
              accuracy$FNR <- confusion[2,1]/sum(confusion[,1])
              accuracy$F1 <- confusion[2,2]*2/(confusion[2,2]*2 + confusion[2,1] + confusion[1,2])
              accuracy$F1.0 <- confusion[1,1]*2/(confusion[1,1]*2 + confusion[2,1] + confusion[1,2])
              colnames(accuracy) <- paste0( colnames(accuracy), "_train")
              rm(confusion)
              
              accuracy.t <- assess.glmnet(
                object=model.loop,
                newx = data.matrix(test),
                newy=Y.test.loop,
                s = model.loop.cv$lambda.min,
                family = c("binomial"))
              accuracy.t <- data.frame(bind_cols(accuracy.t))
              confusion <- confusion.glmnet(
                object=model.loop,
                newx = data.matrix(test),
                newy=Y.test.loop,
                s = model.loop.cv$lambda.min,
                family = c("binomial"))
              confusion <- data.frame(confusion);confusion$True <- as.character(confusion$True);confusion$Predicted <- as.character(confusion$Predicted)
              confusion <- reshape2::dcast(confusion, Predicted~True, value.var = "Freq")
              if(nrow(confusion)<2){
                if(confusion$Predicted==0){
                  confusion <- rbind(confusion[,2:3], 
                                     c(0,0))
                }else{
                  confusion <- rbind(c(0,0),
                                     confusion[,2:3])
                }
              }
              confusion <- data.matrix(confusion)
              accuracy.t$accuracy <- sum(diag(confusion))/sum(confusion)
              accuracy.t$TPR <- confusion[2,2]/sum(confusion[,2])
              accuracy.t$TNR <- confusion[1,1]/sum(confusion[,1])
              accuracy.t$FPR <- confusion[1,2]/sum(confusion[,2])
              accuracy.t$FNR <- confusion[2,1]/sum(confusion[,1])
              accuracy.t$F1 <- confusion[2,2]*2/(confusion[2,2]*2 + confusion[2,1] + confusion[1,2])
              accuracy.t$F1.0 <- confusion[1,1]*2/(confusion[1,1]*2 + confusion[2,1] + confusion[1,2])
              colnames(accuracy.t) <- paste0( colnames(accuracy.t), "_test")
              rm(confusion)
              
              rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
              list_C_index[[b]] <-cbind(accuracy,accuracy.t)
            }
            gc()
            
            df_C_index <- do.call(rbind, list_C_index)
            df_C_index <- data.frame(df_C_index)
            df_C_index$Iteration <- 1:kfold
            df_C_index$ID <- i
            
            df_accuracy.alldata <- rbind(df_accuracy.alldata, df_C_index)
            
            cat("\nAnalysis completed for: ", i,"\n")
            rm(X,Y,  df_C_index, list_C_index, opts, pb)
            
          }
        }
      }
    }
  }
  return(list(Prediction_results = df_accuracy.alldata))
}

Outcome_prediction_res_kfold <- list()
for(d in colnames(df_Y)){
  y <- df_Y[,d]
  remove <- !is.na(y)
  y <- y[!is.na(y)]
  
  list_res_combination1 <- list()
  list_res_clinical1 <- list()
  list_res_drug1 <- list()
  for(clin in c(1,1.1,1.2,2,2.1)){
    if(clin == 1){
      X2 <- df_clinical1
    }else if(clin == 1.1){
      X2 <- df_clinical1.1
    }else if(clin == 1.2){
      X2 <- df_clinical1.2
    }else if(clin == 2){
      X2 <- df_clinical2
    }else if(clin == 2.1){
      X2 <- df_clinical2.1
    }
    
    X2 <- X2[remove,]
    list_res_clinical1[[paste0("clinical_",clin)]] <- Binomial_forecasting_glmnet_kfold(X_data=X2,
                                                                                        y_data=y,
                                                                                        alpha=c(0,0.4,1),
                                                                                        lambda=c(exp(seq(-8,6, 0.1))),
                                                                                        kfold=5,
                                                                                        log_AUC=2,
                                                                                        Patient.Z=2,
                                                                                        Drug.Z =2,
                                                                                        RCPC=0)
    for(m in colnames(df_scores)[c(7,12)]){
      X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
        X1 <- X[,SDs>=tr]
        X1 <- X1[remove,]
        
        list_res_combination1[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=cbind(X2,X1),
                                                                                                                                         y_data=y,
                                                                                                                                         alpha=c(0,0.4,1),
                                                                                                                                         lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                                         kfold=5,
                                                                                                                                         log_AUC=2,
                                                                                                                                         Patient.Z=2,
                                                                                                                                         Drug.Z =2,
                                                                                                                                         RCPC=0)
        
        cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
      }
    }
  }
  
  for(m in colnames(df_scores)[c(7,12)]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y), rownames(X)),]
    
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      X1 <- X1[remove,]
      list_res_drug1[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=X1,
                                                                                                y_data=y,
                                                                                                alpha=c(0,0.4,1),
                                                                                                lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                kfold=5,
                                                                                                log_AUC=2,
                                                                                                Patient.Z=2,
                                                                                                Drug.Z =2,
                                                                                                RCPC=0)
      
      
      cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
  
  Outcome_prediction_res_kfold[[d]] <- list(list_res_combination=list_res_combination1,
                                            list_res_clinical=list_res_clinical1,
                                            list_res_drug=list_res_drug1)
}

save(Outcome_prediction_res_kfold,
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold.RData"))

load(file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold.RData"))

df_prediction_outcomes <- list()
for(d in names(Outcome_prediction_res_kfold)){
  df_C.index.alldata <- c(Outcome_prediction_res_kfold[[d]]$list_res_combination,
                          Outcome_prediction_res_kfold[[d]]$list_res_clinical,
                          Outcome_prediction_res_kfold[[d]]$list_res_drug)
  for(m in names(df_C.index.alldata)){
    df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$Prediction_results
  }
  df_prediction_outcomes[[d]] <- bind_rows(df_C.index.alldata, .id="Metric")
}
df_prediction_outcomes <- bind_rows(df_prediction_outcomes, .id="Outcome")

df_prediction_outcomes$Penalty <- gsub(".*Penalty_","", df_prediction_outcomes$ID)
df_prediction_outcomes$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Pre_selection_N[grepl("clinical",df_prediction_outcomes$Pre_selection_N)] <- ""
df_prediction_outcomes$Data <- gsub(" cutoff.*|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric[grepl("clinical",df_prediction_outcomes$Metric)] <- ""
df_prediction_outcomes$Metric <- factor(df_prediction_outcomes$Metric, levels=unique(df_prediction_outcomes$Metric))
df_prediction_outcomes$Model <- ifelse(df_prediction_outcomes$Penalty==1, "Lasso", "Elastic net")
df_prediction_outcomes$Model[df_prediction_outcomes$Penalty==0] <- "Ridge"
df_prediction_outcomes$Model <- factor(df_prediction_outcomes$Model, levels=unique(df_prediction_outcomes$Model))
df_prediction_outcomes$Pre_selection_N <- factor(df_prediction_outcomes$Pre_selection_N, levels=unique(df_prediction_outcomes$Pre_selection_N)[c(1,6:2)])
df_prediction_outcomes$Outcome <- factor(df_prediction_outcomes$Outcome, levels=unique(df_prediction_outcomes$Outcome))

df_prediction_outcomes$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1.1", "Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1", "Genetics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1_drug_", "Genetics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_prediction_outcomes$Data)




df_prediction_outcomes$Data <- factor(df_prediction_outcomes$Data, 
                                  levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                           "Genetics/Diagnostics",
                                           "Genetics/Diagnostics+Clinical",
                                           
                                           "rAUC_log2","Genetics/rAUC_log2",
                                           "Diagnostics/rAUC_log2","Diagnostics+Clinical/rAUC_log2",
                                           "Genetics/Diagnostics/rAUC_log2","Genetics/Diagnostics+Clinical/rAUC_log2",
                                           
                                           "DSS3","Genetics/DSS3",
                                           "Diagnostics/DSS3","Diagnostics+Clinical/DSS3",
                                           "Genetics/Diagnostics/DSS3","Genetics/Diagnostics+Clinical/DSS3"))

colnames(df_prediction_outcomes)

df_prediction_outcomes <- df_prediction_outcomes[!grepl("Clinical",df_prediction_outcomes$Data),]
df_prediction.stat <- melt(df_prediction_outcomes, measure.vars = colnames(df_prediction_outcomes)[c(3:26)]) %>%
  group_by(Pre_selection_N, Model, Data, Outcome, variable) %>%
  summarise(Mean=mean(value),SD=sd(value)/sqrt(length(value)))
df_prediction.stat$Set <- gsub(".*_", "", df_prediction.stat$variable)
df_prediction.stat$variable <- gsub("_.*", "", df_prediction.stat$variable)

ggexport(ggplot(df_prediction.stat %>% 
                  mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                  subset(Model=="Ridge") %>%
                  subset(grepl("leukemia|Relapse|Status|one_year",Outcome)) %>%
                  subset(grepl("auc",variable)), 
                aes(x=Data, y=Mean, col=Pre_selection_N))+
           ylim(0,1)+
           geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
           geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
           geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
           facet_grid(variable+Set~Outcome)+
           scale_color_manual(values = pal_jama()(7)[c(1,4)])+
           theme_bw()+
           labs(x="", y="ROC AUC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),         
         width = 12, height=5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_1.pdf"))

ggexport(ggplot(df_prediction.stat %>% 
                  mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                  subset(Model=="Ridge") %>%
                  subset(grepl("leukemia|Relapse|Status|one_year",Outcome)) %>%
                  subset(grepl("deviance",variable)), 
                aes(x=Data, y=Mean, col=Pre_selection_N))+
           geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
           facet_grid(variable+Set~Outcome)+
           scale_color_manual(values = pal_jama()(7)[c(1,4)])+
           theme_bw()+
           labs(x="", y="Deviance")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),         
         width = 12, height=5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_2.pdf"))


ggexport(ggplot(df_prediction.stat %>% 
                  mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                  subset(Model=="Ridge") %>%
                  subset(!grepl("leukemia|Relapse|Status|one_year",Outcome)) %>%
                  subset(grepl("auc",variable)), 
                aes(x=Data, y=Mean, col=Pre_selection_N))+
           ylim(0,1)+
           geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
           geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
           geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
           facet_grid(variable+Set~Outcome)+
           scale_color_manual(values = pal_jama()(7)[c(1,4)])+
           theme_bw()+
           labs(x="", y="ROC AUC")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),         
         width = 12, height=5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_1b.pdf"))

ggexport(ggplot(df_prediction.stat %>% 
                  mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                  subset(Model=="Ridge") %>%
                  subset(!grepl("leukemia|Relapse|Status|one_year",Outcome)) %>%
                  subset(grepl("deviance",variable)), 
                aes(x=Data, y=Mean, col=Pre_selection_N))+
           geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD),alpha=0.5)+
           facet_grid(variable+Set~Outcome)+
           scale_color_manual(values = pal_jama()(7)[c(1,4)])+
           theme_bw()+
           labs(x="", y="Deviance")+
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank()),         
         width = 12, height=5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_2b.pdf"))



df_prediction.stat_combi <- bind_rows(df_C.index.stat %>% 
                                        subset(Pre_selection_N %in% c(349,40, "")) %>%
                                        mutate(Outcome = "Survival (CoxPH)",
                                               Pre_selection_N = gsub("349","",Pre_selection_N)),
                                      df_prediction.stat  %>% 
                                       # subset(grepl("leukemia|Relapse|Status|one_year",Outcome)) %>%
                                        subset(grepl("auc",variable) & Set == "test") %>% 
                                        mutate(Outcome = gsub("Status", "Survival (class)", Outcome),
                                               Pre_selection_N = gsub("349","",Pre_selection_N)))
df_prediction.stat_combi$Outcome <- factor(df_prediction.stat_combi$Outcome, levels = rev(unique(df_prediction.stat_combi$Outcome)[c(1,4:9,2:3)]))



ggexport(ggplot(df_prediction.stat_combi, aes(x = Data, y = Outcome, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
           labs(x="", y="", fill="Mean C-index/ROC AUC")+
           facet_grid(Model~Pre_selection_N)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),         
         width = 12, height=6.5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_summary heatmap_1.pdf"))



ggexport(ggplot(df_prediction.stat_combi, aes(x = Data, y = Outcome, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
           labs(x="", y="", fill="Mean C-index/ROC AUC")+
           facet_grid(Pre_selection_N~Model)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),         
         width = 14, height=6.5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_summary heatmap_2.pdf"))


ggexport(ggplot(df_prediction.stat_combi %>% subset(Model == "Ridge"), aes(x = Data, y = Outcome, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
           labs(x="", y="", fill="Mean C-index/ROC AUC")+
           facet_grid(Pre_selection_N~.)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),         
         width = 7.5, height=5.5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_Ridge heatmap.pdf"))


ggexport(ggplot(df_prediction.stat_combi %>% subset(Model == "Lasso"), aes(x = Data, y = Outcome, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
           labs(x="", y="", fill="Mean C-index/ROC AUC")+
           facet_grid(Pre_selection_N~.)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),         
         width = 7.5, height=5.5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_Lasso heatmap.pdf"))


ggexport(ggplot(df_prediction.stat_combi %>% subset(Model == "Elastic net"), aes(x = Data, y = Outcome, fill = Mean)) +
           geom_tile(color = "black") +
           scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
           labs(x="", y="", fill="Mean C-index/ROC AUC")+
           facet_grid(Pre_selection_N~.)+
           coord_fixed()+
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 strip.background = element_blank(),
                 legend.position = "top"),         
         width = 7.5, height=5.5, 
         filename = paste0(path,"/Results_survival predictions/Combination data testing_clinical outcomes_Elastic net heatmap.pdf"))

# Testing for other outcomes - all metrics ----
load(file = paste0(path,"/Drug sensitivity metrics and dose response data_2023-07.RData"))

df_scores <- df_scores %>% subset(!grepl("re", Patient.ID))

X_rAUC <- dcast(data = df_scores, Patient.num ~ drug , value.var="rAUC")
rownames(X_rAUC) <- as.character(X_rAUC$Patient.num); X_rAUC$Patient.num <- NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]


df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$rAUC <- 1-df_scores$rAUC

Outcome_prediction_res_kfold_all <- list()
for(d in colnames(df_Y)){
  y <- df_Y[,d]
  remove <- !is.na(y)
  y <- y[!is.na(y)]
  
  list_res_combination1 <- list()
  list_res_clinical1 <- list()
  list_res_drug1 <- list()
  for(clin in c(1,1.1,1.2,2,2.1)){
    if(clin == 1){
      X2 <- df_clinical1
    }else if(clin == 1.1){
      X2 <- df_clinical1.1
    }else if(clin == 1.2){
      X2 <- df_clinical1.2
    }else if(clin == 2){
      X2 <- df_clinical2
    }else if(clin == 2.1){
      X2 <- df_clinical2.1
    }
    
    X2 <- X2[remove,]
    list_res_clinical1[[paste0("clinical_",clin)]] <- Binomial_forecasting_glmnet_kfold(X_data=X2,
                                                                                        y_data=y,
                                                                                        alpha=c(0,0.4,1),
                                                                                        lambda=c(exp(seq(-8,6, 0.1))),
                                                                                        kfold=5,
                                                                                        log_AUC=2,
                                                                                        Patient.Z=2,
                                                                                        Drug.Z =2,
                                                                                        RCPC=0)
    for(m in colnames(df_scores)[c(6:14)]){
      X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      
      if(grepl("EC50",m)){
        X <- log10(X)
        X[,] <- apply(X,2, function(x) x - mean(x))
      }
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
        X1 <- X[,SDs>=tr]
        X1 <- X1[remove,]
        
        list_res_combination1[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=cbind(X2,X1),
                                                                                                                                         y_data=y,
                                                                                                                                         alpha=c(0,0.4,1),
                                                                                                                                         lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                                         kfold=5,
                                                                                                                                         log_AUC=2,
                                                                                                                                         Patient.Z=2,
                                                                                                                                         Drug.Z =2,
                                                                                                                                         RCPC=0)
        
        cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
      }
    }
  }
  
  for(m in colnames(df_scores)[c(6:14)]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y), rownames(X)),]
    if(grepl("EC50",m)){
      X <- log10(X)
      X[,] <- apply(X,2, function(x) x - mean(x))
    }
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      X1 <- X1[remove,]
      list_res_drug1[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=X1,
                                                                                                y_data=y,
                                                                                                alpha=c(0,0.4,1),
                                                                                                lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                kfold=5,
                                                                                                log_AUC=2,
                                                                                                Patient.Z=2,
                                                                                                Drug.Z =2,
                                                                                                RCPC=0)
      
      
      cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
  
  Outcome_prediction_res_kfold_all[[d]] <- list(list_res_combination=list_res_combination1,
                                            list_res_clinical=list_res_clinical1,
                                            list_res_drug=list_res_drug1)
}

save(Outcome_prediction_res_kfold_all,
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold_all scores.RData"))


load(file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))
df.hillAUC <- df.hillAUC %>% subset(!grepl("re", Patient.ID))

df.hillAUC$rAUC <- 1-df.hillAUC$rAUC
df.hillAUC$hill_AUC <- 1-df.hillAUC$hill_AUC
df.hillAUC$hill0_AUC <- 1-df.hillAUC$hill0_AUC
df.hillAUC$hillB_AUC <- 1-df.hillAUC$hillB_AUC

Outcome_prediction_res_kfold_hill <- list()
for(d in colnames(df_Y)){
  y <- df_Y[,d]
  remove <- !is.na(y)
  y <- y[!is.na(y)]
  
  list_res_combination1 <- list()
  list_res_clinical1 <- list()
  list_res_drug1 <- list()
  for(clin in c(1,1.1,1.2,2,2.1)){
    if(clin == 1){
      X2 <- df_clinical1
    }else if(clin == 1.1){
      X2 <- df_clinical1.1
    }else if(clin == 1.2){
      X2 <- df_clinical1.2
    }else if(clin == 2){
      X2 <- df_clinical2
    }else if(clin == 2.1){
      X2 <- df_clinical2.1
    }
    
    X2 <- X2[remove,]
    list_res_clinical1[[paste0("clinical_",clin)]] <- Binomial_forecasting_glmnet_kfold(X_data=X2,
                                                                                        y_data=y,
                                                                                        alpha=c(0,0.4,1),
                                                                                        lambda=c(exp(seq(-8,6, 0.1))),
                                                                                        kfold=5,
                                                                                        log_AUC=2,
                                                                                        Patient.Z=2,
                                                                                        Drug.Z =2,
                                                                                        RCPC=0)
    for(m in colnames(df.hillAUC)[c(6:8,10:12)]){
      X <- dcast(data = df.hillAUC, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
        X1 <- X[,SDs>=tr]
        X1 <- X1[remove,]
        
        list_res_combination1[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=cbind(X2,X1),
                                                                                                                                         y_data=y,
                                                                                                                                         alpha=c(0,0.4,1),
                                                                                                                                         lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                                         kfold=5,
                                                                                                                                         log_AUC=2,
                                                                                                                                         Patient.Z=2,
                                                                                                                                         Drug.Z =2,
                                                                                                                                         RCPC=0)
        
        cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
      }
    }
  }
  
  for(m in colnames(df.hillAUC)[c(6:8,10:12)]){
    X <- dcast(data = df.hillAUC, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y), rownames(X)),]
    
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      X1 <- X1[remove,]
      list_res_drug1[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Binomial_forecasting_glmnet_kfold(X_data=X1,
                                                                                                y_data=y,
                                                                                                alpha=c(0,0.4,1),
                                                                                                lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                kfold=5,
                                                                                                log_AUC=2,
                                                                                                Patient.Z=2,
                                                                                                Drug.Z =2,
                                                                                                RCPC=0)
      
      
      cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
  
  Outcome_prediction_res_kfold_hill[[d]] <- list(list_res_combination=list_res_combination1,
                                                list_res_clinical=list_res_clinical1,
                                                list_res_drug=list_res_drug1)
}

save(Outcome_prediction_res_kfold_all,####
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold_hill fit.RData"))

load(file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold_hill fit.RData"))
load(file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_various clinical outcomes prediction_kfold_all scores.RData"))

df_prediction_outcomes <- list()
for(d in names(Outcome_prediction_res_kfold_hill)){
  df_C.index.alldata <- c(Outcome_prediction_res_kfold_hill[[d]]$list_res_combination,
                          Outcome_prediction_res_kfold_hill[[d]]$list_res_clinical,
                          Outcome_prediction_res_kfold_hill[[d]]$list_res_drug)
  for(m in names(df_C.index.alldata)){
    df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$Prediction_results
  }
  df_prediction_outcomes[[d]] <- bind_rows(df_C.index.alldata, .id="Metric")
}
df_prediction_outcomes <- bind_rows(df_prediction_outcomes, .id="Outcome")

df_prediction_outcomes$Penalty <- gsub(".*Penalty_","", df_prediction_outcomes$ID)
df_prediction_outcomes$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Pre_selection_N[grepl("clinical",df_prediction_outcomes$Pre_selection_N)] <- ""
df_prediction_outcomes$Data <- gsub(" cutoff.*|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric[grepl("clinical",df_prediction_outcomes$Metric)] <- ""
df_prediction_outcomes$Metric <- factor(df_prediction_outcomes$Metric, levels=unique(df_prediction_outcomes$Metric))
df_prediction_outcomes$Model <- ifelse(df_prediction_outcomes$Penalty==1, "Lasso", "Elastic net")
df_prediction_outcomes$Model[df_prediction_outcomes$Penalty==0] <- "Ridge"
df_prediction_outcomes$Model <- factor(df_prediction_outcomes$Model, levels=unique(df_prediction_outcomes$Model))
df_prediction_outcomes$Pre_selection_N <- factor(df_prediction_outcomes$Pre_selection_N, levels=unique(df_prediction_outcomes$Pre_selection_N))
df_prediction_outcomes$Outcome <- factor(df_prediction_outcomes$Outcome, levels=unique(df_prediction_outcomes$Outcome))

df_prediction_outcomes$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1.1", "Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1", "Genetics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1_drug_", "Genetics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_prediction_outcomes$Data)


comb <- expand.grid(clin=c("","Genetics",
                      "Diagnostics","Diagnostics+Clinical",
                      "Genetics/Diagnostics","Genetics/Diagnostics+Clinical"),
                    metric=unique(df_prediction_outcomes$Metric))
comb$data <- paste(comb$clin, comb$metric, sep="/")
comb$data <- gsub("^\\/|.*\\/$", "", comb$data)

df_prediction_outcomes$Data <- factor(df_prediction_outcomes$Data, 
                                      levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                               "Genetics/Diagnostics",
                                               "Genetics/Diagnostics+Clinical",
                                               comb$data[comb$data != ""] ))


df_prediction_outcomes <- df_prediction_outcomes[!grepl("Clinical",df_prediction_outcomes$Data),]
df_prediction.stat <- melt(df_prediction_outcomes, measure.vars = colnames(df_prediction_outcomes)[c(3:26)]) %>%
  group_by(Pre_selection_N, Model, Data, Outcome, variable) %>%
  summarise(Mean=mean(value),SD=sd(value)/sqrt(length(value)))
df_prediction.stat$Set <- gsub(".*_", "", df_prediction.stat$variable)
df_prediction.stat$variable <- gsub("_.*", "", df_prediction.stat$variable)


df_prediction_outcomes_hill <- df_prediction_outcomes
df_prediction.stat_hill <- df_prediction.stat



df_prediction_outcomes <- list()
for(d in names(Outcome_prediction_res_kfold_all)){
  df_C.index.alldata <- c(Outcome_prediction_res_kfold_all[[d]]$list_res_combination,
                          Outcome_prediction_res_kfold_all[[d]]$list_res_clinical,
                          Outcome_prediction_res_kfold_all[[d]]$list_res_drug)
  for(m in names(df_C.index.alldata)){
    df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$Prediction_results
  }
  df_prediction_outcomes[[d]] <- bind_rows(df_C.index.alldata, .id="Metric")
}
df_prediction_outcomes <- bind_rows(df_prediction_outcomes, .id="Outcome")

df_prediction_outcomes$Penalty <- gsub(".*Penalty_","", df_prediction_outcomes$ID)
df_prediction_outcomes$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Pre_selection_N[grepl("clinical",df_prediction_outcomes$Pre_selection_N)] <- ""
df_prediction_outcomes$Data <- gsub(" cutoff.*|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_prediction_outcomes$Metric)
df_prediction_outcomes$Metric[grepl("clinical",df_prediction_outcomes$Metric)] <- ""
df_prediction_outcomes$Metric <- factor(df_prediction_outcomes$Metric, levels=unique(df_prediction_outcomes$Metric))
df_prediction_outcomes$Model <- ifelse(df_prediction_outcomes$Penalty==1, "Lasso", "Elastic net")
df_prediction_outcomes$Model[df_prediction_outcomes$Penalty==0] <- "Ridge"
df_prediction_outcomes$Model <- factor(df_prediction_outcomes$Model, levels=unique(df_prediction_outcomes$Model))
df_prediction_outcomes$Pre_selection_N <- factor(df_prediction_outcomes$Pre_selection_N, levels=unique(df_prediction_outcomes$Pre_selection_N))
df_prediction_outcomes$Outcome <- factor(df_prediction_outcomes$Outcome, levels=unique(df_prediction_outcomes$Outcome))


df_prediction_outcomes$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1.1", "Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_1", "Genetics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_1_drug_", "Genetics/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_prediction_outcomes$Data)
df_prediction_outcomes$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_prediction_outcomes$Data)


comb <- expand.grid(clin=c("","Genetics",
                           "Diagnostics","Diagnostics+Clinical",
                           "Genetics/Diagnostics","Genetics/Diagnostics+Clinical"),
                    metric=unique(df_prediction_outcomes$Metric))
comb$data <- paste(comb$clin, comb$metric, sep="/")
comb$data <- gsub("^\\/|.*\\/$", "", comb$data)

df_prediction_outcomes$Data <- factor(df_prediction_outcomes$Data, 
                                      levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                               "Genetics/Diagnostics",
                                               "Genetics/Diagnostics+Clinical",
                                               comb$data[comb$data != ""] ))
colnames(df_prediction_outcomes)

df_prediction_outcomes <- df_prediction_outcomes[!grepl("Clinical",df_prediction_outcomes$Data),]
df_prediction.stat <- melt(df_prediction_outcomes, measure.vars = colnames(df_prediction_outcomes)[c(3:26)]) %>%
  group_by(Pre_selection_N, Model, Data, Outcome, variable) %>%
  summarise(Mean=mean(value),SD=sd(value)/sqrt(length(value)))
df_prediction.stat$Set <- gsub(".*_", "", df_prediction.stat$variable)
df_prediction.stat$variable <- gsub("_.*", "", df_prediction.stat$variable)



df_prediction_outcomes_hill
df_prediction.stat$Pre_selection_N


ggplot(df_prediction.stat %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)), 
       aes(x=Data, y=Mean))+
  ylim(0,1)+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_violin(fill="gray", width=1.25)+
  geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD, size=factor(Pre_selection_N), col=Model), pch=16,
                  alpha=0.5, position = position_dodge2(width = 0.5, preserve = "total"))+
  facet_wrap(~Outcome)+
  scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
  scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue"))+
  scale_size_manual(values=c(0.2,0.4,0.6,0.8,1)/2)+
  theme_bw()+
  labs(x="", y="ROC AUC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(width = 10, height=7, 
         filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_1.pdf"))
ggplot(df_prediction.stat %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)), 
       aes(x=Data, y=Mean))+
  ylim(0,1)+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_violin(fill="gray", width=1.25)+
  geom_point(aes(size=factor(Pre_selection_N), col=Model), pch=16,
                  alpha=0.5)+
  facet_wrap(~Outcome)+
  scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
  scale_size_manual(values=c(0.2,0.4,0.6,0.8,1)*3)+
  theme_bw()+
  labs(x="", y="ROC AUC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(width = 10, height=7, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_2.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status|Relaps", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)), 
       aes(x=Data, y=Mean))+
  ylim(0,1)+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_violin(fill="gray", width=1.25)+
  geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD, size=factor(Pre_selection_N), col=Model),
                  alpha=0.5, position = position_dodge2(width = 0.5, preserve = "total"))+
  facet_wrap(~Outcome)+
  scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
  scale_size_manual(values=c(0.2,0.4,0.6,0.8,1)/2)+
  theme_bw()+
  labs(x="", y="ROC AUC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(width = 12, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_3.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status|Relaps", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)), 
       aes(x=Data, y=Mean))+
  ylim(0,1)+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_violin(fill="gray", width=1.25)+
  geom_point(aes(size=factor(Pre_selection_N), col=Model),
             alpha=0.5)+
  facet_wrap(~Outcome)+
  scale_color_manual(values = c("darkorange", "darkgreen", "darkblue"))+
  scale_size_manual(values=c(0.2,0.4,0.6,0.8,1)*3)+
  theme_bw()+
  labs(x="", y="ROC AUC")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(width = 12, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_4.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  facet_grid(Pre_selection_N~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 10, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_5.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  #facet_grid(Pre_selection_N~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 6.5, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_6.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  facet_grid(Model+Pre_selection_N~. )+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 13, height=10, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_7.pdf"))

ggplot(rbind(df_prediction.stat,df_prediction.stat_hill %>% subset(grepl("hill",Data))) %>% 
         subset(!grepl("Status", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  #facet_grid(Model+Pre_selection_N~. )+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 13, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_8.pdf"))


ggplot(rbind(df_prediction.stat) %>% 
         subset(!grepl("Status|four", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  #facet_grid(Model+Pre_selection_N~. )+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 11, height=5, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_9.pdf"))

ggplot(rbind(df_prediction.stat) %>% 
         subset(!grepl("Gen|Diag",Data) & !grepl("Status|four", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  #facet_grid(Pre_selection_N~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 6.5, height=4, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_10.pdf"))


ggplot(rbind(df_prediction.stat) %>% 
         subset(!grepl("\\/",Data) & !grepl("Status|four", Outcome)) %>%
         subset(grepl("auc",variable) & !grepl("train",Set)) %>%
         mutate(Pre_selection_N= as.numeric(ifelse(as.character(Pre_selection_N)=="","349", as.character(Pre_selection_N)))), aes(x = Data, y = Outcome, fill = Mean)) +
  geom_tile(color = "black") +
  #scale_fill_gradientn(colors = c(scales::muted(pal_jco()(2)[1]), "white", scales::muted(pal_jco()(2)[2]))) +
  scale_fill_gradient2(low=scales::muted(pal_jco()(2)[1]), mid = "white", high=scales::muted(pal_jco()(2)[2]), midpoint=0.5) +
  labs(x="", y="", fill="Mean C-index/ROC AUC")+
  #facet_grid(Pre_selection_N~Model)+
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "top")
ggsave(width = 8, height=4, 
       filename = paste0(path,"/Results_survival predictions/Different clinical outcomes testing_11.pdf"))





# Combination testing with clinical variables - all metrics and sts----
load(file = paste0(path,"/Drug sensitivity metrics and dose response data_2023-07.RData"))

df_scores <- df_scores %>% subset(!grepl("re", Patient.ID))

X_rAUC <- dcast(data = df_scores, Patient.num ~ drug , value.var="rAUC")
rownames(X_rAUC) <- as.character(X_rAUC$Patient.num); X_rAUC$Patient.num <- NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]


df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$rAUC <- 1-df_scores$rAUC

df_clinical1 <- df_genetics[match(rownames(Y), rownames(df_genetics)),]
df_clinical1.1 <- df_prognostics_dummy[match(rownames(Y), rownames(df_prognostics_dummy)),grep("ELN|Age|Sex",colnames(df_prognostics_dummy))]
df_clinical1.2 <- df_prognostics_dummy[match(rownames(Y), rownames(df_prognostics_dummy)),grep("ELN|Age|Sex|Primary",colnames(df_prognostics_dummy))]
df_clinical2 <- cbind(df_clinical1, df_clinical1.1)
df_clinical2.1 <- cbind(df_clinical1, df_clinical1.2)

list_res_combination <- list()
list_res_clinical <- list()
list_res_drug <- list()
for(clin in c(1,1.1,1.2,2,2.1)){
  if(clin == 1){
    X2 <- df_clinical1
  }else if(clin == 1.1){
    X2 <- df_clinical1.1
  }else if(clin == 1.2){
    X2 <- df_clinical1.2
  }else if(clin == 2){
    X2 <- df_clinical2
  }else if(clin == 2.1){
    X2 <- df_clinical2.1
  }
  
  list_res_clinical[[paste0("clinical_",clin)]] <- Cox_forecasting_glmnet(X_data=X2,
                                                                          y_data=data.matrix(Y),
                                                                          alpha=c(0,0.4,1),
                                                                          lambda=c(exp(seq(-8,6, 0.1))),
                                                                          free_cores = 2,
                                                                          test.n= c(5,5),
                                                                          nfolds = nrow(Y),
                                                                          iter=200,
                                                                          log_AUC=2,
                                                                          Patient.Z=2,
                                                                          Drug.Z =2,
                                                                          RCPC=0)
  for(m in colnames(df_scores)[c(6:12)]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y), rownames(X)),]
    
    if(grepl("EC50",m)){
      X <- log10(X)
      X[,] <- apply(X,2, function(x) x - mean(x))
    }
    
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    
    for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      list_res_combination[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=cbind(X2,X1),
                                                                                                                           y_data=data.matrix(Y),
                                                                                                                           alpha=c(0,0.4,1),
                                                                                                                           lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                           free_cores = 2,
                                                                                                                           test.n= c(5,5),
                                                                                                                           nfolds = nrow(Y),
                                                                                                                           iter=200,
                                                                                                                           log_AUC=2,
                                                                                                                           Patient.Z=2,
                                                                                                                           Drug.Z =2,
                                                                                                                           RCPC=0)
      
      cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
}

for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
  X <- X[match(rownames(Y), rownames(X)),]
  
  if(grepl("EC50",m)){
    X <- log10(X)
    X[,] <- apply(X,2, function(x) x - mean(x))
  }
  
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_res_drug[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=X1,
                                                                                  y_data=data.matrix(Y),
                                                                                  alpha=c(0,0.4,1),
                                                                                  lambda=c(exp(seq(-8,6, 0.1))),
                                                                                  free_cores = 2,
                                                                                  test.n= c(5,5),
                                                                                  nfolds = nrow(Y),
                                                                                  iter=200,
                                                                                  log_AUC=2,
                                                                                  Patient.Z=2,
                                                                                  Drug.Z =2,
                                                                                  RCPC=0)
    
    
    cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
  }
}

save(list_res_combination, list_res_clinical, list_res_drug,
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction_all metrics.RData"))

load( file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction_all metrics.RData"))

df_C.index.alldata <- c(list_res_clinical,list_res_drug,list_res_combination)
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_C.index.alldata$Metric)
df_C.index.alldata$Pre_selection_N[grepl("clinical",df_C.index.alldata$Pre_selection_N)] <- ""
df_C.index.alldata$Data <- gsub(" cutoff.*|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric[grepl("clinical",df_C.index.alldata$Metric)] <- ""
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))
df_C.index.alldata$Pre_selection_N <- factor(df_C.index.alldata$Pre_selection_N, levels=unique(df_C.index.alldata$Pre_selection_N)[c(1,6:2)])

df_C.index.alldata$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1.1", "Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1", "Genetics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1_drug_", "Genetics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_C.index.alldata$Data)


comb <- expand.grid(clin=c("","Genetics",
                           "Diagnostics","Diagnostics+Clinical",
                           "Genetics/Diagnostics","Genetics/Diagnostics+Clinical"),
                    metric=unique(df_C.index.alldata$Metric))
comb$data <- paste(comb$clin, comb$metric, sep="/")
comb$data <- gsub("^\\/|.*\\/$", "", comb$data)

df_C.index.alldata$Data <- factor(df_C.index.alldata$Data, 
                                      levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                               "Genetics/Diagnostics",
                                               "Genetics/Diagnostics+Clinical",
                                               comb$data[comb$data != ""] ))

df_C.index.alldata <- df_C.index.alldata[!grepl("Clinical", df_C.index.alldata$Data),]
df_C.index.stat <- df_C.index.alldata %>% group_by(Model, Metric, Data, Pre_selection_N) %>% summarise(Median=median(C_index_test),
                                                                                                       Mean=mean(C_index_test))


ggexport(ggarrange(ggplot(df_C.index.alldata %>% 
                            mutate(Pre_selection_N = ifelse(Pre_selection_N==349,"", (as.character(Pre_selection_N)))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5, outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~Metric, scales = "free", space='free')+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=1, ncol=1),
         width = 10, height=6, filename = paste0(path,"/Results_survival predictions/Combination data testing_1_all metrics.pdf"))


ggexport(ggarrange(ggplot(df_C.index.alldata %>% 
                            mutate(Pre_selection_N = ifelse(Pre_selection_N==349,"", (as.character(Pre_selection_N)))) , 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5, outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~Metric, scales = "free", space='free')+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=1, ncol=1),
         width = 10, height=6, filename = paste0(path,"/Results_survival predictions/Combination data testing_2_all metrics.pdf"))




ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()),
                   ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
         width = 4.5, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_3.pdf"))


ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS_AUC_log2","Treatment/DSS_AUC_log2","Genetics/DSS_AUC_log2","Diagnostics/DSS_AUC_log2",
                                                "Genetics/Diagnostics/DSS_AUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()),
                   ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS_AUC_log2","Treatment/DSS_AUC_log2","Genetics/DSS_AUC_log2","Diagnostics/DSS_AUC_log2",
                                                "Genetics/Diagnostics/DSS_AUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
         width = 6, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_4.pdf"))


# sts
Y1 <- Y
Y1$status[Y1$time>365*2] <- 0
Y1$time[Y1$time>365*2] <- 365*2

list_res_combination <- list()
list_res_clinical <- list()
list_res_drug <- list()
for(clin in c(1,1.1,1.2,2,2.1)){
  if(clin == 1){
    X2 <- df_clinical1
  }else if(clin == 1.1){
    X2 <- df_clinical1.1
  }else if(clin == 1.2){
    X2 <- df_clinical1.2
  }else if(clin == 2){
    X2 <- df_clinical2
  }else if(clin == 2.1){
    X2 <- df_clinical2.1
  }
  
  list_res_clinical[[paste0("clinical_",clin)]] <- Cox_forecasting_glmnet(X_data=X2,
                                                                          y_data=data.matrix(Y1),
                                                                          alpha=c(0,0.4,1),
                                                                          lambda=c(exp(seq(-8,6, 0.1))),
                                                                          free_cores = 2,
                                                                          test.n= c(6,4),
                                                                          nfolds = nrow(Y),
                                                                          iter=200,
                                                                          log_AUC=2,
                                                                          Patient.Z=2,
                                                                          Drug.Z =2,
                                                                          RCPC=0)
  for(m in colnames(df_scores)[c(6:12)]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
    X <- X[match(rownames(Y1), rownames(X)),]
    
    if(grepl("EC50",m)){
      X <- log10(X)
      X[,] <- apply(X,2, function(x) x - mean(x))
    }
    
    SDs <- matrixStats::colSds(data.matrix(X_rAUC))
    X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    
    for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
      X1 <- X[,SDs>=tr]
      list_res_combination[[paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=cbind(X2,X1),
                                                                                                                           y_data=data.matrix(Y1),
                                                                                                                           alpha=c(0,0.4,1),
                                                                                                                           lambda=c(exp(seq(-8,6, 0.1))),
                                                                                                                           free_cores = 2,
                                                                                                                           test.n= c(6,4),
                                                                                                                           nfolds = nrow(Y),
                                                                                                                           iter=200,
                                                                                                                           log_AUC=2,
                                                                                                                           Patient.Z=2,
                                                                                                                           Drug.Z =2,
                                                                                                                           RCPC=0)
      
      cat("\n",paste0("combination_",clin,"_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
    }
  }
}

for(m in colnames(df_scores)[c(6:12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  
  if(grepl("EC50",m)){
    X <- log10(X)
    X[,] <- apply(X,2, function(x) x - mean(x))
  }
  
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40,100)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_res_drug[[paste(m,"cutoff",ncol(X[,SDs>tr]))]] <- Cox_forecasting_glmnet(X_data=X1,
                                                                                  y_data=data.matrix(Y1),
                                                                                  alpha=c(0,0.4,1),
                                                                                  lambda=c(exp(seq(-8,6, 0.1))),
                                                                                  free_cores = 2,
                                                                                  test.n= c(6,4),
                                                                                  nfolds = nrow(Y),
                                                                                  iter=200,
                                                                                  log_AUC=2,
                                                                                  Patient.Z=2,
                                                                                  Drug.Z =2,
                                                                                  RCPC=0)
    
    
    cat("\n",paste0("_drug_",m,"_cutoff_",ncol(X[,SDs>tr])), "completed...","\n2")
  }
}

save(list_res_combination, list_res_clinical, list_res_drug,
     file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction_sts_all metrics.RData"))


load( file=paste0(path,"/Results_survival predictions/Drug screen-clinical data model comparisons_survival prediction_sts_all metrics.RData"))

df_C.index.alldata <- c(list_res_clinical,list_res_drug,list_res_combination)
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata$Penalty <- gsub(".*Penalty_","", df_C.index.alldata$ID)
df_C.index.alldata$Pre_selection_N <- gsub(".*cutoff_|.*cutoff ","",df_C.index.alldata$Metric)
df_C.index.alldata$Pre_selection_N[grepl("clinical",df_C.index.alldata$Pre_selection_N)] <- ""
df_C.index.alldata$Data <- gsub(" cutoff.*|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric <- gsub(" .*|.*drug_|_cutoff.*","",df_C.index.alldata$Metric)
df_C.index.alldata$Metric[grepl("clinical",df_C.index.alldata$Metric)] <- ""
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.alldata$Model <- ifelse(df_C.index.alldata$Penalty==1, "Lasso", "Elastic net")
df_C.index.alldata$Model[df_C.index.alldata$Penalty==0] <- "Ridge"
df_C.index.alldata$Model <- factor(df_C.index.alldata$Model, levels=unique(df_C.index.alldata$Model))
df_C.index.alldata$Pre_selection_N <- factor(df_C.index.alldata$Pre_selection_N, levels=unique(df_C.index.alldata$Pre_selection_N)[c(1,6:2)])

df_C.index.alldata$Data <- gsub("clinical_1.2", "Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1.1", "Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_1", "Genetics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2.1", "Genetics/Diagnostics+Clinical",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("clinical_2", "Genetics/Diagnostics",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.2_drug_", "Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1.1_drug_", "Diagnostics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_1_drug_", "Genetics/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2.1_drug_", "Genetics/Diagnostics+Clinical/",df_C.index.alldata$Data)
df_C.index.alldata$Data <- gsub("combination_2_drug_", "Genetics/Diagnostics/",df_C.index.alldata$Data)


comb <- expand.grid(clin=c("","Genetics",
                           "Diagnostics","Diagnostics+Clinical",
                           "Genetics/Diagnostics","Genetics/Diagnostics+Clinical"),
                    metric=unique(df_C.index.alldata$Metric))
comb$data <- paste(comb$clin, comb$metric, sep="/")
comb$data <- gsub("^\\/|.*\\/$", "", comb$data)

df_C.index.alldata$Data <- factor(df_C.index.alldata$Data, 
                                  levels=c("Genetics","Diagnostics","Diagnostics+Clinical",
                                           "Genetics/Diagnostics",
                                           "Genetics/Diagnostics+Clinical",
                                           comb$data[comb$data != ""] ))

df_C.index.alldata <- df_C.index.alldata[!grepl("Clinical", df_C.index.alldata$Data),]
df_C.index.stat <- df_C.index.alldata %>% group_by(Model, Metric, Data, Pre_selection_N) %>% summarise(Median=median(C_index_test),
                                                                                                       Mean=mean(C_index_test))


ggexport(ggarrange(ggplot(df_C.index.alldata %>% 
                            mutate(Pre_selection_N = ifelse(Pre_selection_N==349,"", (as.character(Pre_selection_N)))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5, outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~Metric, scales = "free", space='free')+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=1, ncol=1),
         width = 10, height=6, filename = paste0(path,"/Results_survival predictions/Combination data testing_1_sts_all metrics.pdf"))


ggexport(ggarrange(ggplot(df_C.index.alldata %>% 
                            mutate(Pre_selection_N = ifelse(Pre_selection_N==349,"", (as.character(Pre_selection_N)))) , 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5, outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~Metric, scales = "free", space='free')+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=1, ncol=1),
         width = 10, height=6, filename = paste0(path,"/Results_survival predictions/Combination data testing_2_sts_all metrics.pdf"))



ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()),
                   ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
         width = 4.5, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_3_sts.pdf"))



ggexport(ggarrange(ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS_AUC_log2","Treatment/DSS_AUC_log2","Genetics/DSS_AUC_log2","Diagnostics/DSS_AUC_log2",
                                                "Genetics/Diagnostics/DSS_AUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_test, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()),
                   ggplot(df_C.index.alldata %>% mutate(Pre_selection_N = gsub("349","",Pre_selection_N)) %>% 
                            subset(!(Pre_selection_N %in% c(100,80,60))) %>%
                            subset((Data %in% c("Treatment","Genetics","Diagnostics",
                                                "Genetics/Diagnostics",
                                                "Genetics/Diagnostics/Treatment",
                                                "rAUC_log2","Treatment/rAUC_log2","Genetics/rAUC_log2","Diagnostics/rAUC_log2",
                                                "Genetics/Diagnostics/rAUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS_AUC_log2","Treatment/DSS_AUC_log2","Genetics/DSS_AUC_log2","Diagnostics/DSS_AUC_log2",
                                                "Genetics/Diagnostics/DSS_AUC_log2",
                                                "Treatment/Genetics/Diagnostics/rAUC_log2",
                                                "DSS3","Treatment/DSS3","Genetics/DSS3","Diagnostics/DSS3",
                                                "Genetics/Diagnostics/DSS3",
                                                "Treatment/Genetics/Diagnostics/DSS3"))), 
                          aes(x=Data, y=C_index_train, fill=Model))+
                     geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
                     geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
                     geom_boxplot(alpha=0.5,  outlier.alpha = 0)+
                     facet_grid(Pre_selection_N~.)+
                     scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
                     ylim(0,1)+
                     theme_bw()+
                     labs(x="", y="C-index")+
                     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           strip.background = element_blank()), common.legend = TRUE, nrow=2, ncol=1),
         width = 6, height=9, filename = paste0(path,"/Results_survival predictions/Combination data testing_4_sts.pdf"))


######

# Variable selection and outcome associations ----
path <- "C:/Users/EnserinkLab2019/Desktop/Clinical forecasting - revision 2023-07"
load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(paste0(path, "/Survival, clinical features and ELN2022 classifications_2023-07.RData"))
df_scores <- df_scores %>% subset(!grepl("re", Patient.ID))
X_rAUC <- dcast(data = df_scores, Patient.num ~ drug , value.var="rAUC")
X <- X_rAUC[match(rownames(Y),X_rAUC$Patient.num),]; X $Patient.num <- NULL
rownames(X_rAUC) <-X_rAUC$Patient.num;X_rAUC$Patient.num<-NULL
X_rAUC <- X_rAUC[match(rownames(Y), rownames(X_rAUC)),]
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$rAUC <- 1 - df_scores$rAUC


rownames(df_survival) <- df_survival$ID
df_survival <- df_survival[match(rownames(Y), rownames(df_survival)),]
df_outcomes <- df_outcomes[match(rownames(Y), rownames(df_outcomes)),]


df_Y <- cbind(df_survival[,c("Persist._leukemia_post_first_ind._treatment", "Relapse")],
              df_outcomes[,c("Status", "one_year_survival","two_year_survival","three_year_survival","four_year_survival","five_year_survival")])

df_Y$Persist._leukemia_post_first_ind._treatment[df_Y$Persist._leukemia_post_first_ind._treatment=="Yes"] <- 1
df_Y$Persist._leukemia_post_first_ind._treatment[df_Y$Persist._leukemia_post_first_ind._treatment=="No"] <- 0
df_Y$Persist._leukemia_post_first_ind._treatment <- as.numeric(df_Y$Persist._leukemia_post_first_ind._treatment)
df_Y$Relapse[df_Y$Relapse=="Yes"] <- 1
df_Y$Relapse[df_Y$Relapse=="No"] <- 0
df_Y$Relapse <- as.numeric(df_Y$Relapse)

# sts
Y1 <- Y
Y1$status[Y1$time>365*2] <- 0
Y1$time[Y1$time>365*2] <- 365*2


list_coefficients <- list()
for(d in c("Survival (CoxPH) - long","Survival (CoxPH) - short",colnames(df_Y)[c(1:2,4:8)])){
  if(grepl("CoxPH",d)){
    if(grepl("long",d)){
      y <- Y
    }else{
      y <- Y1
    }
  }else{
    y <- df_Y[,d]
    remove <- !is.na(y)
    y <- y[!is.na(y)]
  }
  

  list_res_drug1 <- list()
  for(a in c(0,0.4,1)){
    
    for(m in colnames(df_scores)[c(7,12)]){
      X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){#
        X1 <- X[,SDs>=tr]
        if(grepl("CoxPH",d)){
          model.loop <- glmnet(data.matrix(X1), data.matrix(y), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), data.matrix(y),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }else{
          X1 <- X1[remove,]
          model.loop <- glmnet(data.matrix(X1), y, family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), y,family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }
      }
    }
    
  }
  list_coefficients[[d]] <- list(list_res_drug=list_res_drug1)
  cat("\n",d, "completed...","\n")
}

save(list_coefficients,
     file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients.RData"))


library("pheatmap")
library("grid")

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
load(file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients.RData"))

df_coefficients <- list_coefficients
for(i in names(df_coefficients)){
  for(j in names(df_coefficients[[i]])){
    for(k in names(df_coefficients[[i]][[j]])){
      df_coefficients[[i]][[j]][[k]] <- data.frame(data.matrix(df_coefficients[[i]][[j]][[k]]), Covariate=rownames(df_coefficients[[i]][[j]][[k]]))
    }
    df_coefficients[[i]][[j]] <- bind_rows(df_coefficients[[i]][[j]], .id="Model")
  }
  df_coefficients[[i]] <- bind_rows(df_coefficients[[i]], .id="Data")
}
df_coefficients <- bind_rows(df_coefficients, .id="Outcome")
df_coefficients$Coef <- df_coefficients$X1
df_coefficients$Coef[is.na(df_coefficients$Coef)] <- df_coefficients$s1[is.na(df_coefficients$Coef)] 

# Outcome association heatmaps ----
df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log|DSS3", df_coefficients$Model) &  !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 2.5,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))
myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
heat

pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_1.pdf"),
    width = 5, height = 7.5)
add.flag(heat,
         kept.labels = rownames(df)[matrixStats::rowMaxs(abs(data.matrix(df[,-1])))>3],
         repel.degree = 1)
dev.off()
# ----
df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log|DSS3", df_coefficients$Model) &  !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 2,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))

myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
heat

pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_2.pdf"),
    width = 5, height = 9.5)
add.flag(heat,
         kept.labels = rownames(df)[matrixStats::rowMaxs(abs(data.matrix(df[,-1])))>2.5],
         repel.degree = 1)
dev.off()

# ----

df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log|DSS3", df_coefficients$Model) &  !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 3,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))
myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
heat

pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_3.pdf"),
    width = 5, height = 7)
heat
dev.off()


# Variable selection and outcome associations - all metrics ----

list_coefficients <- list()
for(d in c("Survival (CoxPH) - long","Survival (CoxPH) - short",colnames(df_Y)[c(1:2,4:8)])){
  if(grepl("CoxPH",d)){
    if(grepl("long",d)){
      y <- Y
    }else{
      y <- Y1
    }
  }else{
    y <- df_Y[,d]
    remove <- !is.na(y)
    y <- y[!is.na(y)]
  }
  
  
  list_res_drug1 <- list()
  for(a in c(0)){#0.4,1
    
    for(m in colnames(df_scores)[c(6:12)]){
      X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      if(grepl("EC50",m)){
        X <- log10(X)
        X[,] <- apply(X,2, function(x) x - mean(x))
      }
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(0)){#quantile(SDs, 1-c(40)/ncol(X_rAUC)),
        X1 <- X[,SDs>=tr]
        if(grepl("CoxPH",d)){
          model.loop <- glmnet(data.matrix(X1), data.matrix(y), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), data.matrix(y),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }else{
          X1 <- X1[remove,]
          model.loop <- glmnet(data.matrix(X1), y, family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), y,family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }
      }
    }
    
  }
  list_coefficients[[d]] <- list(list_res_drug=list_res_drug1)
  cat("\n",d, "completed...","\n")
}

save(list_coefficients,
     file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients_all metrics.RData"))

load(file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))
df.hillAUC <- df.hillAUC %>% subset(!grepl("re", Patient.ID))

df.hillAUC$rAUC <- 1-df.hillAUC$rAUC
df.hillAUC$hill_AUC <- 1-df.hillAUC$hill_AUC
df.hillAUC$hill0_AUC <- 1-df.hillAUC$hill0_AUC
df.hillAUC$hillB_AUC <- 1-df.hillAUC$hillB_AUC


list_coefficients1 <- list()
for(d in c("Survival (CoxPH) - long","Survival (CoxPH) - short",colnames(df_Y)[c(1:2,4:8)])){
  if(grepl("CoxPH",d)){
    if(grepl("long",d)){
      y <- Y
    }else{
      y <- Y1
    }
  }else{
    y <- df_Y[,d]
    remove <- !is.na(y)
    y <- y[!is.na(y)]
  }
  
  
  list_res_drug1 <- list()
  for(a in c(0)){#0.4,1
    
    for(m in colnames(df.hillAUC)[c(6:8,10:12)]){
      X <- dcast(data = df.hillAUC, Patient.num ~ drug , value.var=m)
      rownames(X) <- as.character(X$Patient.num); X$Patient.num<-NULL
      X <- X[match(rownames(Y), rownames(X)),]
      if(grepl("EC50",m)){
        X <- log10(X)
        X[,] <- apply(X,2, function(x) x - mean(x))
      }
      SDs <- matrixStats::colSds(data.matrix(X_rAUC))
      X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
      
      for(tr in c(0)){#quantile(SDs, 1-c(40)/ncol(X_rAUC)),
        X1 <- X[,SDs>=tr]
        if(grepl("CoxPH",d)){
          model.loop <- glmnet(data.matrix(X1), data.matrix(y), family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), data.matrix(y),family = "cox", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }else{
          X1 <- X1[remove,]
          model.loop <- glmnet(data.matrix(X1), y, family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance")
          nfolds=nrow(X1)
          model.loop.cv <- cv.glmnet(data.matrix(X1), y,family = "binomial", alpha = a, standardize = FALSE,lambda=exp(seq(-8,6, 0.1)),  type.measure = "deviance", nfolds=nfolds)
          list_res_drug1[[paste0("drug_",m,"_alpha", a,"_cutoff",ncol(X[,SDs>tr]))]] <- coef.glmnet(model.loop, s = model.loop.cv$lambda.min)
        }
      }
    }
    
  }
  list_coefficients1[[d]] <- list(list_res_drug=list_res_drug1)
  cat("\n",d, "completed...","\n")
}

save(list_coefficients,list_coefficients1,
     file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients_all metrics.RData"))


load(file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients_all metrics.RData"))

df_coefficients <- c(list_coefficients)
for(i in names(df_coefficients)){
  for(j in names(df_coefficients[[i]])){
    for(k in names(df_coefficients[[i]][[j]])){
      df_coefficients[[i]][[j]][[k]] <- data.frame(data.matrix(df_coefficients[[i]][[j]][[k]]), Covariate=rownames(df_coefficients[[i]][[j]][[k]]))
    }
    df_coefficients[[i]][[j]] <- bind_rows(df_coefficients[[i]][[j]], .id="Model")
  }
  df_coefficients[[i]] <- bind_rows(df_coefficients[[i]], .id="Data")
}
df_coefficients <- bind_rows(df_coefficients, .id="Outcome")
df_coefficients$Coef <- df_coefficients$X1
df_coefficients$Coef[is.na(df_coefficients$Coef)] <- df_coefficients$s1[is.na(df_coefficients$Coef)] 

df_coefficients1 <- c(list_coefficients1)
for(i in names(df_coefficients1)){
  for(j in names(df_coefficients1[[i]])){
    for(k in names(df_coefficients1[[i]][[j]])){
      df_coefficients1[[i]][[j]][[k]] <- data.frame(data.matrix(df_coefficients1[[i]][[j]][[k]]), Covariate=rownames(df_coefficients1[[i]][[j]][[k]]))
    }
    df_coefficients1[[i]][[j]] <- bind_rows(df_coefficients1[[i]][[j]], .id="Model")
  }
  df_coefficients1[[i]] <- bind_rows(df_coefficients1[[i]], .id="Data")
}
df_coefficients1 <- bind_rows(df_coefficients1, .id="Outcome")
df_coefficients1$Coef <- df_coefficients1$X1
df_coefficients1$Coef[is.na(df_coefficients1$Coef)] <- df_coefficients1$s1[is.na(df_coefficients1$Coef)] 
#df_coefficients <- rbind(df_coefficients, df_coefficients1)
# Outcome association heatmaps ----
df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & !grepl("EC50", df_coefficients$Model) & !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 3.5,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))
myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")

# annotations <- data.frame(Metric = factor(gsub(" .*","",colnames(df[,-1])), levels = colnames(df_scores)[6:12]))
# rownames(annotations) <- colnames(df[,-1])
# anno_colors <- list(Metric = pal_jco()(7))
# names(anno_colors$Metric) <-  factor(colnames(df_scores)[6:12])

heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 #annotation_col = annotations, annotation_colors = anno_colors,
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)


pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_all metrics_1.pdf"),
    width = 15, height = 8)
heat
dev.off()
# ----

df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & !grepl("EC50", df_coefficients$Model) & !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 2,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))
myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
heat

pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_all metrics_2.pdf"),
    width = 15, height = 8)
add.flag(heat,
         kept.labels = rownames(df)[matrixStats::rowMaxs(abs(data.matrix(df[,-1])))>3.5],
         repel.degree = 1)
dev.off()

# ----

df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & !grepl("EC50", df_coefficients$Model) & !grepl("four",df_coefficients$Outcome) &
                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
            Covariate~Model+Outcome, value.var="Coef")
df <- df[!grepl("Intercept", df$Covariate),]
df[,-1] <- apply(df[,-1], 2, function(x) (x)/sd(x, na.rm = T))
df[,-1] <- df[,-1][,colSums(is.na(df[,-1])) != nrow(df)]
df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
df <- df[matrixStats::rowMaxs(abs(data.matrix(df[,-1])), na.rm = T) > 3,]
rownames(df) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df$Covariate) 
colnames(df)[-1] <- paste(gsub("drug_","",gsub("_alpha.*","",colnames(df)[-1])),
                          gsub(".*cutoff349_","",colnames(df)[-1]),
                          sep=" -> ")
colnames(df) <- gsub("DSS_AUC","loess_AUC",colnames(df))
colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
colnames(df) <- gsub("_","-",colnames(df))
myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
myBreaks <- c(seq(min(df[,-1]), 0, length.out=ceiling(100/2) + 1), 
              seq(max(df[,-1])/100, max(df[,-1]), length.out=floor(100/2)))

HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
heat <- pheatmap(df[,-1],
                 cluster_rows=HC_r,
                 cluster_cols=HC_c,
                 color = myColor,#bluered(100),
                 breaks=myBreaks,
                 #cutree_cols = 6,
                 #cutree_rows = 8,
                 
                 border_color="black",
                 fontsize_row=8,
                 fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
heat

pdf(file=paste0(path,"/Results_survival predictions/Clinical outcome associations_heatmap_all metrics_3.pdf"),
    width = 15, height = 8)
heat
dev.off()



# Variable enrichment ----
load(file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients.RData"))

df_coefficients <- list_coefficients
for(i in names(df_coefficients)){
  for(j in names(df_coefficients[[i]])){
    for(k in names(df_coefficients[[i]][[j]])){
      df_coefficients[[i]][[j]][[k]] <- data.frame(data.matrix(df_coefficients[[i]][[j]][[k]]), Covariate=rownames(df_coefficients[[i]][[j]][[k]]))
    }
    df_coefficients[[i]][[j]] <- bind_rows(df_coefficients[[i]][[j]], .id="Model")
  }
  df_coefficients[[i]] <- bind_rows(df_coefficients[[i]], .id="Data")
}
df_coefficients <- bind_rows(df_coefficients, .id="Outcome")


library(readxl)
library(viridis)
df_drugclasses <- read_excel(paste0(path,"/Selleck_Drug_classes_clean.xlsx") )

Class2Drug <- data.frame(Class=df_drugclasses$Group, Drug=df_drugclasses$`Item Name`)

Target2Drug <- data.frame(Class=df_drugclasses$Target, Drug=df_drugclasses$`Item Name`)
Target2Drug <- splitstackshape::cSplit(Target2Drug, "Class", "/", "long")


list_enrichment <- list()
df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Lasso"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))

df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Elastic net"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))


df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["DSS3 Lasso"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))



df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))
list_enrichment[["DSS3 Elastic net"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))


df_class_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_class_enrichment$Rank <- stringr::str_pad(df_class_enrichment$Rank*10, 3, pad = "0")

ggexport(ggplot(df_class_enrichment)+
           geom_bar(aes(x=paste(Rank, ID),y=Enrichment, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
           geom_text(aes(x=paste(Rank, ID),y=Enrichment, label=ifelse(pvalue<0.05,"*","")), size=7)+
           facet_wrap(~Data, scales="free_y")+
           labs(x="")+
           scale_fill_viridis(option="magma")+
           theme_classic()+
           theme(strip.background = element_blank())+
           coord_flip(),
         width=8,height=4.75,
         filename = paste0(path,"/Results_survival predictions/Enrichment test_variable selection_drug class.pdf"))




list_enrichment <- list()
df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Lasso"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))

df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Elastic net"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))


df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["DSS3 Lasso"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))



df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))
list_enrichment[["DSS3 Elastic net"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))


df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")




ggexport(ggplot(df_target_enrichment)+
           geom_bar(aes(x=paste(Rank, ID),y=Enrichment, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
           geom_text(aes(x=paste(Rank, ID),y=Enrichment, label=ifelse(pvalue<0.05,"*","")), size=7)+
           facet_wrap(~Data, scales="free_y")+
           labs(x="")+
           scale_fill_viridis(option="magma")+
           theme_classic()+
           theme(strip.background = element_blank())+
           coord_flip(),
         width=8,height=4.75,
         filename = paste0(path,"/Results_survival predictions/Enrichment test_variable selection_drug target.pdf"))


Class2Drug <- data.frame(Class=df_drugclasses$Group, Drug=df_drugclasses$`Item Name`)



list_enrichment <- list()
df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Lasso"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))

df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Elastic net"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))


df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["DSS3 Lasso"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))



df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Class2Drug)
df_class_enrichment <- Class_enrichment@result
df_class_enrichment$Enrichment <- log2((df_class_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_class_enrichment$BgRatio))/nrow(df_drugclasses)))
list_enrichment[["DSS3 Elastic net"]] <- df_class_enrichment %>% mutate(Rank=rank(-log10(pvalue)))


df_class_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_class_enrichment$Rank <- stringr::str_pad(df_class_enrichment$Rank*10, 3, pad = "0")

ggexport(ggplot(df_class_enrichment)+
           geom_bar(aes(x=paste(Rank, ID),y=Enrichment, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
           geom_text(aes(x=paste(Rank, ID),y=Enrichment, label=ifelse(pvalue<0.05,"*","")), size=7)+
           facet_wrap(~Data, scales="free_y")+
           labs(x="")+
           scale_fill_viridis(option="magma")+
           theme_classic()+
           theme(strip.background = element_blank())+
           coord_flip(),
         width=8,height=4.75,
         filename = paste0(path,"/Results_survival predictions/Enrichment test_variable selection_drug class_sts.pdf"))



Target2Drug <- data.frame(Class=df_drugclasses$Target, Drug=df_drugclasses$`Item Name`)
Target2Drug <- splitstackshape::cSplit(Target2Drug, "Class", "/", "long")

list_enrichment <- list()
df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Lasso"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))

df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["rAUC Elastic net"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))


df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha1_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))

list_enrichment[["DSS3 Lasso"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))



df_drugclasses_reduced <- df_drugclasses[df_drugclasses$`Item Name` %in% 
                                           df_coefficients$Covariate[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                                                       grepl("cutoff349", df_coefficients$Model) & grepl("alpha0.4_", df_coefficients$Model) & 
                                                                       grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0],]

set.seed(1)
Class_enrichment <- clusterProfiler::enricher(df_drugclasses_reduced$`Item Name`,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              universe=df_drugclasses$`Item Name`,
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              qvalueCutoff = 1,
                                              TERM2GENE=Target2Drug)
df_target_enrichment <- Class_enrichment@result
df_target_enrichment$Enrichment <- log2((df_target_enrichment$Count/nrow(df_drugclasses_reduced))/(as.numeric(gsub("\\/.*","",df_target_enrichment$BgRatio))/nrow(df_drugclasses)))
list_enrichment[["DSS3 Elastic net"]] <- df_target_enrichment %>% mutate(Rank=rank(-log10(pvalue), ties.method = "random"))


df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")




ggexport(ggplot(df_target_enrichment)+
           geom_bar(aes(x=paste(Rank, ID),y=Enrichment, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
           geom_text(aes(x=paste(Rank, ID),y=Enrichment, label=ifelse(pvalue<0.05,"*","")), size=7)+
           facet_wrap(~Data, scales="free_y")+
           labs(x="")+
           scale_fill_viridis(option="magma")+
           theme_classic()+
           theme(strip.background = element_blank())+
           coord_flip(),
         width=8,height=4.75,
         filename = paste0(path,"/Results_survival predictions/Enrichment test_variable selection_drug target_sts.pdf"))






# Variable enrichment GSEA ----
Class2Drug1 <- Class2Drug
Class2Drug1$Drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", Class2Drug1$Drug)

list_enrichment <- list()
df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Class2Drug1)
df_class_enrichment <- Class_enrichment@result


clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 7)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "rAUC_log2 Ridge")
ggsave(width=10, height=9,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_rAUC_drug class.pdf"))

list_enrichment[["rAUC Ridge"]] <- df_class_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Class2Drug1)
df_class_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 7)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "DSS3 Ridge")
ggsave(width=10, height=9,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_DSS3_drug class.pdf"))

list_enrichment[["DSS3 Ridge"]] <- df_class_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_class_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_class_enrichment$Rank <- stringr::str_pad(df_class_enrichment$Rank, 2, pad = "0")

ggplot(df_class_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y")+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()
ggsave(width=9, height=3.25,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug class.pdf"))



Target2Drug1 <- Target2Drug
Target2Drug1$Drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", Target2Drug1$Drug)

list_enrichment <- list()
df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Target2Drug1)
df_target_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 10)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "rAUC_log2 Ridge")
ggsave(width=8, height=6.75,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_rAUC_drug target.pdf"))

list_enrichment[["rAUC Ridge"]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))

df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <-gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Target2Drug1)
df_target_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 10)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "DSS3 Ridge")
ggsave(width=8, height=6.75,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_DSS3_drug target.pdf"))

list_enrichment[["DSS3 Ridge"]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y")+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()
ggsave(width=10, height=5,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug target.pdf"))


Class2Drug1 <- Class2Drug
Class2Drug1$Drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", Class2Drug1$Drug)

list_enrichment <- list()
df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Class2Drug1)
df_class_enrichment <- Class_enrichment@result


clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 7)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "rAUC_log2 Ridge")
ggsave(width=10, height=9,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_rAUC_drug class_sts.pdf"))

list_enrichment[["rAUC Ridge"]] <- df_class_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Class2Drug1)
df_class_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 7)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "DSS3 Ridge")
ggsave(width=10, height=9,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_DSS3_drug class_sts.pdf"))

list_enrichment[["DSS3 Ridge"]] <- df_class_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_class_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_class_enrichment$Rank <- stringr::str_pad(df_class_enrichment$Rank, 2, pad = "0")

ggplot(df_class_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y")+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()
ggsave(width=9, height=3.25,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug class_sts.pdf"))



Target2Drug1 <- Target2Drug
Target2Drug1$Drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", Target2Drug1$Drug)

list_enrichment <- list()
df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("rAUC_log2", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Target2Drug1)
df_target_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 10)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "rAUC_log2 Ridge")
ggsave(width=8, height=6.75,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_rAUC_drug target_sts.pdf"))

list_enrichment[["rAUC Ridge"]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))

df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl("DSS3", df_coefficients$Model) & 
                                            grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                            grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
ranked_list <- df_drugclasses_reduced$X1
names(ranked_list) <-gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
ranked_list <- sort(ranked_list, decreasing = TRUE)
set.seed(1)
Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                          pvalueCutoff = 1,
                                          pAdjustMethod = "BH",
                                          minGSSize = 3,
                                          maxGSSize = 500,
                                          TERM2GENE=Target2Drug1)
df_target_enrichment <- Class_enrichment@result
clusterProfiler::cnetplot(Class_enrichment, categorySize="pvalue", foldChange=ranked_list, showCategory = 10)+
  #scale_color_gradientn(colours=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100)))+
  scale_colour_gradient2(low=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[1],
                         mid=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[50],
                         high=rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))[100],
                         midpoint = 0)+
  labs(color="Coefficient", title = "DSS3 Ridge")
ggsave(width=8, height=6.75,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_DSS3_drug target_sts.pdf"))

list_enrichment[["DSS3 Ridge"]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))


df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y")+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()
ggsave(width=10, height=5,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug target_sts.pdf"))


# Enrichment all metrics----
load(file=paste0(path,"/Results_survival predictions/Clinical outcome_model coefficients_all metrics.RData"))

df_coefficients <- list_coefficients
for(i in names(df_coefficients)){
  for(j in names(df_coefficients[[i]])){
    for(k in names(df_coefficients[[i]][[j]])){
      df_coefficients[[i]][[j]][[k]] <- data.frame(data.matrix(df_coefficients[[i]][[j]][[k]]), Covariate=rownames(df_coefficients[[i]][[j]][[k]]))
    }
    df_coefficients[[i]][[j]] <- bind_rows(df_coefficients[[i]][[j]], .id="Model")
  }
  df_coefficients[[i]] <- bind_rows(df_coefficients[[i]], .id="Data")
}
df_coefficients <- bind_rows(df_coefficients, .id="Outcome")


list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  
  df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl(m, df_coefficients$Model) & 
                                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                              grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
  ranked_list <- df_drugclasses_reduced$X1
  names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  set.seed(1)
  Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            minGSSize = 3,
                                            maxGSSize = 500,
                                            TERM2GENE=Target2Drug1)
  df_target_enrichment <- Class_enrichment@result
  
  list_enrichment[[m]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
  
}

df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")
df_target_enrichment$Data <- factor(df_target_enrichment$Data, levels=unique(df_target_enrichment$Data))

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y", ncol=4)+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()

ggsave(width=20, height=10,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug targets_all metrics.pdf"))

list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  
  df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl(m, df_coefficients$Model) & 
                                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                              grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
  ranked_list <- df_drugclasses_reduced$X1
  names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  set.seed(1)
  Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            minGSSize = 3,
                                            maxGSSize = 500,
                                            TERM2GENE=Target2Drug1)
  df_target_enrichment <- Class_enrichment@result
  
  list_enrichment[[m]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
  
}

df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")
df_target_enrichment$Data <- factor(df_target_enrichment$Data, levels=unique(df_target_enrichment$Data))

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y", ncol=4)+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()

ggsave(width=20, height=10,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug targets_all metrics_sts.pdf"))



list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  
  df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl(m, df_coefficients$Model) & 
                                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                              grepl("long", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
  ranked_list <- df_drugclasses_reduced$X1
  names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  set.seed(1)
  Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            minGSSize = 3,
                                            maxGSSize = 500,
                                            TERM2GENE=Class2Drug1)
  df_target_enrichment <- Class_enrichment@result
  
  list_enrichment[[m]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
  
}

df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")
df_target_enrichment$Data <- factor(df_target_enrichment$Data, levels=unique(df_target_enrichment$Data))

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y", ncol=4)+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()

ggsave(width=17, height=6,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug classes_all metrics.pdf"))

list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  
  df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & grepl(m, df_coefficients$Model) & 
                                              grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                              grepl("short", df_coefficients$Outcome) & df_coefficients$X1 != 0,]
  ranked_list <- df_drugclasses_reduced$X1
  names(ranked_list) <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  set.seed(1)
  Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            minGSSize = 3,
                                            maxGSSize = 500,
                                            TERM2GENE=Class2Drug1)
  df_target_enrichment <- Class_enrichment@result
  
  list_enrichment[[m]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
  
}

df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")
df_target_enrichment$Data <- factor(df_target_enrichment$Data, levels=unique(df_target_enrichment$Data))

ggplot(df_target_enrichment)+
  geom_hline(yintercept = c(-1,1), lty=2)+
  geom_bar(aes(x=paste(Rank, ID),y=NES, fill=-log10(pvalue)), stat = "identity", col="black", size=0.1)+
  geom_text(aes(x=paste(Rank, ID),y=NES, label=ifelse(pvalue<0.05,"*","")), size=7)+
  facet_wrap(~Data, scales="free_y", ncol=4)+
  labs(x="")+
  scale_fill_viridis(option="magma")+
  theme_classic()+
  theme(strip.background = element_blank())+
  coord_flip()

ggsave(width=17, height=6,filename =  paste0(path,"/Results_survival predictions/Enrichment test_GSEA_drug classes_all metrics_sts.pdf"))



# Enrichment over all outcomes----
df_coefficients$beta <- df_coefficients$X1
df_coefficients$beta[is.na(df_coefficients$beta)] <- df_coefficients$s1[is.na(df_coefficients$beta)]

list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  for(j in unique(df_coefficients$Outcome)){
    df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & gsub("drug_|_alpha.*","",df_coefficients$Model) == m & 
                                                grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                                df_coefficients$Outcome == j,]
    ranked_list <- df_drugclasses_reduced$beta
    names(ranked_list) <- df_drugclasses_reduced$Covariate#gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
    ranked_list <- sort(ranked_list, decreasing = TRUE)
    set.seed(1)
    Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              TERM2GENE=Class2Drug)
    df_target_enrichment <- Class_enrichment@result
    
    list_enrichment[[paste(m,":",j)]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
    
  }
}
df_class_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_class_enrichment$Rank <- stringr::str_pad(df_class_enrichment$Rank, 2, pad = "0")
df_class_enrichment$Data <- factor(df_class_enrichment$Data, levels=unique(df_class_enrichment$Data))

df_class_enrichment$Metric <- gsub(" \\:.*","", df_class_enrichment$Data)
df_class_enrichment$Metric <- factor(df_class_enrichment$Metric, levels=unique(df_class_enrichment$Metric))
df_class_enrichment$Outcome <- gsub(".* \\:","", df_class_enrichment$Data)

list_enrichment <- list()
for(m in colnames(df_scores)[6:12]){
  for(j in unique(df_coefficients$Outcome)){
    df_drugclasses_reduced <- df_coefficients[!grepl("combination", df_coefficients$Model) & gsub("drug_|_alpha.*","",df_coefficients$Model) == m & 
                                                grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model) & 
                                                df_coefficients$Outcome == j,]
    ranked_list <- df_drugclasses_reduced$beta
    names(ranked_list) <- df_drugclasses_reduced$Covariate#gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_drugclasses_reduced$Covariate)
    ranked_list <- sort(ranked_list, decreasing = TRUE)
    set.seed(1)
    Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                              pvalueCutoff = 1,
                                              pAdjustMethod = "BH",
                                              minGSSize = 3,
                                              maxGSSize = 500,
                                              TERM2GENE=Target2Drug)
    df_target_enrichment <- Class_enrichment@result
    
    list_enrichment[[paste(m,":",j)]] <- df_target_enrichment %>% mutate(Rank=rank(NES, ties.method = "random"))
    
  }
}

library(pheatmap)
df_target_enrichment <- bind_rows(list_enrichment, .id = "Data")
df_target_enrichment$Rank <- stringr::str_pad(df_target_enrichment$Rank, 2, pad = "0")
df_target_enrichment$Data <- factor(df_target_enrichment$Data, levels=unique(df_target_enrichment$Data))

df_target_enrichment$Metric <- gsub(" \\:.*","", df_target_enrichment$Data)
df_target_enrichment$Metric <- factor(df_target_enrichment$Metric, levels=unique(df_target_enrichment$Metric))
df_target_enrichment$Outcome <- gsub(".* \\:","", df_target_enrichment$Data)


m <- colnames(df_scores)[6:12][2]
for(m in  colnames(df_scores)[6:12]){
  df <- dcast(df_class_enrichment[df_class_enrichment$Metric == m & !grepl("four", df_class_enrichment$Data),] %>% mutate(value = -log10(pvalue)*ifelse(NES<0,-1,1)),
              Description~Outcome, value.var = "value")
  df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
  #df <- df[matrixStats::rowMaxs(data.matrix(abs(df[,-1])))>-log10(0.1),]
  colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
  colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
  colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
  colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
  colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
  rownames(df) <- df$Description
  myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
  myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1), 
                seq(max(df[,-1])/100, 2, length.out=floor(100/2)))
  
  HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
  HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
  heat <- pheatmap(df[,-1],
                   cluster_rows=HC_r,
                   cluster_cols=HC_c,
                   color = myColor,
                   breaks=myBreaks,
                   #cutree_cols = 6,
                   #cutree_rows = 8,
                   
                   cellwidth = 8, cellheight = 8, filename = paste0(path,"/Results_survival predictions/Enrichment GSEA_all outcomes_drug classes_",m,".pdf"),
                   
                   border_color="black",
                   fontsize_row=8,
                   fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
  
  # pdf(file=paste0(path,"/Results_survival predictions/Enrichment GSEA_all outcomes_drug classes_",m,".pdf"),
  #     width = 4, height = 4.5)
  # heat
  # dev.off()
  
  df <- dcast(df_target_enrichment[df_target_enrichment$Metric == m & !grepl("four", df_target_enrichment$Data),] %>% mutate(value = -log10(pvalue)*ifelse(NES<0,-1,1)),
              Description~Outcome, value.var = "value")
  df[,grepl("year", colnames(df))] <- -1*df[,grepl("year", colnames(df))]
  df <- df[df$Description != "ATR",]; df$Description <- gsub("ATM", "ATM/ATR", df$Description)   #redundant
  df <- df[df$Description != "progestogen Receptor",] #redundant
  #df <- df[matrixStats::rowMaxs(data.matrix(abs(df[,-1])))>-log10(0.1),]
  colnames(df) <- gsub("Persist._leukemia_post_first_ind._treatment","Response to induction",colnames(df))
  colnames(df) <- gsub("one_year_survival","One year survival",colnames(df))
  colnames(df) <- gsub("two_year_survival","Two year survival",colnames(df))
  colnames(df) <- gsub("three_year_survival","Three year survival",colnames(df))
  colnames(df) <- gsub("five_year_survival","Five year survival",colnames(df))
  rownames(df) <- df$Description
  myColor = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"PuOr"))(100))
  myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1), 
                seq(max(df[,-1])/100, 2, length.out=floor(100/2)))
  
  HC_r <- hclust(dist((df[,-1]), method="euclidean"), method="ward.D2")
  HC_c <- hclust(dist(t(df[,-1]), method="euclidean"), method="ward.D2")
  heat <- pheatmap(df[,-1],
                   cluster_rows=HC_r,
                   cluster_cols=HC_c,
                   color = myColor,
                   breaks=myBreaks,
                   #cutree_cols = 6,
                   #cutree_rows = 8,
                   cellwidth = 8, cellheight = 8, filename = paste0(path,"/Results_survival predictions/Enrichment GSEA_all outcomes_drug targets_",m,".pdf"),
                   
                   border_color="black",
                   fontsize_row=8,
                   fontsize_col=8, legend = TRUE, annotation_legend = TRUE)
  
  # pdf(file=paste0(path,"/Results_survival predictions/Enrichment GSEA_all outcomes_drug targets_",m,".pdf"),
  #     width = 4, height = 10)
  # heat
  # dev.off()
}

# Relapse sensitivity analysis ----
library(ggrepel)


load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
df_scores$DSS3 <- df_scores$DSS3/100
df_scores$DSS2 <- df_scores$DSS2/100
df_scores$DSS1 <- df_scores$DSS1/100
df_scores$DSS_AUC <- df_scores$DSS_AUC/100
df_scores$rAUC <- 1 - df_scores$rAUC
df_scores_re <- df_scores %>% subset(Patient.num %in% df_scores$Patient.num[grepl("re", Patient.ID)])

m <- "rAUC_log2"
pt_st <- TRUE

for(m in colnames(df_scores)[6:12]){
  X_relapse <- dcast(data = df_scores_re, Patient.ID ~ drug , value.var=m)
  rownames(X_relapse) <- X_relapse$Patient.ID; X_relapse$Patient.ID<-NULL
  X_relapse <- X_relapse[-grep("44", rownames(X_relapse)),]
  if(pt_st){X_relapse <- (X_relapse- rowMeans(X_relapse))/matrixStats::rowSds(data.matrix(X_relapse))}
  
  
  df <- dcast(df_coefficients[!grepl("combination", df_coefficients$Model) & gsub("drug_|_alpha.*","",df_coefficients$Model) == m & !grepl("four",df_coefficients$Outcome) &
                                grepl("cutoff349", df_coefficients$Model) & grepl("alpha0_", df_coefficients$Model),],
              Covariate~Outcome, value.var="beta")
  df <- df[!grepl("Intercept", df$Covariate),]
  #rownames(df) <- gsub(" \\(.*","",df$Covariate) 
  
  
  # Relapse association
  df_relapse_t.test <- cbind(Mean=rowMeans(t(X_relapse[grepl("re", rownames(X_relapse)),]) - t(X_relapse[!grepl("re", rownames(X_relapse)),])),
                             SD=matrixStats::rowSds(t(X_relapse[grepl("re", rownames(X_relapse)),]) - t(X_relapse[!grepl("re", rownames(X_relapse)),])))
  df_relapse_t.test <- data.frame(df_relapse_t.test)
  df_relapse_t.test$SE <- df_relapse_t.test$SD/sqrt(8)
  df_relapse_t.test$t.stat <- df_relapse_t.test$Mean/df_relapse_t.test$SE
  df_relapse_t.test$pval <- pt(abs(df_relapse_t.test$t.stat), 8, lower.tail = FALSE)
  df_relapse_t.test$Drug <- gsub(" \\(.*","",rownames(df_relapse_t.test))
  
  
  df_relapse_t.test$Coef_survival <- df$`Survival (CoxPH) - long`[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$Coef_survival_sts <- df$`Survival (CoxPH) - short`[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$Coef_relapse <- df$Relapse[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$Coef_induction <- df$Persist._leukemia_post_first_ind._treatment[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$five_year_survival <- df$five_year_survival[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$three_year_survival <- df$three_year_survival[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$two_year_survival <- df$two_year_survival[match(rownames(df_relapse_t.test),  df$Covariate)]
  df_relapse_t.test$one_year_survival <- df$one_year_survival[match(rownames(df_relapse_t.test),  df$Covariate)]
  
  ggarrange(ggplot(df_relapse_t.test,aes(Coef_survival,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Survival association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(Coef_survival_sts,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Survival association (short)",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(Coef_relapse,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Relapse association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(Coef_induction,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Induction response association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(-one_year_survival,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="One year survival association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(-two_year_survival,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Two year survival association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(-three_year_survival,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Three year survival association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(-five_year_survival,Mean))+
              geom_point(aes(), alpha=0.75, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              labs(x="Five year survival association",y="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            common.legend = T)
  ggsave(filename = paste0(path, "/Results_relapse analysis/Relapse differential vs model coefficients_",m,".pdf"),
         width=9, height=7.75)
  
  
  
  ggarrange(ggplot(df_relapse_t.test,aes(Mean, -log10(pval)))+
              geom_point(aes(col=pval<0.05), size=1)+
              geom_text_repel(data=df_relapse_t.test[df_relapse_t.test$pval < 0.03 | 
                                                       df_relapse_t.test$pval < 0.1 & abs(df_relapse_t.test$Mean) > 0.5,],
                              aes(label=Drug),
                              nudge_y = 0.1, force=1, nudge_x = 0,
                              max.overlaps = Inf,min.segment.length = 0,
                              size=3)+
              labs(x="Differential drug sensitivity\nRelapse - Treatment naive", y="-log10(p-value)",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw()+
              theme(legend.position = "none"),
            ggplot(df_relapse_t.test,aes(y=Coef_survival,x=Mean))+
              geom_point(aes(), alpha=1, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              geom_text_repel(data=df_relapse_t.test[abs(df_relapse_t.test$Mean)>0.1 &
                                                       abs(df_relapse_t.test$Coef_survival)/sd(df_relapse_t.test$Coef_survival)>2,],
                              aes(label=Drug),
                              nudge_y = 0, force=1,min.segment.length = 0,
                              max.overlaps = Inf,
                              size=3)+
              labs(y="Survival association",x="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(y=Coef_relapse,x=Mean))+
              geom_point(aes(), alpha=1, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              geom_text_repel(data=df_relapse_t.test[abs(df_relapse_t.test$Mean)>0.1 &
                                                       abs(df_relapse_t.test$Coef_relapse)/sd(df_relapse_t.test$Coef_relapse)>2,],
                              aes(label=Drug),
                              nudge_y = 0, force=1,min.segment.length = 0,
                              max.overlaps = Inf,
                              size=3)+
              labs(y="Relapse association",x="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(),
            ggplot(df_relapse_t.test,aes(y=Coef_induction,x=Mean))+
              geom_point(aes(), alpha=1, size=1, col=pal_jco()(3)[3])+
              stat_cor()+
              stat_smooth(method = "lm", col = pal_jco()(1)) +
              geom_text_repel(data=df_relapse_t.test[abs(df_relapse_t.test$Mean)>0.1 &
                                                       abs(df_relapse_t.test$Coef_induction)/sd(df_relapse_t.test$Coef_induction)>2,],
                              aes(label=Drug),
                              nudge_y = 0, force=1,min.segment.length = 0,
                              max.overlaps = Inf,
                              size=3)+
              labs(y="Induction response association",x="Differential drug sensitivity\nRelapse - Treatment naive",
                   col="Significant")+
              scale_color_manual(values=pal_jco()(3)[c(3,1)])+
              theme_bw(), align = "v")
  ggsave(filename = paste0(path, "/Results_relapse analysis/Relapse differential vs model coefficients_",m,"_2.pdf"),
         width=7, height=6)
  
  
  Class2Drug <- data.frame(Class=df_drugclasses$Group, Drug=df_drugclasses$`Item Name`)
  #Class2Drug$Drug <- gsub(" \\(.*", "", Class2Drug$Drug)
  
  Target2Drug <- data.frame(Class=df_drugclasses$Target, Drug=df_drugclasses$`Item Name`)
  Target2Drug <- splitstackshape::cSplit(Target2Drug, "Class", "/", "long")
  #Target2Drug$Drug <- gsub(" \\(.*", "", Target2Drug$Drug)
  
  
  ranked_list <- df_relapse_t.test$Mean
  names(ranked_list) <- rownames(df_relapse_t.test)
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  set.seed(1)
  Class_enrichment <- clusterProfiler::GSEA(ranked_list,
                                            pvalueCutoff = 1,
                                            pAdjustMethod = "BH",
                                            minGSSize = 3,
                                            maxGSSize = 500,
                                            TERM2GENE=Class2Drug)
  df_class_enrichment <- Class_enrichment@result
  
  set.seed(1)
  Target_enrichment <- clusterProfiler::GSEA(ranked_list,
                                             pvalueCutoff = 1,
                                             pAdjustMethod = "BH",
                                             minGSSize = 3,
                                             maxGSSize = 500,
                                             TERM2GENE=Target2Drug)
  df_target_enrichment <- Target_enrichment@result
  
  
  df_target_enrichment <- df_target_enrichment[df_target_enrichment$Description != "ATR",]; df_target_enrichment$Description <- gsub("ATM", "ATM/ATR", df_target_enrichment$Description)   #redundant
  df_target_enrichment <- df_target_enrichment[df_target_enrichment$Description != "progestogen Receptor",] #redundant
  
  
  ggplot(rbind(df_target_enrichment %>% mutate(Set="Drug target"),
               df_class_enrichment %>% mutate(Set="Drug class")), aes(enrichmentScore, -log10(pvalue)))+
    geom_point(aes(fill=-log10(pvalue), size=stringr::str_count(core_enrichment, "/")+1), pch=21)+
    geom_text_repel(data=rbind(df_target_enrichment %>% mutate(Set="Drug target"),
                               df_class_enrichment %>% mutate(Set="Drug class")) %>% 
                      subset(pvalue<0.05),aes(label=ID),
                    max.overlaps = Inf,
                    size=3)+
    facet_wrap(~Set)+
    scale_fill_viridis(option="magma")+
    scale_size_continuous(guide = "none")+
    theme_bw()+
    theme(legend.position = "top",
          strip.background = element_blank())
  ggsave(filename = paste0(path, "/Results_relapse analysis/GSEA relapse differential_",m,".pdf"),
         width=4.75, height=3.25)
  
}







load(paste0(path,"/Drug sensitivity QC results_2022-11.RData"))
df_clinical <- read.csv(paste0(path,"/Clinical_data_221122.csv"), 
                        comment.char="#", stringsAsFactors=FALSE)
df_relapse_time <- df_Z.prime_avg[df_Z.prime_avg$Patient.ID %in% rownames(X_relapse),]
df_relapse_time$Sample <- ifelse(grepl("re", df_relapse_time$Patient.ID), "Relapse", "Treatment_naive")
df_relapse_time <- reshape2::dcast(df_relapse_time, Patient.num~Sample, value.var = "Date_time")
df_relapse_time$Sample_time <- difftime(df_relapse_time$Relapse, df_relapse_time$Treatment_naive, units = "days")
attr(df_relapse_time$Sample_time, "units") <- NULL
df_relapse_time$Sample_time <- round(df_relapse_time$Sample_time)

df_relapse_time$Date_diagnose <- df_clinical$Date_diagnose[match(df_relapse_time$Patient.num, df_clinical$PatientCode_JElist)]
df_relapse_time$Relapse_date <- df_clinical$Residiv_date[match(df_relapse_time$Patient.num, df_clinical$PatientCode_JElist)]

df_relapse_time$Relapse_time <- difftime(df_relapse_time$Relapse_date, df_relapse_time$Date_diagnose, units = "days")
attr(df_relapse_time$Relapse_time, "units") <- NULL
df_relapse_time$Relapse_time <- round(df_relapse_time$Relapse_time)


df_relapse_time <- cbind(df_relapse_time, Y[match(df_relapse_time$Patient.num, rownames(Y)),])



ggplot(df_relapse_time)+
  geom_segment(aes(rank(-time), time, 
                   xend=rank(-time), yend=0))+
  geom_point(aes(rank(-time), time, pch=factor(status), fill=factor(status)), size=2)+
  geom_point(aes(rank(-time), Relapse_time),pch=15, size=2, fill=pal_jco()(2)[2])+
  geom_point(aes(rank(-time), Sample_time),pch=16, size=2, fill=pal_jco()(2)[2])+
  scale_shape_manual(values=c(23,21), name="Status")+
  scale_fill_jama()+
  theme_classic()+
  coord_flip()


#----

