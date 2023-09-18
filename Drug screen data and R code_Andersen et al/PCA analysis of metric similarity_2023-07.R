library(ggsci)
library(viridis)
library(ggrepel)

path <- "C:/Users/EnserinkLab2019/Desktop/Clinical forecasting - revision 2023-07"
load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(paste0(path, "/Survival, clinical features and ELN2022 classifications_2023-07.RData"))

df_Y <- data.frame(Y, x="")
fit <- survfit(formula = Surv(time, status) ~ x , data = df_Y)

ggexport(ggsurvplot(fit, data = df_Y,
                    surv.median.line = "hv", # Add medians survival
                    pval = TRUE,
                    
                    conf.int = TRUE,
                    # Add risk table
                    risk.table = TRUE,
                    tables.height = 0.2,
                    tables.theme = theme_cleantable(),
                    palette = pal_jco()(2)[1],
                    ggtheme = theme_classic()),
         filename = paste0(path, "/Survival plot.pdf"), width=3, height=4)

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

load(file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))
df.hillAUC <- df.hillAUC %>% subset(!grepl("re", Patient.ID))

df.hillAUC$rAUC <- 1-df.hillAUC$rAUC
df.hillAUC$hill_AUC <- 1-df.hillAUC$hill_AUC
df.hillAUC$hill0_AUC <- 1-df.hillAUC$hill0_AUC
df.hillAUC$hillB_AUC <- 1-df.hillAUC$hillB_AUC


# Drug sensitivity metric similarity ----
df1 <- cbind(df_scores[,6:12], 
             df.hillAUC[match(paste(df_scores$Patient.ID, df_scores$drug), paste(df.hillAUC$Patient.ID, df.hillAUC$drug)),
                        c(6:8,10:12)])
df1[,] <- apply(df1, 2, function(x) (x)/sd(x, na.rm = T))
pca <- prcomp(t(df1[,]), scale. = F)
df_pca <- data.frame(pca$x)
df_pca$Metric <- gsub(" .*","",rownames(df_pca))
df_pca$Metric <- factor(df_pca$Metric, levels = c(colnames(df_scores)[6:12], colnames(df.hillAUC)[c(6:8,10:12)]))
df_pca$Transformed <- "" 
df_pca$Transformed[grepl("hill", df_pca$Metric)] <- "hill"
df_pca$Transformed[grepl("log", df_pca$Metric) & grepl("hill", df_pca$Metric)] <- "hill + log2"
pca_stat <- summary(pca)

ggplot(df_pca,aes(PC1,PC2))+
  geom_point(aes(fill=Metric, alpha=Transformed, pch=Transformed),size=5)+
  geom_text_repel(aes(label=Metric), size=3)+
  labs(x=paste0("PC1 (",pca_stat$importance[2,1]*100,"%)"),
       y=paste0("PC2 (",pca_stat$importance[2,2]*100,"%)"))+
  scale_shape_manual(values=c(21,24,23))+
  scale_alpha_manual(values=c(1,1,1))+
  scale_fill_manual(values = c(pal_jco()(7), pal_jama()(3), pal_jama()(3)))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(filename = paste0(path, "/Results_metric comparison/PCA_metric comparison.pdf"),
       width=3, height=2.75)

df2 <- df1

df1 <- cbind(df_scores[,1:12], 
             df.hillAUC[match(paste(df_scores$Patient.ID, df_scores$drug), paste(df.hillAUC$Patient.ID, df.hillAUC$drug)),
                        c(6:8,10:12)])
df1 <- df1 %>% group_by(Patient.ID) %>% mutate_at(c(colnames(df_scores)[6:12], colnames(df.hillAUC)[c(6:8,10:12)]), function(x) (x-mean(x))/sd(x))
df1[,-c(1:5)] <- apply(df1[,-c(1:5)], 2, function(x) (x)/sd(x, na.rm = T))
pca <- prcomp(t(df1[,-c(1:5)]), scale. = F)
df_pca <- data.frame(pca$x)
df_pca$Metric <- gsub(" .*","",rownames(df_pca))
df_pca$Metric <- factor(df_pca$Metric, levels = c(colnames(df_scores)[6:12], colnames(df.hillAUC)[c(6:8,10:12)]))
df_pca$Transformed <- "" 
df_pca$Transformed[grepl("hill", df_pca$Metric)] <- "hill"
df_pca$Transformed[grepl("log", df_pca$Metric) & grepl("hill", df_pca$Metric)] <- "hill + log2"
pca_stat <- summary(pca)

ggplot(df_pca,aes(PC1,PC2))+
  geom_point(aes(fill=Metric, alpha=Transformed, pch=Transformed),size=5)+
  geom_text_repel(aes(label=Metric), size=3)+
  labs(x=paste0("PC1 (",pca_stat$importance[2,1]*100,"%)"),
       y=paste0("PC2 (",pca_stat$importance[2,2]*100,"%)"))+
  scale_shape_manual(values=c(21,24,23))+
  scale_alpha_manual(values=c(1,1,1))+
  scale_fill_manual(values = c(pal_jco()(7), pal_jama()(3), pal_jama()(3)))+
  theme_classic()+
  theme(legend.position = "none")
ggsave(filename = paste0(path, "/Results_metric comparison/PCA_metric comparison_z-scores.pdf"),
       width=3, height=2.75)


ggplot(df.hillAUC, aes(rAUC, hill_AUC))+
  geom_hline(yintercept=0, lty=2)+
  #geom_point()+
  geom_hex(bins = 50) +
  #stat_cor()+
  theme_bw()+
  scale_fill_viridis(discrete=F,option="magma") +
  theme(panel.grid = element_blank())
ggsave(filename = paste0(path, "/Results_metric comparison/Hill AUC cutoff.pdf"),
       width=4, height=3)


path <- "C:/Users/EnserinkLab2019/Desktop/Clinical forecasting - revision 2023-07"
load(paste0(path, "/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(paste0(path, "/Survival, clinical features and ELN2022 classifications_2023-07.RData"))
load(paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))

df_pred_Hill$ID <- NULL
df_scores$ID <- NULL
df.hillAUC$ID <- NULL

write.csv(df_pred_Hill, file = paste0(path,"/Dose response data and Hill curve fits_2023-07.csv"), 
          row.names = FALSE)
write.csv(df_scores, file = paste0(path,"/Drug sensitivity metrics_2023-07.csv"), 
          row.names = FALSE)
write.csv(df.hillAUC, file = paste0(path,"/Drug sensitivity metrics_Hill rAUCs_2023-07.csv"), 
          row.names = FALSE)


