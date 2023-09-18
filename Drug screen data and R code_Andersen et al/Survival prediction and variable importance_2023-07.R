library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(drc)
library(data.table)
library(CoxTools)
library(doSNOW)
library(survival)
library(survminer)
library(ggrepel)
library(pheatmap)
library(superheat)



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
df_scores$DSS_AUC_log2 <- NULL

# Variable association statistics - bootstrapping ----
list_variable_association <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                          y_data=data.matrix(Y),
                                                                          alpha=0,
                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                          free_cores = 5,
                                                                          pre.CV = F,
                                                                          lambda_opt = 0,
                                                                          iter=200,
                                                                          log_AUC=2,
                                                                          Patient.Z=2,
                                                                          Drug.Z =2,
                                                                          RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}

list_variable_association_elnet <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association_elnet[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                                y_data=data.matrix(Y),
                                                                                alpha=0.4,
                                                                                lambda=exp(seq(-8,6, 0.1)),
                                                                                free_cores = 5,
                                                                                pre.CV = F,
                                                                                lambda_opt = 0,
                                                                                iter=200,
                                                                                log_AUC=2,
                                                                                Patient.Z=2,
                                                                                Drug.Z =2,
                                                                                RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}

list_variable_association_lasso <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association_lasso[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                                y_data=data.matrix(Y),
                                                                                alpha=1,
                                                                                lambda=exp(seq(-8,6, 0.1)),
                                                                                free_cores = 5,
                                                                                pre.CV = F,
                                                                                lambda_opt = 0,
                                                                                iter=200,
                                                                                log_AUC=2,
                                                                                Patient.Z=2,
                                                                                Drug.Z =2,
                                                                                RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}


save(list_variable_association, list_variable_association_elnet, list_variable_association_lasso,
     file = paste0(path, "/Results_survival predictions/Variable association results.RData"))


# Variable withdrawal testing ----
Cox_forecasting_drug_withdrawal <- function(X_data,
                                            y_data,
                                            Reduce=c("None", "Naive"),
                                            alpha=0,
                                            lambda=c(exp(seq(-4,6, 0.1))),
                                            free_cores = 2,
                                            test.n= c(5,5),
                                            iter=200,
                                            log_AUC=2,
                                            Patient.Z=2,
                                            Drug.Z =2,
                                            RCPC=0,
                                            path=path){
  # Data preparation
  a = alpha
  Transform = c("log2(AUC)", "AUC")[log_AUC[1]]
  pt.st = c(TRUE, FALSE)[Patient.Z[1]]
  cpc = RCPC[1]
  drug.st = c(TRUE, FALSE)[Drug.Z[1]]
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
  colnames(X) <- colnames(X_data)
  X_data <- X
  
  Drug_WD_test <- function(X_data,
                           y_data,
                           alpha=a,
                           lambda=c(exp(seq(-4,6, 0.1))),
                           free_cores = 2,
                           test.n= c(6,4),
                           iter=200,
                           log_AUC=c(1:2),
                           Patient.Z=c(1:2),
                           Drug.Z =c(1:2),
                           RCPC=c(0,1,2,3,4),
                           i=i){
    a=alpha
    # Prediction function
    # Pred_cox <- function(X,
    #                      Y,
    #                      test.n=test.n,
    #                      a=a,lambda=lambda, iter=iter, i=i){
    #   
    #   cores <- min(c(parallel::detectCores()-free_cores, iter))
    #   cluster.cores<-makeCluster(cores)
    #   registerDoSNOW(cluster.cores)
    #   pb <- txtProgressBar(max=iter, style=3)
    #   progress <- function(n) setTxtProgressBar(pb, n)
    #   opts <- list(progress=progress)
    #   list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
    #     set.seed(b)
    #     setTxtProgressBar(pb,b)
    #     X.0 <- X[Y[,2]==0,]
    #     X.1 <- X[Y[,2]==1,]
    #     Y.0 <- Y[Y[,2]==0,]
    #     Y.1 <- Y[Y[,2]==1,]
    #     
    #     ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
    #     ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)
    #     
    #     train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
    #     test <- rbind(X.0[ind.0,], X.1[ind.1,])
    #     
    #     Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
    #     Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])
    #     
    #     model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
    #     nfolds=nrow(train)
    #     model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
    #     
    #     c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
    #     ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))
    #     
    #     rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
    #     
    #     return(c(c,ct))
    #   }
    #   stopCluster(cluster.cores)
    #   closeAllConnections()
    #   rm(cluster.cores)
    #   rm(cores)
    #   gc()
    #   
    #   df_C_index <- do.call(rbind, list_C_index)
    #   df_C_index <- data.frame(df_C_index)
    #   colnames(df_C_index) <- c("C_index_test", "C_index_train")
    #   df_C_index$ID <- i
    #   
    #   rm(X,Y, list_C_index, opts, pb)
    #   
    #   
    #   return(df_C_index)
    # }
    # 
    
    
    # Run on full data
    # df_parent <- Pred_cox(X_data, y_data, test.n, a,lambda, iter, i)
    
    # Prepare drug withdrawal data
    mat.C.index.test <- data.frame(array(data=NA, dim=c(iter, ncol(X_data))))
    colnames(mat.C.index.test) <- colnames(X_data)
    mat.C.index.train <- mat.C.index.test
    
    cores <- min(c(parallel::detectCores()-free_cores, ncol(X_data)))
    cluster.cores<-makeCluster(cores)
    registerDoSNOW(cluster.cores)
    pb <- txtProgressBar(max=ncol(X_data), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_C_index <-foreach(d = 0:ncol(X_data), .packages = c("glmnet"), .options.snow=opts) %dopar%{
      
      
      if(d != 0){
        dr <- colnames(X_data)[d]
        X_data1 <- X_data[,colnames(X_data) != dr]
      }else{
        X_data1 <- X_data
        
      }
      
      Y <- data.matrix(y_data)
      
      list_C_index <- list()
      for(b in 1:iter){
        set.seed(b)
        X.0 <- X_data1[Y[,2]==0,]
        X.1 <- X_data1[Y[,2]==1,]
        Y.0 <- Y[Y[,2]==0,]
        Y.1 <- Y[Y[,2]==1,]
        
        ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
        ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)
        
        train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
        test <- rbind(X.0[ind.0,], X.1[ind.1,])
        
        Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
        Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])
        
        model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
        nfolds=nrow(train)
        model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
        
        c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
        ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))
        
        #rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
        list_C_index[[b]] <- c(c,ct)
      }
      
      df_C_index <- do.call(rbind, list_C_index)
      df_C_index <- data.frame(df_C_index)
      colnames(df_C_index) <- c("C_index_test", "C_index_train")
      df_C_index$ID <- i
      
      gc()
      
      df_wd <- df_C_index
      
      #rm(X_data1)
      return(df_wd)
    }
    stopCluster(cluster.cores)
    closeAllConnections()
    rm(cluster.cores)
    rm(cores)
    gc()
    
    cat("\n", "WD test completed...", "\n")
    
    names(list_C_index) <- c("Complete",colnames(X_data))
    for(dr in  names(list_C_index)){
      mat.C.index.test[,dr] <- list_C_index[[dr]]$C_index_test
      mat.C.index.train[,dr] <- list_C_index[[dr]]$C_index_train
    }
    
    C_loss_test <- colMeans((mat.C.index.test - list_C_index[["Complete"]]$C_index_test))
    C_loss_train <- colMeans((mat.C.index.train - list_C_index[["Complete"]]$C_index_train))
    
    df_loss <- data.frame(C_loss_test=C_loss_test,
                          C_loss_train=C_loss_train,
                          Var=names(C_loss_train))
    
    
    return(list(df_loss=df_loss,
                mat.C.index.test=mat.C.index.test,
                mat.C.index.train=mat.C.index.train))
  }
  
  list_res_initial <- Drug_WD_test(X_data,
                                   y_data,
                                   alpha=alpha,
                                   lambda=lambda,
                                   free_cores = free_cores,
                                   test.n= test.n,
                                   iter=iter,
                                   log_AUC=log_AUC,
                                   Patient.Z=Patient.Z,
                                   Drug.Z =Drug.Z,
                                   RCPC=RCPC,
                                   i=i)
  return(list_res_initial)
}


list_variable_withdrawal <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                    y_data=data.matrix(Y),
                                                                                    Reduce="Naive",
                                                                                    alpha=0,
                                                                                    lambda=exp(seq(-8,6, 0.1)),
                                                                                    free_cores = 5,
                                                                                    test.n= c(5,5),
                                                                                    iter=50,
                                                                                    log_AUC=2,
                                                                                    Patient.Z=2,
                                                                                    Drug.Z =2,
                                                                                    RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}

list_variable_withdrawal_elnet <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal_elnet[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                          y_data=data.matrix(Y),
                                                                                          Reduce="Naive",
                                                                                          alpha=0.4,
                                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                                          free_cores = 5,
                                                                                          test.n= c(5,5),
                                                                                          iter=50,
                                                                                          log_AUC=2,
                                                                                          Patient.Z=2,
                                                                                          Drug.Z =2,
                                                                                          RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}

list_variable_withdrawal_lasso <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal_lasso[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                          y_data=data.matrix(Y),
                                                                                          Reduce="Naive",
                                                                                          alpha=1,
                                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                                          free_cores = 5,
                                                                                          test.n= c(5,5),
                                                                                          iter=50,
                                                                                          log_AUC=2,
                                                                                          Patient.Z=2,
                                                                                          Drug.Z =2,
                                                                                          RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}
save(list_variable_withdrawal,list_variable_withdrawal_elnet,list_variable_withdrawal_lasso,
     file = paste0(path, "/Results_survival predictions/Variable withdrawal results.RData"))



# Variable association statistics - bootstrapping - sts ----
Y1 <- Y
Y1$status[Y1$time>365*2] <- 0
Y1$time[Y1$time>365*2] <- 365*2
list_variable_association <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                          y_data=data.matrix(Y1),
                                                                          alpha=0,
                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                          free_cores = 5,
                                                                          pre.CV = F,
                                                                          lambda_opt = 0,
                                                                          iter=200,
                                                                          log_AUC=2,
                                                                          Patient.Z=2,
                                                                          Drug.Z =2,
                                                                          RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}

list_variable_association_elnet <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association_elnet[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                                y_data=data.matrix(Y1),
                                                                                alpha=0.4,
                                                                                lambda=exp(seq(-8,6, 0.1)),
                                                                                free_cores = 5,
                                                                                pre.CV = F,
                                                                                lambda_opt = 0,
                                                                                iter=200,
                                                                                log_AUC=2,
                                                                                Patient.Z=2,
                                                                                Drug.Z =2,
                                                                                RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}

list_variable_association_lasso <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  SDs <- matrixStats::colSds(data.matrix(X_rAUC))
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  for(tr in c(quantile(SDs, 1-c(40)/ncol(X_rAUC)),0)){
    X1 <- X[,SDs>=tr]
    list_variable_association_lasso[[paste(m, ncol(X1))]] <-  Cox_bootstrapping(X_data=X1,
                                                                                y_data=data.matrix(Y1),
                                                                                alpha=1,
                                                                                lambda=exp(seq(-8,6, 0.1)),
                                                                                free_cores = 5,
                                                                                pre.CV = F,
                                                                                lambda_opt = 0,
                                                                                iter=200,
                                                                                log_AUC=2,
                                                                                Patient.Z=2,
                                                                                Drug.Z =2,
                                                                                RCPC=0)
    cat("\n", m, ncol(X1), "\n")
  }
}


save(list_variable_association, list_variable_association_elnet, list_variable_association_lasso,
     file = paste0(path, "/Results_survival predictions/Variable association results_sts.RData"))


# Variable withdrawal testing  - sts ----
Cox_forecasting_drug_withdrawal <- function(X_data,
                                            y_data,
                                            Reduce=c("None", "Naive"),
                                            alpha=0,
                                            lambda=c(exp(seq(-4,6, 0.1))),
                                            free_cores = 2,
                                            test.n= c(5,5),
                                            iter=200,
                                            log_AUC=2,
                                            Patient.Z=2,
                                            Drug.Z =2,
                                            RCPC=0,
                                            path=path){
  # Data preparation
  a = alpha
  Transform = c("log2(AUC)", "AUC")[log_AUC[1]]
  pt.st = c(TRUE, FALSE)[Patient.Z[1]]
  cpc = RCPC[1]
  drug.st = c(TRUE, FALSE)[Drug.Z[1]]
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
  colnames(X) <- colnames(X_data)
  X_data <- X
  
  Drug_WD_test <- function(X_data,
                           y_data,
                           alpha=a,
                           lambda=c(exp(seq(-4,6, 0.1))),
                           free_cores = 2,
                           test.n= c(6,4),
                           iter=200,
                           log_AUC=c(1:2),
                           Patient.Z=c(1:2),
                           Drug.Z =c(1:2),
                           RCPC=c(0,1,2,3,4),
                           i=i){
    a=alpha
    # Prediction function
    # Pred_cox <- function(X,
    #                      Y,
    #                      test.n=test.n,
    #                      a=a,lambda=lambda, iter=iter, i=i){
    #   
    #   cores <- min(c(parallel::detectCores()-free_cores, iter))
    #   cluster.cores<-makeCluster(cores)
    #   registerDoSNOW(cluster.cores)
    #   pb <- txtProgressBar(max=iter, style=3)
    #   progress <- function(n) setTxtProgressBar(pb, n)
    #   opts <- list(progress=progress)
    #   list_C_index <-foreach(b = 1:iter, .packages = c("glmnet"), .options.snow=opts) %dopar%{
    #     set.seed(b)
    #     setTxtProgressBar(pb,b)
    #     X.0 <- X[Y[,2]==0,]
    #     X.1 <- X[Y[,2]==1,]
    #     Y.0 <- Y[Y[,2]==0,]
    #     Y.1 <- Y[Y[,2]==1,]
    #     
    #     ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
    #     ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)
    #     
    #     train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
    #     test <- rbind(X.0[ind.0,], X.1[ind.1,])
    #     
    #     Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
    #     Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])
    #     
    #     model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
    #     nfolds=nrow(train)
    #     model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
    #     
    #     c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
    #     ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))
    #     
    #     rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
    #     
    #     return(c(c,ct))
    #   }
    #   stopCluster(cluster.cores)
    #   closeAllConnections()
    #   rm(cluster.cores)
    #   rm(cores)
    #   gc()
    #   
    #   df_C_index <- do.call(rbind, list_C_index)
    #   df_C_index <- data.frame(df_C_index)
    #   colnames(df_C_index) <- c("C_index_test", "C_index_train")
    #   df_C_index$ID <- i
    #   
    #   rm(X,Y, list_C_index, opts, pb)
    #   
    #   
    #   return(df_C_index)
    # }
    # 
    
    
    # Run on full data
    # df_parent <- Pred_cox(X_data, y_data, test.n, a,lambda, iter, i)
    
    # Prepare drug withdrawal data
    mat.C.index.test <- data.frame(array(data=NA, dim=c(iter, ncol(X_data))))
    colnames(mat.C.index.test) <- colnames(X_data)
    mat.C.index.train <- mat.C.index.test
    
    cores <- min(c(parallel::detectCores()-free_cores, ncol(X_data)))
    cluster.cores<-makeCluster(cores)
    registerDoSNOW(cluster.cores)
    pb <- txtProgressBar(max=ncol(X_data), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_C_index <-foreach(d = 0:ncol(X_data), .packages = c("glmnet"), .options.snow=opts) %dopar%{
      
      
      if(d != 0){
        dr <- colnames(X_data)[d]
        X_data1 <- X_data[,colnames(X_data) != dr]
      }else{
        X_data1 <- X_data
        
      }
      
      Y <- data.matrix(y_data)
      
      list_C_index <- list()
      for(b in 1:iter){
        set.seed(b)
        X.0 <- X_data1[Y[,2]==0,]
        X.1 <- X_data1[Y[,2]==1,]
        Y.0 <- Y[Y[,2]==0,]
        Y.1 <- Y[Y[,2]==1,]
        
        ind.0 <- sample(seq_len(nrow(Y.0)), size = test.n[1] , replace=FALSE)
        ind.1 <- sample(seq_len(nrow(Y.1)), size = test.n[2] , replace=FALSE)
        
        train <- rbind(X.0[-ind.0,], X.1[-ind.1,])
        test <- rbind(X.0[ind.0,], X.1[ind.1,])
        
        Y.train.loop <- rbind(Y.0[-ind.0,], Y.1[-ind.1,])
        Y.test.loop <- rbind(Y.0[ind.0,], Y.1[ind.1,])
        
        model.loop <- glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance")
        nfolds=nrow(train)
        model.loop.cv <- cv.glmnet(train, Y.train.loop, family = "cox", alpha = a, standardize = FALSE,lambda=lambda,  type.measure = "deviance", nfolds=nfolds)
        
        c=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(test)), y=data.matrix(Y.test.loop))
        ct=Cindex(predict(model.loop, s = model.loop.cv$lambda.min, newx= data.matrix(train)), y=data.matrix(Y.train.loop))
        
        #rm(X.0, X.1, Y.1, Y.0, train, test, Y.train.loop, Y.test.loop, model.loop, model.loop.cv, ind.0, ind.1)
        list_C_index[[b]] <- c(c,ct)
      }
      
      df_C_index <- do.call(rbind, list_C_index)
      df_C_index <- data.frame(df_C_index)
      colnames(df_C_index) <- c("C_index_test", "C_index_train")
      df_C_index$ID <- i
      
      gc()
      
      df_wd <- df_C_index
      
      #rm(X_data1)
      return(df_wd)
    }
    stopCluster(cluster.cores)
    closeAllConnections()
    rm(cluster.cores)
    rm(cores)
    gc()
    
    cat("\n", "WD test completed...", "\n")
    
    names(list_C_index) <- c("Complete",colnames(X_data))
    for(dr in  names(list_C_index)){
      mat.C.index.test[,dr] <- list_C_index[[dr]]$C_index_test
      mat.C.index.train[,dr] <- list_C_index[[dr]]$C_index_train
    }
    
    C_loss_test <- colMeans((mat.C.index.test - list_C_index[["Complete"]]$C_index_test))
    C_loss_train <- colMeans((mat.C.index.train - list_C_index[["Complete"]]$C_index_train))
    
    df_loss <- data.frame(C_loss_test=C_loss_test,
                          C_loss_train=C_loss_train,
                          Var=names(C_loss_train))
    
    
    return(list(df_loss=df_loss,
                mat.C.index.test=mat.C.index.test,
                mat.C.index.train=mat.C.index.train))
  }
  
  list_res_initial <- Drug_WD_test(X_data,
                                   y_data,
                                   alpha=alpha,
                                   lambda=lambda,
                                   free_cores = free_cores,
                                   test.n= test.n,
                                   iter=iter,
                                   log_AUC=log_AUC,
                                   Patient.Z=Patient.Z,
                                   Drug.Z =Drug.Z,
                                   RCPC=RCPC,
                                   i=i)
  return(list_res_initial)
}


list_variable_withdrawal <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                    y_data=data.matrix(Y1),
                                                                                    Reduce="Naive",
                                                                                    alpha=0,
                                                                                    lambda=exp(seq(-8,6, 0.1)),
                                                                                    free_cores = 5,
                                                                                    test.n= c(6,4),
                                                                                    iter=50,
                                                                                    log_AUC=2,
                                                                                    Patient.Z=2,
                                                                                    Drug.Z =2,
                                                                                    RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}

list_variable_withdrawal_elnet <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal_elnet[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                          y_data=data.matrix(Y1),
                                                                                          Reduce="Naive",
                                                                                          alpha=0.4,
                                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                                          free_cores = 5,
                                                                                          test.n= c(6,4),
                                                                                          iter=50,
                                                                                          log_AUC=2,
                                                                                          Patient.Z=2,
                                                                                          Drug.Z =2,
                                                                                          RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}

list_variable_withdrawal_lasso <- list()
for(m in colnames(df_scores)[c(7,12)]){
  X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
  rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
  X <- X[match(rownames(Y1), rownames(X)),]
  X <- (X - rowMeans(X))/matrixStats::rowSds(data.matrix(X))
  
  list_variable_withdrawal_lasso[[paste(m, ncol(X))]] <-  Cox_forecasting_drug_withdrawal(X_data=X,
                                                                                          y_data=data.matrix(Y1),
                                                                                          Reduce="Naive",
                                                                                          alpha=1,
                                                                                          lambda=exp(seq(-8,6, 0.1)),
                                                                                          free_cores = 5,
                                                                                          test.n= c(6,4),
                                                                                          iter=50,
                                                                                          log_AUC=2,
                                                                                          Patient.Z=2,
                                                                                          Drug.Z =2,
                                                                                          RCPC=0)
  cat("\n", m, ncol(X), "\n")
  
}
save(list_variable_withdrawal,list_variable_withdrawal_elnet,list_variable_withdrawal_lasso,
     file = paste0(path, "/Results_survival predictions/Variable withdrawal results_sts.RData"))


# Variable association statistics - results ----

load(file = paste0(path, "/Results_survival predictions/Variable association results.RData"))



for(m in names(list_variable_association)){
  list_variable_association[[m]]$Betas_statistics <-cbind(list_variable_association[[m]]$Betas_statistics,
                                                          data.frame(list_variable_association[[m]]$Betas %>% 
                                                                       group_by(drug) %>% 
                                                                       summarise(Q025 = quantile(Beta, 0.025),
                                                                                 Q5 = quantile(Beta, 0.5),
                                                                                 Q975 = quantile(Beta, 0.975)))[,2:4])
  list_variable_association_elnet[[m]]$Betas_statistics <-cbind(list_variable_association_elnet[[m]]$Betas_statistics,
                                                                data.frame(list_variable_association_elnet[[m]]$Betas %>% 
                                                                             group_by(drug) %>% 
                                                                             summarise(Q025 = quantile(Beta, 0.025),
                                                                                       Q5 = quantile(Beta, 0.5),
                                                                                       Q975 = quantile(Beta, 0.975)))[,2:4])
  list_variable_association_lasso[[m]]$Betas_statistics <-cbind(list_variable_association_lasso[[m]]$Betas_statistics,
                                                                data.frame(list_variable_association_lasso[[m]]$Betas %>% 
                                                                             group_by(drug) %>% 
                                                                             summarise(Q025 = quantile(Beta, 0.025),
                                                                                       Q5 = quantile(Beta, 0.5),
                                                                                       Q975 = quantile(Beta, 0.975)))[,2:4])
}




df_drugs <- data.frame(Drug=colnames(X_rAUC), Label= NA)
df_drugs$Label[df_drugs$Drug %in% c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")] <- "Standard treatment"
df_drugs$Label[df_drugs$Drug %in% c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
                                    "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")] <- "Alternate treament"
df_drugs <- na.omit(df_drugs)


df_betas <- rbind(list_variable_association$`rAUC_log2 349`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)

save(df_betas_lab, file = paste0(path, "/Results_survival predictions/Signficiant variable associations.RData"))



ggexport(ggplot(df_betas, aes(x=Rank, y=Mean))+
           geom_errorbar(aes(ymin=Q025 , ymax=Q975),size=0.1,alpha=0.25,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=df_betas_lab, aes(x=Rank, y=Mean,col=Label),size=0.9, alpha=0.99)+
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
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "Model coefficient (95% CI)", x = "t-statistic rank", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 12, height=8.5, 
         filename = paste0(path,"/Results_survival predictions/Model coefficients_full model.pdf"))



df_betas_modelcomp <- dcast(df_betas, drug+Model ~ Data, value.var = "Mean")
df_betas_modelcomp$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_modelcomp$drug)
df_betas_modelcomp$Label <- df_betas_lab$Label[match(df_betas_modelcomp$drug, df_betas_lab$drug)]

ggexport(ggplot(df_betas_modelcomp, aes(x=DSS3, y=rAUC_log2))+
           geom_abline(lty=2)+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15 | Label != ""), 
                      aes(col=Label),size=0.9, alpha=0.99)+
           stat_cor(label.sep=',')+
           geom_text_repel(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15) ,
                           aes(x=DSS3,y=rAUC_log2, label=drug),
                           nudge_y      = 0.02,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(~Model, scales="free")+
           labs(y= "Model coefficients\nrAUC_log2", x = "Model coefficients\nDSS3", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 8, height=3.85, 
         filename = paste0(path,"/Results_survival predictions/Model coefficient correlation_full model.pdf"))





df_betas <- rbind(list_variable_association$`rAUC_log2 40`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 40`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 40`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)


df_betas_modelcomp <- dcast(df_betas, drug+Model ~ Data, value.var = "Mean")
df_betas_modelcomp$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_modelcomp$drug)
df_betas_modelcomp$Label <- df_betas_lab$Label[match(df_betas_modelcomp$drug, df_betas_lab$drug)]

ggexport(ggplot(df_betas_modelcomp, aes(x=DSS3, y=rAUC_log2))+
           geom_abline(lty=2)+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15 | Label != ""), 
                      aes(col=Label),size=0.9, alpha=0.99)+
           stat_cor(label.sep=',')+
           geom_text_repel(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15) ,
                           aes(x=DSS3,y=rAUC_log2, label=drug),
                           nudge_y      = 0.02,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(~Model, scales="free")+
           labs(y= "Model coefficients\nrAUC_log2", x = "Model coefficients\nDSS3", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 8, height=3.85, 
         filename = paste0(path,"/Results_survival predictions/Model coefficient correlation_reduced model (40).pdf"))
# Variable withdrawal results ----

load(file = paste0(path, "/Results_survival predictions/Variable withdrawal results.RData"))

for(m in names(list_variable_withdrawal)){
  list_variable_withdrawal[[m]]$df_loss <-cbind(list_variable_withdrawal[[m]]$df_loss[match(colnames(list_variable_withdrawal[[m]]$mat.C.index.test),list_variable_withdrawal[[m]]$df_loss$Var),],
                                                Mean=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                SD=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                Q025=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                Q5=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                Q975=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                C_index_mean=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                C_index_sd=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                C_index_Q025=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                C_index_Q5=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                C_index_Q975=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
  list_variable_withdrawal_elnet[[m]]$df_loss <-cbind(list_variable_withdrawal_elnet[[m]]$df_loss[match(colnames(list_variable_withdrawal_elnet[[m]]$mat.C.index.test),list_variable_withdrawal_elnet[[m]]$df_loss$Var),],
                                                      Mean=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                      SD=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                      Q025=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                      Q5=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                      Q975=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                      C_index_mean=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                      C_index_sd=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                      C_index_Q025=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                      C_index_Q5=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                      C_index_Q975=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
  list_variable_withdrawal_lasso[[m]]$df_loss <-cbind(list_variable_withdrawal_lasso[[m]]$df_loss[match(colnames(list_variable_withdrawal_lasso[[m]]$mat.C.index.test),list_variable_withdrawal_lasso[[m]]$df_loss$Var),],
                                                      Mean=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                      SD=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                      Q025=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                      Q5=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                      Q975=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                      C_index_mean=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                      C_index_sd=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                      C_index_Q025=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                      C_index_Q5=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                      C_index_Q975=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
}


df_drugs <- data.frame(Drug=colnames(X_rAUC), Label= NA)
df_drugs$Label[df_drugs$Drug %in% c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")] <- "Standard treatment"
df_drugs$Label[df_drugs$Drug %in% c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
                                    "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")] <- "Alternate treament"
df_drugs <- na.omit(df_drugs)

df_betas <- rbind(list_variable_association$`rAUC_log2 349`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)

df_loss <- rbind(list_variable_withdrawal$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Ridge"),
                 list_variable_withdrawal$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Ridge"),
                 list_variable_withdrawal_elnet$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Elastic net"),
                 list_variable_withdrawal_elnet$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Elastic net"),
                 list_variable_withdrawal_lasso$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Lasso"),
                 list_variable_withdrawal_lasso$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Lasso"))
df_loss$Model <- factor(df_loss$Model, levels=unique(df_loss$Model))
df_loss <- na.omit(df_loss)
df_loss$SD <- df_loss$SD#/sqrt(50)
df_loss$Beta_mean <- df_betas$Mean[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                         paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss$Beta_sd <- df_betas$StErr[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                        paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss <- df_loss %>% group_by(Model, Data) %>% #
  mutate(Rank = rank(C_loss_test,ties.method="random")) %>% as.data.frame()#
df_loss_lab <- subset(df_loss, Mean < -0.0025 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse(df_loss_lab$Q025 > 0 & df_loss_lab$Q975 > 0 | df_loss_lab$Q025 < 0 & df_loss_lab$Q975 < 0, "*", "")
df_loss_lab$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss_lab$Var)
df_loss$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss$Var)

df_betas_lab$Mean_loss <- df_loss$Mean[match(paste(as.character(df_betas_lab$Model), df_betas_lab$Data, df_betas_lab$drug),
                                             paste(as.character(df_loss$Model), df_loss$Data, df_loss$Var))]

ggexport(ggplot(df_loss, aes(x=abs(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=abs(Beta_mean)-Beta_sd , xmax=abs(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=abs(Mean), y=Mean_loss, col=Label),size=0.9)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=5.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_beta relation_1.pdf"))


ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=(Mean), y=Mean_loss, col=Label),size=0.9, alpha=0.99)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=5.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_beta relation_2.pdf"))



ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=(Mean), y=Mean_loss, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_betas_lab ,
                           
                           aes(x=(Mean), y=Mean_loss, label=paste0(drug, Sign), col=Label),
                           force      =10,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           hjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 9.5, height=7.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_beta relation_3.pdf"))


ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_loss_lab, aes(x=(Beta_mean), y=Mean, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_loss_lab ,
                           
                           aes(x=(Beta_mean), y=Mean, label=paste0(Var, Sign), col=Label),
                           force      =10,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           hjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 9.5, height=7.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_beta relation_4.pdf"))







df_loss <- rbind(list_variable_withdrawal$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Ridge"),
                 list_variable_withdrawal$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Ridge"),
                 list_variable_withdrawal_elnet$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Elastic net"),
                 list_variable_withdrawal_elnet$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Elastic net"),
                 list_variable_withdrawal_lasso$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Lasso"),
                 list_variable_withdrawal_lasso$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Lasso"))
df_loss$Model <- factor(df_loss$Model, levels=unique(df_loss$Model))
df_loss <- na.omit(df_loss)
df_loss$SD <- df_loss$SD#/sqrt(50)
df_loss$Beta_mean <- df_betas$Mean[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                         paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss$Beta_sd <- df_betas$StErr[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                        paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss <- df_loss %>% group_by(Model, Data) %>% subset(Mean != 0) %>%
  mutate(Rank = rank(C_loss_test,ties.method="random")) %>% as.data.frame()#
df_loss_lab <- subset(df_loss, Mean < -0.0025 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse(df_loss_lab$Q025 > 0 & df_loss_lab$Q975 > 0 | df_loss_lab$Q025 < 0 & df_loss_lab$Q975 < 0, "*", "")
df_loss_lab$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss_lab$Var)
df_loss$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss$Var)



ggexport(ggplot(df_loss, aes(x=Rank, y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),size=0.1,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_loss_lab, aes(x=Rank, y=Mean, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_loss_lab ,
                           
                           aes(x=Rank, y=Mean, label=paste0(Var, Sign), col=Label),
                           nudge_y      =-0.02,
                           nudge_x=50,
                           direction    = "y",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Rank", title="Variable importance")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 12, height=8.5, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_rank_1.pdf"))


# Variable association statistics - results -sts ----

load(file = paste0(path, "/Results_survival predictions/Variable association results_sts.RData"))



for(m in names(list_variable_association)){
  list_variable_association[[m]]$Betas_statistics <-cbind(list_variable_association[[m]]$Betas_statistics,
                                                          data.frame(list_variable_association[[m]]$Betas %>% 
                                                                       group_by(drug) %>% 
                                                                       summarise(Q025 = quantile(Beta, 0.025),
                                                                                 Q5 = quantile(Beta, 0.5),
                                                                                 Q975 = quantile(Beta, 0.975)))[,2:4])
  list_variable_association_elnet[[m]]$Betas_statistics <-cbind(list_variable_association_elnet[[m]]$Betas_statistics,
                                                                data.frame(list_variable_association_elnet[[m]]$Betas %>% 
                                                                             group_by(drug) %>% 
                                                                             summarise(Q025 = quantile(Beta, 0.025),
                                                                                       Q5 = quantile(Beta, 0.5),
                                                                                       Q975 = quantile(Beta, 0.975)))[,2:4])
  list_variable_association_lasso[[m]]$Betas_statistics <-cbind(list_variable_association_lasso[[m]]$Betas_statistics,
                                                                data.frame(list_variable_association_lasso[[m]]$Betas %>% 
                                                                             group_by(drug) %>% 
                                                                             summarise(Q025 = quantile(Beta, 0.025),
                                                                                       Q5 = quantile(Beta, 0.5),
                                                                                       Q975 = quantile(Beta, 0.975)))[,2:4])
}




df_drugs <- data.frame(Drug=colnames(X_rAUC), Label= NA)
df_drugs$Label[df_drugs$Drug %in% c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")] <- "Standard treatment"
df_drugs$Label[df_drugs$Drug %in% c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
                                    "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")] <- "Alternate treament"
df_drugs <- na.omit(df_drugs)


df_betas <- rbind(list_variable_association$`rAUC_log2 349`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)

save(df_betas_lab, file = paste0(path, "/Results_survival predictions/Signficiant variable associations_sts.RData"))



ggexport(ggplot(df_betas, aes(x=Rank, y=Mean))+
           geom_errorbar(aes(ymin=Q025 , ymax=Q975),size=0.1,alpha=0.25,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=df_betas_lab, aes(x=Rank, y=Mean,col=Label),size=0.9, alpha=0.99)+
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
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "Model coefficient (95% CI)", x = "t-statistic rank", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 12, height=8.5, 
         filename = paste0(path,"/Results_survival predictions/Model coefficients_full model_sts.pdf"))



df_betas_modelcomp <- dcast(df_betas, drug+Model ~ Data, value.var = "Mean")
df_betas_modelcomp$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_modelcomp$drug)
df_betas_modelcomp$Label <- df_betas_lab$Label[match(df_betas_modelcomp$drug, df_betas_lab$drug)]

ggexport(ggplot(df_betas_modelcomp, aes(x=DSS3, y=rAUC_log2))+
           geom_abline(lty=2)+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15 | Label != ""), 
                      aes(col=Label),size=0.9, alpha=0.99)+
           stat_cor(label.sep=',')+
           geom_text_repel(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15) ,
                           aes(x=DSS3,y=rAUC_log2, label=drug),
                           nudge_y      = 0.02,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(~Model, scales="free")+
           labs(y= "Model coefficients\nrAUC_log2", x = "Model coefficients\nDSS3", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 8, height=3.85, 
         filename = paste0(path,"/Results_survival predictions/Model coefficient correlation_full model_sts.pdf"))





df_betas <- rbind(list_variable_association$`rAUC_log2 40`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 40`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 40`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 40`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)


df_betas_modelcomp <- dcast(df_betas, drug+Model ~ Data, value.var = "Mean")
df_betas_modelcomp$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_modelcomp$drug)
df_betas_modelcomp$Label <- df_betas_lab$Label[match(df_betas_modelcomp$drug, df_betas_lab$drug)]

ggexport(ggplot(df_betas_modelcomp, aes(x=DSS3, y=rAUC_log2))+
           geom_abline(lty=2)+
           geom_point(size=0.9,col=pal_jco()(3)[3], alpha=0.99, pch=16)+
           geom_point(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15 | Label != ""), 
                      aes(col=Label),size=0.9, alpha=0.99)+
           stat_cor(label.sep=',')+
           geom_text_repel(data=subset(df_betas_modelcomp, sqrt(DSS3*rAUC_log2)>0.15) ,
                           aes(x=DSS3,y=rAUC_log2, label=drug),
                           nudge_y      = 0.02,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(~Model, scales="free")+
           labs(y= "Model coefficients\nrAUC_log2", x = "Model coefficients\nDSS3", title="log(HR) association")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 8, height=3.85, 
         filename = paste0(path,"/Results_survival predictions/Model coefficient correlation_reduced model (40)_sts.pdf"))

# Variable withdrawal results -sts ----

load(file = paste0(path, "/Results_survival predictions/Variable withdrawal results_sts.RData"))

for(m in names(list_variable_withdrawal)){
  list_variable_withdrawal[[m]]$df_loss <-cbind(list_variable_withdrawal[[m]]$df_loss[match(colnames(list_variable_withdrawal[[m]]$mat.C.index.test),list_variable_withdrawal[[m]]$df_loss$Var),],
                                                Mean=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                SD=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                Q025=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                Q5=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                Q975=apply((list_variable_withdrawal[[m]]$mat.C.index.test - list_variable_withdrawal[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                C_index_mean=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                C_index_sd=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                C_index_Q025=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                C_index_Q5=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                C_index_Q975=apply(list_variable_withdrawal[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
  list_variable_withdrawal_elnet[[m]]$df_loss <-cbind(list_variable_withdrawal_elnet[[m]]$df_loss[match(colnames(list_variable_withdrawal_elnet[[m]]$mat.C.index.test),list_variable_withdrawal_elnet[[m]]$df_loss$Var),],
                                                      Mean=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                      SD=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                      Q025=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                      Q5=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                      Q975=apply((list_variable_withdrawal_elnet[[m]]$mat.C.index.test - list_variable_withdrawal_elnet[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                      C_index_mean=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                      C_index_sd=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                      C_index_Q025=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                      C_index_Q5=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                      C_index_Q975=apply(list_variable_withdrawal_elnet[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
  list_variable_withdrawal_lasso[[m]]$df_loss <-cbind(list_variable_withdrawal_lasso[[m]]$df_loss[match(colnames(list_variable_withdrawal_lasso[[m]]$mat.C.index.test),list_variable_withdrawal_lasso[[m]]$df_loss$Var),],
                                                      Mean=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  mean(x)),
                                                      SD=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  sd(x)),
                                                      Q025=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.025)),
                                                      Q5=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.5)),
                                                      Q975=apply((list_variable_withdrawal_lasso[[m]]$mat.C.index.test - list_variable_withdrawal_lasso[[m]]$mat.C.index.test$Complete), 2, function(x)  quantile(x, 0.975)),
                                                      C_index_mean=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  mean(x)),
                                                      C_index_sd=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  sd(x)),
                                                      C_index_Q025=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.025)),
                                                      C_index_Q5=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.5)),
                                                      C_index_Q975=apply(list_variable_withdrawal_lasso[[m]]$mat.C.index.test, 2, function(x)  quantile(x, 0.975)))
}


df_drugs <- data.frame(Drug=colnames(X_rAUC), Label= NA)
df_drugs$Label[df_drugs$Drug %in% c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")] <- "Standard treatment"
df_drugs$Label[df_drugs$Drug %in% c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
                                    "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")] <- "Alternate treament"
df_drugs <- na.omit(df_drugs)

df_betas <- rbind(list_variable_association$`rAUC_log2 349`$Betas_statistics %>% mutate(Data="rAUC_log2", Model="Ridge"),
                  list_variable_association$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Ridge"),
                  list_variable_association_elnet$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Elastic net"),
                  list_variable_association_elnet$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Elastic net"),
                  list_variable_association_lasso$`rAUC_log2 349`$Betas_statistics%>% mutate(Data="rAUC_log2", Model="Lasso"),
                  list_variable_association_lasso$`DSS3 349`$Betas_statistics%>% mutate(Data="DSS3", Model="Lasso"))
df_betas$Model <- factor(df_betas$Model, levels=unique(df_betas$Model))
df_betas <- na.omit(df_betas)
df_betas_lab <- subset(df_betas, Q025 >= 0 & Q5 > 0 | Q975 <= 0 & Q5 < 0 | drug %in% df_drugs$Drug)
df_betas_lab$Label <- df_drugs$Label[match(df_betas_lab$drug, df_drugs$Drug)]
df_betas_lab$Label[is.na(df_betas_lab$Label)] <- ""
df_betas_lab$Sign <- ifelse(df_betas_lab$Q025 > 0 & df_betas_lab$Q975 > 0 | df_betas_lab$Q025 < 0 & df_betas_lab$Q975 < 0, "*", "")
df_betas_lab$drug <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_betas_lab$drug)

df_loss <- rbind(list_variable_withdrawal$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Ridge"),
                 list_variable_withdrawal$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Ridge"),
                 list_variable_withdrawal_elnet$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Elastic net"),
                 list_variable_withdrawal_elnet$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Elastic net"),
                 list_variable_withdrawal_lasso$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Lasso"),
                 list_variable_withdrawal_lasso$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Lasso"))
df_loss$Model <- factor(df_loss$Model, levels=unique(df_loss$Model))
df_loss <- na.omit(df_loss)
df_loss$SD <- df_loss$SD#/sqrt(50)
df_loss$Beta_mean <- df_betas$Mean[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                         paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss$Beta_sd <- df_betas$StErr[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                        paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss <- df_loss %>% group_by(Model, Data) %>% #
  mutate(Rank = rank(C_loss_test,ties.method="random")) %>% as.data.frame()#
df_loss_lab <- subset(df_loss, Mean < -0.0025 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse(df_loss_lab$Q025 > 0 & df_loss_lab$Q975 > 0 | df_loss_lab$Q025 < 0 & df_loss_lab$Q975 < 0, "*", "")
df_loss_lab$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss_lab$Var)
df_loss$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss$Var)

df_betas_lab$Mean_loss <- df_loss$Mean[match(paste(as.character(df_betas_lab$Model), df_betas_lab$Data, df_betas_lab$drug),
                                             paste(as.character(df_loss$Model), df_loss$Data, df_loss$Var))]

ggexport(ggplot(df_loss, aes(x=abs(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=abs(Beta_mean)-Beta_sd , xmax=abs(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=abs(Mean), y=Mean_loss, col=Label),size=0.9, alpha=0.99)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=5.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_sts_beta relation_1.pdf"))


ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=(Mean), y=Mean_loss, col=Label),size=0.9, alpha=0.99)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 7, height=5.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_sts_beta relation_2.pdf"))



ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_betas_lab, aes(x=(Mean), y=Mean_loss, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_betas_lab ,
                           
                           aes(x=(Mean), y=Mean_loss, label=paste0(drug, Sign), col=Label),
                           force      =10,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           hjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 9.5, height=7.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_sts_beta relation_3.pdf"))


ggexport(ggplot(df_loss, aes(x=(Beta_mean), y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),
                         size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_errorbarh(aes(xmin=(Beta_mean)-Beta_sd , xmax=(Beta_mean)+Beta_sd),
                          size=0.2,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_loss_lab, aes(x=(Beta_mean), y=Mean, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_loss_lab ,
                           
                           aes(x=(Beta_mean), y=Mean, label=paste0(Var, Sign), col=Label),
                           force      =10,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           hjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Model coefficient")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 9.5, height=7.6, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_sts_beta relation_4.pdf"))







df_loss <- rbind(list_variable_withdrawal$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Ridge"),
                 list_variable_withdrawal$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Ridge"),
                 list_variable_withdrawal_elnet$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Elastic net"),
                 list_variable_withdrawal_elnet$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Elastic net"),
                 list_variable_withdrawal_lasso$`rAUC_log2 349`$df_loss %>% mutate(Data="rAUC_log2", Model="Lasso"),
                 list_variable_withdrawal_lasso$`DSS3 349`$df_loss %>% mutate(Data="DSS3", Model="Lasso"))
df_loss$Model <- factor(df_loss$Model, levels=unique(df_loss$Model))
df_loss <- na.omit(df_loss)
df_loss$SD <- df_loss$SD#/sqrt(50)
df_loss$Beta_mean <- df_betas$Mean[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                         paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss$Beta_sd <- df_betas$StErr[match(paste(df_loss$Model, df_loss$Data, df_loss$Var),
                                        paste(df_betas$Model, df_betas$Data, df_betas$drug))]
df_loss <- df_loss %>% group_by(Model, Data) %>% subset(Mean != 0) %>%
  mutate(Rank = rank(C_loss_test,ties.method="random")) %>% as.data.frame()#
df_loss_lab <- subset(df_loss, Mean < -0.0025 | Var %in% df_drugs$Drug)
df_loss_lab$Label <- df_drugs$Label[match(df_loss_lab$Var, df_drugs$Drug)]
df_loss_lab$Label[is.na(df_loss_lab$Label)] <- ""
df_loss_lab$Sign <- ifelse(df_loss_lab$Q025 > 0 & df_loss_lab$Q975 > 0 | df_loss_lab$Q025 < 0 & df_loss_lab$Q975 < 0, "*", "")
df_loss_lab$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss_lab$Var)
df_loss$Var <- gsub(" \\(.*| HCl| 2HCl| sodium salt", "", df_loss$Var)



ggexport(ggplot(df_loss, aes(x=Rank, y=Mean))+
           geom_errorbar(aes(ymin=Mean-SD , ymax=Mean+SD),size=0.1,alpha=0.5,col=pal_jco()(3)[3])+
           geom_point(size=0.9,col=pal_jco()(3)[3], pch=16, alpha=0.99)+
           geom_point(data=df_loss_lab, aes(x=Rank, y=Mean, col=Label),size=0.9, alpha=0.99)+
           geom_text_repel(data=df_loss_lab ,
                           
                           aes(x=Rank, y=Mean, label=paste0(Var, Sign), col=Label),
                           nudge_y      =-0.02,
                           nudge_x=50,
                           direction    = "y",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           angle        = 0,
                           vjust        = 0.2,
                           segment.size = 0.2,
                           size=2.5)+
           scale_color_manual(values=c("black", "#00618D", pal_jco()(9)[9]))+
           facet_wrap(Data~Model, scales="free")+
           labs(y= "C-index change", x = "Rank", title="Variable importance")+
           theme_bw()+
           theme(strip.background = element_blank(),
                 legend.position = "top"),
         width = 12, height=8.5, 
         filename = paste0(path,"/Results_survival predictions/Variable withdrawal results_sts_rank_1.pdf"))





# Prediction testing with treatment drugs only ----

Std_treat <-  c("Cytarabine", "Daunorubicin HCl (Daunomycin HCl)", "Idarubicin HCl")
Alt_treat <- c("ABT-263 (Navitoclax)", "Azacitidine (Vidaza)", "Fludarabine (Fludara)","Fludarabine Phosphate (Fludara)",
               "Etoposide (VP-16)", "Mitoxantrone HCl","Decitabine")

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

table(Y$status)/nrow(Y)*10

list_res_treatment_drugs <- list()
for(tr in c("Standard", "Alternative", "Combined")){
  for(m in colnames(df_scores)[6:14]){
    X <- dcast(data = df_scores, Patient.num ~ drug , value.var=m)
    rownames(X) <- as.character(X$Patient.num); X$Patient.num <- NULL
    X <- X[match(rownames(Y), rownames(X)),]
    if(grepl("EC50",m)){
      X <- log10(X)
      X[,] <- apply(X,2, function(x) x - mean(x))
    }
    X <- (X- rowMeans(X))/matrixStats::rowSds(data.matrix(X))
    
    if(tr == "Standard"){
      X <- X[,which(colnames(X_rAUC) %in% Std_treat)]
    }else if(tr == "Alternative"){
      X <- X[,which(colnames(X_rAUC) %in% Alt_treat)]
    }else{
      X <- X[,which(colnames(X_rAUC) %in% c(Std_treat,Alt_treat))]
    }
    
    list_res_treatment_drugs[[paste(m,tr)]] <-  Cox_forecasting_glmnet(X_data=X,
                                                        y_data=data.matrix(Y),
                                                        alpha=c(0,0.4,1),
                                                        lambda=c(exp(seq(-8,6, 0.1))),
                                                        free_cores = 5,
                                                        test.n= c(5,5),
                                                        nfolds = nrow(Y),
                                                        iter=200,
                                                        log_AUC=2,
                                                        Patient.Z=2,
                                                        Drug.Z =2,
                                                        RCPC=0)
  }
}

save(list_res_treatment_drugs, file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_treatment drugs.RData"))

load(file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_Ridge.RData"))
load(file=paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_elnet and lasso.RData"))

df_C.index.alldata <- list_res_DS_metrics_L1
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- data.frame(df_C.index.alldata[[m]]$C_index_results)
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Penalty <- as.factor(as.numeric(gsub(".*Penalty_","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}
df_C.index.alldata1 <-  bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata <- list_res_DS_metrics
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- data.frame(df_C.index.alldata[[m]]$C_index_results)
  df_C.index.alldata[[m]]$RCPC <- as.factor(as.numeric(gsub(".*RCPC_|\\/Penalty.*","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Penalty <- as.factor(as.numeric(gsub(".*Penalty_","",df_C.index.alldata[[m]]$ID)))
  df_C.index.alldata[[m]]$Data_preparation <- as.factor(gsub("\\/RCPC.*","",df_C.index.alldata[[m]]$ID))
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Metric")
df_C.index.alldata <- rbind(df_C.index.alldata,df_C.index.alldata1)
df_C.index.alldata$Patient_stdz <- ifelse(grepl("Patient", df_C.index.alldata$ID), "Patient standardization", "")
df_C.index.alldata$Drug_stdz <- ifelse(grepl("Drug", df_C.index.alldata$ID), "Drug standardization", "")
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat <- df_C.index.alldata %>% subset(Patient_stdz != "" & Drug_stdz == "") %>% 
  group_by(ID, Metric, Penalty) %>% 
  summarise(Median=median(C_index_test),
            Mean=mean(C_index_test))


df_C.index.alldata <- list_res_treatment_drugs
for(m in names(df_C.index.alldata)){
  df_C.index.alldata[[m]] <- df_C.index.alldata[[m]]$C_index_results
}
df_C.index.alldata <- bind_rows(df_C.index.alldata, .id="Data")
df_C.index.alldata$Penalty <- gsub(".*_","",df_C.index.alldata$ID)
df_C.index.alldata$Drugs <- gsub(".* ","",df_C.index.alldata$Data)
df_C.index.alldata$Drugs <- factor(df_C.index.alldata$Drugs, levels=unique(df_C.index.alldata$Drugs))
df_C.index.alldata$Metric <- gsub(" .*","",df_C.index.alldata$Data)
df_C.index.alldata$Metric <- factor(df_C.index.alldata$Metric, levels=unique(df_C.index.alldata$Metric))
df_C.index.stat1 <- df_C.index.alldata %>% group_by(ID, Metric,Drugs,Penalty) %>% summarise(Median=median(C_index_test),
                                                                             Mean=mean(C_index_test))


ggplot(df_C.index.alldata %>% subset(!grepl("EC50", Metric)), aes(x=Drugs, y=C_index_test, fill=Penalty))+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_hline(data=df_C.index.stat %>% subset(!grepl("EC50", Metric)), aes(yintercept=Median, col=Penalty), alpha=0.5)+
  geom_boxplot(alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
  facet_grid(~Metric)+
  scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
  scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
  labs(x="", y="C-index")+
  ylim(0,1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_treatment drugs_1.pdf"), width=7.5, height=2.5)

ggplot(df_C.index.alldata, aes(x=Drugs, y=C_index_test, fill=Penalty))+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_hline(data=df_C.index.stat, aes(yintercept=Median, col=Penalty), alpha=0.5)+
  geom_boxplot(alpha=0.5, outlier.size = NULL, outlier.alpha = 0)+
  facet_grid(~Metric)+
  scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
  scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
  labs(x="", y="C-index")+
  ylim(0,1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_treatment drugs_2.pdf"), width=8.5, height=2.5)


ggplot(df_C.index.stat1 , aes(x=Drugs, y=Mean, fill=Penalty))+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_hline(data=df_C.index.stat , aes(yintercept=Mean, col=Penalty), alpha=0.5)+
  geom_bar(alpha=0.5, outlier.size = NULL, outlier.alpha = 0, stat = "identity", position = position_dodge2())+
  facet_grid(~Metric)+
  scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
  scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
  labs(x="", y="C-index")+
  #ylim(0,1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_treatment drugs_3.pdf"), width=8.5, height=2.5)

ggplot(df_C.index.stat1, aes(x=Drugs, y=Median, fill=Penalty))+
  geom_hline(yintercept = 0.5, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0.75, lty=2, alpha=0.5)+
  geom_hline(data=df_C.index.stat, aes(yintercept=Median, col=Penalty), alpha=0.5)+
  geom_bar(alpha=0.5, outlier.size = NULL, outlier.alpha = 0, stat = "identity", position = position_dodge2())+
  facet_grid(~Metric)+
  scale_fill_manual(values = c("darkorange","darkgreen","darkblue"))+
  scale_color_manual(values = c("darkorange","darkgreen","darkblue"))+
  labs(x="", y="C-index")+
  #ylim(0,1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_survival predictions/Drug sensitivity score comparison_treatment drugs_4.pdf"), width=8.5, height=2.5)



