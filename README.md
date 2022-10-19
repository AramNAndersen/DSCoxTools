# DSCoxTools
Repository of tools for testing and analysis of survival predictions from ex vivo drug screens using regularized Cox regression based on the Glmnet framework.
For usage and reference see https://doi.org/10.1101/2022.10.11.509866.

The package requires dplyr, reshape, doSNOW, glmnet and glmnetUtils.

# Installation:
``` r
devtools::install_github("Enserink-lab/DSCoxTools",
                         ref="main")
``` 
# Example:
``` r
library(glmnet)
library(glmnetUtils)
library(reshape)
library(dplyr)
library(doSNOW)
library(CoxTools)

# Model prediction testing: drug sensitivity data
list_Cox_testing <- Cox_forecasting_glmnet_CVA(X_data=X.AUC, 
                                               y_data=Y, 
                                               alpha=c(0,1), 
                                               lambda=c(exp(seq(-4,6, 0.1))),
                                               free_cores = 2,
                                               test.n= c(6,4), 
                                               nfolds = nrow(y.train.red),
                                               iter=200,
                                               log_AUC=1:2,
                                               Patient.Z=1:2,
                                               Drug.Z =1:2,
                                               RCPC=0:4)
                                                              
# Drug importance scoring by variable withdrawal
Results_Cox_drug_WD <- Cox_forecasting_drug_withdrawal(X.AUC, 
                                                       Y, 
                                                       Reduce="None",
                                                       alpha=0, 
                                                       lambda=c(exp(seq(-4,6, 0.1))),
                                                       free_cores = 2,
                                                       test.n= c(6,4), 
                                                       iter=200,
                                                       log_AUC=1,
                                                       Patient.Z=1,
                                                       Drug.Z =2,
                                                       RCPC=0)

# Drug survival association statistics
Results_Cox_betas <- Cox_bootstrapping(X.AUC, 
                                       Y, 
                                       alpha=0, 
                                       lambda=c(exp(seq(-4,6, 0.1))),
                                       pre.CV=FALSE,
                                       lambda_opt = 0,
                                       free_cores = 2,
                                       iter=200,
                                       log_AUC=1,
                                       Patient.Z=1,
                                       Drug.Z =2,
                                       RCPC=0)
                                       
# Prediction testing: combination data 
list_combined_data_Cox_testing <- Cox_forecasting_glmnet_combination(X.AUC, 
                                                                     Y, 
                                                                     X.clinical,
                                                                     alpha=0, 
                                                                     lambda=c(exp(seq(-4,6, 0.1))),
                                                                     free_cores = 2,
                                                                     test.n= c(6,4), 
                                                                     iter=200,
                                                                     log_AUC=c(1:2),
                                                                     Patient.Z=c(1:2),
                                                                     Drug.Z =c(1:2),
                                                                     RCPC=c(0,2))
``` 
