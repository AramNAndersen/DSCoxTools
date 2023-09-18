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

# Path ----
path <- "/mnt/CommonStorageRAI/Aram/Drug screen - clinical forecasting/Clinical forecasting manuscript - revision 2023-07"

# Load ----
load(file = paste0(path,"/Drug sensitivity metrics and dose response data_2023-07.RData"))
load(file = paste0(path,"/Drug sensitivity QC results and batch data_2023-07.RData"))
load(file = paste0(path,"/Survival, clinical features and ELN2022 classifications_2023-07.RData"))


#Raw AUC functions
log10.AUC <- function(conc, resp){
  resp = resp[order(conc)]
  conc = conc[order(conc)]
  conc = log10(conc)
  
  A = 0
  for(j in 2:length(conc)){
    a <- (resp[j]) *(conc[j]-conc[j-1])/(max(conc)-min(conc))
    A = A + a
    rm(a,j)
  }
  A
}

#DSS/EC50 functions adapted from BREEZE (https://github.com/potdarswapnil/Breeze)
dss<-function(ic50,slope,max,min.conc.tested,max.conc.tested,y=10,DSS.type=2,concn_scale=1e-9){
  #rdata should be in in format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration
  
  a=as.numeric(unname(max))
  
  b=as.numeric(unname(slope))
  d=0 # min response
  ic50 = as.numeric(unname(ic50))
  min.conc.tested = as.numeric(unname(min.conc.tested))
  max.conc.tested = as.numeric(unname(max.conc.tested))
  Min.Conc<- log10(min.conc.tested*concn_scale) #
  Max.Conc<- max.conc.tested
  x2<-log10(Max.Conc*concn_scale)
  
  
  if(is.na(ic50)||is.na(b)||is.na(a)||is.na(Min.Conc)||is.na(Max.Conc)){
    dss<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    dss<-0
  }
  else if(isTRUE(b==0)){
    dss<-0
  }
  else{
    if(a>100){ a<-100  }
    if(isTRUE(b<0)){ b<--b  }
    c<-log10(ic50*concn_scale)
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))
        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
        else if(isTRUE(x1 > x2)){x1<-x2}
      }
      else {x1<-Min.Conc}
      
      # This is a logistic function used in Dotmatics.com
      # y = d+(a-d)/(1+10^(b*(c-x)))
      #inverse function
      # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
      int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
      
      total_area<-(x2-Min.Conc)*(100-y)
      
      if(DSS.type==1){
        norm_area<-((int_y/total_area)*100)#DSS1
      }
      if(DSS.type==2){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)/log10(a)#DSS2 #AUC1
        if(isTRUE(norm_area > 50)){ norm_area <- 0}
      }
      if(DSS.type==3){
        #	if(a>100){a<-100}
        norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3 #AUC5
      }
      if(isTRUE(norm_area < 0|norm_area > 100)){
        dss<-0
      }else{
        dss<-round(norm_area,digits=4)}
    } else {dss<-0}
  }
  return (dss)
}
#hist(100-df_dose_responses$response*100, breaks=100, xlim=c(-500,100))

CALC_IC50_EC50_DSS <- compiler::cmpfun(function(i, drug_wells_, xpr_tbl, DSS_typ, readoutCTX = F, path){
  
  tryCatch({
    #gc(T);
    TEC50 = ifelse(readoutCTX, "TC50", "EC50"); drug_wells = drug_wells_[i,];
    
    #find indices of wells with drugs
    idx_filt <- xpr_tbl$ID %in% drug_wells$ID #& xpr_tbl$ProductName %in% drug_wells$ProductName
    #extract inhib. and viab. for wells with drugs in current plate
    inhibition = inhibition2 <- xpr_tbl$inhibition_percent[idx_filt]; viability2 = 100 - inhibition2; # with 2 used for ploting of real values.
    # if there are identical values in inhibition, add a bit noise
    if(all(inhibition <= 0)) inhibition <- rep(0, length(inhibition))
    if(any(duplicated(inhibition))) inhibition <- seq(from = 0, length.out = length(inhibition), by = 0.01) + inhibition;
    #if(any(duplicated(inhibition))) inhibition <- rnorm(length(inhibition),0,0.01) + inhibition
    
    inhibition3 = inhibition2
    if(any(duplicated(inhibition3))) inhibition3 <- rnorm(length(inhibition3),0,0.01) + inhibition3
    
    viability = 100-inhibition; believe_ = T;
    
    # extract concentrations, unique drug names and product ids for wells with drugs in current plate
    dose <- as.numeric(xpr_tbl$Concentration[idx_filt])
    drug_name <- unique(as.character(xpr_tbl$ID)[idx_filt])
    product_id <- unique(as.character(xpr_tbl$ID)[idx_filt])
    
    #combine the data and sort by dose.
    mat_tbl <- data.frame(inhibition,dose,logconc = log10(dose),viability, inhibition2, viability2,inhibition3)
    mat_tbl <- mat_tbl[order(mat_tbl[,"dose"]),]
    
    
    if(DSS_typ == "AUC"){
      
      mat_tbl$indexx = 1:nrow(mat_tbl)
      model = approx(x = mat_tbl$indexx, y = mat_tbl$inhibition2, xout = seq(1,nrow(mat_tbl),length.out = 100), method="linear")
      loess_fit <- loess(y ~ x, model)
      model$y = predict(loess_fit)
      
      y_pred <- predict(loess_fit, newdata = mat_tbl %>% mutate(x=indexx))
      pred_err <- mat_tbl$inhibition2 - y_pred
      y_pred <- data.frame(y_pred=y_pred, pred_err=pred_err, dose=mat_tbl$dose)
      
      AUC <- round(sum(diff(model$x) * (head(model$y,-1)+tail(model$y,-1)))/2 / 5, 2)
   
      list(AUC=AUC,
           y_pred=y_pred)
      
    }else if(nrow(mat_tbl) <= 3 || (length(unique(mat_tbl$dose)) <= 3)){
      
      print("Less than 3 rows... skipping...")
      NULL
    } else {
      
      #############################
      #############    IC50
      
      estimate_param <- tryCatch({drm(inhibition ~ logconc, data = mat_tbl, fct = LL.4(fixed = c(NA, NA, NA,NA),names = c("SLOPE","MIN","MAX","IC50")),logDose=10,control = drmc(errorm = F))},
                                 warning=function(w){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)},
                                 error=function(e){drm(inhibition ~ logconc, data = mat_tbl, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10)})
      # (extract and name coefficients)
      coef_estim <- coef(estimate_param); names(coef_estim) <- c("SLOPE","MIN","MAX","IC50")
      # see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696819/
      coef_estim["SLOPE"] <- coef_estim["SLOPE"]*-1
      
      # if curve decreases or IC50 is higher than max (i.e. IC50 is "outlier"), set IC50 to max conc.
      coef_estim["IC50"] <- ifelse(coef_estim["MAX"]<=coef_estim["MIN"] | coef_estim["IC50"]>max(mat_tbl$dose,na.rm=T), max(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # if IC50 is less than 0 set it to min. conc. and if even min. conc. < 0, then set IC50 to mean of all conc.
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,min(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<0,mean(mat_tbl$dose,na.rm=T),coef_estim["IC50"])
      # similar to previous step but now compare log10(IC50) with log(min. conc.).
      coef_estim["IC50"] <- log10(coef_estim["IC50"])
      coef_estim["IC50"] <- ifelse(coef_estim["IC50"]<min(mat_tbl$logconc),max(mat_tbl$logconc),coef_estim["IC50"])
      # if all inhib. < 0 set IC50 to max. log. conc !!!!! not obvious why!
      coef_estim["IC50"] <- ifelse(all(mat_tbl$inhibition<0),max(mat_tbl$logconc,na.rm=T),coef_estim["IC50"])
      #(Trying to fix curves that need outlier kickout)
      coef_estim["MIN"] <- 0; coef_estim["MAX"] <- max(mat_tbl$inhibition,na.rm=T)
      #(Fix off minimums) Find lowest inhibition value. If it is not in (0:100), fix it whether to 0 or 99.
      min_lower <- ifelse(min(mat_tbl$inhibition,na.rm=T) > 0,min(mat_tbl$inhibition,na.rm=T),0)
      min_lower <- ifelse(min_lower >= 100,99,min_lower)
      #similar to previous step but for MAX
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]>100,100,coef_estim["MAX"])
      coef_estim["MAX"] <- ifelse(coef_estim["MAX"]<0,100,coef_estim["MAX"])
      #max_lower and max_upper - lower and upper bounds for 'nl2sol' algorithm in nonlinear least-squares
      max_lower <- ifelse(max(mat_tbl$inhibition,na.rm=T)>100,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,coef_estim["MAX"],max(mat_tbl$inhibition,na.rm=T))
      max_lower <- ifelse(max_lower < 0,0,max_lower)
      max_lower <- ifelse(max_lower > 100,100,max_lower)
      #(Fix upper maximum for negative slopes)
      run_avg <- caTools::runmean(mat_tbl$inhibition, 10)
      max_upper <- ifelse(any(run_avg[-nrow(mat_tbl)]>run_avg[nrow(mat_tbl)]),max(mat_tbl$inhibition[run_avg>run_avg[nrow(mat_tbl)]]),coef_estim["MAX"])
      max_upper <- ifelse(any(mat_tbl$inhibition > max_upper),mean(mat_tbl$inhibition[mat_tbl$inhibition > max_upper])+5,max_upper)
      max_upper <- ifelse(max_upper < 0,coef_estim["MAX"],max_upper)
      max_upper <- ifelse(max_upper > 100,100,max_upper) #coef_estim["MAX"]
      max_upper <- ifelse(max_lower > max_upper,coef_estim["MAX"],max_upper)
      # left it as it was, just rewritten a bit (ALEKS). not clear how values 25, 60 and 5 are chosen.
      mean_inh_last = mean(tail(mat_tbl$inhibition,2),na.rm=T)
      if(mean_inh_last < 60) {
        if(mean_inh_last > 25) coef_estim["IC50"] <- mean(mat_tbl$logconc,na.rm=T)
        else if(mean_inh_last < 25) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)}
      if(mean(mat_tbl$inhibition[1:3],na.rm=T)<5) coef_estim["IC50"] <- max(mat_tbl$logconc,na.rm=T)
      
      #add a bit of positive noise to MAX if it is the same as MIN.
      if(unname(coef_estim["MIN"]) == unname(coef_estim["MAX"])) coef_estim["MAX"] <- coef_estim["MAX"] + 0.001
      
      #adaptive nonlinear Least-Squares algorithm NL2SOL to estimate parameters.
      nls_result_ic50_old <- function(){
        tryCatch({
          nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)), upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {
          
          # allows higher residual sum-of-squares
          minpack.lm::nlsLM(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                            start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                            lower=c(SLOPE=0.1, MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),
                            upper=c(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)))
        })
      }
      
      # IC50 first
      nls_result_ic50 <- nls_result_ic50_old();
      
      # IC50 second
      nls_result_ic50_2 <- tryCatch({
        # allows higher residual sum-of-squares
        nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",  start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim["MAX"][[1]], IC50=median(mat_tbl$logconc)),lower=list(SLOPE=0.1,MIN=0,MAX=max_lower, IC50=min(mat_tbl$logconc)),upper=list(SLOPE=2.5,MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {
        nls_result_ic50_old()
      },error = function(e) {
        nls_result_ic50_old()
      })
      
      #element (4, 4) is zero, so the inverse cannot be computed
      nls_result_ic50 = tryCatch({summary(nls_result_ic50); nls_result_ic50},error=function(e){nls_result_ic50_2})
      
      #Calculate the standard error scores
      sumIC50 = list(summary(nls_result_ic50), summary(nls_result_ic50_2))
      
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2))
      nls_result_ic50 = list(nls_result_ic50, nls_result_ic50_2)[[switch_]]
      
      #if SLOPE <= 0.2, decrease IC50, change lower bound for SLOPE to 0.1 and repeat.
      if(coef(nls_result_ic50)["SLOPE"] <= 0.2)
      {
        if(mean_inh_last > 60)
          coef_estim["IC50"] <- min(mat_tbl$logconc,na.rm=T)
        nls_result_ic50 <- nls(inhibition ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",
                               start=list(SLOPE=1, MIN=unname(coef_estim["MIN"]),MAX=unname(coef_estim["MAX"]),IC50=unname(coef_estim["IC50"])),
                               lower=list(SLOPE=0.1,MIN=min_lower,MAX=max_lower,IC50=min(mat_tbl$logconc)),
                               upper=list(SLOPE=2.5, MIN=0,MAX=max_upper, IC50=max(mat_tbl$logconc)),control=list(warnOnly=T,minFactor = 1/2048))
      }
      
      ########
      
      
      ######## Add bi-directional fit - Aram N. Andersen 2023
      coef_estim1 <- coef(estimate_param); names(coef_estim1) <- c("SLOPE","MIN","MAX","IC50")
      coef_estim1["MAX"] <- mat_tbl$inhibition[which.max(abs(mat_tbl$inhibition))]
      coef_estim1["MAX"] <- ifelse(coef_estim1["MAX"]>100,100,coef_estim1["MAX"])
      coef_estim1["MAX"] <- ifelse(coef_estim1["MAX"]< -100,-100,coef_estim1["MAX"])
      
      nls_result_ic50_old.1 <- function(){
        tryCatch({
          nls(inhibition3 ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port", 
              start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim1["MAX"][[1]],IC50=coef_estim["IC50"][[1]]), 
              lower=list(SLOPE=0.1,MIN=0,MAX=-100, IC50=min(mat_tbl$logconc)), 
              upper=list(SLOPE=2.5,MIN=0,MAX=100, IC50=max(mat_tbl$logconc)),
              control=list(warnOnly=T,minFactor = 1/2048))
        }, error = function(e) {
          
          # allows higher residual sum-of-squares
          minpack.lm::nlsLM(inhibition3 ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl,
                            start=list(SLOPE=1, MIN=coef_estim["MIN"][[1]],MAX=coef_estim1["MAX"][[1]],IC50=coef_estim["IC50"][[1]]),
                            lower=c(SLOPE=0.1, MIN=0,MAX=-100, IC50=min(mat_tbl$logconc)),
                            upper=c(SLOPE=2.5, MIN=0,MAX=100, IC50=max(mat_tbl$logconc)))
        })
      }
      
      # IC50 first
      nls_result_ic50.1 <- nls_result_ic50_old.1();
      
      # IC50 second
      nls_result_ic50_2.1 <- tryCatch({
        # allows higher residual sum-of-squares
        nls(inhibition3 ~ MIN + (MAX - MIN)/ (1 + (10^(SLOPE * (IC50 - logconc)))), data=mat_tbl, algorithm="port",
            start=list(SLOPE=1,MIN=coef_estim["MIN"][[1]],MAX=coef_estim1["MAX"][[1]], IC50=median(mat_tbl$logconc)),
            lower=list(SLOPE=0.1,MIN=0,MAX=-100, IC50=min(mat_tbl$logconc)),
            upper=list(SLOPE=2.5,MIN=0,MAX=100, IC50=max(mat_tbl$logconc)),
            control=list(warnOnly=T,minFactor = 1/2048))
      },warning = function(w) {
        nls_result_ic50_old()
      },error = function(e) {
        nls_result_ic50_old()
      })
      
      #element (4, 4) is zero, so the inverse cannot be computed
      nls_result_ic50.1 = tryCatch({summary(nls_result_ic50.1); nls_result_ic50.1},error=function(e){nls_result_ic50_2.1})
      
      #Calculate the standard error scores
      sumIC50 = list(summary(nls_result_ic50.1), summary(nls_result_ic50_2.1))
      
      ic50std_resid <- round(sqrt(sum((sumIC50[[1]]$residuals)^2)/(length(sumIC50[[1]]$residuals)-1)),1);
      ic50std_resid2 <- round(sqrt(sum((sumIC50[[2]]$residuals)^2)/(length(sumIC50[[2]]$residuals)-1)),1);
      
      # continue with the best
      switch_ = which.min(c(ic50std_resid, ic50std_resid2))
      nls_result_ic50.1 = list(nls_result_ic50.1, nls_result_ic50_2.1)[[switch_]]
      
      sumIC50 = summary(nls_result_ic50.1)
      resid_free = as.numeric(sumIC50$residuals)
      coef_ec50 <- coef(nls_result_ic50.1)[c("IC50", "SLOPE","MAX","MIN")]; coef_ec50["IC50"] <- 10^coef_ec50["IC50"]
      
      ##########
      
      
      
      #Calculate the standard error scores
      sumIC50 = summary(nls_result_ic50);
      resid_constr = as.numeric(sumIC50$residuals)
      
      ic50std_Error <- sumIC50$coefficients["IC50","Std. Error"]; #tec50std_Error <- sumTEC50$coefficients["TEC50","Std. Error"]
      ic50std_resid <- round(sqrt(sum((sumIC50$residuals)^2)/(length(sumIC50$residuals)-1)),1);
      max_signal <- max(mat_tbl$dose,na.rm=T); min_signal <- min(mat_tbl$dose,na.rm=T)
      
      #############################
      #############   Final modification & STD error
      
      #prepare final data and convert IC50 back from log scale (inverse)
      coef_ic50 <- coef(nls_result_ic50)[c("IC50", "SLOPE","MAX","MIN")]; coef_ic50["IC50"] <- 10^coef_ic50["IC50"]
      #(Fix ic50 for curves in wrong direction)
      coef_ic50["IC50"] <- ifelse(coef_ic50["SLOPE"]<0,max_signal,coef_ic50["IC50"])
      #(Fix based on MAX)
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<0,max_signal,coef_ic50["IC50"])
      coef_ic50["IC50"] <- ifelse(coef_ic50["MAX"]<10,max_signal,coef_ic50["IC50"])
      coef_ic50["MAX"] <- ifelse(coef_ic50["MAX"]<0,0,coef_ic50["MAX"])
      #(Fix over sensitive drugs)
      coef_ic50["IC50"] <- ifelse(all(c(max(mat_tbl$inhibition,na.rm=T),min(mat_tbl$inhibition,na.rm=T))>50),min_signal,coef_ic50["IC50"])
      
      # for ploting
      x <- seq(min(mat_tbl$logconc),max(mat_tbl$logconc), length=100)
      yic <- predict(nls_result_ic50, data.frame(logconc=x))
      auc <- MESS::auc(x,yic)
      
      ##average replicates
      mat_tblCp <- mat_tbl[, c("inhibition", "dose")]
      cols_ <- colnames(mat_tblCp)[!grepl("inhibition", colnames(mat_tblCp))] # columns which should be equal to average PI
      X <- as.data.table(mat_tblCp)
      mat_tblCp <- as.data.frame(X[,list(inhibition = mean(inhibition)),cols_], stringAsFactors = !1)
      
      
      perInh <- t(matrix(mat_tblCp[,"inhibition"],dimnames=
                           list(paste0(rep("D", length(mat_tblCp[,"inhibition"])), 1:length(mat_tblCp[,"inhibition"])))))
      
      coef_tec50 = coef_ic50;
      coef_tec50["IC50"] <- ifelse(coef_tec50["MAX"] > 25, coef_tec50["IC50"], max(mat_tbl$dose,na.rm=T))
      if(readoutCTX){
        names(coef_tec50) <- c("TC50","SLOPE","MAX","MIN"); ytec <- yic; perViaTox <- perInh;
      } else{
        names(coef_tec50) <- c("EC50","SLOPE","MAX","MIN");
        coef_tec50["SLOPE"] = -1 * coef_tec50["SLOPE"]; # min - 0, max - 77 in ec50 it is max - 100, min - 23
        tmp = coef_tec50["MAX"]; coef_tec50["MAX"] = 100 - coef_tec50["MIN"]; coef_tec50["MIN"] = 100 - tmp; ytec <- 100 - yic;
        perViaTox <- 100 - perInh;
      }
      
      
      #############################
      #############    DSS
      
      #dss_score <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=as.integer(DSS_typ))),1);
      dss_score1 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=1)),1);
      dss_score2 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=2)),1);
      dss_score3 <- round(as.numeric(dss(coef_ic50["IC50"],coef_ic50["SLOPE"],coef_ic50["MAX"],min_signal,max_signal, DSS.type=3)),1);
      coef_ic50 <- c(coef_ic50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,IC50_std_error=ic50std_Error)
      coef_tec50 <- c(coef_tec50,Min.Conc.tested=min_signal,Max.Conc.tested=max_signal,TEC50_std_error=ic50std_Error)
      
      ####
      Hill_fun = function(x,parm){parm[3]*x^parm[2]/(x^parm[2]+parm[1]^parm[2])}
      y_pred.breeze = Hill_fun(mat_tbl$dose, coef_ic50)
      
      
      y_true <- mat_tbl$inhibition2
      y_fit = predict(nls_result_ic50.1, mat_tbl$dose)
      #y_fit.adj = predict(nls_result_ic50, mat_tbl$dose)
      y_fit.adj <- ifelse(y_fit.adj<0, 0, y_fit)
      
      y_pred <- data.frame(y_true=y_true,
                           y_fit=y_fit,
                           y_fit.adj=y_fit.adj,
                           y_fit.breeze=y_pred.breeze,
                           pred_err=y_true-y_fit, 
                           pred_err.adj=y_true-y_fit.adj,
                           pred_err.breeze=y_true-y_pred.breeze,
                           dose=mat_tbl$dose)
      
      
      list(DSS1=dss_score1,
           DSS2=dss_score2,
           DSS3=dss_score3,
           coef_ic50 = coef_ic50,
           coef_tec50 = coef_tec50,
           coef_ec50 = coef_ec50,
           y_pred=y_pred)
    }
  }, error = function(e) {
    print(paste0("error in ", product_id, ", in ", xpr_tbl$screen_id[idx_filt][[1]]));
    print(e);
  })
})

# Compute raw AUC
df.rAUC <- df_dose_responses %>%
  group_by(ID,Patient.ID, Patient.num, drug, platenumber) %>%
  summarise(rAUC = log10.AUC(concentration, response)) %>%
  as.data.frame()
df.rAUC$rAUC_log2 <- -log2(df.rAUC$rAUC + min(df.rAUC$rAUC[df.rAUC$rAUC != 0]))

# Compute BREEZE drug sensitivity metrics
set.seed(1)
xpr_tbl <- df_dose_responses
xpr_tbl$Concentration <- xpr_tbl$concentration
drug_wells_ <- df.rAUC
iter <- nrow(drug_wells_)

list.pred_Hill <- list()
list.summary_Hill <- list()
for(it in 1:iter){
  
  df.DSS.ctx <- CALC_IC50_EC50_DSS(it, drug_wells_, xpr_tbl, DSS_typ="1", readoutCTX = T, path)
  
  list.pred_Hill[[it]] <- df.DSS.ctx$y_pred
  list.summary_Hill[[it]] <- df.DSS.ctx
  
  rm(df.DSS.ctx)
  
  svMisc::progress(it/iter*100)
}
rm(xpr_tbl,iter, drug_wells_)

for(it in 1:length(list.pred_Hill)){
  list.pred_Hill[[it]]$y_fit.adj <- ifelse(list.pred_Hill[[it]]$y_fit<0, 0, list.pred_Hill[[it]]$y_fit)
  list.pred_Hill[[it]]$pred_err.adj <- list.pred_Hill[[it]]$y_true - list.pred_Hill[[it]]$y_fit.adj
  list.summary_Hill[[it]]$y_pred$y_fit.adj <- ifelse(list.summary_Hill[[it]]$y_pred$y_fit<0, 0, list.summary_Hill[[it]]$y_pred$y_fit)
  list.summary_Hill[[it]]$y_pred$pred_err.adj <- list.summary_Hill[[it]]$y_pred$y_true - list.summary_Hill[[it]]$y_pred$y_fit.adj
}


df_pred_Hill <- bind_rows(list.pred_Hill, .id = "ID")


df_pred_Hill$Patient.ID <- df.rAUC$Patient.ID[match(df_pred_Hill$ID, 1:nrow(df.rAUC))]
df_pred_Hill$Patient.num <- df.rAUC$Patient.num[match(df_pred_Hill$ID, 1:nrow(df.rAUC))]
df_pred_Hill$drug <- df.rAUC$drug[match(df_pred_Hill$ID, 1:nrow(df.rAUC))]


ggarrange(ggplot(df_pred_Hill)+
            geom_point(aes(y_true, y_fit, col=factor(dose)), size=0.5, pch=16)+
            scale_color_viridis_d()+
            labs(x="Inhibition (%)", y="Inhibition curvefit (%)", col="Dose")+
            ylim(-100,100)+
            xlim(-100,100)+
            theme_classic(),
          ggplot(df_pred_Hill)+
            geom_point(aes(y_true, y_fit.adj, col=factor(dose)), size=0.5, pch=16)+
            scale_color_viridis_d()+
            labs(x="Inhibition (%)", y="Inhibition curvefit - adjusted (%)")+
            ylim(0,100)+
            xlim(-100,100)+
            theme_classic(),
          ggplot(df_pred_Hill)+
            geom_point(aes(y_true, y_fit.breeze, col=factor(dose)), size=0.5, pch=16)+
            scale_color_viridis_d()+
            labs(x="Inhibition (%)", y="Inhibition curvefit - Breeze (%)")+
            ylim(0,100)+
            xlim(-100,100)+
            theme_classic(), common.legend = T)
ggsave(width=6, height=5.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/Hill curve fit.tiff"))


df.hillAUC <- df_pred_Hill %>%
  group_by(ID,Patient.ID, Patient.num, drug) %>%
  summarise(rAUC = log10.AUC(dose, 1-y_true/100),
            hill_AUC = log10.AUC(dose, 1-y_fit/100),
            hill0_AUC = log10.AUC(dose, 1-y_fit.adj/100),
            hillB_AUC = log10.AUC(dose, 1-y_fit.breeze/100)) %>%
  as.data.frame()
df.hillAUC$rAUC_log2 <- -log2(df.hillAUC$rAUC + min(df.hillAUC$rAUC[df.hillAUC$rAUC != 0]))
df.hillAUC$hill_AUC_log2 <- -log2(df.hillAUC$hill_AUC + min(df.hillAUC$hill_AUC[df.hillAUC$hill_AUC != 0]))
df.hillAUC$hill0_AUC_log2 <- -log2(df.hillAUC$hill0_AUC + min(df.hillAUC$hill0_AUC[df.hillAUC$hill0_AUC != 0]))
df.hillAUC$hillB_AUC_log2 <- -log2(df.hillAUC$hillB_AUC + min(df.hillAUC$hillB_AUC[df.hillAUC$hillB_AUC != 0]))

save(df_pred_Hill, list.pred_Hill, df.hillAUC,
     file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))



ggarrange(ggplot(df.hillAUC, aes(rAUC, hill_AUC))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC", y="Hill AUC")+
            theme_classic(),
          ggplot(df.hillAUC, aes(rAUC, hill0_AUC))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC", y="Hill AUC - adjusted")+
            theme_classic(),
          ggplot(df.hillAUC, aes(rAUC, hillB_AUC))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC", y="Hill AUC - Breeze")+
            theme_classic(),
          ggplot(df.hillAUC, aes(rAUC_log2, hill_AUC_log2))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC-log2", y="Hill AUC-log2")+
            theme_classic(),
          ggplot(df.hillAUC, aes(rAUC_log2, hill0_AUC_log2))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC-log2", y="Hill AUC-log2 - adjusted")+
            theme_classic(),
          ggplot(df.hillAUC, aes(rAUC_log2, hillB_AUC_log2))+
            geom_point(size=0.5, pch=16, col="darkgray", alpha=0.7)+
            stat_cor()+
            scale_color_viridis_d()+
            labs(x="rAUC-log2", y="Hill AUC-log2 - Breeze")+
            theme_classic(),common.legend = T, ncol=3, nrow=2)
ggsave(width=7.5, height=4.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/Hill AUC comparison.tiff"))
ggsave(width=7.5, height=4.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/Hill AUC comparison.pdf"))


# Analysis of confounding components ----

load(file = paste0(path,"/Drug sensitivity curvefit predictions and AUC_2023-07.RData"))

mat_rAUC <- reshape2::dcast(df.hillAUC %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = "rAUC")
rownames(mat_rAUC) <- as.character(mat_rAUC$Patient.num); mat_rAUC$Patient.num <- NULL
mat_rAUC <- mat_rAUC[match(rownames(Y), rownames(mat_rAUC)),]
SDs <- matrixStats::colSds(data.matrix(mat_rAUC))


X_biological <- cbind(df_genetics, df_prognostics_dummy[,grep("Primary|Age|Sex|FAB|ELN",colnames(df_prognostics_dummy))])


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


df.hillAUC$rAUC <- 1-df.hillAUC$rAUC
df.hillAUC$hill_AUC <- 1-df.hillAUC$hill_AUC
df.hillAUC$hill0_AUC <- 1-df.hillAUC$hill0_AUC
df.hillAUC$hillB_AUC <- 1-df.hillAUC$hillB_AUC
# df.hillAUC$r0AUC <- ifelse(df.hillAUC$rAUC <0,0,df.hillAUC$rAUC)
# df.hillAUC$r0AUC_log2 <- ifelse(df.hillAUC$rAUC_log2 <0,0,df.hillAUC$rAUC_log2)

df_r2 <- c()
df_PC <- c()
for(t in c(quantile(SDs, 1-c(100)/ncol(mat_rAUC)),0)){
  for(pt.st in c(T, F)){
    for(metric in colnames(df.hillAUC)[5:12]){
      
      X <- reshape2::dcast(df.hillAUC %>% subset(!grepl("re", Patient.ID)), Patient.num~drug ,value.var = metric)
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

df_r2$Percent_of_variance <- df_PC$Percent_of_variance[match(paste(df_r2$Metric,df_r2$Pre_selection_N, df_r2$Patient_stdz,as.character(df_r2$PC)),
                                                             paste(df_PC$Metric,df_PC$Pre_selection_N, df_PC$Patient_stdz,as.character(df_PC$PC)))]
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

df_r2$Metric3 <- gsub("_log2", "",df_r2$Metric)
df_r2$Metric3 <- factor(df_r2$Metric3, levels=unique(df_r2$Metric3))
df_r2$log2 <- ifelse(grepl("log", df_r2$Metric), "log2", "")

df_PC$Metric3 <- gsub("_log2", "",df_PC$Metric)
df_PC$Metric3 <- factor(df_PC$Metric3, levels=unique(df_PC$Metric3))
df_PC$log2 <- ifelse(grepl("log", df_PC$Metric), "log2", "")

df_PC$Patient_stdz <- ifelse(df_PC$Patient_stdz == "","", "z-score")
df_r2$Patient_stdz <- ifelse(df_r2$Patient_stdz == "","", "z-score")


ggplot(df_r2 %>% subset(grepl("Batch|Biol",Variable) & Metric3!="rAUC"))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric3, lty=Patient_stdz), alpha=1)+
  theme_bw()+
  scale_color_manual(values=pal_jama()(7))+
  facet_grid(Variable~log2+Pre_selection_N, scales="free")+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=5.75, height=2.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_cumulative variance explained_1.pdf"))

ggplot(df_r2 %>% subset(grepl("Batch|Biol",Variable) & Metric3!="rAUC"))+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric3), alpha=1)+
  geom_line(data=df_r2 %>% subset(grepl("Batch|Biol",Variable) & Metric3=="rAUC"),
            aes(y=Cumulative_variance_explained, x=as.numeric(PC)), alpha=1, lty=2)+
  theme_bw()+
  scale_color_manual(values=pal_jama()(7))+
  facet_grid(Variable+Pre_selection_N~Patient_stdz+log2, scales="free")+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=5.75, height=4, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_cumulative variance explained_2.pdf"))


ggplot(df_PC %>% subset(Pre_selection_N == 349 & as.numeric(PC) %in% 1:10 & Metric3!="rAUC"))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  geom_point(data=df_PC %>% subset(Pre_selection_N == 349 & as.numeric(PC) %in% 1:10 & Metric3=="rAUC"),
             aes(y=Percent_of_variance*100, x=as.numeric(PC)), 
             alpha=1)+
  theme_bw()+
  facet_grid(Pre_selection_N~Patient_stdz+log2, scales = "free")+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(width=9, height=2.5, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_PC variance_349_1.pdf"))


ggplot(df_PC %>% subset(Pre_selection_N != 349 & as.numeric(PC) %in% 1:10 & Metric3!="rAUC"))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  geom_point(data=df_PC %>% subset(Pre_selection_N != 349 & as.numeric(PC) %in% 1:10 & Metric3=="rAUC"),
             aes(y=Percent_of_variance*100, x=as.numeric(PC)), 
             alpha=1)+
  theme_bw()+
  facet_grid(Pre_selection_N~Patient_stdz+log2, scales = "free")+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  labs(x="Principal components", y = "Cumulative variance explained")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())

ggsave(width=9, height=2.5, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_PC variance_100_1.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N == 349 & Metric3!="rAUC") %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  geom_point(data=df_r2 %>% subset(Pre_selection_N == 349 & as.numeric(PC) %in% 1:10 & Metric3=="rAUC"),
             aes(y=adj.r.squared, x=as.numeric(PC)), 
             alpha=1)+
  theme_bw()+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  facet_grid(Variable~Patient_stdz+log2, scales = "free")+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=9, height=5.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_PC variance explained_349_1.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N != 349 & Metric3!="rAUC") %>%
         subset(as.numeric(PC) %in% 1:10))+
  ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  geom_point(data=df_r2 %>% subset(Pre_selection_N != 349 & as.numeric(PC) %in% 1:10 & Metric3=="rAUC"),
             aes(y=adj.r.squared, x=as.numeric(PC)), 
             alpha=1)+
  theme_bw()+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  facet_grid(Variable~Patient_stdz+log2, scales = "free")+
  labs(x="", y = "adjusted r-squared")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=9, height=5.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_PC variance explained_100_1.pdf"))


df_r2 <- df_r2 %>% group_by(Patient_stdz, log2, Variable,Pre_selection_N, PC) %>% mutate(adj.r.squared_ref = mean(ifelse(Metric3=="rAUC", adj.r.squared, NA), na.rm=T),
                                                                                         r.squared_ref = mean(ifelse(Metric3=="rAUC", r.squared, NA), na.rm=T),
                                                                                         Variance_explained_ref = mean(ifelse(Metric3=="rAUC", Variance_explained, NA), na.rm=T))

ggplot(df_r2 %>% subset(Pre_selection_N == 349 & Metric3!="rAUC") %>%
         subset(as.numeric(PC) %in% 1:10))+
  #ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared-adj.r.squared_ref, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  facet_grid(Variable~Patient_stdz+log2, scales = "free")+
  labs(x="", y = "adjusted r-squared change")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=9, height=5.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_diff. PC variance explained_349_1.pdf"))


ggplot(df_r2 %>% subset(Pre_selection_N != 349 & Metric3!="rAUC") %>%
         subset(as.numeric(PC) %in% 1:10))+
  #ylim(0,1)+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  geom_bar(aes(y=adj.r.squared-adj.r.squared_ref, x=as.numeric(PC), fill=Metric3), alpha=1, stat="identity", position = position_dodge(width=0.85))+
  theme_bw()+
  scale_fill_manual(values=pal_jama()(7)[c(1:7)])+
  facet_grid(Variable~Patient_stdz+log2, scales = "free")+
  labs(x="", y = "adjusted r-squared change")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(width=9, height=5.75, 
       filename = paste0(path, "/Results_Hill fit and PCA/PCA of confounding factors_Hill fit_diff. PC variance explained_100_1.pdf"))


# PCA analysis of DR data directly
mat_y.pt <- reshape2::dcast(df_pred_Hill %>% subset(!grepl("re", Patient.ID)), Patient.num~drug+dose, value.var = "y_true")
rownames(mat_y.pt) <- mat_y.pt$Patient.num; mat_y.pt$Patient.num <- NULL
mat_error.pt <- reshape2::dcast(df_pred_Hill%>% subset(!grepl("re", Patient.ID)), Patient.num~drug+dose, value.var = "pred_err")
rownames(mat_error.pt) <- mat_error.pt$Patient.num; mat_error.pt$Patient.num <- NULL
mat_y_fit.pt <- reshape2::dcast(df_pred_Hill%>% subset(!grepl("re", Patient.ID)), Patient.num~drug+dose, value.var = "y_fit")
rownames(mat_y_fit.pt) <- mat_y_fit.pt$Patient.num; mat_y_fit.pt$Patient.num <- NULL
mat_y_fit.adj.pt <- reshape2::dcast(df_pred_Hill%>% subset(!grepl("re", Patient.ID)), Patient.num~drug+dose, value.var = "y_fit.adj")
rownames(mat_y_fit.adj.pt) <- mat_y_fit.adj.pt$Patient.num; mat_y_fit.adj.pt$Patient.num <- NULL
mat_y_fit.breeze.pt <- reshape2::dcast(df_pred_Hill%>% subset(!grepl("re", Patient.ID)), Patient.num~drug+dose, value.var = "y_fit.breeze")
rownames(mat_y_fit.breeze.pt) <- mat_y_fit.breeze.pt$Patient.num; mat_y_fit.breeze.pt$Patient.num <- NULL

df_r2 <- c()
df_PC <- c()
for(st in c(T, F)){
  for(metric in c("y","y_Hill","y_Hill0","y_HillB","error")){
    
    if(metric=="y"){
      X <- data.matrix(mat_y.pt)
    }else if(metric=="y_Hill"){
      X <- data.matrix(mat_y_fit.pt)
    }else if(metric=="y_Hill0"){
      X <- data.matrix(mat_y_fit.adj.pt)
    }else if(metric=="y_HillB"){
      X <- data.matrix(mat_y_fit.breeze.pt)
    }else{
      X <- data.matrix(mat_error.pt)
    }
    
    if(st){
      X <- (X - rowMeans(X))/matrixStats::rowSds(X)
    }
    
    X <- t((t(X) - colMeans(X))/matrixStats::colSds(X))
    
    
    if(length(which(is.na(X)))){X <- replace(X,is.na(X),0)}
    svdz <- svd(t(X))
    df_PCA <- data.frame(svdz$v[,1:50]); colnames(df_PCA) <- paste0("PC", 1:50); rownames(df_PCA) <- rownames(mat_y.pt)
    
    m <- data.frame(Percent_of_variance=((svdz$d^2)/sum(svdz$d^2))[1:50], PC = colnames(df_PCA))
    m$Metric <- metric
    m$Patient_stdz <- ifelse(st,"Standardized", "")
    m$Pre_selection_N <- ncol(X)
    df_PC <- rbind(df_PC, m)
    
    for(i in colnames(df_PCA)[1:50]){
      
      
      dat <- cbind(df_PCA[,i],df_batch[match(rownames(df_PCA), rownames(df_batch)),])
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Batch covariates")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(st,"Standardized", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
      
      dat <- cbind(df_PCA[,i],
                   X_biological[match(rownames(df_PCA), rownames(X_biological)),])
      
      colnames(dat)[1] <- "PC"
      lm.fit <- summary(lm(PC~.,dat))
      m <- data.frame(r.squared=lm.fit$r.squared, adj.r.squared=lm.fit$adj.r.squared,rmse=sqrt(mean(lm.fit$residuals^2)), PC=i, Variable="Biological/clinical covariates")
      m$Metric <- metric
      m$Patient_stdz <- ifelse(st,"Standardized", "")
      m$Pre_selection_N <- ncol(X)
      df_r2 <- rbind(df_r2, m); rm(lm.fit, m, dat)
      
    }
  }
}

df_r2$PC <- factor(df_r2$PC, levels = unique(df_r2$PC))
df_PC$PC <- factor(df_PC$PC, levels = unique(df_PC$PC))

df_r2$Percent_of_variance <- df_PC$Percent_of_variance[match(paste(df_r2$Metric, df_r2$Patient_stdz,as.character(df_r2$PC)),
                                                             paste(df_PC$Metric, df_PC$Patient_stdz,as.character(df_PC$PC)))]
df_r2 <- df_r2 %>%  mutate(Variance_explained=adj.r.squared*Percent_of_variance)
df_r2 <- df_r2 %>% group_by(Variable,Metric,Patient_stdz) %>% mutate(Cumulative_variance_explained=cumsum(Variance_explained))
df_r2$Metric <- factor(df_r2$Metric, levels=unique(df_r2$Metric))
df_PC$Metric <- factor(df_PC$Metric, levels=unique(df_PC$Metric))






ggplot(df_r2)+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric, lty=ifelse(Metric=="y", "y", "")), alpha=1)+
  theme_bw()+
  scale_color_manual(values=c("black",pal_jama()(7)[c(1:7)]))+
  scale_linetype_manual(values = c(1,2))+
  facet_grid(Patient_stdz~Variable)+
  labs(x="Principal components", y = "Cumulative variance explained", lty="Reference")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(filename = paste0(path,"/Results_Hill fit and PCA/Curve fit data_cumulative variance explained_1.pdf"), 
       width = 4.75, height=3)

ggplot(df_r2)+
  geom_line(aes(y=Cumulative_variance_explained, x=as.numeric(PC), col=Metric, lty=ifelse(Metric=="y", "y", "")), alpha=1)+
  theme_bw()+
  scale_color_manual(values=c("black",pal_jama()(7)[c(1:7)]))+
  scale_linetype_manual(values = c(1,2))+
  facet_grid(Variable~Patient_stdz, scales = "free")+
  labs(x="Principal components", y = "Cumulative variance explained", lty="Reference")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())
ggsave(filename = paste0(path,"/Results_Hill fit and PCA/Curve fit data_cumulative variance explained_2.pdf"), 
       width = 4.75, height=3)

ggplot(df_PC %>% subset(as.numeric(PC) %in% 1:10))+
  geom_bar(aes(y=Percent_of_variance*100, x=as.numeric(PC), fill=Metric), 
           alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  facet_grid(~Patient_stdz)+
  scale_fill_manual(values=c("gray",pal_jama()(7)[c(1:7)]))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="Principal components", y = "Variance explained (%)", lty="Reference")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_Hill fit and PCA/Curve fit data_PC variance.pdf"), 
       width = 7, height=2.25)

ggplot(df_r2 %>% subset(as.numeric(PC) %in% 1:10))+
  geom_bar(aes(y=adj.r.squared, x=as.numeric(PC), fill=Metric), alpha=1, stat="identity", position = position_dodge(width=0.75))+
  theme_bw()+
  facet_grid(Variable~Patient_stdz)+
  ylim(0,1)+
  scale_fill_manual(values=c("gray",pal_jama()(7)[c(1:7)]))+
  scale_x_continuous(breaks=c(2,4,6,8,10))+
  labs(x="Principal components", y = "Adjusted r-squared", lty="Reference")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_blank())
ggsave(filename = paste0(path,"/Results_Hill fit and PCA/Curve fit data_PC variance explained.pdf"), 
       width = 7, height=3.25)



