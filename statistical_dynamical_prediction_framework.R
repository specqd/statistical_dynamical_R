### FRAMEWORK FOR STATISTICAL-DYNAMICAL PREDICTION

### The following parameters are assumed fixed : longitude,latitude,lead time and start date
### A set of reforecast start dates (rfcst_start_dates) is considered for training
### All precipitation values are weekly averages over the selected lead time

### The following data is assumed available :

### Forecast data :
### pr_mod_test : vector of length nmb1 (ensemble size with nmb1 members)
### bridging_preds_mod_test : array of dimension nmb1 x N for N bridging predictors
### pr_mod_train : matrix of dimension nmb2 x length(rfcst_start_dates) (reforecasts with ensemble size nmb2)
### bridging_preds_mod_train : array of dimension nmb2 x length(rfcst_start_dates) x N

### Observation data :
### pr_obs_train : vector of length(rfcst_start_dates)
### clim_obs_train : matrix of dimension 3 x length(rfcst_start_dates), adding the previous and next week to the observed climatology
### bridging_preds_obs_train : array of dimension length(rfcst_start_dates) x N for N bridging predictors


###

library(EnvStats)

source("functions/precip_to_gaussian.R")
source("functions/gaussian_to_precip.R")
source("functions/climatological_bridging.R")
source("functions/uncertainty_predictors.R")
source("statistical_dynamical_prediction.R")

### Step 1 : Transforming forecast and observed precipitation into a gaussian variable

dim(pr_mod_train)=nmb2*length(rfcst_start_dates)
pr_mod_test=precip_to_gaussian(pr_mod_test,as.numeric(pr_mod_train))
pr_mod_train=precip_to_gaussian(as.numeric(pr_mod_train),as.numeric(pr_mod_train))
dim(pr_mod_train)=c(nmb2,length(rfcst_start_dates))

dim(clim_obs_train)=3*length(rfcst_start_dates)
pr_obs_train=precip_to_gaussian(pr_obs_train,as.numeric(clim_obs_train))

### Step 2 : Taking the ensemble mean of forecast data

pr_mod_test=mean(pr_mod_test)
pr_mod_train=apply(pr_mod_train,2,mean)
if (N>0){
  bridging_preds_mod_test=apply(bridging_preds_mod_test,2,mean)
  bridging_preds_mod_train=apply(bridging_preds_mod_train,2:3,mean)
}

### Step 3 : Learning the climatological bridging relationships


if (N>0){
  
  a0=c()
  b0=c()
  std0=c()
  for (i in 1:N){
    
    tmp=climatological_bridging(pr_obs_train,bridging_preds_obs_train[,i])
    a0=c(a0,tmp[[1]])
    b0=c(b0,tmp[[2]])
    std0=c(std0,tmp[[3]])
    rm(tmp)
    
  }
  
}

### Step 4 : Learning the uncertainty about predictors

tmp=uncertainty_predictors(pr_obs_train,pr_mod_train) ### Uncertainty in precipitation forecasts
af=tmp[[1]]
bf=tmp[[2]]
stdf=tmp[[3]]
rm(tmp)

if (N>0){
  
  for (i in 1:N){
    
    tmp=uncertainty_predictors(bridging_preds_obs_train[,i],bridging_preds_mod_train[,i]) ### Uncertainty in bridging predictors forecasts
    af=c(af,tmp[[1]])
    bf=c(bf,tmp[[2]])
    stdf=c(stdf,tmp[[3]])
    rm(tmp)
    
  }
  
}

### Step 5 : Performing the statistical-dynamical prediction

tmp=statistical_dynamical_prediction(pr_mod_test,bridging_preds_mod_test,a0,b0,std0,af,bf,stdf)
pr_mod_statdyn=tmp[[1]] ### pr_mod_statdyn is also a gaussian variable
std_statdyn=tmp[[2]]
rm(tmp)

### Application 1 : Computing the probability of an event (e.g exceedance of the 80th percentile)

threshold=qnorm(0.8,mean=0,sd=1)
prob_statdyn=pnorm(threshold,mean=pr_mod_statdyn,sd=std_statdyn,lower.tail=F)

### Application 2 : Transforming back into a precipitation value

pr_mod_statdyn=gaussian_to_precip(pr_mod_statdyn,clim_obs_train)