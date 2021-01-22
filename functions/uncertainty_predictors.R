uncertainty_predictors <-function(obs_pred,fcst_pred){
  
  if (length(obs_pred)!=length(fcst_pred)){
    stop("Observed and forecast values must have the same length.")
  }
  
  fitf=lm(fcst_pred~obs_pred)
  af=fitf[["coefficients"]][1]
  bf=fitf[["coefficients"]][2]
  stdf=summary(fitf)[["sigma"]]
  
  return(list(af,bf,stdf))
  
}