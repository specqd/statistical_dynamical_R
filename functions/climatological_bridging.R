climatological_bridging <-function(obs,bridging_pred){
  
  if (length(obs)!=length(bridging_pred)){
    stop("Predictor and predictand must have the same length.")
  }
  
  fit0=lm(bridging_pred~obs)
  a0=fit0[["coefficients"]][1]
  b0=fit0[["coefficients"]][2]
  std0=summary(fit0)[["sigma"]]
  
  return(list(a0,b0,std0))
  
}