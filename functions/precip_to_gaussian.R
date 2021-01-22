precip_to_gaussian<-function(precip,precip_clim){
  
  require(EnvStats)
  
  precip_p=pemp(precip,precip_clim,discrete=F)
  normalized_precip=qnorm(precip_p,mean=0,sd=1)
  
  return(normalized_precip)
  
}