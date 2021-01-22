gaussian_to_precip<-function(gaussian,precip_clim){
  
  require(EnvStats)
  
  precip_p=pnorm(gaussian,mean=0,sd=1)
  precip=qemp(precip_p,precip_clim,discrete=F)
  
  return(precip)
  
}