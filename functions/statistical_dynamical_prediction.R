statistical_dynamical_prediction <- function(calib_pred=NULL,bridging_preds=NULL,a0=NULL,b0=NULL,std0=NULL,af,bf,stdf){

	### calib_pred : the forecast predictand value (length 1)
	### bridging_preds : the forecast bridging predictors (e.g ENSO or MJO indices, length N)
	### a0 : intercepts of the climatological linear regressions of bridging predictors against the predictand (length N), computed by function climatological_bridging
	### b0 : slopes of the climatological linear regressions of bridging predictors against the predictand (length N), computed by function climatological_bridging
	### std0 : standard deviations of the climatological linear regressions of bridging predictors against the predictand (length N), computed by function climatological_bridging
	### af : intercepts of the linear regressions of all observed predictors against all forecast predictors (length N+1), computed by function uncertainty_predictors
	### bf : slopes of the linear regressions of all observed predictors against all forecast predictors (length N+1), computed by function uncertainty_predictors
	### stdf : standard deviations of the linear regressions of all observed predictors against all forecast predictors (length N+1), computed by function uncertainty_predictors

	### (a0, b0, std0) and (af, bf, stdf) are determined on the reforecast training period

	if (is.null(calib_pred) & is.null(bridging_preds)){
		stop("No predictor provided")
	}

	### Putting all predictors together (length N+1)

	preds=c(calib_pred,bridging_preds)
	
	### Adding the "climatological" coefficients for the calibration predictor (i.e predicting an observed variable with itself)
	
	if(!is.null(calib_pred)){
		a0=c(0,a0)
		b0=c(1,b0)
		std0=c(0,std0)
	}
	
	### Checking data consistency
	
	if (length(a0)!=length(preds) | length(b0)!=length(preds) | length(std0)!=length(preds)){
		stop("The number of climatological coefficients must be the same as the number of predictors")
	}
	
	if (length(a0)!=length(preds) | length(b0)!=length(preds) | length(std0)!=length(preds)){
		stop("The number of climatological coefficients must be the same as the number of predictors")
	}
	
	### Computing the resulting coefficients
	
	al=af+bf*a0
  bl=bf*b0
	stdl=sqrt(bf^2*std0^2+stdf^2)
	
	###
	
	mean_prior=0
	std_prior=1
	
	inv_var=1/std_prior^2+sum(bl^2/stdl^2)
	statdyn_fcst=(1/inv_var)*((mean_prior/std_prior^2)+sum((bl/stdl)^2*(preds-al)/bl))
	std=sqrt(1/inv_var)
	
	return(list(statdyn_fcst,std))

}