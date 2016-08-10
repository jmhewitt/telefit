#' Basic evaluation of fit
#'
#' Provides basic measures for evalutating the fit.  Includes Brier skill score
#' against the climatology, MSPE, PPL, overall correlation, and a computation
#' of the coverage probabilities for confidence intervals
#'
#' @param clim the climatology for the location in Y
#'
#'
#' @export
#'
#' 


stEval = function(forecast, Y, clim) {
  
  nt = ncol(Y)
  if(is.null(nt)) {
    nt = 1
    Y = matrix(Y, ncol=1)
  }
  
  if(class(forecast)=='stPredict')
    forecast = list(forecast)
  
  clim.resid = c()
  
  for(t in 1:nt) {
    fcst = forecast[[t]]
    
    fcst$err = list()
    
    fcst$pred$resid = Y[,t] - fcst$pred$Y
    fcst$pred$covered = ifelse( Y[,t]>=fcst$pred$Y.lwr,
                                  ifelse(Y[,t]<=fcst$pred$Y.upr, T, F),
                                  F )
    
    # standard fit measurements
    fcst$err$mspe = mean(fcst$pred$resid^2)
    fcst$err$ppl = sum(fcst$pred$resid^2)/2 + sum(fcst$pred$se^2)
    fcst$err$r2 = 1 - var(fcst$pred$resid)/var(Y[,t])
    fcst$err$cor = cor(fcst$pred$Y, Y[,t])
    fcst$err$coverage = mean(fcst$pred$covered, na.rm = T)
    
    # brier skill score if climatology was provided for reference forecast
    if(!is.null(clim)) {
      resid = Y[,t] - clim
      clim.resid = c(clim.resid, resid)
      mspe.ref = mean(resid^2)
      fcst$err$bss = 1 - fcst$err$mspe / mspe.ref
    }
    
    forecast[[t]] = fcst
  }
 
  # collect all predictions and residuals
  Y = as.numeric(Y)
  Y.hat = as.numeric(sapply(forecast, function(f) { f$pred$Y }))
  resid = as.numeric(sapply(forecast, function(f) { f$pred$resid }))
  se = as.numeric(sapply(forecast, function(f) { f$pred$se }))
  coverages = as.numeric(sapply(forecast, function(f) { f$pred$covered }))
  
  # add overall error evaluations
  attr(forecast, 'err.mspe') = mean(resid^2)
  attr(forecast, 'err.ppl') = sum(resid^2)/2 + sum(se^2)
  attr(forecast, 'err.r2') = 1 - var(resid)/var(Y)
  attr(forecast, 'err.cor') = cor(Y.hat, Y)
  attr(forecast, 'err.coverage') = mean(coverages, na.rm = T)
  attr(forecast, 'err.bss') = 1 - mean(resid^2) / mean(clim.resid^2)
  
  
  if(nt==1)
    forecast[[1]]
  else
    forecast
}