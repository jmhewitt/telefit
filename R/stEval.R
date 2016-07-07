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
      mspe.ref = mean((Y[,t] - clim)^2)
      fcst$err$bss = 1 - fcst$err$mspe / mspe.ref
    }
    
    
    
    forecast[[t]] = fcst
  }
 
  forecast
}