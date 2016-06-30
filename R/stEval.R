#' Basic evaluation of fit
#'
#'
#'
#' @export
#' 
#' 
#' 


stEval = function(forecast, Y) {
  
  nt = ncol(Y)
  
  for(t in 1:nt) {
    fcst = forecast[[t]]
    
    fcst$err = list()
    
    fcst$pred$resid = Y[,t] - fcst$pred$y
    fcst$err$mspe = mean(fcst$pred$resid^2)
    fcst$err$ppl = sum(fcst$pred$resid^2)/2 + sum(fcst$pred$se^2)
    
    forecast[[t]] = fcst
  }
  
  forecast
}