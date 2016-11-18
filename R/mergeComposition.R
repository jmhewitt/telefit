#' Combine results from composition sampler
#' 
#' @param xfull Raw output from one run of the Rcpp/Armadillo composition sampler
#' @param yfull Raw output from another run of the Rcpp/Armadillo composition sampler
#'  
#' 

mergeComposition = function(xfull, yfull) {
  if(is.null(yfull)) {
    xfull
  } else {
    
    # merge forecasts
    if(!is.null(xfull$forecast)) {
      # extract forecast parts
      x = xfull$forecast
      y = yfull$forecast
      
      # merge composition samples
      x$forecast = abind3(x$forecast, y$forecast)
      x$local = abind3(x$local, y$local)
      x$remote = abind3(x$remote, y$remote)
      
      xfull$forecast = x
    }
    
    
    # merge teleconnection effects
    if(!is.null(xfull$alpha_knots)) {
      # extract alpha parts
      x = xfull$alpha_knots
      y = yfull$alpha_knots
      
      z = list(
        est = mergeMean(x$est, y$est, x$nSamples, y$nSamples),
        sd = sqrt(mergeVar(x$sd^2, y$sd^2, x$est, y$est, x$nSamples, y$nSamples)),
        nSamples = x$nSamples + y$nSamples,
        samples = rbind(x$samples, y$samples)
      )
      
      if(!is.null(x$cov)) {
        z$cov = mergeCovmat(x$cov, y$cov, x$est, x$est, y$est, y$est, x$nSamples, y$nSamples)
      }
      
      xfull$alpha_knots = z
    }
    
    # merge full-field teleconnection effects
    if(!is.null(xfull$alpha)) {
      # extract alpha parts
      x = xfull$alpha
      y = yfull$alpha
      
      z = list(
        est = mergeMean(x$est, y$est, x$nSamples, y$nSamples),
        sd = sqrt(mergeVar(x$sd^2, y$sd^2, x$est, y$est, x$nSamples, y$nSamples)),
        nSamples = x$nSamples + y$nSamples,
        samples = rbind(x$samples, y$samples)
      )
      
      xfull$alpha = z
    }
    
    xfull
  }
}