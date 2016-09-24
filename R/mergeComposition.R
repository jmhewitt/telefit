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
      x$forecasts = abind3(x$forecasts, y$forecasts)
      x$local = abind3(x$local, y$local)
      x$remote = abind3(x$remote, y$remote)
      x$noise = abind3(x$noise, y$noise)
      
      xfull$forecast = x
    }
    
    
    # merge teleconnection effects
    if(!is.null(xfull$alpha)) {
      # extract alpha parts
      x = xfull$alpha
      y = yfull$alpha
      
      z = list(
        est = mergeMean(x$est, y$est, x$nsamples, y$nsamples),
        sd = sqrt(mergeVar(x$sd^2, y$sd^2, x$est, y$est, x$nsamples, y$nsamples)),
        nsamples = x$nsamples + y$nsamples,
        samples = rbind(x$samples, y$samples)
      )
      
      xfull$alpha = z
    }
    
    xfull
  }
}