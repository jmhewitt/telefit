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
      
      # use parallel algorithms for merging means, variances, etc.
      ZVar = x$ZVar
      r = nrow(ZVar)
      for(i in 1:dim(ZVar)[3]) {
        rStart = (i-1) * r + 1
        rEnd =  rStart + r - 1
        ZVar[,,i] = mergeCovmat(x$ZVar[,,i], y$ZVar[,,i],
                                (x$est[,rStart:rEnd]), (x$est[,rStart:rEnd]),
                                (y$est[,rStart:rEnd]), (y$est[,rStart:rEnd]),
                                x$nsamples, y$nsamples)
      }
      
      z = list(
        est = mergeMean(x$est, y$est, x$nsamples, y$nsamples),
        sd = sqrt(mergeVar(x$sd^2, y$sd^2, x$est, y$est, x$nsamples, y$nsamples)),
        covBetaAlpha = mergeCovmat(x$covBetaAlpha, y$covBetaAlpha,
                                   t(x$beta), t(x$est), t(y$beta), t(y$est),
                                   x$nsamples, y$nsamples),
        ZVar = ZVar,
        beta = mergeMean(x$beta, y$beta, x$nsamples, y$nsamples),
        nsamples = x$nsamples + y$nsamples,
        samples = rbind(x$samples, y$samples)
      )
      
      xfull$alpha = z
    }
    
    xfull
  }
}