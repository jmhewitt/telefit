#' Compute forecasts based on posterior samples
#'
#'
#'
#' @export
#' 
#' @importFrom doMC registerDoMC
#' @import doRNG
#' @importFrom foreach foreach
#' @importFrom fields rdist.earth
#' @importFrom mvtnorm rmvnorm
#' 
#'
#' 
#' 

stPredict = function(stFit, coords.local, X, Z, yearLabs, ncores=1) {
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  n = nrow(coords.local)
  nt = ncol(Z)
  
  I.ns = diag(n)
  
  
  # process each timepoint
  Y = foreach(t = 1:nt, .combine='c') %do% {
    
    # extract covariates for this timepoint
    x = X[,,t]
    z = Z[,t]
    
    # predict mean
    y.local = x %*% colMeans(stFit$parameters$samples$beta)
    y.remote = dgemkmm(I.ns, t(z), stFit$alpha$summary$alpha)
    
    # compute standard errors
    se = sqrt(mean(stFit$parameters$samples$sigmasq_y * 
                     (1+stFit$parameters$samples$sigmasq_eps)) + 
              2*diag(dgemkmm(I.ns, t(z), t(stFit$alpha$covBetaAlpha)) %*% t(x)) + 
              apply(stFit$alpha$ZVar, 3, function(A) { t(z) %*% A %*% z}))
    
    # return results
    list(list(
      pred = data.frame(
        y = y.local + y.remote,
        y.local = y.local,
        y.remote = y.remote,
        se = se,
        lon = coords.local[,1],
        lat = coords.local[,2]
      ),
      yrLab = yearLabs[t]
    ))
  }
  
  
 Y 
}