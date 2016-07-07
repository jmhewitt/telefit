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

stPredict = function(stFit, coords.local, X, Z, yearLabs, ncores=1, 
                     localOnly=F, conf = .95) {
  
  # get critical value for forming approximate normal confidence intervals
  zcrit = qnorm((1-conf)/2, lower.tail = F)
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  n = nrow(coords.local)
  nt = dim(X)[3]
  
  I.ns = diag(n)
  
  # process each timepoint
  Y = foreach(t = 1:nt, .combine='c') %do% {
    
    # extract covariates for this timepoint
    x = X[,,t]
    if(!localOnly)
      z = Z[,t]
    
    # predict mean
    y.local = x %*% colMeans(stFit$parameters$samples$beta)
    if(!localOnly)
      y.remote = dgemkmm(I.ns, t(z), stFit$alpha$summary$alpha)
    
    # compute standard errors
    if(localOnly) {
      se = sqrt(mean(stFit$parameters$samples$sigmasq_y * 
                       (1+stFit$parameters$samples$sigmasq_eps)) 
                  )
    } else {
      se = sqrt(mean(stFit$parameters$samples$sigmasq_y * 
                       (1+stFit$parameters$samples$sigmasq_eps)) + 
                  2*diag(dgemkmm(I.ns, t(z), t(stFit$alpha$covBetaAlpha)) %*% t(x)) +
                  apply(stFit$alpha$ZVar, 3, function(A) { t(z) %*% A %*% z}))
    }
    
    # package results
    pred = 
    # package results
    if(localOnly) {
      pred = data.frame(
        Y = y.local,
        se = se,
        Y.lwr = y.local - zcrit * se,
        Y.upr = y.local + zcrit * se,
        lon = coords.local[,1],
        lat = coords.local[,2]
      )
    } else {
      pred = data.frame(
        Y = y.local + y.remote,
        Y.local = y.local,
        Y.remote = y.remote,
        se = se,
        Y.lwr = y.local + y.remote - zcrit * se,
        Y.upr = y.local + y.remote + zcrit * se,
        lon = coords.local[,1],
        lat = coords.local[,2]
      )
    }
    
    
    
    # return results
    list(list(
      pred = pred,
      yrLab = yearLabs[t]
    ))
  }
  
  
 Y 
}