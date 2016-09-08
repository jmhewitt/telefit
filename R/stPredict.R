#' Compute forecasts based on posterior samples
#'
#'
#'
#' @export
#' 
#' @importFrom doMC registerDoMC
#' @import doRNG
#' @import foreach
#' @import Matrix
#' @importFrom fields rdist.earth
#' @importFrom mvtnorm rmvnorm
#' 
#'
#' @param stFit object of class stFit that includes estimated teleconnection 
#'  effects.  This information will be used to make predictions
#' @param stData object of class stData that includes information needed for 
#'  making forecasts.  If response data is included, this function will 
#'  automatically run stEval using the empirical climatology as the reference
#'  forecast
#' @param burn Number of posterior samples to ignore
#' 

stPredict = function( stFit, stData, ncores = 1, conf = .95,
                      coords.s = stData$coords.s, X = stData$X, 
                      Z = stData$Z, tLabs = stData$tLabs,
                      burn = 1) {
  
  # extract localOnly and varying
  localOnly = stFit$localOnly
  varying = stFit$varying
    
  # get critical value for forming approximate normal confidence intervals
  zcrit = qnorm((1-conf)/2, lower.tail = F)
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  n = nrow(coords.s)
  nt = dim(X)[3]
  nt = ifelse(is.na(nt), 1, nt)
  
  I.ns = diag(n)
  
  # mean beta estimate
  betaH = colMeans(stFit$parameters$samples$beta[-(1:burn),])
  
  # process each timepoint
  Y = foreach(t = 1:nt) %dopar% {

    # extract covariates for this timepoint
    if(nt==1) {
      x = X[,]
    } else {
      x = X[,,t]
    }
    
    if(varying) {
      # convert to sparse spatially varying format
      x = t(sparseMatrix(i = 1:prod(dim(x)),
                         j = rep(1:nrow(x), rep(ncol(x),nrow(x))),
                         x = as.numeric(t(x))))
    }
      
    if(!localOnly) {
      if(nt==1) {
        z = Z
      } else {
        z = Z[,t]
      }
    }

    # predict mean
    y.local = as.numeric(x %*% betaH)
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
    if(localOnly) {
      pred = data.frame(
        Y = y.local,
        se = se,
        Y.lwr = y.local - zcrit * se,
        Y.upr = y.local + zcrit * se
      )
    } else {
      pred = data.frame(
        Y = y.local + y.remote,
        Y.local = y.local,
        Y.remote = y.remote,
        se = se,
        Y.lwr = y.local + y.remote - zcrit * se,
        Y.upr = y.local + y.remote + zcrit * se
      )
    }

    # package and return results
    r = list(
      pred = pred,
      yrLab = tLabs[t]
    )

    r
  }
  

  ret = list(
    pred = Y,
    coord.s = coords.s,
    localOnly = localOnly,
    varying = varying,
    tLabs = tLabs,
    Y.lab = stData$Y.lab
  )
  class(ret) = 'stPredict'
  
  # evaluate performance if response data is given
  if(!is.null(stData$Y)) {
    if(is.null(ncol(stData$Y)))
      ret = stEval(ret, stData$Y, stData$Y)
    else 
      ret = stEval(ret, stData$Y, rowMeans(stData$Y))
  }
  
  ret
}