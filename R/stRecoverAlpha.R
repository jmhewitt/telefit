#' Estimate teleconnection effects
#' 
#' This method uses composition sampling from the posterior parameter estimates
#' and can be executed in parallel.
#'
#' @export
#' 
#' @importFrom doMC registerDoMC
#' @import doRNG
#' @importFrom itertools ichunk
#' @importFrom foreach foreach "%dopar%"
#' @importFrom fields rdist.earth
#' 
#' @useDynLib telefit
#'
#' @param summaryOnly Only return alpha estimates and key variances; do not 
#'        return all posterior samples.
#' @param stFit Object with class 'stFit' containing posterior parameter samples
#'  needed to composition sample the teleconnection effects. The data needed 
#'  from stFit need only be manually entered if not using a stData object.
#' @param stData Object with class 'stData' containing data needed to fit this 
#'  model. The data need only be manually entered if not using a stData object.
#' @param burn number of posterior samples to burn before drawing composition
#'  samples
#' @param prob confidence level for approximate confidence intervals of 
#'  teleconnection effects
#' @param X [ns, p, nt] array of design matrices with local covariates
#' @param Y [ns, nt] matrix with response data
#' @param Z [nr, nt] matrix with remote covariates
#' @param coords.s matrix with coordinates where responses were 
#'  observed (lon, lat)
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param ncores Since the teleconnection effects can be sampled in parallel,
#'  this parameter lets users specify the number of cores to use to sample 
#'  teleconnection effects
#'  
#' 

stRecoverAlpha = function( stFit, stData, burn, ncores=1, prob=.95, 
                           summaryOnly=T, X = stData$X, Y = stData$Y, 
                           Z = stData$Z, coords.s = stData$coords.s, 
                           coords.r = stData$coords.r ) {
  
  maxIt = length(stFit$parameters$samples$ll)
  nAlphas = length(burn:maxIt)
  
  n = nrow(coords.s)
  r = nrow(coords.r)
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.s, miles=stFit$miles)
  Dz = rdist.earth(coords.r, miles=stFit$miles)
  
  Xl = as.matrix(arrayToLong(X, coords.s, 1)[,-(1:3)])
  Yl = matrix(as.numeric(Y), ncol=1)
  Z = as.matrix(Z)
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  # initialize composition sample object
  alpha = list()
  
  
  mergeAlpha = function(x, y) {
    if(is.null(y)) {
      x
    } else {
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
      
      z
    }
  }
  
  
  # estimate chunksize that will minimize number of function calls
  chunkSize = ceiling((maxIt-burn+1)/ncores)
  
  if(!summaryOnly) {
    # correct the estimate so that each rep will return <2GB of data since that
    # is currently the limit of data that mclapply can transfer. also subtract
    # a few samples to account for the summary data that gets returned too
    chunkSize = max(min(chunkSize, floor(2 / (nrow(Y)*nrow(Z)*8/1024^3) - 5)), 1)
  }
  
  # make sure that no chunks have fewer than 5 elements
  
  # draw composition samples
  alpha = foreach( inds = ichunk(burn:maxIt, chunkSize = chunkSize),
                           .combine = mergeAlpha,
                           .options.multicore=mcoptions ) %dorng% {
    inds = unlist(inds)            
    
    if(stFit$varying) {
      .Call("_stvCompositionAlpha", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, 
            Dz, stFit$priors$cov.s$smoothness, stFit$priors$cov.r$smoothness, 
            stFit$parameters$samples$beta[inds,], 
            stFit$parameters$samples$sigmasq_y[inds], 
            stFit$parameters$samples$rho_y[inds], 
            stFit$parameters$samples$rho_r[inds], 
            stFit$parameters$samples$sigmasq_r[inds], 
            stFit$parameters$samples$sigmasq_eps[inds], 0, summaryOnly)
    } else {
      .Call("_compositionAlpha", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, 
            Dz, stFit$priors$cov.s$smoothness, stFit$priors$cov.r$smoothness, 
            stFit$parameters$samples$beta[inds,], 
            stFit$parameters$samples$sigmasq_y[inds], 
            stFit$parameters$samples$rho_y[inds], 
            stFit$parameters$samples$rho_r[inds], 
            stFit$parameters$samples$sigmasq_r[inds], 
            stFit$parameters$samples$sigmasq_eps[inds], 0, summaryOnly)
    }
    
  }
  
  # remove unwanted information
  attr(alpha, 'rng') = NULL
  alpha$beta = NULL
  if(summaryOnly)
    alpha$samples = NULL
  
  # convert results away from unnecessary matrix format
  alpha$est = as.numeric(alpha$est)
  alpha$sd = as.numeric(alpha$sd)
  
  # compute approximate intervals, etc.
  alpha$summary = summariseAlpha(alpha, burn, prob, coords.s, coords.r)
  
  # remove information redundant with the summary
  alpha$est = NULL
  alpha$sd = NULL
  
  stFit$alpha = alpha
  stFit   
}