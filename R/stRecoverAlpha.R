#' Matern covariance
#'
#' This function evaluates the Matern covariance function for the elements of 
#' a spatial distance matrix, d---i.e., a symmetric matrix with 0's along its
#' main diagonal.
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
#'
#' 
#' 

stRecoverAlpha = function( stFit, X, Y, Z, coords.local, coords.remote, nu_y, 
                           nu_r, burn, prob=.95, miles=T, ncores=1, 
                           summaryOnly=T  ) {
  
  maxIt = length(stFit$ll)
  nAlphas = length(burn:maxIt)
  
  n = nrow(coords.local)
  r = nrow(coords.remote)
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.local, miles=miles)
  Dz = rdist.earth(coords.remote, miles=miles)
  
  Xl = as.matrix(arrayToLong(X, coords.local, 1)[,-(1:3)])
  Yl = matrix(as.numeric(Y), ncol=1)
  Z = as.matrix(Z)
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  # initialize composition sample object
  alpha = list()
  
  
  mergeAlpha = function(x, y) {
    # use paralell algorithms for merging means, variances, etc.
    
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
  
  
  # estimate chunksize that will minimize number of function calls
  chunkSize = ceiling((maxIt-burn)/ncores)
  
  if(!summaryOnly) {
    # correct the estimate so that each rep will return <2GB of data since that
    # is currently the limit of data that mclapply can transfer. also subtract
    # a few samples to account for the summary data that gets returned too
    chunkSize = min(chunkSize, floor(2 / (nrow(Y)*nrow(Z)*8/1024^3) - 5))
  }
  
  
  # draw composition samples
  alpha = foreach( inds = ichunk(burn:maxIt, chunkSize = chunkSize),
                           .combine = mergeAlpha,
                           .options.multicore=mcoptions ) %dorng% {
    inds = unlist(inds)            
    
    .Call("_compositionAlpha", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, 
          Dz, nu_y, nu_r, stFit$beta[inds,], stFit$sigmasq_y[inds], 
          stFit$rho_y[inds], stFit$rho_r[inds], stFit$sigmasq_r[inds], 
          stFit$sigmasq_eps[inds], 0, summaryOnly)
  }
  
  # remove unwanted information
  attr(alpha, 'rng') = NULL
  alpha$beta = NULL
  
  # convert results away from unnecessary matrix format
  alpha$est = as.numeric(alpha$est)
  alpha$sd = as.numeric(alpha$sd)
  
  # compute approximate intervals, etc.
  alpha$summary = summariseAlpha(alpha, burn, prob, coords.local, coords.remote)
  
  # remove information redundant with the summary
  alpha$est = NULL
  alpha$sd = NULL
  
  stFit$alpha = alpha
  stFit   
}