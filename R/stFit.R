#' Fits the spatial teleconnection model
#'
#' Fits the teleconnection model with L = Sigma
#'
#' @export
#'
#' @importFrom fields rdist.earth
#' @useDynLib telefit
#'
#' @param localOnly TRUE to fit the model without the teleconnection effects
#'  (typically for evaluating impact of teleconnection effects)
#' @param stData Object with class 'stData' containing data needed to fit this 
#'  model. The data need only be manually entered if not using a stData object.
#' @param X [ns, p, nt] array of design matrices with local covariates
#' @param Y [ns, nt] matrix with response data
#' @param Z [nr, nt] matrix with remote covariates
#' @param coords.s matrix with coordinates where responses were 
#'  observed (lon, lat)
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param priors A list containing parameters for the prior distributions. The
#'  list needs to contain the following values:
#'    \describe{
#'      \item{beta} { list(Lambda=matrix) specifying the prior covariance matrix
#'        for the local effects if varying==F; otherwise 
#'        list(Psi=matrix, nu=double) specifying the Inverse wishart prior 
#'        distribution for the spatially varying coefficient process if 
#'        varying==T.   }
#'      
#'      \item{cov.s}{ list(smoothness=double, range=c(min, max), 
#'        variance=c(shape, rate), nugget=c(shape, rate)) }
#'        
#'      \item{cov.r}{ list(smoothness=double, range=c(min, max), 
#'        variance=c(shape, rate)) }
#'    }
#' @param rw.initsd A list containing initial standard deviation parameters for
#'  the MCMC parameters requiring random walk updates:
#'    \describe{
#'      \item{cov.s}{ list(range=double, nugget=double) }
#'      \item{cov.r}{ list(range=double, variance=double) }
#'    }
#' @param maxIt number of iterations to run the MCMC chain for
#' @param returnll TRUE to compute the model log-likelihood at each iteration
#' @param miles TRUE if covariance matrix distances should be in miles, FALSE 
#'  for kilometers
#' @param C scaling factor used in adapting random walk proposal variances.
#' @param alpha target acceptance rate for random walk proposals.
#' @param varying TRUE to fit the model with spatially varying local coefficients
#' 


stFit = function( stData = NULL, priors, maxIt, X = stData$X, Y = stData$Y, 
                  Z = stData$Z, coords.s = stData$coords.s, 
                  coords.r = stData$coords.r, rw.initsd = NULL, 
                  returnll = T, miles = T, C=1, alpha=.44, localOnly = F,
                  varying = T ) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.s, miles=miles)
  
  Xl = as.matrix(arrayToLong(X, coords.s, 1)[,-(1:3)])
  Yl = matrix(as.numeric(Y), ncol=1)
  
  # default random walk proposal standard deviations
  if(is.null(rw.initsd))
    rw.initsd = list(
      cov.s = list(range = .07, nugget = .09),
      cov.r = list(range = .15, variance = .1)
    )
  
  # fit model
  if(localOnly) {
    
    ptm = proc.time()
    
    res = .Call("_sfit", PACKAGE = 'telefit', p, n, t, Xl, Yl, Dy, 
                priors$cov.s$smoothness, priors$cov.s$variance[1], 
                priors$cov.s$variance[2], priors$cov.s$range[1], 
                priors$cov.s$range[2], priors$cov.s$nugget[1], 
                priors$cov.s$nugget[2], priors$beta$Lambda, rw.initsd$cov.s$range, 
                rw.initsd$cov.s$nugget, maxIt, returnll, errDump, C, alpha)
    
    ptm = proc.time() - ptm
  } else {
    
    r = nrow(coords.r)
    Dz = rdist.earth(coords.r, miles=miles)
    Z = as.matrix(Z)
    
    ptm = proc.time()
    
    if(varying) {
      res = .Call("_stvfit", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, Dz, 
                  priors$cov.s$smoothness, priors$cov.r$smoothness, 
                  priors$cov.s$variance[1], priors$cov.s$variance[2], 
                  priors$cov.r$variance[1], priors$cov.r$variance[2], 
                  priors$cov.s$range[1], priors$cov.s$range[2],
                  priors$cov.r$range[1], priors$cov.r$range[2],
                  priors$cov.s$nugget[1], priors$cov.s$nugget[2], priors$beta$Psi,
                  rw.initsd$cov.s$range, rw.initsd$cov.r$range, 
                  rw.initsd$cov.s$nugget, rw.initsd$cov.r$variance, maxIt, 
                  returnll, errDump, C, alpha, priors$beta$nu)
      
    } else {
      res = .Call("_stfit", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, Dz, 
                  priors$cov.s$smoothness, priors$cov.r$smoothness, 
                  priors$cov.s$variance[1], priors$cov.s$variance[2], 
                  priors$cov.r$variance[1], priors$cov.r$variance[2], 
                  priors$cov.s$range[1], priors$cov.s$range[2],
                  priors$cov.r$range[1], priors$cov.r$range[2],
                  priors$cov.s$nugget[1], priors$cov.s$nugget[2], priors$beta$Lambda,
                  rw.initsd$cov.s$range, rw.initsd$cov.r$range, 
                  rw.initsd$cov.s$nugget, rw.initsd$cov.r$variance, maxIt, 
                  returnll, errDump, C, alpha)
    }
    
    ptm = proc.time() - ptm
  }
    
    
  message('Total time (min): ', signif(ptm[3]/60, 3))
  message('Samples per second: ', signif(maxIt/ptm[3],3))
  
  reslist = list(
    parameters = list(samples = res),
    priors = priors,
    miles = miles,
    localOnly = localOnly,
    varying = varying
  )
  
  class(reslist) = 'stFit'
  
  reslist
}