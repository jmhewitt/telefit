#' Simulate the spatio-temporal teleconnection model
#'
#' This function will let you provide as much information as you'd like, but 
#' will fill in the gaps as needed
#'
#'
#' @export
#' 
#' @importFrom foreach foreach "%do%"
#' @importFrom fields rdist.earth
#' @importFrom mvtnorm rmvnorm
#' 
#' 

stSimulate = function( X=NULL, Z=NULL, coords.s=NULL, coords.r=NULL, cov.s=NULL,
                       cov.r=NULL, alpha=NULL, beta=NULL, nt, ns, nr, p, 
                       miles=T ) {
  
  
  # default covariance parameters
  
  if(is.null(cov.s))
    cov.s = list()
  if(is.null(cov.s$var))
    cov.s$var = 1
  if(is.null(cov.s$range))
    cov.s$range = 1
  if(is.null(cov.s$smoothness))
    cov.s$smoothness = 2
  if(is.null(cov.s$nugget))
    cov.s$nugget = 0
    
  if(is.null(cov.r))
    cov.r = list()
  if(is.null(cov.r$var))
    cov.r$var = 1
  if(is.null(cov.r$range))
    cov.r$range = 1
  if(is.null(cov.r$smoothness))
    cov.r$smoothness = 2
  
  
  # simulate locations
  
  if(is.null(coords.s)) {
    coords.s = matrix(runif(2*ns), ncol=2)
  } else 
    ns = nrow(coords.s)
  
  if(is.null(coords.r)) { 
    coords.r = matrix(-runif(2*nr), ncol=2)
  } else
    nr = nrow(coords.r)
  
  
  # build covariance matrices
  
  Dy = rdist.earth(coords.s, miles=miles)
  Dz = rdist.earth(coords.r, miles=miles)
  
  Sigma = maternCov( Dy, scale = cov.s$var, range = cov.s$range, 
                     smoothness = cov.s$smoothness, 
                     nugget = cov.s$smoothness * cov.s$nugget )
  
  R = maternCov( Dz, scale = cov.r$var, range = cov.r$range, 
                     smoothness = cov.r$smoothness, nugget = 0 )
  
  Sigma.chol = chol(Sigma, pivot=T)
  Sigma.chol = t(Sigma.chol[, order(attr(Sigma.chol, 'pivot'))])
  
  
  # simulate covariate data
  
  if(is.null(Z)) {
    Z = matrix(rnorm(nr*nt), ncol=nt)
  } else
    nt = ncol(Z)
  
  if(is.null(X)) {
    X = foreach(t = 1:nt, .combine='abind3') %do% { 
      cbind( 1, matrix( rnorm(ns*(p-1)), ncol=p-1 ) )
    }
  } else
    p = ncol(X)

  
  # simulate parameters
  
  if(is.null(alpha))
    alpha = rmatnorm(1, R, Sigma)
  
  if(is.null(beta))
    beta = rnorm(p)
  
  
  # simulate response data
  
  Y = foreach(t = 1:nt, .combine='cbind') %do% {
    Zta = apply(alpha, 2, function(a) { t(Z[,t]) %*% a })
    X[,,t] %*% beta + Zta + Sigma.chol %*% rnorm(ns)
  }
  
  
  # package and return data
  
  res = list(
    tLabs = 1:nt,
    coords.s = coords.s,
    coords.r = coords.r,
    X = X,
    Y = Y,
    Z = Z,
    X.lab = 'X',
    Y.lab = 'Y',
    Z.lab = 'Z',
    beta = beta,
    alpha = as.numeric(alpha)
  )
  
  class(res) = 'stData'
  
  res
}