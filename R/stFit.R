#' Fits the spatial teleconnection model
#'
#' Fits the teleconnection model with L = Sigma
#'
#' @export
#'
#' @importFrom fields rdist.earth
#' @useDynLib telefit
#'
#' @param nugget Spatial covariance nugget.  
#' @param rho_y_sd Initial sd for rho_y RW updates
#'
#' 
#' 

stFit = function(X, Y, Z, coords.local, coords.remote, Lambda, nu_y, nu_r, ay, 
                 by, ar, br, ary, bry, arr, brr, aeps=2, beps=1, rho_y_sd=.07, 
                 rho_r_sd=.15, eps_sd=.09, sigmasq_r_sd=.1, 
                 maxIt, returnll=T, miles=T, C=1, RWrate=.44) {
  
  
  n = nrow(coords.local)
  r = nrow(coords.remote)
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.local, miles=miles)
  Dz = rdist.earth(coords.remote, miles=miles)
  
  Xl = as.matrix(arrayToLong(X, coords.local, 1)[,-(1:3)])
  Yl = matrix(as.numeric(Y), ncol=1)
  Z = as.matrix(Z)
  
  
  ptm = proc.time()
  
  res = .Call("_stfit", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, Dz, 
              nu_y, nu_r, ay, by, ar, br, ary, bry, arr, brr, aeps, beps, 
              Lambda, rho_y_sd, rho_r_sd, eps_sd, sigmasq_r_sd, maxIt, returnll, 
              errDump, C, RWrate)
  
  ptm = proc.time() - ptm
     
  message('Total time (min): ', signif(ptm[3]/60, 3))
  message('Samples per second: ', signif(maxIt/ptm[3],3))
  
  reslist = list(
    parameters = list(samples = res)
  )
  class(reslist) = 'stFit'
  
  reslist
}