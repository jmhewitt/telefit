#' Make predictions using a fitted varying coefficient model
#'
#'
#' @export
#'
#' @importFrom fields rdist.earth
#' @useDynLib telefit
#'
#' 
#' @example examples/svcMod.R


svcPredict = function(fit, Xn, Zn, burn=0) {
  
  D = rdist.earth(fit$coords, miles=fit$miles)
  
  if(burn>0) {
    
  }
    
  res = .Call("_svcpredict", PACKAGE = 'telefit', fit$parameters$samples, Xn, 
              Zn, D, fit$priors$cov$nu)
  
  # allow fit evaluations through a separate function, perhaps stEval
  
  reslist = list(
    samples = res
  )
  
  class(reslist) = 'svcPredict'
  
  reslist
}