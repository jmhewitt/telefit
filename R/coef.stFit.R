#' Compute point estimates for parameters from posterior samples
#'
#' @export
#'
#' 
#' @param stFit stFit object containing posterior samples for model
#' @param burn number of posterior samples to reject before computing estimates
#' @param fun function for computing point estimates
#'
#'  

coef.stFit = function(stFit, burn = 1, fun = mean) {
  lapply(stFit$parameters$samples, function(s) {
    apply(as.matrix(s[-(1:burn),]), 2, fun)
  })
}