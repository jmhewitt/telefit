#' Compute Highest posterior density intervals from posterior samples
#'
#' @export
#'
#' @importFrom coda HPDinterval mcmc
#' 
#' @param stFit stFit object containing posterior samples for model
#' @param burn number of posterior samples to reject before computing estimates
#' @param prob The target probability content of the intervals
#'  

HPDinterval.stFit = function(stFit, burn = 1, prob = .95) {
  res = lapply(stFit$parameters$samples, function(s) {
    HPDinterval(mcmc(as.matrix(s[-(1:burn),])))
  })
  
  # add names for the betas
  if(!is.null(stFit$parameters$beta.names)) {
    rownames(res$beta) = stFit$parameters$beta.names
  }
  
  res
}