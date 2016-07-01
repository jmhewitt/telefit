#' Creates a pairwise plot of the sampled parameters
#'
#' @export
#'
#' 
#' @param stFit return object from telefit::stFit
#' @param burn number of samples to discard before plotting
#' 
#' 

pairsplot = function(stFit, burn=1) {
  
  # extract posterior samples
  res.df = data.frame(stFit$parameters$samples)
  maxIt = nrow(res.df)
  
  # discard burned samples
  res.df = res.df[burn:maxIt,]
  
  plot(res.df)
}