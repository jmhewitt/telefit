#' Creates a density plot of the sampled parameters
#'
#' @export
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' 
#' @param stFit return object from telefit::stFit
#' @param burn number of samples to discard before plotting
#' 
#' @return a ggplot object with the density plots
#' 

densityplot = function(stFit, burn=1, text.size=36, axis.text.size=24) {
  
  # extract posterior samples
  res.df = data.frame(stFit$parameters$samples)
  maxIt = nrow(res.df)
  
  # discard burned samples
  res.df = res.df[burn:maxIt,]
  res.df$Iteration = burn:maxIt
  
  # build plot
  ggplot(melt(res.df, id.vars = 'Iteration', variable.name = 'param', 
              value.name = 'Value'),
         aes(x=Value)) +
    geom_density() +
    facet_wrap(~param, scales='free') +
    theme( text = element_text(size=text.size),
           axis.text = element_text(size=axis.text.size) ) +
    ylab('Posterior density')
}