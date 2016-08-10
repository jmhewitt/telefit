#' Plot stData objects
#'
#' This function provides basic plotting for telefit package data.
#' 
#'
#' @export
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' 
#' @param boxsize size of grid boxes plotted
#' @param map name of map provided by the maps package. These include county, 
#'  france, italy, nz, state, usa, world, world2.  By default, all stData plots
#'  will include us state outlines.
#' @param region name of subregions to include. Defaults to . which includes 
#'  all subregions. See documentation for map for more details.
#' @param type Either 'traceplot', 'density', 'pairs', or 'teleconnection' 
#'  to specify which part of stFit to plot. Note that the value for type can
#'  be an abbreviation since partial matching is used during plotting.
#' @param stFit Object of class stFit to plot.
#' @param coord.s if plot type is 'teleconnection', specifies the longitude and 
#'  latitude of local coordinate for which to plot estimated teleconnection 
#'  effects. if NULL, the middle local coordinate will be plotted.
#' @param zlim c(min, max) vector that specifies the colorscale limits
#' @param text.size number specifying the size of text labels
#' @param axis.text.size number specifying the size of axis text labels
#' @param burn number of observations to exclude from graph
#' @param stData Object of class stData to provide coordinate and related
#'  information for plotting estimated teleconnection effects
#' 
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stFit = function( stFit, type='density', boxsize=NULL, stData=NULL,
                       map='world', region='.', coord.s=NULL, zlim=NULL,
                       text.size=NULL, axis.text.size=NULL, burn = 1) {

  # determine which type of plot is requested
  match.opts = c('traceplot', 'density', 'pairs', 'teleconnection')
  type = match.opts[pmatch(type, match.opts)]
  
  # extract posterior samples if necessary
  if( type %in% c('traceplot', 'density', 'pairs') ) {
    # extract posterior samples
    res.df = data.frame(stFit$parameters$samples)
    maxIt = nrow(res.df)
    
    # discard burned samples
    res.df = res.df[burn:maxIt,]
    res.df$Iteration = burn:maxIt
    
    # coerce to plottable form
    res.plottable = melt(res.df, id.vars = 'Iteration', variable.name = 'param', 
                         value.name = 'Value')
  }
  
  
  # build plots
  if( type=='traceplot' ) {
    ret = ggplot(res.plottable, aes(x=Iteration, y=Value)) +
      geom_line() +
      facet_wrap(~param, scales='free')
    
  } else if( type=='density' ) {
    ret = ggplot(res.plottable, aes(x=Value)) +
      geom_density() +
      facet_wrap(~param, scales='free') +
      ylab('Posterior density')
    
  } else if( type=='pairs' ) {
    plot(res.df)
    
  } else if( type=='teleconnection' ) {
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    
    stData$alpha = stFit$alpha$summary$alpha
    ret = plot.stData(stData, 'tele', boxsize = boxsize, map = map, region = region,
               coord.s = coord.s, zlim = zlim, 
               lab.teleconnection = expression(hat(alpha))) + 
      ggtitle('Estimated teleconnection effects')
  }
  
  
  # modify text sizes if requested
  if(!is.null(text.size))
    ret = ret + theme( text = element_text(size=text.size))
  if(!is.null(axis.text.size))
    ret = ret + theme( axis.text = element_text(size=axis.text.size))

  # return plot
  if(!is.null(ret))
    ret
}