#' Plot stFit objects
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
#' @param type Either 'traceplot', 'density', 'pairs', 'teleconnection',
#'  'teleconnection_knot', or 'beta' to specify which part of stFit to plot. 
#'  Note that the value for
#'  type can be an abbreviation since partial matching is used during plotting.
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
#' @param signif.telecon if TRUE, will highlight significant teleconnection
#'  effects when type=='teleconnection'
#' @param p If stFit was fit with spatially varying coefficients, p specifies 
#'  the index of the spatially varying coefficient to plot
#' 
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stFit = function( stFit, type='density', boxsize=NULL, stData=NULL,
                       map='world', region='.', coord.s=NULL, zlim=NULL,
                       text.size=NULL, axis.text.size=NULL, burn = 1,
                       signif.telecon = F, p = 1 ) {

  # determine which type of plot is requested
  match.opts = c('traceplot', 'density', 'pairs', 'teleconnection', 'beta',
                 'teleconnection_knot')
  type = match.opts[pmatch(type, match.opts)]
  
  # extract posterior samples if necessary
  if( type %in% c('traceplot', 'density', 'pairs') ) {
    # extract posterior samples
    if(stFit$varying) {
      listInd = which(!(names(stFit$parameters$samples) %in% c('beta', 'T')))
      res.df = data.frame(stFit$parameters$samples[listInd])
    } else {
      res.df = data.frame(stFit$parameters$samples)
    }
    maxIt = nrow(res.df)
    
    # discard burned samples
    res.df = res.df[burn:maxIt,]
    res.df$Iteration = burn:maxIt
    
    # add names for the betas
    if(!is.null(stFit$parameters$beta.names)) {
      colnames(res.df)[
        1:length(stFit$parameters$beta.names)] = stFit$parameters$beta.names
    }
    
    # coerce to plottable form
    res.plottable = melt(res.df, id.vars = 'Iteration', variable.name = 'param', 
                         value.name = 'Value')
  }
  
  
  # build plots
  ret = NULL
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
    
    coord.s = unlist(coord.s)
    
    if(signif.telecon) {
      teleCor = list(
        cor = matrix(stFit$alpha$summary$alpha, nrow = nrow(stData$coords.s),
                     byrow = T),
        coords.s = stData$coords.s,
        coords.r = stData$coords.r,
        signif = matrix(stFit$alpha$summary$signif, nrow = nrow(stData$coords.s),
                        byrow = T)
      )
      class(teleCor) = 'teleCor'
      
      ret = plot.teleCor(teleCor, signif = T, coord.s = coord.s, 
                         boxsize = boxsize, map = map, region = region, 
                         zlim = zlim )
      
    } else {
      stData$alpha = stFit$alpha$summary$alpha
      ret = plot.stData(stData, 'tele', boxsize = boxsize, map = map, 
                        region = region, coord.s = coord.s, zlim = zlim, 
                        lab.teleconnection = 'alpha') + 
        ggtitle('Estimated teleconnection effects')
    }
    
  } else if( type=='teleconnection_knot' ) {
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    
    coord.s = unlist(coord.s)
    
    stData$alpha_knots = stFit$alpha_knots$summary$alpha
    stData$alpha_knots_signif = stFit$alpha_knots$summary$signif
    stData$coords.knots = stFit$coords.knots
    
    ret = plot.stData(stData, 'teleconnection_knot', boxsize = boxsize, map = map, 
                      region = region, coord.s = coord.s, zlim = zlim, 
                      lab.teleconnection = 'alpha', 
                      signif.telecon = signif.telecon) + 
      ggtitle('Estimated teleconnection effects')
    
  } else if( type=='beta' ) {
    if(is.null(stData)) {
      stop('stData object required for plotting estimated spatially varying coefficients.')
    }
    
    # extract coefficient estimates
    betaH = colMeans(stFit$parameters$samples$beta[-(1:burn),
        seq(from = p, to = prod(dim(stData$X)[1:2]), by = ncol(stData$X) )])
    
    # put estimates into plottable structure
    stData$Y[,1] = betaH
    stData$Y.lab = 'Coefficient'
    
    # build plot
    ret = plot.stData(stData, 'response') + 
      ggtitle('Estimated spatially varying coefficient')
    
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