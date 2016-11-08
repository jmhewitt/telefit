#' Plot stPredict objects
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
#' @param type Either 'prediction', 'residual', 'observed', 'standard_error' (or 'se'), 
#'  'local', 'remote', 'correlation', 'teleconnection', 'teleconnection_knot',
#'  or 'cat.prediction'
#'  to specify which part of stPredict to plot. Note that the value for type can
#'  be an abbreviation since partial matching is used during plotting. 'truth'
#'  and 'residual' plots are only available if the model has been evaluted and
#'  the predictions have been compared to another response dataset.
#' @param stPredict Object of class stPredict to plot.
#' @param t timepoint to plot.  Will automatically plot the first timepoint if
#'  t=NULL.
#' @param stData Object of class stData to provide coordinate and related
#'  information for plotting estimated teleconnection effects
#' @param stFit Object of class stFit to provide related
#'  information and structures for plotting estimated teleconnection effects
#' @param ... additional arguments to be passed to lower-level plotting functions
#'  
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stPredict = function( stPredict, type='prediction', t=NULL, stFit=NULL, 
                           stData=NULL, ... ) {

  # determine which type of plot is requested
  match.opts = c('prediction', 'residual', 'observed', 'standard_error', 'se', 
                 'local', 'remote', 'correlation', 'teleconnection', 
                 'cat.prediction', 'teleconnection_knot')
  type = match.opts[pmatch(type, match.opts)]
  
  # create a basic stData object for plotting
  if(is.null(stData)) {
    stData = list(
      coords.s = stPredict$coords.s
    )
    class(stData) = 'stData'
  }
  
  # get timepoint to plot
  if(is.null(t))
    t=stPredict$tLabs[1]
  
  # extract data for the right timepoint to plot
  if(length(stPredict$tLabs)==1) {
    pred = stPredict$pred
    stData$tLabs = stPredict$tLabs
  } else {
    pred = stPredict$pred[[match(t, stPredict$tLabs)]]
    stData$tLabs = t
  }
  
  # build plot
  if( type=='prediction' ) {
    stData$Y = pred$pred$Y
    stData$Y.lab = paste('Predicted', stPredict$Y.lab)
    ret = plot.stData(stData, ...)
  } else if( type=='residual' ) {
    stData$Y = pred$pred$resid
    stData$Y.lab = 'Residual'
    ret = plot.stData(stData, ...)
  } else if( type=='observed' ) {
    stData$Y = pred$pred$Y + pred$pred$resid
    stData$Y.lab = paste('Observed', stPredict$Y.lab)
    ret = plot.stData(stData, ...)
  } else if( type=='standard_error' || type=='se' ) {
    stData$Y = pred$pred$se
    stData$Y.lab = 'SE'
    ret = plot.stData(stData, ...)
  } else if( type=='local' ) {
    stData$Y = pred$pred$Y.local
    stData$Y.lab = paste('Local contribution to', stPredict$Y.lab)
    ret = plot.stData(stData, ...)
  } else if( type=='remote' ) {
    stData$Y = pred$pred$Y.remote
    stData$Y.lab = paste('Remote contribution to', stPredict$Y.lab)
    ret = plot.stData(stData, ...)
  } else if( type=='correlation' ) {
    ret = ggplot(pred$pred, aes(y=Y, x=Y+resid)) +
      geom_abline(slope=1, intercept=0) +
      geom_errorbar(aes(x=Y+resid, ymax = Y.upr, ymin = Y.lwr), 
                    alpha=.2) +
      geom_point(alpha=.5) +
      ylab(paste('Predicted', stPredict$Y.lab)) +
      xlab(paste('Observed', stPredict$Y.lab))
  } else if( type=='teleconnection' ) {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    stFit$alpha$summary = stPredict$alpha
    
    ret = plot.stFit(stFit = stFit, stData = stData, type='teleconnection', ...)
  } else if( type=='teleconnection_knot') {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    stFit$alpha_knots$summary = stPredict$alpha_knots
    
    ret = plot.stFit(stFit = stFit, stData = stData, type='teleconnection_knot', ...)
  } else if( type=='cat.prediction' ) {
    stData$Y.cat = factor(pred$pred$Y.cat)
    stData$Y.lab = paste('Predicted', stPredict$Y.lab)
    ret = plot.stData(stData, type='cat.response', ...)
  } 
  
  ret
}