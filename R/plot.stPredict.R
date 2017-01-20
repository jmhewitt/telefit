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
#' @importFrom fields rdist.earth
#' 
#' @param type Either 'prediction', 'residual', 'observed', 'standard_error' 
#'  (or 'se'), 'eof_alpha',
#'  'local', 'remote', 'correlation', 'teleconnection', 'teleconnection_knot',
#'  'teleconnection_knot_transect', 'errors', or 'cat.prediction'
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
#' @param err.comparison data.frame with Year column  and a column for a variable
#'  that will be used to plot annual errors against
#' @param err.var name of variable in err.comparison for plotting against
#' @param err.lab label for name of variable in err.comparison for plotting against
#' @param pattern if type=='eof_alpha', this specified which eof the remote 
#'  coefficients should be mapped onto and then plotted over the local domain
#' @param burn number of observations to exclude from graph
#' @param ... additional arguments to be passed to lower-level plotting functions
#'  
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stPredict = function( stPredict, type='prediction', t=NULL, stFit=NULL, 
                           stData=NULL, err.comparison=NULL, err.var=NULL,
                           err.lab=err.var, pattern=1, dots=NULL, burn=1, ... ) {

  # merge unique list of dots
    dots = c(dots, list(...))
    dots = dots[!duplicated(dots)]
  # overwrite arguments to function if they exist in dots
    for(x in setdiff(names(formals(eval(match.call()[[1]]))), c('dots', '...'))) {
      if(x %in% names(dots)) {
        assign(eval(x), dots[[x]])
      }
    }
  
  # determine which type of plot is requested
  match.opts = c('prediction', 'residual', 'observed', 'standard_error', 'se', 
                 'local', 'remote', 'correlation', 'teleconnection', 
                 'cat.prediction', 'teleconnection_knot', 'eof_alpha',
                 'teleconnection_knot_transect', 'errors', 'teleconnection_knot_local')
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
    # stData$Y.lab = paste('Predicted', stPredict$Y.lab)
    stData$Y.lab = stPredict$Y.lab
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='residual' ) {
    stData$Y = pred$pred$resid
    stData$Y.lab = 'Residual'
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='observed' ) {
    stData$Y = pred$pred$Y + pred$pred$resid
    stData$Y.lab = paste('Observed', stPredict$Y.lab)
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='standard_error' || type=='se' ) {
    stData$Y = pred$pred$se
    stData$Y.lab = 'SE'
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='local' ) {
    stData$Y = pred$pred$Y.local
    stData$Y.lab = paste('Local contribution to', stPredict$Y.lab)
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='remote' ) {
    stData$Y = pred$pred$Y.remote
    stData$Y.lab = paste('Remote contribution to', stPredict$Y.lab)
    ret = plot.stData(stData, dots=dots, ...)
  } else if( type=='correlation' ) {
    ret = ggplot(pred$pred, aes(y=Y, x=Y+resid)) +
      geom_abline(slope=1, intercept=0) +
      geom_errorbar(aes(x=Y+resid, ymax = Y.upr, ymin = Y.lwr), 
                    alpha=.2) +
      geom_point(alpha=.5) +
      ylab(paste('Predicted', stPredict$Y.lab)) +
      xlab(paste('Observed', stPredict$Y.lab))
  } else if( type=='eof_alpha' ) {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    # compute eofs
    eof = prcomp(stData$Z, center = F)
    W = -eof$x
    Tmat = -eof$rotation
    
    # compute correlation matrices
    Dz_knots = rdist.earth(stFit$coords.knots, miles=stFit$miles)
    Dz_to_knots = rdist.earth(stData$coords.r, stFit$coords.knots, miles=stFit$miles)
    Rst = maternCov( Dz_knots, smoothness = stFit$priors$cov.r$smoothness,
                     scale = mean(stFit$parameters$samples$sigmasq_r[-(1:burn)]),
                     range = mean(stFit$parameters$samples$rho_r[-(1:burn)]) )
    cst = maternCov( Dz_to_knots, smoothness = stFit$priors$cov.r$smoothness,
                     scale = mean(stFit$parameters$samples$sigmasq_r[-(1:burn)]),
                    range = mean(stFit$parameters$samples$rho_r[-(1:burn)]) )
    
    # compute eof coefficients for all locations in "local" domain
    A=matrix(stPredict$alpha_knots$alpha, nrow = nrow(stFit$coords.knots))
    Ast = t(W) %*% cst %*% solve(Rst) %*% A
    
    # extract and plot data
    stData$Y = t(Ast)
    stData$tLabs = 1:nrow(Ast)
    stData$Y.lab = 'Coef.'
    ret = plot.stData(stData, type='response', t=pattern) + 
      ggtitle(bquote(alpha*"'"[.(pattern)]))
    
  } else if( type=='teleconnection' ) {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    stFit$alpha$summary = stPredict$alpha
    
    ret = plot.stFit(stFit = stFit, stData = stData, type='teleconnection', 
                     dots=dots, ...)
  } else if( type=='teleconnection_knot') {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    stFit$alpha_knots$summary = stPredict$alpha_knots
    
    ret = plot.stFit(stFit = stFit, stData = stData, type='teleconnection_knot', 
                     dots=dots, ...)
  } else if( type=='teleconnection_knot_local') {
    
    if(is.null(stData)) {
      stop('stData object required for plotting estimated teleconnection effects.')
    }
    if(is.null(stFit)) {
      stop('stFit object required for plotting estimated teleconnection effects.')
    }
    
    stFit$alpha_knots$summary = stPredict$alpha_knots
    
    ret = plot.stFit(stFit = stFit, stData = stData, type='teleconnection_knot_local', 
                     dots=dots, ...)
    
  } else if( type=='teleconnection_knot_transect') {
    
    stFit$alpha_knots$summary = stPredict$alpha_knots
    
    ret = plot.stFit(stFit = stFit, type='teleconnection_knot_transect', 
                     dots=dots, ...)
    
  } else if( type=='cat.prediction' ) {
    stData$Y.cat = factor(pred$pred$Y.cat)
    # stData$Y.lab = paste('Predicted', stPredict$Y.lab)
    stData$Y.lab = stPredict$Y.lab
    ret = plot.stData(stData, type='cat.response', 
                      category.breaks = stPredict$category.breaks, 
                      dots=dots, ...)
  } else if( type=='errors' ) {
    
    # extract yearly errors
    errsum = foreach(f=stPredict$pred, .combine='rbind') %do% {
      data.frame(Year=f$yrLab, f$err)
    } %>% mutate(RMSPE = sqrt(mspe), mspe = NULL,
                 'Brier skill' = bss, bss = NULL,
                 'Categorical Accuracy' = cat.correct, cat.correct = NULL,
                 'Heidke skill' = cat.heidke, cat.heidke = NULL,
                 'Response Correlation' = cor, cor = NULL,
                 'R sq.' = r2, r2 = NULL,
                 'PPL' = ppl, ppl=NULL,
                 'Coverage' = coverage, coverage=NULL
    )
    
    errsum = cbind(errsum, err.comparison)
    
    errsum.plottable = melt(errsum, id.vars = c('Year', err.var))
    
    ret = ggplot(errsum.plottable, aes_string(x=err.var, y='value')) +
      facet_wrap(~variable, scales='free') +
      geom_smooth(span=1) +
      geom_label(aes(label=Year), alpha=.8) +
      ylab('') +
      xlab(err.lab) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  }
  
  ret
}