#' Basic evaluation of fit
#'
#' Provides basic measures for evalutating the fit.  Includes Brier skill score
#' against the climatology, MSPE, PPL, overall correlation, and a computation
#' of the coverage probabilities for confidence intervals
#'
#' @param clim the climatology for the location in Y
#'
#'
#' @export
#'
#' 


stEval = function(forecast, Y, clim) {
  
  # ensure Y is in matrix format
  nt = ncol(Y)
  if(is.null(nt)) {
    nt = 1
    Y = matrix(Y, ncol=1)
  }
  
  # categorize Y according to empirical breakpoints
  if(!is.null(forecast$cat.probs)) {
    Y.cat = Y
    for(s in 1:nrow(Y)) {
      Y.cat[s,] = 1 + findInterval(Y[s,], forecast$category.breaks[s,])
    }
  }
  
  
  # evaluation for categorical predictions
  eval.cat = function(pred, obs) {
    
    # compute probability of correct prediction
    pc = mean(pred==obs)
    
    # compute probability of correct prediction (reference forecast)
    pe = 0
    # loop over all categories seen in this timepoint's preds. and obs.
    for(i in union(unique(pred), unique(obs)) ) {
      # add up the probability that a random prediction is correct
      pe = pe + mean(pred==i) * mean(obs==i)
    }
    
    # format return
    data.frame(
      pct.correct = pc,
      heidke.skill = (pc - pe)/(1 - pe)
    )
  }
  
  # initialize overall residual holder
  clim.resid = c()
  
  # compute errors for each timepoint
  for(t in 1:nt) {
    # pull out prediction data frame
    fcst = forecast$pred[[t]]
    
    # evaluate categorical errors
    if(!is.null(forecast$cat.probs)) {
      fcst.cat.eval = eval.cat(fcst$pred$Y.cat, Y.cat[,t])
    }
    
    # compute residual
    fcst$pred$resid = Y[,t] - fcst$pred$Y
    # determine if CI covers observation
    fcst$pred$covered = ifelse( Y[,t]>=fcst$pred$Y.lwr,
                                  ifelse(Y[,t]<=fcst$pred$Y.upr, T, F),
                                  F )
    
    # standard fit measurements
    fcst$err = list(
      mspe = mean(fcst$pred$resid^2),
      ppl = sum(fcst$pred$resid^2)/2 + sum(fcst$pred$se^2),
      r2 = 1 - var(fcst$pred$resid)/var(Y[,t]),
      cor = cor(fcst$pred$Y, Y[,t]),
      coverage = mean(fcst$pred$covered, na.rm = T)
    )
    
    if(!is.null(forecast$cat.probs)) {
      fcst$err$cat.correct = fcst.cat.eval$pct.correct
      fcst$err$cat.heidke = fcst.cat.eval$heidke.skill
    }
    
    # brier skill score if climatology was provided for reference forecast
    if(!is.null(clim)) {
      resid = Y[,t] - clim
      clim.resid = c(clim.resid, resid)
      mspe.ref = mean(resid^2)
      fcst$err$bss = 1 - fcst$err$mspe / mspe.ref
    }
    
    forecast$pred[[t]] = fcst
  }
 
  # collect all predictions and residuals
  Y = as.numeric(Y)
  Y.cat = as.numeric(Y.cat)
  Y.hat = as.numeric(sapply(forecast$pred, function(f) { f$pred$Y }))
  Y.cat.hat = as.numeric(sapply(forecast$pred, function(f) { f$pred$Y.cat }))
  resid = as.numeric(sapply(forecast$pred, function(f) { f$pred$resid }))
  se = as.numeric(sapply(forecast$pred, function(f) { f$pred$se }))
  coverages = as.numeric(sapply(forecast$pred, function(f) { f$pred$covered }))
  
  # add overall error evaluations
  attr(forecast, 'err.mspe') = mean(resid^2)
  attr(forecast, 'err.ppl') = sum(resid^2)/2 + sum(se^2)
  attr(forecast, 'err.r2') = 1 - var(resid)/var(Y)
  attr(forecast, 'err.cor') = cor(Y.hat, Y)
  attr(forecast, 'err.coverage') = mean(coverages, na.rm = T)
  attr(forecast, 'err.bss') = 1 - mean(resid^2) / mean(clim.resid^2)
  if(!is.null(forecast$cat.probs)) {
    fcst.cat.eval = eval.cat(pred = Y.cat.hat, obs = Y.cat)
    attr(forecast, 'err.cat') = fcst.cat.eval$pct.correct
    attr(forecast, 'err.heidke') = fcst.cat.eval$heidke.skill
  }
  
  
  forecast
}