#' Compute forecasts based on posterior samples
#'
#' Predict response at new timepoints by drawing samples of the response from
#' the posterior predictive distribution.  Since this requires sampling 
#' teleconnection effects, this method can return estimates of the 
#' teleconnection effects as a by-product.
#'
#'
#' @export
#' 
#' @importFrom coda mcmc HPDinterval
#' @importFrom doMC registerDoMC
#' @import doRNG
#' @importFrom itertools ichunk
#' @import foreach
#' @importFrom fields rdist.earth
#' 
#' @useDynLib telefit
#' 
#' 
#' @param stFit Object with class 'stFit' containing posterior parameter samples
#'  needed to composition sample the teleconnection effects and generate 
#'  posterior predictions. The data needed from stFit need only be manually 
#'  entered if not using a stData object.
#' @param stData Object with class 'stData' containing data needed to fit this 
#'  model. The data need only be manually entered if not using a stData object.
#' @param stDataNew object of class stData that includes information needed for 
#'  making forecasts.  If response data is included, this function will 
#'  automatically run stEval using the empirical climatology as the reference
#'  forecast
#' @param burn number of posterior samples to burn before drawing composition
#'  samples
#' @param prob confidence level for approximate confidence intervals of 
#'  teleconnection effects (only needed if returnAlphas==TRUE)
#' @param ncores Since the teleconnection effects and posterior predictions can 
#'  be sampled in parallel, this parameter lets users specify the number of 
#'  cores to use to draw teleconnection and prediction samples
#' @param returnAlphas TRUE to return the teleconnection effects sampled 
#'  alongside the forecasts.  Note that only basic summary information about the 
#'  teleconnection effects will be returned.
#' @param conf Parameter specifying the HPD level to compute for posterior 
#'  predictive samples
#' @param tLabs Forecast timepoint labels
#' @param X [ns, p, nt] array of design matrices with local covariates
#' @param Y [ns, nt] matrix with response data
#' @param Z [nr, nt] matrix with remote covariates
#' @param Xnew [ns, p, nt0] array of design matrices with local covariates 
#'  at forecast timepoints
#' @param Znew [nr, nt0] matrix with remote covariates at forecast timepoints
#' @param coords.s matrix with coordinates where responses were 
#'  observed (lon, lat)
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param cat.probs vector of probabilities for also returning categorical 
#'  predictions from the posterior prediction samples; NULL otherwise
#'  
  
stPredict = function( stFit, stData, stDataNew, burn = 1, prob = .95, 
                      ncores = 1, conf = .95, tLabs = stDataNew$tLabs,
                      X = stData$X, Y = stData$Y, Z = stData$Z, 
                      Xnew = stDataNew$X, Znew = stDataNew$Z,
                      coords.s = stData$coords.s, coords.r = stData$coords.r,
                      returnAlphas = T, cat.probs = c(1/3, 2/3) ) {
  
  # TODO: add support for localOnly and/or non-varying models
  
  # extract localOnly and varying
  localOnly = stFit$localOnly
  varying = stFit$varying
  
  maxIt = length(stFit$parameters$samples$ll)
  nSamples = length(burn:maxIt)
  
  n = nrow(coords.s)
  r = nrow(coords.r)
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.s, miles=stFit$miles)
  Dz = rdist.earth(coords.r, miles=stFit$miles)
  
  Xl = as.matrix(arrayToLong(X, coords.s, 1)[,-(1:3)])
  Yl = matrix(as.numeric(Y), ncol=1)
  Z = as.matrix(Z)
  
  Xlnew = as.matrix(arrayToLong(Xnew, coords.s, 1)[,-(1:3)])
  Znew = as.matrix(Znew)
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  # estimate chunksize that will minimize number of function calls
  chunkSize = ceiling((maxIt-burn+1)/ncores)

  # draw composition samples
  composition = foreach( inds = ichunk(burn:maxIt, chunkSize = chunkSize),
                         .combine = mergeComposition,
                         .options.multicore=mcoptions ) %dorng% {
    inds = unlist(inds)            
                           
    if(stFit$varying) {
      .Call("_stvcomposition", PACKAGE = 'telefit', p, r, n, t, Xl, Z, Yl, Dy, 
            Dz, stFit$priors$cov.s$smoothness, stFit$priors$cov.r$smoothness, 
            stFit$parameters$samples$beta[inds,], 
            stFit$parameters$samples$sigmasq_y[inds], 
            stFit$parameters$samples$rho_y[inds], 
            stFit$parameters$samples$rho_r[inds], 
            stFit$parameters$samples$sigmasq_r[inds], 
            stFit$parameters$samples$sigmasq_eps[inds], 0, T, 
            Xlnew, Znew, returnAlphas, T)
    } else {
      
    }
  }
  
  # remove unwanted information
  attr(composition, 'rng') = NULL
  
  # post process teleconnection effects
  if(returnAlphas) {
    # remove unwanted information
    composition$alpha$beta = NULL
    composition$alpha$samples = NULL
    
    # convert results away from unnecessary matrix format
    composition$alpha$est = as.numeric(composition$alpha$est)
    composition$alpha$sd = as.numeric(composition$alpha$sd)
    
    # compute approximate intervals, etc.
    composition$alpha$summary = summariseAlpha(composition$alpha, burn, prob, 
                                               coords.s, coords.r)
    
    # remove information redundant with the summary
    composition$alpha$est = NULL
    composition$alpha$sd = NULL
  }
  
  # compute empirical breakpoints at each location to define forecast categories
  if(!is.null(cat.probs)) {
    category.breaks = t(apply(stData$Y, 1, 
                              function(r) { quantile(r, probs = cat.probs)}))
  }
  
  # package results
  nt0 = ncol(composition$forecast$forecasts)
  Y = foreach(t= 1:nt0) %dopar% {
    
    # TODO: move this code somewhere where it can be called outside of sampling
    
    # generate HPD intervals
    forecast.mcmc = mcmc(t(composition$forecast$forecasts[,t,]))
    forecast.hpd = HPDinterval(forecast.mcmc, prob = conf)
    
    if(!is.null(cat.probs)) {
      # build categorical predictions (process by location)
      Y.cat = foreach(s = 1:nrow(composition$forecast$forecasts), 
                      .combine='rbind') %do% {
        # extract posterior samples for specified location and timepoint
        y = composition$forecast$forecasts[s,t,]
        # return label for posterior mode of categories
        1 + as.numeric(names(which.max(table(findInterval(y,category.breaks[s,])))))
      }
    }
    
    pred = data.frame(
      Y = colMeans(forecast.mcmc),
      Y.local = colMeans(t(composition$forecast$local[,t,])),
      Y.remote = colMeans(t(composition$forecast$remote[,t,])),
      se = apply(forecast.mcmc, 2, sd),
      Y.lwr = forecast.hpd[,1],
      Y.upr = forecast.hpd[,2],
      Y.cat = Y.cat
    )
    
    r = list(
      pred = pred,
      yrLab = tLabs[t]
    )
    r
  }
  

  # format return
  ret = list(
    pred = Y,
    samples = composition$forecast,
    alpha = composition$alpha$summary,
    coords.s = coords.s,
    localOnly = localOnly,
    varying = varying,
    tLabs = tLabs,
    Y.lab = stData$Y.lab,
    cat.probs = cat.probs,
    category.breaks = category.breaks
  )
  class(ret) = 'stPredict'
  
  # evaluate performance if response data is given
  if(!is.null(stDataNew$Y)) {
    if(is.null(ncol(stDataNew$Y)))
      ret = stEval(ret, stDataNew$Y, stDataNew$Y)
    else
      ret = stEval(ret, stDataNew$Y, rowMeans(stDataNew$Y))
  }
  
  ret
}