#' Compute log likelihood for model
#'
#'
#' @export
#' 
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
#' @param returnFullAlphas TRUE to return the teleconnection effects sampled 
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
  
stLL = function( stData, stFit, beta, sigmasq_y, sigmasq_r, sigmasq_eps, rho_y, rho_r,
                 X = stData$X, Y = stData$Y, Z = stData$Z, 
                 coords.s = stData$coords.s, coords.r = stData$coords.r,
                 coords.knots = stFit$coords.knots, miles=T, sigmasq_r_eps) {
  
  n = nrow(coords.s)
  r = nrow(coords.r)
  r_knots = nrow(coords.knots)
  p = dim(X)[2]
  t = dim(X)[3]
  
  Dy = rdist.earth(coords.s, miles=miles)
  Dz_knots = rdist.earth(coords.knots, miles=miles)
  Dz_to_knots = rdist.earth(coords.r, coords.knots, miles=miles)
  
  Z = as.matrix(Z)
  
  # format data
  Yl = matrix(as.numeric(Y), ncol=1)
  
  # format design matrix
  Xl = as.matrix(arrayToLong(X, coords.s, 1)[,-(1:3)])
  
  .Call("_ll", PACKAGE = 'telefit', matrix(Xl, ncol=p), Z, Yl, 
        Dy, Dz_knots, Dz_to_knots, p, n, r, r_knots, t, 
        stFit$priors$cov.s$smoothness, stFit$priors$cov.r$smoothness,
        matrix(beta, ncol=p), sigmasq_y, sigmasq_r, sigmasq_eps, rho_y, rho_r,
        sigmasq_r_eps)
}