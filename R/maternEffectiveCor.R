#' Compute effective range for Matern correlation to drop to a specified level
#'
#'
#' @export
#'
#' @param cor Effective correlation to check for
#' @param range Matern range parameter.  Controls the decay of pointwise 
#'        correlations as a function of distance.
#' @param smoothness Matern smoothness parameter.  Controls the number of 
#'        process derivatives.
#'
#' 
#' 

maternEffectiveCor = function(cor = .05, range = 1,  smoothness = .5) {
  optim(1, function(d) {
    abs(cor - maternArray(d, range = range, smoothness = smoothness))
  }, method='BFGS')$par
}