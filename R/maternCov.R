#' Matern covariance
#'
#' This function evaluates the Matern covariance function for the elements of 
#' a (potentially non-square) spatial distance matrix
#'
#' @useDynLib telefit
#' 
#' @export
#'
#' @param d A numeric vector or matrix of distances at which the Matern 
#'        correlation function should be evaluated.
#' @param scale Scales correlations to covariances.
#' @param range Matern range parameter.  Controls the decay of pointwise 
#'        correlations as a function of distance.
#' @param smoothness Matern smoothness parameter.  Controls the number of 
#'        process derivatives.
#' @param nugget Spatial covariance nugget.  
#' 

maternCov = function(d, scale = 1, range = 1, smoothness = .5, nugget = 0) {
  .Call("_maternCov", PACKAGE = 'telefit', as.matrix(d), 
        scale, range, smoothness, nugget)
}