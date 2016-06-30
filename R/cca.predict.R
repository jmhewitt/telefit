#' Make predictions using canonical correlation analysis results
#'
#' @export
#'
#' @param obj \code{cca} object containing results of a cca analysis.
#' 
#' @param X An \code{(nvars x nobs)} data frame or matrix in which each column contains 
#'  all observations of measured (predictor) variables for a given timepoint or
#'  sample.  For example, if X represents a spatial variable that was recorded 
#'  at several timepoints, then each row of X should contain the variable's 
#'  measurement for all timepoints at a single location.
#'  
#'  
#' @return Y Predicted values for Y.
#'
#' @seealso \link{cca} For details about canonical correlation analyses.
#' 
#' 
#' @references Cook, E.R., Briffa, K.R., and Jones, P.D., 1994, Spatial regression methods in dendroclimatology: A review and comparison of two techniques: International Journal of Climatology, v. 14, p. 379–402.
#' 
#' 
#' 

cca.predict = function(obj, X) {
  
  if(class(obj)!='cca')
    stop('obj is not of class "cca".  
         cca output is required to make predictions.')
  
  if(obj$X.k < obj$Y.k)
    stop('Y.k must be less than X.k in order to make predictions.')
  
  # data prep
  if(class(X)=='data.frame')
    X = data.matrix(X)
  
  # standardize data since we perform cca on a standardized scale
  X = sweep(sweep(X, 1, obj$scaling$X.center, '-'),
            1, obj$scaling$X.scale, '/')
  
  # predict (standardized) Y following Cook, et. al. (1994)  
  Y.std = t(t(X) %*% obj$X.eof$Phi[,1:obj$X.k] %*% obj$A %*% diag(obj$cor^2) %*% 
    solve(obj$B) %*% t(obj$Y.eof$Phi[,1:obj$Y.k]))
  
  # return Y to original data scale
  sweep(sweep(Y.std, 1, obj$scaling$Y.scale, '*'),
        1, obj$scaling$Y.center, '+')
}