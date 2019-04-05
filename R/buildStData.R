#' Assemble stData object from raw components
#'
#' \code{telefit} fitting and prediction objects rely on \code{stData} objects,
#' which are \code{list} objects with some common structure.  \code{buildStData}
#' is a convenience function for creating the list object that also adds some 
#' basic validation that input objects have the right dimensions, etc.
#' 
#' @export
#' 
#' @param X 3-dimensional array with \eqn{n} rows, \eqn{p} columns, and \eqn{t}
#'   slices. \code{X[i,j,k]} represents the \code{j}th covariate observed at 
#'   time \code{k} at location \code{i}.
#' @param Y matrix with \eqn{n} rows and \eqn{t} columns. \code{Y[i,j]} 
#'   represents the \code{j}th response observed at location \code{i}.
#' @param tLabs vector of labels for each timepoint
#' @param coords.s (optional) matrix with coordinates where responses were 
#'  observed (lon, lat).  If \code{coords.s} is not provided, then \code{Q}
#'  must be provided
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param Q if spatial locations for \code{X} and \code{Y} represent areal 
#'  units, then \code{Q} represents the adjacency matrix for the spatial 
#'  configuration
#' @param X.lab name for X data (optional)
#' @param Y.lab name for Y data (optional)
#' @param Z.lab name for Z data (optional)
#'  
#' @example examples/assemble.R

buildStData = function(X, Y, Z, tLabs, coords.r, coords.s=NULL, Q=NULL,
                       X.lab='', Y.lab='', Z.lab='') {
    
  if(!is.array(X)) { stop('X must be an array') }
  if(length(dim(X)) != 3) { stop('X array must have 3 array dimensions.') }
  
  # extract dimensions that other structures should match
  n = nrow(X)[1]
  t = dim(X)[3]
  
  if(!is.matrix(Y)) { stop('Y must be a matrix') }
  if(nrow(Y) != n) { stop('X and Y must have the same number of rows') }
  if(ncol(Y) != t) { stop('X and Y must have the same number timepoints') }
  
  if(!is.matrix(Z)) { stop('Z must be a matrix') }
  if(ncol(Z) != t) { stop('X and Z must have the same number timepoints') }
  
  # extract more dimensions that other structures should share
  nz = nrow(Z)
  
  if(length(tLabs) != t) { 
    stop('tLabs must describe same number of timepoints as X') 
  }
  
  if(!is.matrix(coords.r)) { stop('coords.r must be a matrix') }
  if(nrow(coords.r) != nz) { 
    stop('coords.r and Z must have same number of rows')
  }
  
  
  if(all(is.null(coords.s), is.null(Q))) { 
    stop('Spatial information about observation locations for X,Y must be 
         provided via coords.s or Q')
  }
  
  if(!is.null(coords.s)) {
    if(!is.matrix(coords.s)) { stop('coords.s must be a matrix') }
    if(nrow(coords.s) != n) { 
      stop('coords.s and X must have same number of rows')
    }
  }
  
  if(!is.null(Q)) {
    if(!is.matrix(Q)) { stop('Q must be a matrix') }
    if(!all(dim(Q)==n)) { stop('Q must be a square matrix with same number of 
                               rows as X')}
  }
    
  res = list(tLabs = tLabs, coords.s = coords.s, coords.r = coords.r, Q = Q,
    X = X, Y = Y, Z = Z, X.lab = X.lab, Y.lab = Y.lab, Z.lab = Z.lab
  )
  
  if(is.null(Q)) { res$Q = NULL }
  if(is.null(coords.s)) { res$coords.s = NULL }
  
  class(res) = 'stData'
  
  res
}
