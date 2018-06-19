#' Performs an EOF decomposition of the data
#' 
#' Uses the stats::prcomp function to implement EOF decompositions of data
#' 
#' @export
#' 
#' @importFrom stats prcomp
#' 
#' @param X [variable x observation] matrix of data for which to compute EOFs
#' @param center TRUE/FALSE to center columns of X in call to prcomp
#' @param scale TRUE/FALSE to scale columns of X in call to prcomp
#' 
#' @return A list containing EOF patterns as columns, and their scores 
#' 
#' 

eof = function(X, center = F, scale = F) {
  e = prcomp(X, center = center, scale. = scale)
  
  dimnames(e$rotation)[[2]] = 1:ncol(X)
  
  list(patterns = -e$x,
       scores = -e$rotation,
       sd = e$sdev)
}