#' Performs an EOF decomposition of the data
#' 
#' @export
#' 
#' @param X [variable x observation] matrix of data for which to compute EOFs
#' 
#' @return A list containing EOF patterns as columns, and their scores 
#' 

eof = function(X, center = F, scale = F) {
  e = prcomp(X, center = center, scale. = scale)
  
  dimnames(e$rotation)[[2]] = 1:ncol(X)
  
  list(patterns = -e$x,
       scores = -e$rotation,
       sd = e$sdev)
}