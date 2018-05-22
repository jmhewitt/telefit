#' Evaluate kron(A,B) * C without storing kron(A,B)
#' 
#' Evaluate kron(A,B) * C without storing kron(A,B)
#' 
#' @useDynLib telefit
#' 
#' @param A (m x n) matrix
#' @param B (p x q) matrix
#' @param C (nq x r) matrix
#'
#' 

dgemkmm = function(A, B, C) {
  
  .Call("_dgemkmm", PACKAGE = 'telefit', 
        as.matrix(A), as.matrix(B), as.matrix(C))
}