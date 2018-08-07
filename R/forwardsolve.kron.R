#' Solves a triangular system with a Kronecker product structure
#'
#' Solves \eqn{kron(A, B) x = y} where \eqn{A} and \eqn{B} are lower triangular
#' matrices.
#'
#' @export
#' 
#' @param A an \eqn{m x n} matrix
#' @param B an \eqn{p x q} matrix
#' @param y an \eqn{mp x 1} column vector
#' 
#' @return x
#' 
#' @example examples/kronSolve.R

forwardsolve.kron = function(A, B, y) {
  Y = matrix(y, nrow = nrow(B))
  Xp = forwardsolve(B, Y)
  as.numeric(t(forwardsolve(A, t(Xp))))
}