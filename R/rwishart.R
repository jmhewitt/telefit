#' Random wishart matrix
#'
#'
#' @useDynLib telefit
#'
#' @param V symmetric positive definite p x p scale matrix
#' @param n degrees of freedom (greater than p-1)
#'
#' 

rwishart = function(V, n) {
  .Call("_rwishart", PACKAGE = 'telefit', V, n)
}