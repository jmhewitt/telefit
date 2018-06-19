#' Compute pointwise correlations for an exploratory teleconnection analysis
#'
#'
#' @export
#'
#' @import foreach
#' @importFrom doMC registerDoMC
#' @importFrom itertools ichunk
#' @importFrom stats cor
#' 
#' @param Y [ny x nt]
#' @param Z [nz x nt]
#' @param coords.s coordinates of locations in Y
#' @param coords.r coordinates of locations in Z
#' @param stData stData object containing data to analyze
#' 
#' @return list with a matrix 'cor' containing correlations.  The columns index
#'  the remote coordinates, while the rows index the local coordinates

teleCor = function( stData = NULL, Y = stData$Y, Z = stData$Z, 
                    coords.s = stData$coords.s, coords.r = stData$coords.r) {
  
  # coerce to a more compact form
  res = list(
    cor = cor(t(Y), t(Z)),
    coords.r = coords.r,
    coords.s = coords.s
  )
  
  class(res) = 'teleCor'
  
  res
}
