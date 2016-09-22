#' Compute pointwise correlations for an exploratory teleconnection analysis
#'
#'
#' @export
#'
#' @import foreach
#' @importFrom doMC registerDoMC
#' 
#' @param Y [ny x nt]
#' @param Z [nz x nt]
#' @param coords.s coordinates of locations in Y
#' @param coords.r coordinates of locations in Z
#' 
#' @return list with a matrix 'cor' containing correlations.  The columns index
#'  the remote coordinates, while the rows index the local coordinates

teleCor = function( stData = NULL, Y = stData$Y, Z = stData$Z, 
                    coords.s = stData$coords.s, coords.r = stData$coords.r, 
                    ncores = 1 ) {
  
  ny = nrow(Y)
  nz = nrow(Z)
  
  # replace with less memory intensive backend
  registerDoMC(ncores)
  
  # compute pointwise correlation of local variate with each remote location
  res = foreach(y = 1:ny, .combine='rbind') %dopar% {
    foreach(z = 1:nz, .combine='rbind') %do% {
      data.frame( lon.Y = coords.s[y,1], lat.Y = coords.s[y,2],
                  lon.Z = coords.r[z,1], lat.Z = coords.r[z,2],
                  cor = cor(Y[y,], Z[z,]) )
    }
  }
  
  # coerce to a more compact form
  res = list(
    cor = matrix(res$cor, nrow=ny, byrow = T),
    coords.r = coords.r,
    coords.s = coords.s
  )
  
  class(res) = 'teleCor'
  
  res
}