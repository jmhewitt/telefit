#' Compute pointwise correlations for an exploratory teleconnection analysis
#'
#'
#' @export
#'
#' @importFrom foreach foreach "%dopar%" "%do%"
#' @importFrom doMC registerDoMC
#' 
#' @param Y [ny x nt]
#' @param Z [nz x nt]
#' @param coords.local coordinates of locations in Y
#' @param coords.remote coordinates of locations in Z
#' @param wide if TRUE, then the results will be coerced into a matrix
#' 

cor.tel = function( Y, Z, coords.local, coords.remote, ncores=1, 
                    wide=F, wide.label.signif=5 ) {
  
  ny = nrow(Y)
  nz = nrow(Z)
  
  # replace with less memory intensive backend
  registerDoMC(ncores)
  
  # compute pointwise correlation of local variate with each remote location
  res = foreach(y = 1:ny, .combine='rbind') %dopar% {
    foreach(z = 1:nz, .combine='rbind') %do% {
      data.frame( lon.Y = coords.local[y,1], lat.Y = coords.local[y,2],
                  lon.Z = coords.remote[z,1], lat.Z = coords.remote[z,2],
                  cor = cor(Y[y,], Z[z,]) )
    }
  }
  
  # coerce to a matrix, if requested
  if(wide) {
    res = matrix(res$cor, nrow=ny, byrow = T)
    colnames(res) = apply(signif(coords.remote, wide.label.signif), 1, 
                               paste, collapse=', ')
    rownames(res) = apply(signif(coords.local, wide.label.signif), 1, 
                               paste, collapse=', ')
  }
  
  res
}