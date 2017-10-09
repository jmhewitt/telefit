#' Compute pointwise correlations for an exploratory teleconnection analysis
#'
#'
#' @export
#'
#' @import foreach
#' @importFrom doMC registerDoMC
#' @importFrom itertools ichunk
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
  
  # set up basic parallel backend if none is registered
  if(!getDoParRegistered()) {
    registerDoMC(ncores)
  }
  
  
  chunkSize = ceiling(nt0/ncores)
  Y = foreach(inds = ichunk(1:nt0, chunkSize = chunkSize), .combine='c',
              .export = c('cat.probs', 'category.breaks', 'composition', 
                          'tLabs')) %dopar% {
                            foreach(t = unlist(inds)) %do% {
                              
                              
  # compute pointwise correlation of local variate with each remote location
  # (use chunkSize to exactly split tasks across ncores)
  chunkSize = ceiling(ny/ncores)
  res = foreach(inds = ichunk(1:ny, chunkSize = chunkSize), .combine = 'rbind',
                .export = c('coords.s', 'coords.r', 'Y', 'Z')) %dopar% {
    foreach(y = unlist(inds), .combine='rbind') %do% {
    foreach(z = 1:nz, .combine='rbind') %do% {
      data.frame( lon.Y = coords.s[y,1], lat.Y = coords.s[y,2],
                  lon.Z = coords.r[z,1], lat.Z = coords.r[z,2],
                  cor = cor(Y[y,], Z[z,]) )
    }
  }}
  
  # coerce to a more compact form
  res = list(
    cor = matrix(res$cor, nrow=ny, byrow = T),
    coords.r = coords.r,
    coords.s = coords.s
  )
  
  class(res) = 'teleCor'
  
  res
}