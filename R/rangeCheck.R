#' Determine basic interval over which a covariance range parameter makes numerically
#' invertible matrices
#' 
#' Note that the remote range numbers don't necessarily mean anything since the
#' full model is quite generous.  All the same, it's informative.
#'
#'
#' @export
#'
#' @importFrom fields rdist.earth
#' @useDynLib telefit
#'
#'
#' 
#' 

rangeCheck = function( coords.local, coords.remote, Z, miles = T,
                       smoothness.local = 2, smoothness.remote = 1.76,
                       nugget.local = 0, scale.remote = 1,
                       maxIt = 1000, tol=1e-10 ) {
  
  # extend to also check for numerical symmetry?
  
  # fix this so that it is accurate for the remote range
  
  Dy = rdist.earth(coords.local, miles = miles)
  Dz = rdist.earth(coords.remote, miles = miles)
  
  nt = ncol(Z)
  
  min.offdiag = function(d) {
    # find the minimum non-diagonal entry in a matrix, d
    diag(d) = max(d)
    min(d)
  }
    
  rC = function(d, local=T, minR, maxR) {
    # initialize bisection search
    a = minR
    b = maxR
    
    # find 'a', the largest possible range parameter for which the key matrices
    # are invertible
    converged = F
    for(i in 1:maxIt) {
      # get next point to consider moving to
      nxt = mean(c(a, b))
      
      # compute matrix condition numbers
      if(local) {
        rnum = rcond(maternCov(d, 1, nxt, smoothness.local, nugget.local))
      } else {
        rnum = rcond(t(Z) %*% 
          maternCov(d, scale.remote, nxt, smoothness.remote, 0) %*% Z + 
            diag(1, nt))
      }
      
      # compute update
      if(rnum <= .Machine$double.eps)
        b = nxt
      else
        a = nxt
      
      # check for convergence
      if(b-a < tol) {
        converged = T
        break
      }
      
    }
    
    data.frame(range=ifelse(local, 'local','remote'),
      dist.max = maxR,
      dist.maxinv = a,
      dist.halfmax = maxR/2,
      converged = converged)
  }
  
  rbind( rC(Dy, T, min.offdiag(Dy), max(Dy)),
         rC(Dz, F, min.offdiag(Dz), max(Dz)) )
}