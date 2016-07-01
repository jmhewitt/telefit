#' Convenience function for stacking matrices into an array.
#'
#' This function extends the abind function from the abind package.
#'
#' @export
#'
#' @useDynLib telefit
#'
#' 
#' 

test = function() {
  
  r = 12
  n = 100
  
  L = matrix(rnorm(n^2), ncol=n)
  L = L + t(L) + n*diag(n)
  
  print(tail(eigen(L)$values))
  
  R = matrix(rnorm(r^2), ncol=r)
  R = R + t(R) + r*diag(r)
  
  print(tail(eigen(R)$values))
  
  set.seed(2012)
  
  A = .Call("_mvnShort", PACKAGE = 'telefit', L, R)
  
  set.seed(2012)
  
  B = .Call("_mvnLong", PACKAGE = 'telefit', L, R)
    
    
  print(str(A))
  print(str(B))
  
  max(abs(A-B))
  
}