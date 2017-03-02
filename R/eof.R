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


# #' Performs an EOF decomposition of the data
# #'
# #' Based on Wikle's text
# #'
# #' Each column in the amplitude contains the time series for one EOF pattern
# #' Each column in the pattern contains the spatial loadings for one EOF
# #'
# #' @export
# #'
# #' @param Z A spatio-temporal observation matrix of a field for which spatial
# #'  EOFs will be computed.  Each column should have all spatial observations
# #'  for a given timepoint.
# #'
# #' @return A list with the following components:
# #'
# #'  patterns (columns of Phi)
# #'  amplitudes (sc)
# #'  cumvar (cumulative proportion of variance explained by each pattern)
# #'  mu (climatology)
# #'
# 
# eof = function(Z, center=T, scale=T) {
# 
#   # set number of "non-degenerate" eofs to return
#   t = min(dim(Z)) - 1
# 
#   # center and scale data
#   if(center | scale)
#     Z = t(scale(t(Z), center = center, scale = scale))
# 
#   # compute eofs using computational efficiency fact:
#   #  - non-zero eigenvalues of (Z^t)Z are the same as for ZZ^t
#   #  - eigenvectors for either matrix can be computed from the other
#   #  - faster to compute eigen decomposition for smaller of the two matrices
#   if(nrow(Z) < ncol(Z)) {
#     # compute the numerically simpler eigen decomposition
#     Z.e = eigen(Z %*% t(Z))
# 
#     # extract eofs and variances
#     Phi = Z.e$vectors[,1:t]
#     v = Z.e$values[1:t] / (ncol(Z) - 1)
# 
#   } else {
#     # compute the numerically simpler eigen decomposition
#     A = t(Z) %*% Z
#     A.e = eigen(A)
# 
#     # recover eofs and variances
#     Phi = Z %*% A.e$vectors[,1:t]
#     for(i in 1:ncol(Phi))
#       Phi[,i] = Phi[,i] / sqrt( t(A.e$vectors[,i]) %*% A %*% A.e$vectors[,i] )
#     v = A.e$values[1:t] / (ncol(Z) - 1)
#   }
# 
#   # compute principal component time series (a.k.a. amplitudes or scores)
#   sc = t(t(Phi) %*% Z)
#   attr(sc, 'dimnames') = NULL
# 
#   # standardize coefficients according to von storch
#   v.sqrt = sqrt(v)
#   Phi = sweep(Phi, 2, v.sqrt, '*')
#   sc = sweep(sc, 2, v.sqrt, '/')
# 
#   res = list(
#     Phi = Phi,
#     sc = sc,
#     sd = v.sqrt,
#     cumvar = cumsum(v)/sum(v),
#     mu = attr(Z, 'scaled:center'),
#     sigma = attr(Z, 'scaled:scale')
#   )
#   class(res) = 'eof'
# 
#   res
# }