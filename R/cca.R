#' Compute canonical correlations between multivariate datasets using EOF decompositions
#'
#' @export
#' 
#' @description Based on Glahn (1968) and Cook, et. al. (1994).
#' 
#' @details This will center and scale data for efficiency.
#'
#'
#' @param X An \code{(nvars x nobs)} data frame or matrix in which each column contains 
#'  all observations of measured (predictor) variables for a given timepoint or
#'  sample.  For example, if X represents a spatial variable that was recorded 
#'  at several timepoints, then each row of X should contain the variable's 
#'  measurement for all timepoints at a single location.
#'
#' @param Y An \code{(nresponse x nobs)} data frame or matrix in which each  
#'  column contains all observations of measured (response) variables for a 
#'  given timepoint or sample.  For example, if Y represents a spatial variable 
#'  that  was recorded  at several timepoints, then each row of Y should contain 
#'  the variable's  measurement for all timepoints at a single location.
#'  
#' @param X.k Number of empirical orthogonal functions (i.e., principal 
#'  components) from the X data to use in the analysis.
#'  
#' @param Y.k Number of empirical orthogonal functions (i.e., principal 
#'  components) from the Y data to use in the analysis.
#'  
#'  
#' @return A list containing the following components:
#'  \describe{
#'    \item{\code{X.eof}}{EOF decomposition of X.}
#'    \item{\code{Y.eof}}{EOF decomposition of Y.}
#'    \item{\code{X.k}}{Number of EOFs of X used in analysis.}
#'    \item{\code{Y.k}}{Number of EOFs of Y used in analysis.}
#'    \item{\code{cor}}{Correlations of canonical variables.}
#'    \item{\code{A}}{Canonical patterns for X in EOF-space. Each column 
#'      represents one pattern.}
#'    \item{\code{B}}{Canonical patterns for Y in EOF-space. Each column 
#'      represents one pattern.}
#'    \item{\code{C}}{Time series of canonical amplitudes (a.k.a. scores) for X. 
#'      Each column represents the amplitude of one pattern over time.}
#'    \item{\code{D}}{Time series of canonical amplitudes (a.k.a. scores) for Y. 
#'      Each column represents the amplitude of one pattern over time.}
#'    \item{\code{F.X}}{Canonical patterns for X in the same coordinate space 
#'      as X. Each column represents one pattern.}
#'    \item{\code{F.Y}}{Canonical patterns for Y in the same coordinate space 
#'      as Y. Each column represents one pattern.}
#'  } 
#'  
#'  See Cook, et. al. (1994) for details about objects A, B, C, and D; see 
#'  von Storch and Zwiers (1999) for details about F.X and F.Y.
#'  
#' @seealso \link{eof} For details about EOF decompositions.
#' 
#' 
#' @references Cook, E.R., Briffa, K.R., and Jones, P.D., 1994, Spatial regression methods in dendroclimatology: A review and comparison of two techniques: International Journal of Climatology, v. 14, p. 379–402.
#' 
#' @references Glahn, H.R., 1968, Canonical Correlation and Its Relationship to Discriminant Analysis and Multiple Regression: Journal of the Atmospheric Sciences, v. 25, p. 23–31, doi: 10.1175/1520-0469(1968)025<0023:CCAIRT>2.0.CO;2.
#' 
#' @references Storch, H. Von, and Zwiers, F.W., 1999, Statistical Analysis in Climate Research: Cambridge, Cambrudge University Press.
#' 
#' 
cca = function(X, Y, X.k, Y.k) {
  
    # cca uses correlations between datasets, thus requiring at least one 
    #   equal dimension
    if(ncol(X)!=ncol(Y))
      stop('X and Y do not have same number of columns.')
  
    # data prep
    if(class(X)=='data.frame')
      X = data.matrix(X)
    if(class(Y)=='data.frame')
      Y = data.matrix(Y)
    
    # standardize data
    X = t(scale(t(X), center = T, scale = T))
    Y = t(scale(t(Y), center = T, scale = T))
    
    # compute eofs
    #  - following von Storch and Zwiers (1999), data is centered and scaled 
    #    so that cca analysis is computationally simplified as U and V 
    #    (defined next) will be orthogonalized datasets
    X.eof = eof(X, center = T, scale = T)
    Y.eof = eof(Y, center = T, scale = T)
    
    # truncate data
    U = X.eof$sc[,1:X.k]
    V = Y.eof$sc[,1:Y.k]
    
    
    #
    # the cca analysis structure follows Cook, et. al. (1994), but uses
    #   computational techniques from Glahn (1968)
    #
    
    # change notation
    p = X.k
    q = Y.k
    
    if(p<q)
      warning('Y.k must be less than X.k in order to make predictions.')
    
    # compute intercorrelation matrix
    pRq = cor(U, V)
    
    # compute canonical patterns using computational efficiency facts:
    #  - Glahn (1968) to use one eigendecomposition to get both sets of patterns
    #  - compute the eigendecomposition for the smaller of p, q
    if(p<q) {
      cc = eigen( pRq %*% t(pRq) )
      cors = sqrt(cc$values)
      A = cc$vectors
      B = t(pRq) %*% A %*% diag(1/cors)
    } else {
      cc = eigen( t(pRq) %*% pRq )
      cors = sqrt(cc$values)
      B = cc$vectors
      A = pRq %*% B %*% diag(1/cors)
    }
    
    # compute canonical time series
    C = U %*% A
    D = V %*% B
    
    # follow von Storch and Zwiers (1999) to recover cca patterns in data-space
    F.X = X.eof$Phi[,1:X.k] %*% diag(1/X.eof$sd[1:X.k]) %*% A
    F.Y = Y.eof$Phi[,1:Y.k] %*% diag(1/Y.eof$sd[1:Y.k]) %*% B
    
    res = list(
      X.eof = X.eof, Y.eof = Y.eof, X.k = X.k, Y.k = Y.k, cor = cors,
      A = A, B = B, C = C, D = D, F.X = F.X, F.Y = F.Y,
      scaling = list(X.center = attr(X, 'scaled:center'), 
                     X.scale = attr(X, 'scaled:scale'),
                     Y.center = attr(Y, 'scaled:center'),
                     Y.scale = attr(Y, 'scaled:scale'))
    )
    class(res) = 'cca'

    res
}