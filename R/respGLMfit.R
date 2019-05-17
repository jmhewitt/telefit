#' Fit the remote effects spatial process (RESP) model
#' 
#' @export
#' 
#' @importFrom fields rdist.earth
#' @importFrom methods as
#' @importFrom Matrix rowSums
#' @importFrom stats poisson
#' @importFrom stats glm.fit
#' 
#' @useDynLib telefit, .registration = TRUE
#'
#' @param stData Object with class 'stData' containing data needed to fit this 
#'  model. The data need only be manually entered if not using a stData object.
#' @param X [ns, p, nt] array of design matrices with local covariates
#' @param Y [ns, nt] matrix with response data
#' @param Z [nr, nt] matrix with remote covariates
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param coords.knots matrix with coordinates where remote teleconnections
#'  will be based (lon, lat)
#' @param miles TRUE if covariance matrix distances should be in miles, FALSE 
#'  for kilometers
#' @param Q  the adjacency/neighborhood matrix for the spatial configuration
#' @param nSamples number of posterior samples to draw
#' @param priors list of prior distributions for model parameters.
#'    \describe{
#'      \item{cov}{ list(nu=double, sigmasq=c(a,b), rho=c(L,U), kappa=c(a, b) )}
#'    }
#' @param inits initial values \code{c(kappa, sigmasq, rho)} for covariance 
#'  parameters
#' @param sds initial proposal standard deviations for covariance parameter RW 
#'  samplers
#' @param C scaling constants used to update proposal standard deviations for 
#'  covariance parameter RW samplers.  Set \code{C=rep(0,3)} to not adapt the 
#'  proposal sd's.
#' @param alpha0 matrix of initial parameters for teleconnection effects.  by 
#'  default will be the MLE fits for a marginal poisson regression.
#' @param beta0 vector of initial parameters for regression coefficients. by 
#'  default will be the MLE fits for a marginal poisson regression.
#' @param family the GLM likelihood to fit.  currently only accepts 'poisson'  
#' @param k number of ocean EOFs to use in fitting. this can induce dimension 
#'  reduction
#' @param eta0 (optional) initial values of loaded teleconnection effects.
#' @param thin amount by which the Gibbs sampler will thin its output
#'  
#' @example examples/respGLMfit.R
#' 
respGLMfit = function( stData = NULL, X = stData$X, Y = stData$Y, Z = stData$Z, 
                       coords.r = stData$coords.r, Q = stData$Q, miles = TRUE, 
                       sds = rep(1,3), C = rep(.1,3), alpha0 = NULL, 
                       beta0 = NULL, family = 'poisson', k = ncol(stData$Z),
                       coords.knots, nSamples, priors, inits, eta0 = NULL,
                       thin = 1) {
  
  # verify family
  if(family!='poisson') { 
    stop("Only 'poisson' likelihoods are currently supported")
  }
  
  # check priors
  if(is.null(priors)) { stop('Must specify prior distributions') } 
  else {
    if(is.null(priors$nu)) { 
      stop('Must specify prior covariance: nu') 
    }
    if(is.null(priors$rho)) { 
      stop('Must specify prior parameters for covariance: rho') 
    }
    if(is.null(priors$kappa)) { 
      stop('Must specify prior parameters for covariance: kappa') 
    }
    if(is.null(priors$sigmasq)) { 
      stop('Must specify prior parameters for covariance: sigmasq') 
    }
    if(is.null(priors$beta)) {
      stop("Must specify prior variances for regression coefficients: beta")
    }
    
    # formatting for passing to C++
    priors$rho_L = priors$rho[1]
    priors$rho_U = priors$rho[2]
    priors$kappa_a = priors$kappa[1]
    priors$kappa_b = priors$kappa[2]
    priors$sigmasq_a = priors$sigmasq[1]
    priors$sigmasq_b = priors$sigmasq[2]
  }
  
  # coerce Q to sparse matrix
  if(class(Q)!='dgCMatrix') { Q = as(Q, 'dgCMatrix') }
  
  # compute distances for remote covariance matrices
  dknots = rdist.earth(coords.knots, miles = TRUE)
  dzknots = rdist.earth(coords.r, coords.knots, miles = TRUE)
  
  # compute EOFs for remote covariates
  Z.eof = eof(Z)
  W = Z.eof$patterns[,1:k]
  A = Z.eof$scores[,1:k]
  
  # compute structural form for precision matrix (Rue and Held, 2005, eq. 3.30)
  Q = -1 * Q
  diag(Q) = -rowSums(Q)
  
  # format response
  Yl = matrix(as.numeric(Y), ncol=1)
  
  # format design matrix
  Xl = matrix(NA, nrow = nrow(X) * dim(X)[3], ncol = ncol(X))
  for(i in 1:dim(X)[3]) {
    Xl[(i-1)*nrow(X) + 1:nrow(X),] = X[,,i]
  }
  
  # check initial regression coefficient values
  if(is.null(beta0)) {
    if(family=='poisson') {
      fit.beta = glm.fit(x = Xl, y = Yl, family = poisson())
      beta0 = fit.beta$coefficients
    }
  } else if(length(beta0) != ncol(Xl)) { 
    stop('Initial beta vector must have ncol(X) parameters')
  } else {
    fit.beta = list(
      linear.predictors = Xl %*% beta0
    )
  }
  
  if(is.null(eta0)) {
    
    # check initial teleconnection coefficient values
    if(is.null(alpha0)) {
      if(family=='poisson') {
        alpha0 = matrix(NA, nrow = nrow(Y), ncol = k)
        offsets = matrix(fit.beta$linear.predictors, nrow = nrow(Y))
        for(i in 1:nrow(Y)) {
          fit.alpha = glm.fit(x = A, y = Y[i,], family = poisson(), 
                              offset = offsets[i,])
          alpha0[i,] = fit.alpha$coefficients
        }
      }
    } else if(any(dim(alpha0) != c(nrow(Y), k))) {
      stop('Initial alpha matrix has wrong dimensions')
    }
  
    # transform initial guesses for teleconnection effects
    eta0 = as.numeric(alpha0 %*% t(A))
  } else if (all(!is.null(alpha0), !is.null(eta0))) { 
    warning('Initial alpha matrix is ignored because initial eta0 vector is
            also provided.')
  }
  
  # compute rank deficiency in GMRF structural component
  df = nrow(Q) - qr(Q)$rank
  
  res = .Call(`_telefit_respglm_fit`, dknots, dzknots, W, as.integer(nSamples), 
              priors, inits, sds, C, Q, eta0, beta0, Yl, Xl, t(A), df, thin)
  
  class(res) = 'respGLMfit'
  res
}