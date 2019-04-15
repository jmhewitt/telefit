#' Compute forecasts based on posterior samples
#'
#' Predict response at new timepoints by using weighted mixtures to approximate
#' the posterior predictive distribution.
#'
#'
#' @export
#'
#' @import foreach
#' @importFrom fields rdist.earth
#' @importFrom stats quantile
#' @importFrom bisque wMix wBuild
#'
#' @useDynLib telefit, .registration = TRUE
#'
#'
#' @param stData Object with class 'stData' containing data needed to fit this
#'  model. The data need only be manually entered if not using a stData object.
#' @param stDataNew object of class stData that includes information needed for
#'  making forecasts.  If response data is included, this function will
#'  automatically run stEval using the empirical climatology as the reference
#'  forecast
#' @param priors A list containing parameters for the prior distributions. The
#'  list needs to contain the following values
#'    \describe{
#'      \item{beta}{ list(Lambda=matrix) specifying the prior covariance matrix
#'        for the local effects if varying==F, otherwise
#'        list(Psi=matrix, nu=double) specifying the Inverse wishart prior
#'        distribution for the spatially varying coefficient process if
#'        varying==T. }
#'
#'      \item{cov.s}{ list(smoothness=double, range=c(min, max),
#'        variance=c(shape, rate), nugget=c(shape, rate)) }
#'
#'      \item{cov.r}{ list(smoothness=double, range=c(min, max),
#'        variance=c(shape, rate), nugget=c(shape, rate)) }
#'    }
#' @param coords.knots matrix with coordinates where remote teleconnections
#'  will be based (lon, lat)
#' @param miles TRUE if covariance matrix distances should be in miles, FALSE
#'  for kilometers

#' @param prob confidence level for approximate confidence intervals of
#'  teleconnection effects (only needed if returnAlphas==TRUE)
#' @param ncores Since the teleconnection effects and posterior predictions can
#'  be sampled in parallel, this parameter lets users specify the number of
#'  cores to use to draw teleconnection and prediction samples
#' @param conf Parameter specifying the prediction interval
#' @param tLabs Forecast timepoint labels
#' @param X [ns, p, nt] array of design matrices with local covariates
#' @param Y [ns, nt] matrix with response data
#' @param Z [nr, nt] matrix with remote covariates
#' @param Xnew [ns, p, nt0] array of design matrices with local covariates
#'  at forecast timepoints
#' @param Znew [nr, nt0] matrix with remote covariates at forecast timepoints
#' @param coords.s matrix with coordinates where responses were
#'  observed (lon, lat)
#' @param coords.r matrix with coordinates where remote covariates
#'  were observed (lon, lat)
#' @param cat.probs vector of probabilities for also returning categorical
#'  predictions from the posterior prediction samples; NULL otherwise
#' @param level set the quadrature integration level used by wtdMix
#' @param w.control list of control parameters for the optimization used to 
#'  build the quadrature integration
#'   
#' @example examples/stWtdMixPredict.R
#'
stWtdMixPredict = function( stData, stDataNew, priors, coords.knots, miles,
                            conf = .95,
                            tLabs = stDataNew$tLabs, X = stData$X, Y = stData$Y,
                            Z = stData$Z, Xnew = stDataNew$X,
                            Znew = stDataNew$Z, coords.s = stData$coords.s,
                            coords.r = stData$coords.r,
                            cat.probs = c(1/3, 2/3), level = 2, ncores = 1,
                            w.control = list(method = 'Nelder', maxit = 5e4) ) {

  n = nrow(coords.s)
  r = nrow(coords.r)
  r_knots = coords.knots
  p = dim(X)[2]
  t = dim(X)[3]
  t0 = length(tLabs)

  Dy = rdist.earth(coords.s, miles=miles)
  Dz_knots = rdist.earth(coords.knots, miles=miles)
  Dz_to_knots = rdist.earth(coords.r, coords.knots, miles=miles)

  Z = as.matrix(Z)
  Znew = as.matrix(Znew)

  # format data
  Yl = matrix(as.numeric(Y), ncol=1)

  # format design matrix
  Xl = as.matrix(arrayToLong(X, coords.s, 1)[,-(1:3)])
  Xlnew = as.matrix(arrayToLong(Xnew, coords.s, 1)[,-(1:3)])

  # compute empirical breakpoints at each location to define forecast categories
  if(!is.null(cat.probs)) {
    category.breaks = t(apply(stData$Y, 1,
                              function(r) { quantile(r, probs = cat.probs)}))
  }

  
  # joint conditional likelihood
  ll = function(Y, sigmasq_eps, sigmasq_w, sigmasq_alpha, rho_w, rho_alpha,
                nu_w, nu_alpha, log = TRUE) {
    # Evaluates marginalized log-likelihood for RESP model given covariance
    # parameter values.
    #
    # Parameters:
    # Y - column vector of observations

    #
    # build spatial covariances
    #

    Sigma = maternCov(d = Dy, scale = sigmasq_w, range = rho_w,
                      smoothness = nu_alpha, nugget = sigmasq_w * sigmasq_eps)

    R = maternCov(d = Dz_knots, scale = sigmasq_alpha, range = rho_alpha,
                  smoothness = nu_alpha, nugget = 0)

    cst = maternCov(d = Dz_to_knots, scale = sigmasq_alpha, range = rho_alpha,
                    smoothness = nu_alpha, nugget = 0)

    # decompose covariance matrices
    L.Sigma = t(chol(Sigma))
    L.R = t(chol(R))

    # sequentially build, then decompose temporal covariance matrix
    Cinv = forwardsolve(L.R, t(cst) %*% Z)
    Cinv = t(Cinv) %*% Cinv
    diag(Cinv) = diag(Cinv) + 1
    L.Cinv = t(chol(Cinv))

    # sequentially build, then decompose quadratic form components
    LinvX = forwardsolve.kron(L.Cinv, L.Sigma, Xl)
    L.Wmat = t(chol(solve(priors$beta$Lambda) + t(LinvX) %*% LinvX))

    #
    # log-determinant components
    #

    ldet.CinvSigma = n * 2 * sum(log(diag(L.Cinv))) +
                     t * 2 * sum(log(diag(L.Sigma)))

    ldet.Delta = sum(log(diag(priors$beta$Lambda)))

    ldet.Covariates = 2 * sum(log(diag(L.Wmat)))

    ldet = ldet.CinvSigma + ldet.Delta + ldet.Covariates

    #
    # quadratic form components
    #

    u = forwardsolve.kron(L.Cinv, L.Sigma, Y)
    v = t(LinvX) %*% u
    x = forwardsolve(l = L.Wmat, x = v)
    qform = t(u) %*% u - t(x) %*% x

    # assemble and return log-likelihood
    res = - .5 * (n*t*log(2*pi) + ldet + qform)
    if(log) { res } else { exp(res) }
  }

  dinvgamma = function(x, a, b, log = TRUE) {
    res = a * log(b) - lgamma(a) - (a+1) * log(x) - b/x
    if(log) { res } else { exp(res) }
  }

  post.joint = function(params, log = TRUE, ...) {
    # Parameters
    #  params - c(sigmasq_eps, sigmasq_w, sigmasq_alpha, rho_w, rho_alpha)

    if(!is.matrix(params)) {
      params = matrix(params, ncol = 5)
    }

    res = apply(params, 1, function(params) {
      # log-likelihood
      ll(Y = Yl, sigmasq_eps = params[1], sigmasq_w = params[2],
         sigmasq_alpha = params[3], rho_w = params[4], rho_alpha = params[5],
         nu_w = priors$cov.s$smoothness, nu_alpha = priors$cov.r$smoothness,
         log = TRUE) -
      # uniform range priors
      log(diff(priors$cov.s$range)) -
      log(diff(priors$cov.r$range)) +
      # inv-gamma scale priors
      dinvgamma(x = params[1], a = priors$cov.s$nugget[1],
                b = priors$cov.s$nugget[2], log = TRUE) +
      dinvgamma(x = params[2], a = priors$cov.s$variance[1],
                b = priors$cov.s$variance[2], log = TRUE) +
      dinvgamma(x = params[3], a = priors$cov.r$variance[1],
                b = priors$cov.r$variance[2], log = TRUE)
    })

    if(log) { res } else { exp(res) }
  }

  y0.mix = function(y0, params, x.ind, t.ind, log = FALSE, ...) {
    # mixture component for the posterior for an individual x0; note that the
    # parameters vector contains the posterior distribution parameters for all
    # y0
    #
    # Parameters:
    #  params - c(mu_cond, sds_cond)
    #  x.ind - index of spatial coordinate
    #  t.ind - index of temporal coordinate

    # extract posterior parameters
    ind = (t.ind-1)*n + x.ind
    mu_cond = params[ind]
    sds_cond = params[n*t0 + ind]

    dnorm(y0, mean = mu_cond, sd = sds_cond, log = log)
  }

  y0.precompute = function(params, ...) {
    # Parameters
    #  params - c(sigmasq_eps, sigmasq_w, sigmasq_alpha, rho_w, rho_alpha)

    # Note, due to potential memory constraints, this method only returns
    # marginal means and variances for f(y0 | y, params).

    # extract parameters
    sigmasq_eps = params[1]
    sigmasq_w = params[2]
    sigmasq_alpha = params[3]
    rho_w = params[4]
    rho_alpha = params[5]
    nu_w = priors$cov.s$smoothness
    nu_alpha = priors$cov.r$smoothness

    #
    # build spatial covariances
    #

    Sigma = maternCov(d = Dy, scale = sigmasq_w, range = rho_w,
                      smoothness = nu_alpha, nugget = sigmasq_w * sigmasq_eps)

    R = maternCov(d = Dz_knots, scale = sigmasq_alpha, range = rho_alpha,
                  smoothness = nu_alpha, nugget = 0)

    cst = maternCov(d = Dz_to_knots, scale = sigmasq_alpha, range = rho_alpha,
                    smoothness = nu_alpha, nugget = 0)

    # decompose covariance matrices
    L.Sigma = t(chol(Sigma))
    L.R = t(chol(R))

    # sequentially build complete temporal covariance matrix
    Cinv = forwardsolve(L.R, t(cst) %*% cbind(Znew, Z))
    Cinv = t(Cinv) %*% Cinv
    diag(Cinv) = diag(Cinv) + 1
    L.C11inv = t(chol(Cinv[-(1:t0), -(1:t0)]))

    # sequentially build, then decompose additional components
    L11invX11 = forwardsolve.kron(L.C11inv, L.Sigma, Xl)
    L.W11mat = t(chol(solve(priors$beta$Lambda) + t(L11invX11) %*% L11invX11))

    # compute conditional means and s.d.'s at each prediction timepoint
    cond.mu = matrix(NA, nrow = n, ncol = t0)
    cond.sd = matrix(NA, nrow = n, ncol = t0)
    for(i in 1:t0) {
      # conditioning components
      C10inv = Cinv[-(1:t0), i]
      M10 = Xl %*% priors$beta$Lambda %*% t(Xlnew[(i-1)*n + 1:n,])
      V = forwardsolve.kron(L.C11inv, L.Sigma, kronecker(C10inv, Sigma) + M10 )
      u = forwardsolve.kron(L.C11inv, L.Sigma, Yl)
      B = forwardsolve(l = L.W11mat, x = t(L11invX11) %*% V)

      # conditional means
      cond.mu[,i] = t(V) %*% u - t(B) %*% forwardsolve(l = L.W11mat,
                                                       x = t(L11invX11) %*% u)

      # conditional s.d.'s
      C00inv = Cinv[i, i]
      M00 = Xlnew[(i-1)*n + 1:n,] %*%
            priors$beta$Lambda %*%
            t(Xlnew[(i-1)*n + 1:n,])
      cond.sd[,i] = sqrt(diag(
        kronecker(C00inv, Sigma) + M00 - (t(V) %*% V - t(B) %*% B)
      ))
    }

    # return concatenated means, s.d.'s
    c(as.numeric(cond.mu), as.numeric(cond.sd))
  }

  # build weighted mixture of posteriors
  
  w.y0 = wBuild(
    f = post.joint, 
    init = c(priors$cov.s$nugget[2] / (priors$cov.s$nugget[2]  + 1),
             priors$cov.s$variance[2] / (priors$cov.s$variance[2]  + 1),
             priors$cov.r$variance[2] / (priors$cov.r$variance[2]  + 1),
             mean(priors$cov.s$range),
             mean(priors$cov.r$range)), 
    approx = 'gaussian', 
    link = c(rep('log', 3), rep('logit', 2)), 
    link.params = list(NA, NA, NA, priors$cov.s$range, priors$cov.r$range), 
    optim.control = w.control
  )
  
  post.y0 = wMix(f1 = y0.mix, f2 = post.joint, w = w.y0, 
                 f1.precompute = y0.precompute, level = level, ncores = ncores)

  y0.F.mix = function(y0, mix, wts, x.ind, t.ind, lower.tail = TRUE) {
    # given a set of mixture objects, compute the marginal posterior CDF
    # for y0(ind) | y.  the mixture parameters provide information for y0 | y
    # at all locations, but this function evaluates the CDF at a single
    # location and timepoint, specified by ind.
    #
    # Parameters:
    #  y0 - vector of values at which to evaluate CDF
    #  mix - matrix of mixture parameters
    #  wts - weight of each mixture component
    #  x.ind - index of spatial coordinate
    #  t.ind - index of temporal coordinate

    # extract posterior parameters
    ind = (t.ind-1)*n + x.ind
    mu_cond = mix[, ind]
    sds_cond = mix[, n*t0 + ind]

    res = numeric(length(y0))
    for(i in 1:length(mu_cond)) {
      res = res + pnorm(q = y0, mean = mu_cond[i], sd = sds_cond[i],
                        lower.tail = lower.tail) * wts[i]
    }
    res
  }

  y0.Q.mix = function(y0, x.ind, t.ind, lower.tail = TRUE) {
    # given a set of mixture objects, compute the marginal posterior quantile
    # for y0(ind) | y.  the mixture parameters provide information for y0 | y
    # at all locations, but this function evaluates the quantile at a single
    # location and timepoint, specified by ind.
    #
    # Parameters:
    #  y0 - vector of probabilities for which quantile should be computed
    #  x.ind - index of spatial coordinate
    #  t.ind - index of temporal coordinate

    if(!lower.tail) { y0 = 1 - y0 }

    sapply(y0, function(y0) {
      o = optim(0, function(q) {
        abs(y0.F.mix(y0 = q, mix = post.y0$mix, wts = post.y0$wts,
          x.ind = x.ind, t.ind = t.ind, lower.tail = TRUE) - y0)
      }, method = 'CG')
      o$par
    })
  }

  # get quantiles for equal-tailed credible interval
  conf.q = (1-conf)/2
  conf.q = c(0, conf) + conf.q

  # package results
  Y = foreach(t = 1:t0, .combine='c',
              .export = c('cat.probs', 'category.breaks', 'tLabs')) %do% {

    # generate 95% posterior intervals
    forecast.intervals = matrix(0, nrow = n, ncol = 2)
    for(i in 1:n) {
      forecast.intervals[i,] = y0.Q.mix(y0 = conf.q,
                                        x.ind = i, t.ind = t)
    }

    if(!is.null(cat.probs)) {
      # build categorical predictive distribution (process by location)
      pred.cat = matrix(0, nrow = n, ncol = length(cat.probs) + 1)
      for(i in 1:n) {
        pred.cat[i,1:2] = y0.F.mix(y0 = category.breaks[i,], mix = post.y0$mix,
                                wts = post.y0$wts, x.ind = i, t.ind = t)
        pred.cat[i,2:3] = c(diff(pred.cat[i,1:2]), 1 - pred.cat[i,2])
      }

      # extract categorical predictions (process by location)
      Y.cat = apply(pred.cat, 1, which.max)
    }

    # posterior means
    cond.mu = colSums(sweep(post.y0$mix[, (t-1)*n + 1:n], 1,
                            post.y0$wts, FUN = '*'))

    # posterior sds
    nt0 = n * t0
    cond.sd = sapply(1:n, function(i){
      sqrt(
        # variance of expectation
        sum((post.y0$mix[, (t-1)*n + i] - cond.mu[i])^2 * post.y0$wts) +
        # expectation of variance
        sum(post.y0$mix[, nt0 + (t-1)*n + i]^2 * post.y0$wts)
      )
    })

    pred = data.frame(
      Y = cond.mu,
      se = cond.sd,
      Y.lwr = forecast.intervals[,1],
      Y.upr = forecast.intervals[,2],
      Y.cat = Y.cat
    )

    r = list(
      pred = pred,
      pred.cat = pred.cat,
      yrLab = tLabs[t]
    )
    list(r)
  }

  # format return
  ret = list(
    pred = Y,
    coords.s = coords.s,
    tLabs = tLabs,
    Y.lab = stData$Y.lab,
    cat.probs = cat.probs,
    category.breaks = category.breaks
  )

  class(ret) = 'stPredict'

  # evaluate performance if response data is given
  if(!is.null(stDataNew$Y)) {
    if(is.null(ncol(stDataNew$Y)))
      ret = stEval(ret, stDataNew$Y, stDataNew$Y)
    else
      ret = stEval(ret, stDataNew$Y, rowMeans(stDataNew$Y))
  }

  ret
}
