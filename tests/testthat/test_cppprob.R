context("Custom C++ probability routines")

test_that("Evaluating probabilities for intervals of normal r.v.s", {
  set.seed(2018)
  
  mu = rnorm(1)
  sigma = runif(1)
  
  breaks = seq(0, 1, length.out = 4)
  
  expect_equal(
    as.numeric(.Call(`_telefit_r_qintnorm`,
          qnorm(breaks[-c(1, length(breaks))], mu, sigma), mu, sigma)),
    diff(breaks),
    tolerance = 1e-3
  )
})

test_that("GLM likelihood Taylor expansions: betas", {
  
  # Goal: Verify GLM likelihood gradient and Hessians are properly coded in C++
  # 
  # Test concept: gradient is 0 at MLE, and Hessian is related to MLE variance.
  #   Assume R's GLM code is properly implemented for betas.  If we can recover 
  #   the MLE variance and output associated with 0 gradient (i.e., MLE 
  #   estimates for betas), then we can be reasonably confident C++ is coded 
  #   correctly.
  
  # run poisson GLM example from stats::glm
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  d.AD <- data.frame(treatment, outcome, counts)
  glm.D93 <- stats::glm(counts ~ outcome + treatment, family = poisson())
  
  # extract MLE components
  bhat = coef(glm.D93)
  bcov = vcov(glm.D93)
  bhess = -solve(bcov)
  
  # compute gradient/hessian using coefficients
  r = test_taylor_beta(beta0 = bhat, eta0 = rep(0,9), y = counts,
                       x = model.matrix(glm.D93), n = 3, t = 3, p = 5)
  
  # compute gradient/hessian using linear predictors
  r2 = test_taylor_beta(beta0 = rep(0, length(bhat)), 
                        eta0 = glm.D93$linear.predictors, y = counts,
                        x = model.matrix(glm.D93), n = 3, t = 3, p = 5)
  
  expect_lt(
    max(abs(solve(bcov) %*% bhat - r$b), # taylor-shift of gradient
        abs(r2$b),                       # gradient
        abs(bhess + r$C),                # -hessian
        abs(bhess + r2$C)),              # -hessian v2
    1e-4
  )
})


test_that("GLM likelihood Taylor expansions: etas", {
  
  # Goal: Verify GLM likelihood gradient and Hessians are properly coded in C++
  # 
  # Test concept: gradient is 0 at MLE, and Hessian is related to MLE variance.
  #   For the etas, these quantities are easy to see analytically.  WLOG, we 
  #   may assume the betas are equal to 0.  The MLE etas are equal to the log 
  #   counts, and the negative hessian is equal to the counts.
  
  # data for poisson GLM example from stats::glm
  counts <- c(18,17,15,20,10,20,25,13,12)
  
  # compute gradient/hessian using linear predictors
  r = test_taylor_eta0(beta = c(0, 0), eta0 = log(counts), y = counts,
                       x = matrix(0, nrow=length(counts), ncol = 2), n = 3, 
                       t = 3, p = 2)
  
  expect_lt(
    max(abs(counts - r$C),                 # -hessian
        abs( counts * log(counts) - r$b)), # taylor-shift of gradient
    1e-10
  )
})


test_that("GLM posterior approximation: betas", {
  
  # Goal: Verify Gaussian approximation to posterior is properly coded in C++
  # 
  # Test concept: We look for agreement of the Gaussian approximation parameters
  #   as computed with our C++ implementation of a Newton-Raphson approach, and 
  #   the BFGS approach implemented in stats::optim.  
  
  # construct design elements for poisson GLM example from stats::glm
  counts <- c(18,17,15,20,10,20,25,13,12)
  outcome <- gl(3,1,9)
  treatment <- gl(3,3)
  x = model.matrix(counts ~ outcome + treatment)
  
  # set dimensions
  n = 3
  t = 3
  p = ncol(x)
  
  # set prior precision for betas; fix eta0
  Q = diag(1/3,p)
  sds = sqrt(1/diag(Q))
  eta0 = rep(0,n*t)
  
  # compute gaussian approx. at the posterior mode
  o = stats::optim(rep(0,p), function(beta) {
    lambda = exp(x %*% beta + eta0)
    sum(dpois(counts, lambda, log = TRUE)) + 
      sum(dnorm(beta, sd = sds, log = TRUE))
  }, control = list(fnscale=-1), hessian = TRUE, method = 'BFGS')

  # test C++ implementation of Newton-Raphson approach for approximation
  r = test_gaussian_approx_beta(beta0 = rep(0,p), Q = diag(Q), p = p, it = 100, 
                                eta0 = eta0, y = counts, x = x, n = n, 
                                t = t)
  
  expect_lt(
    max(abs(o$par - r$mu),
        abs(t(chol(-o$hessian)) - r$L)),
    1e-5
  )
})
