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
  #   Assume R's GLM code is properly implemented.  If we can recover the MLE 
  #   variance and output associated with 0 gradient (i.e., MLE estimates for 
  #   betas), then we can be reasonably confident C++ is coded correctly.
  
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
    max(abs(bhat %*% solve(bcov) - r$b), # taylor-shift of gradient
        abs(r2$b),                       # gradient
        abs(bhess - r$C),                # hessian
        abs(bhess - r2$C)),              # hessian v2
    1e-4
  )
})
