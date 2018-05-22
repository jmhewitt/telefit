context("Custom C++ probability routines")

test_that("Evaluating probabilities for intervals of normal r.v.s", {
  set.seed(2018)
  
  mu = rnorm(1)
  sigma = runif(1)
  
  breaks = seq(0, 1, length.out = 4)
  
  expect_equal(
    as.numeric(.Call('_qintnorm', PACKAGE = 'telefit', 
          qnorm(breaks[-c(1, length(breaks))], mu, sigma), mu, sigma)),
    diff(breaks)
  )
})
