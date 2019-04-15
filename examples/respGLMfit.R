Rcpp::compileAttributes('r/temporal_LSigma/packages/telefit/')
devtools::load_all('r/temporal_LSigma/packages/telefit/')
devtools::document('r/temporal_LSigma/packages/telefit/')
devtools::test('r/temporal_LSigma/packages/telefit/')


data("caprecip.counts")
data("coprecip.fit")

# scale watershed areas for stability
caprecip.counts$X[,3,] = scale(caprecip.counts$X[,3,1])

devtools::load_all('r/temporal_LSigma/packages/telefit/')

# devtools::load_all('.')

# library(telefit)

data("caprecip.counts")
data("coprecip.fit")

# scale watershed areas for stability
caprecip.counts$X[,3,] = scale(caprecip.counts$X[,3,1])
fit = respGLMfit( stData = caprecip.counts, 
                  coords.knots = coprecip.fit$coords.knots,
                  priors = list(
                    nu = .5, sigmasq = c(2,.1), rho = c(1,600), 
                    kappa = c(2,1/30)
                  ),
                  nSamples = 1e2,
                  inits = c(1e1, 1, 500),
                  sds = c(.05,.05,.05),
                  C = rep(0,3)
                )
  
library(coda)

1-rejectionRate(mcmc(fit$kappa))

plot(fit$ll, type='l')
plot(fit$kappa, type='l')
plot(fit$sigmasq, type='l')
plot(fit$rho, type='l')

apply(fit$eta0,1,function(x) (sd(x)))
