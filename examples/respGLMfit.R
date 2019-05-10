# Rcpp::compileAttributes('r/temporal_LSigma/packages/telefit/')
# devtools::load_all('r/temporal_LSigma/packages/telefit/')
# devtools::document('r/temporal_LSigma/packages/telefit/')
# devtools::test('r/temporal_LSigma/packages/telefit/')
# 
# 
# data("caprecip.counts")
# data("coprecip.fit")
# 
# # scale watershed areas for stability
# caprecip.counts$X[,3,] = scale(caprecip.counts$X[,3,1])
# 
# devtools::load_all('r/temporal_LSigma/packages/telefit/')
# 
# # scale watershed areas for stability
# fit = respGLMfit( stData = caprecip.counts, 
#                   coords.knots = coprecip.fit$coords.knots,
#                   priors = list(
#                     nu = .5, sigmasq = c(2,.1), rho = c(1,600), 
#                     kappa = c(2,1/30), beta = rep(3, ncol(caprecip.counts$X))
#                   ),
#                   nSamples = 1e2,
#                   inits = c(1e1, 1, 500),
#                   sds = c(kappa = .5, ss = .25, rho = 1) * 0,
#                   C = rep(0e4,3)
#                 )
#   
# library(coda)
# 
# 1-rejectionRate(mcmc(fit$kappa))
# 
# plot(fit$ll, type='l')
# plot(fit$kappa, type='l')
# plot(fit$sigmasq, type='l')
# plot(fit$rho, type='l')
# plot(fit$beta[,1], type='l')
# plot(fit$beta[,2], type='l')
# plot(fit$beta[,3], type='l')
# 
# apply(fit$eta0,1,function(x) (sd(x)))
