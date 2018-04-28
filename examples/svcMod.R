# demonstration of svcFit and svcPredict methods for spatially varying 
# coefficient models with remote effects.
#
# this example may take several minutes to run.  additionally, the MCMC sample
# paths for the covariance parameters sigmasq, sigmasqeps, and rho have large
# autocorrelation. this example is mainly for code demonstration purposes only
# since this example does not remedy this issue.

library(fields)
library(mvtnorm)

set.seed(2018)


# helper sampling functions
kronSamp = function(A, B) {
  y = matrix(0, nrow = nrow(A) * nrow(B), ncol=1)
  .Call('_mvrnorm_postKron', PACKAGE = 'telefit', y, solve(A), solve(B), 1, T)
}
invWSamp = function(Psi, n) {
  .Call('_mc2_rinvwishart', PACKAGE = 'telefit', Psi, n)
}


# set key parameters
dims = list(N=200, nt=30, k=2, p=2)
params = list(sigmasq=.2, rho=.3, eps=.5, nu=.5)


# generate parameters and data
coords = matrix( runif(2 * dims$N), ncol = 2 )
X = matrix( rnorm(dims$p * dims$N * dims$nt), ncol = dims$p )
beta = c(-1, .5)
z = matrix( rnorm(dims$k * dims$nt), ncol = dims$nt)
H = maternCov(rdist.earth(coords), scale = params$sigmasq, range = params$rho,
              smoothness = params$nu, nugget = params$sigmasq * params$eps)
Hinv = solve(H)
Tm = matrix(c(.5,.2, .2, .5), ncol=2)/2
theta = kronSamp(Hinv, Tm)


# generate response
xb = X %*% beta
zt = as.numeric(apply(z, 2, function(d) { 
  kronecker(diag(dims$N), t(d)) %*% theta }))
w = kronSamp(diag(dims$nt), H)
y =  xb + zt + w


# view response and components

par(mfrow=c(2,2))
hist(y)
hist(xb)
hist(zt)
hist(w)

quilt.plot(y[1:dims$N], coords, main = 'y')
quilt.plot(xb[1:dims$N], coords, main = 'xb')
quilt.plot(zt[1:dims$N], coords, main = 'zt')
quilt.plot(w[1:dims$N], coords, main = 'w')


# fit model
it = 4e3
priors = list(
  T = list(Psi = .1*diag(dims$k), nu = dims$k),
  beta = list(Linv = diag(dims$p) * 1e-2),
  sigmasq = list(a=2, b=1),
  rho = list(L=0, U=1),
  cov = list(nu=.5)
)

fit = svcFit(y=y, X=X, z=z, coords=coords, priors=priors, nSamples=it)


# check convergence
library(coda)
burn = 1e3/2

par(mfrow=c(1,1))

plot(mcmc(fit$parameters$samples$beta[-(1:burn),]))

plot(theta, colMeans(fit$parameters$samples$theta[-(1:burn),]))
cor(theta, colMeans(fit$parameters$samples$theta[-(1:burn),]))
abline(0,1)

plot(mcmc(fit$parameters$samples$T[-(1:burn),1:2]))
cor(as.numeric(Tm), colMeans(fit$parameters$samples$T[-(1:burn),]))

# very slow sampler for covariance scales
plot(mcmc(fit$parameters$samples$sigmasq[-(1:burn),]))
summary(mcmc(fit$parameters$samples$sigmasq[-(1:burn),]))
effectiveSize(mcmc(fit$parameters$samples$sigmasq[-(1:burn),]))
acf(mcmc(fit$parameters$samples$sigmasq[-(1:burn),]), lag.max = 2e2)

plot(mcmc(fit$parameters$samples$rho[-(1:burn),]))
summary(mcmc(fit$parameters$samples$rho[-(1:burn),]))
HPDinterval(mcmc(fit$parameters$samples$rho[-(1:burn),]))
effectiveSize(mcmc(fit$parameters$samples$rho[-(1:burn),]))

plot(mcmc(fit$parameters$samples$sigmasqeps[-(1:burn),]))
summary(mcmc(fit$parameters$samples$sigmasqeps[-(1:burn),]))
effectiveSize(mcmc(fit$parameters$samples$sigmasqeps[-(1:burn),]))


#
# predict at new timepoints
#

# generate parameters and data
nt0 = dims$nt/2
Xn = matrix( rnorm(dims$p * dims$N * nt0), ncol = dims$p )
zn = matrix( rnorm(dims$k * nt0), ncol = nt0)

# generate response
xbn = Xn %*% beta
ztn = as.numeric(apply(zn, 2, function(d) { 
  kronecker(diag(dims$N), t(d)) %*% theta }))
wn = kronSamp(diag(nt0), H)
yn =  xbn + ztn + wn

# predict responses
pred = svcPredict(fit, Xn, zn)
yh = colMeans(mcmc(pred$samples$y[-(1:burn),]))

# R^2
1 - var(yh - yn)/var(yn)

# coverage
yh.hpd = HPDinterval(mcmc(pred$samples$y[-(1:burn),]))
coverage = numeric(length(yn))
for(i in 1:length(coverage))
  coverage[i] = (yh.hpd[i,1] <= yn[i]) && (yh.hpd[i,2] >= yn[i])
mean(coverage)

# example posterior densities
plot(mcmc(pred$samples$y[-(1:burn),1:4]))

# example posterior densities overlapped with response
par(mfrow=c(2,2))
for(i in 1:4) {
  plot(density(pred$samples$y[-(1:burn),i]), main=i)
  abline(v = yn[i], lty=3)
}

# visually compare point predictions
par(mfrow=c(1,2))
quilt.plot(yn[1:dims$N], coords, main = 'ynew')
quilt.plot(yh[1:dims$N], coords, main = 'yhat')