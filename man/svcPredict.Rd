% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svcPredict.R
\name{svcPredict}
\alias{svcPredict}
\title{Make predictions using a fitted varying coefficient model}
\usage{
svcPredict(
  fit,
  Xn = NULL,
  Zn = NULL,
  stData = NULL,
  stDataNew = NULL,
  burn = 0,
  cat.probs = c(1/3, 2/3),
  conf = 0.95
)
}
\arguments{
\item{fit}{svcFit object containing posterior samples}

\item{Xn}{[nr*nt, p] matrix of local covariates at new timepoint}

\item{Zn}{[nr, nt] matrix of remote covariates at new timepoints}

\item{stData}{Object with class 'stData' containing data needed to fit this 
model. The data is used to compute empirical quantiles for making 
categorical predictions.}

\item{stDataNew}{object of class stData that includes information needed for 
making forecasts.}

\item{burn}{number of posterior samples to burn from fit}

\item{cat.probs}{vector of probabilities for also returning categorical 
predictions from the posterior prediction samples; NULL otherwise}

\item{conf}{Parameter specifying the HPD level to compute for posterior 
predictive samples}
}
\description{
Make predictions using a fitted varying coefficient model
}
\examples{
library(fields)
library(mvtnorm)

set.seed(2018)

# set key parameters
dims = list(N=100, nt=3, k=2, p=2)
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
xb = X \%*\% beta
zt = as.numeric(apply(z, 2, function(d) { 
  kronecker(diag(dims$N), t(d)) \%*\% theta }))
w = kronSamp(diag(dims$nt), H)
y =  xb + zt + w


# fit model
it = 100
priors = list(
  T = list(Psi = .1*diag(dims$k), nu = dims$k),
  beta = list(Linv = diag(dims$p) * 1e-2),
  sigmasq = list(a=2, b=1),
  rho = list(L=0, U=1),
  cov = list(nu=.5)
)

fit = svcFit(y=y, X=X, z=z, coords=coords, priors=priors, nSamples=it)


#
# predict at new timepoints
#

# generate parameters and data
nt0 = 3
Xn = matrix( rnorm(dims$p * dims$N * nt0), ncol = dims$p )
zn = matrix( rnorm(dims$k * nt0), ncol = nt0)

# generate response
xbn = Xn \%*\% beta
ztn = as.numeric(apply(zn, 2, function(d) { 
  kronecker(diag(dims$N), t(d)) \%*\% theta }))
wn = kronSamp(diag(nt0), H)
yn =  xbn + ztn + wn

# predict responses
pred = svcPredict(fit, Xn, zn, burn = 50)

}
