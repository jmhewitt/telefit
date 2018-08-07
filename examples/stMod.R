\dontrun{
  
# demonstration of core methods for working with the remote effects 
# spatial proces (RESP) model. 
#  
# while this example may take several minutes to run, it only demonstrates use 
# of the model fitting and evaluation functions.  the demonstration does not 
# draw enough posterior samples to be able to conduct inference.

set.seed(2018)
  
library(dplyr)
library(foreach)
library(itertools)

data("coprecip")
attach(coprecip)

# specify prior distributions
priors = list(
  beta = list( Lambda = diag(10, ncol(X[,,1]))), 
  cov.s = list( smoothness = .5, range = c(1, 600), variance = c(2, 30), 
                nugget = c(2, 1) ),
  cov.r = list( smoothness = .5, range = c(1, 600), variance = c(2,1e-1), 
                nugget = c(2,1) )
)


#
# place nearly equally-spaced knots over the remote covariates
#

n.lon = 7
n.lat = 6

lon.unique = unique(coords.r[,1])
lon.range = (180-min(lon.unique[lon.unique>0])) + 
  max(lon.unique[lon.unique<0]) + 180
lon.step = lon.range/(n.lon-1)

lat.unique = unique(coords.r[,2])
lat.range = max(lat.unique[lat.unique>0]) - min(lat.unique[lat.unique<0])
lat.step = lat.range/(n.lat-1)

coords.ideal = expand.grid(
  lon = seq(from = min(lon.unique[lon.unique>0]),
            by = lon.step, length.out = n.lon),
  lat = seq(from = max(lat.unique[lat.unique>0]),
            by = -lat.step, length.out = n.lat)) %>%
  mutate(lon = ifelse(lon>180, lon - 360, lon))

# map ideal knot coordinates to closest available coordinates in data
coords.knots = foreach(r = iter(coords.ideal, by = 'row'), 
                       .combine = 'rbind') %do% {
  coords.r[which.min(fields::rdist.earth(r, coords.r)),]
}
# remove duplicate knots
coords.knots = unique(coords.knots)
rownames(coords.knots) = NULL


# configure sampler
maxIt = 100
burn = 10

# fit RESP model (run sampler)
coprecip.fit = stFit(stData = coprecip, priors = priors, maxIt = maxIt,
            coords.knots = coords.knots)

# fit is a stFit object
class(coprecip.fit)

# extract posterior means and HPD intervals for model parameters
ests = coef(coprecip.fit, burn = burn)
HPDinterval.stFit(coprecip.fit, burn = burn)

# graph posterior distributions and traceplots
plot(coprecip.fit, burn = burn)
plot(coprecip.fit, type='trace')

# compute a profile likelihood for covariance parameters
sseq = seq(.1, .7, length.out = 15)
rseq = seq(100, 300, length.out = 15)
seqs = expand.grid(sigmasq_y = sseq, rho_y = rseq)
n = nrow(seqs)
ll = stLL(stData = coprecip, stFit = coprecip.fit, 
          beta = matrix(ests$beta, nrow = n, ncol = 2, byrow = TRUE), 
          sigmasq_y = seqs$sigmasq_y, sigmasq_r = rep(ests$sigmasq_r, n), 
          sigmasq_eps = rep(ests$sigmasq_eps, n),
          rho_y = seqs$rho_y, rho_r = rep(ests$rho_r, n), 
          sigmasq_r_eps = rep(0, n))

# plot the profile likelihood
contour(sseq, rseq, matrix(ll, length(sseq), length(rseq)), 
        xlab = expression(sigma[w]^2), ylab = expression(rho[w]))

# compute VIFs
stVIF(stData = coprecip, stFit = coprecip.fit, burn = burn)

# sample posterior predictive distributions AND estimate teleconnection effects
coprecip.precict = stPredict(stFit = coprecip.fit, stData = coprecip, 
                             stDataNew = coprecip, burn = burn, 
                             returnFullAlphas = TRUE)

# fcst is a stPredict object
class(coprecip.precict)

# evaluate posterior predictions (includes Brier skill against climatology)
clim = rowMeans(Y)
coprecip.precict = stEval(coprecip.precict, Y, clim)

# quickly look at fit summaries
summary(coprecip.precict)

# view posterior predictions for 1982
plot(coprecip.precict, burn = burn, t = 1982)
plot(coprecip.precict, burn = burn, t = 1982, type = 'se')
plot(coprecip.precict, burn = burn, t = 1982, type = 'cat.pred')
  
# extract posterior means for estimates of EOF teleconnection effects
alpha.eof = coef(coprecip.precict, stData = coprecip, stFit = coprecip.fit, 
                 burn = burn) 

# recompute approximate 90% credible intervals for teleconnection effects
alpha.90 = summariseAlpha(alpha = coprecip.precict$alpha, prob = .9, 
                          coords.s = coords.s, coords.r = coords.r)
alpha.eof.90 = summariseEOFAlpha(eof_alpha = coprecip.precict$eof_alpha_knots, 
                                 prob = .9, coords.s = coords.s)

# view posterior means for the ENSO effect on Colorado precipitation
fields::quilt.plot(coords.s, alpha.eof[,1])

# alternatively, plot.stPredict can do this as well, but can also use the 
# uncertainty estimates to highlight locations with significant effects
plot(coprecip.precict, type='eof_alpha_knots', stFit = coprecip.fit, 
     stData = coprecip, signif.telecon = TRUE, lwd=0.7, alpha=.6)

# compare with significance using 90% intervals
alpha.eof.95 = coprecip.precict$eof_alpha_knots
coprecip.precict$eof_alpha_knots = alpha.eof.90
plot(coprecip.precict, type='eof_alpha_knots', stFit = coprecip.fit, 
     stData = coprecip, signif.telecon = TRUE, lwd=0.7, alpha=.6)


# simulate new responses using covariates and the estimated model parameters
cosim = stSimulate(coprecip, coprecip, coords.knots, 
                   params = list(
                     beta = ests$beta,
                     cov.s = list(smoothness = priors$cov.s$smoothness,
                                  range = ests$rho_y, 
                                  variance = ests$sigmasq_r,
                                  nugget = ests$sigmasq_r_eps),
                     cov.r = list(smoothness = priors$cov.r$smoothness,
                                  range = ests$rho_r,
                                  variance = ests$sigmasq_r,
                                  nugget = 0)
                   ))

# view some of the simulated data
plot(cosim$dat.sim.train, t = 1982)

# simulation output includes simulated teleconnection effects
str(cosim$params$alphaKnots)
str(cosim$params$alpha)

}
