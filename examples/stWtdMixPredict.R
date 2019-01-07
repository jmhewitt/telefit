# load data
data("coprecip")
data("coprecip.fit")

#
# subset data so that demonstration is fast
#

# subsample spatial coordinates
inds.s = sample(nrow(coprecip$coords.s), 10)
inds.r = sample(nrow(coprecip$coords.r), 100)
inds.knots = sample(nrow(coprecip.fit$coords.knots), 10)

# train model on two timepoints
inds.t = 1:2

coprecip.fit$coords.knots = coprecip.fit$coords.knots[inds.knots,]

coprecip.new = coprecip
coprecip.new$coords.s = coprecip.new$coords.s[inds.s,]
coprecip.new$coords.r = coprecip.new$coords.r[inds.r,]
coprecip.new$X = coprecip.new$X[inds.s,, -inds.t, drop = FALSE]
coprecip.new$Y = coprecip.new$Y[inds.s, -inds.t, drop = FALSE]
coprecip.new$Z = coprecip.new$Z[inds.r, -inds.t, drop = FALSE]
coprecip.new$tLabs = coprecip.new$tLabs[-inds.t]

coprecip$coords.s = coprecip$coords.s[inds.s,]
coprecip$coords.r = coprecip$coords.r[inds.r,]
coprecip$X = coprecip$X[inds.s,, inds.t, drop = FALSE]
coprecip$Y = coprecip$Y[inds.s, inds.t, drop = FALSE]
coprecip$Z = coprecip$Z[inds.r, inds.t, drop = FALSE]
coprecip$tLabs = coprecip$tLabs[inds.t]


#
# build weighted mixture
#

fit.wtdMix = stWtdMixPredict(stData = coprecip, stDataNew = coprecip.new,
                             priors = coprecip.fit$priors, 
                             coords.knots = coprecip.fit$coords.knots, 
                             miles = coprecip.fit$miles,
                             w.control = list(method = 'Nelder', maxit = 1))