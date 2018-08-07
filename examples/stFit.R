library(dplyr)
library(foreach)
library(itertools)

set.seed(2018)

data("coprecip")
data("coprecip.fit")
attach(coprecip)

coprecip.fit = stFit(stData = coprecip, priors = coprecip.fit$priors, 
                     maxIt = 10, coords.knots = coprecip.fit$coords.knots)
