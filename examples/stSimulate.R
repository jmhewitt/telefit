data("coprecip")
data("coprecip.fit")

attach(coprecip.fit)

ests = coef(coprecip.fit, burn = 50)

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
