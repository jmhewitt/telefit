\dontrun{
data("coprecip")
data("coprecip.fit")
attach(coprecip)

# sample posterior predictive distributions AND estimate teleconnection effects
coprecip.precict = stPredict(stFit = coprecip.fit, stData = coprecip, 
                             stDataNew = coprecip, burn = 90, 
                             returnFullAlphas = TRUE)

alpha.90 = summariseAlpha(alpha = coprecip.precict$alpha, prob = .9, 
                          coords.s = coords.s, coords.r = coords.r)
}
