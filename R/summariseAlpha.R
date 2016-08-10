#' Summarize alphas
#'
#' This function computes approximate normal intervals, etc. for fitted alphas.
#'
#' @export
#' 
#' 
#'
#' 
#' 

summariseAlpha = function( alpha, burn, prob=.95, coords.s, coords.r ) {
  
  n = nrow(coords.s)
  r = nrow(coords.r)
  
  # compute approximate normal intervals and significance
  z = qnorm((1-prob)/2, lower.tail = F)
  lower = alpha$est - alpha$sd*z
  upper = alpha$est + alpha$sd*z
  sig = apply(cbind(lower, upper), 1, function(x) {
    ifelse(prod(sign(x))==-1, F, T)
  })
  
  # report basic stats for alphas
  message(sum(sig), ' (', signif(mean(sig)*100, 3), 
          '%) significant teleconnection effects (', prob*100, 
          '% approximate normal intervals)')
  
  # compile and return information in a data frame
  data.frame(
    alpha = alpha$est,
    sd = alpha$sd,
    lower = lower,
    upper = upper,
    signif = sig,
    lon.Z = coords.r[,1],
    lat.Z = coords.r[,2],
    lon.Y = rep(coords.s[,1], rep(r,n)),
    lat.Y = rep(coords.s[,2], rep(r,n))
  )
}