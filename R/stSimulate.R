#' Compute forecasts based on posterior samples
#'
#'
#'
#' @export
#' 
#' @importFrom doMC registerDoMC
#' @import doRNG
#' @importFrom foreach foreach
#' @importFrom fields rdist.earth
#' @importFrom mvtnorm rmvnorm
#' 
#'
#' 
#' 

stSimulate = function(stFit, coords.local, X, Z, yearLabs, nu_y, miles=T, 
                     ncores=1) {
  
  registerDoMC(ncores)
  mcoptions = list(preschedule=FALSE)
  
  n = nrow(coords.local)  
  nt = ncol(Z)
  Dy = rdist.earth(coords.local, miles=miles)
  
  maxIt = length(stFit$ll)
  burn = maxIt - nrow(stFit$alpha$samples) + 1
  
  I.ns = diag(n)
  
  # process each timepoint
  Y = foreach(t=1:nt, .combine='c') %do% {
    
    # extract covariates for this timepoint
    x = X[,,t]
    z = Z[,t]
    
    # draw posterior samples
    y = foreach(it=burn:(maxIt-1), .combine='c', 
                .options.multicore=mcoptions) %dorng% {
      
      # compute posterior sample mean and variance
      xb = x %*% stFit$beta[it,]
      za = kronecker(I.ns, t(z)) %*% stFit$alpha$samples[it-burn+1,]
      
      # create posterior sample object
      r = list(
        local = t(xb),
        remote = t(za),
        y = rmvnorm( 1, mean = xb + za,
                     sigma = maternCov(Dy, stFit$sigmasq_y[it],
                                       stFit$rho_y[it], nu_y, 
                                       stFit$sigmasq_y[it]*stFit$sigmasq_eps[it]) )
      )
      
      # remote unwanted information
      attr(r, 'rng') = NULL
      
      # return posterior sample
      list(r)
    }
    
    # remote unwanted information
    attr(y, 'rng') = NULL
    
    # extract data into data frames
    y.sample = foreach(r=y, .combine='rbind', 
                       .options.multicore=mcoptions) %dopar% { r$y }
    local.sample = foreach(r=y, .combine='rbind',
                           .options.multicore=mcoptions) %dopar% { r$local }
    remote.sample = foreach(r=y, .combine='rbind', 
                            .options.multicore=mcoptions) %dopar% { r$remote }
    
    # compute means and variances; put data into plottable frame
    pred = data.frame(
      y = colMeans(y.sample),
      y.local = colMeans(local.sample),
      y.remote = colMeans(remote.sample),
      se = apply(y.sample, 2, sd),
      lon = coords.local[,1],
      lat = coords.local[,2]
    )
    
    # return posterior samples for this timepoint
    list(list(
      samples = list(
        y = y.sample,
        local = local.sample,
        remote = remote.sample
      ),
      pred = pred,
      yrLab = yearLabs[t]
    ))
  }
  
 Y 
}