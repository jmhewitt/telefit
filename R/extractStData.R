#' Basic extraction of SpatialGridDataFrame data for teleconnection analysis
#'
#' @export
#' 
#' @importFrom sp coordinates
#' @importFrom raster extent
#' @importFrom foreach foreach "%do%"
#' 
#' @param X SpatialGridDataFrame with local covariates.  If X is a list, each
#'  SpatialGridDataFrame will be included as one covariate.
#' @param Y SpatialGridDataFrame with response data
#' @param Z SpatialGridDataFrame with remote covariates. If Z is a list, this 
#'  function assumes each element of the list contains observations for the same
#'  covariate, but from different spatial regions.  If Z is a list, D.r and 
#'  mask.r must
#'  also be lists so that this function can know which regions to extract from
#'  each SpatialGridDataFrame
#' @param t Column indexes in X, Y, and Z from which to extract data
#' @param D.s c(xmin, xmax, ymin, ymax) region from which to extract data from 
#'  X and Y
#' @param D.r c(xmin, xmax, ymin, ymax) region from which to extract data from Z
#' @param intercept If TRUE, an intercept will be added to the design matrix
#' @param mask.s SpatialGridDataFrame to be used as a mask when extracting data
#'  from X and Y.  Locations in mask.s with NA values will be ignored when 
#'  extracting data from X and Y.
#' @param mask.r SpatialGridDataFrame to be used as a mask when extracting data
#'  from Z.  Locations in mask.s with NA values will be ignored when 
#'  extracting data from Z.
#' @param type.s 'response' 'anomaly' or 'std.anomaly' depending on whether
#'  data extracted from X and Y should be the observed data, anomalies, or
#'  standardized anomalies (where the climatology is computed from the 
#'  observations as the pointwise temporal average)
#' @param type.r 'response' 'anomaly' or 'std.anomaly' depending on whether
#'  data extracted from Z should be the observed data, anomalies, or
#'  standardized anomalies (where the climatology is computed from the 
#'  observations as the pointwise temporal average)
#'  @param X.lab name for X data (optional)
#'  @param Y.lab name for Y data (optional)
#'  @param Z.lab name for Z data (optional)
#'  

extractStData = function( X, Y, Z, t, D.s, D.r, mask.s = NULL, mask.r = NULL,
                          aggfact.s = NULL, aggfact.r = NULL, intercept = T,
                          type.s = 'response', type.r = 'response',
                          X.lab = NULL, Y.lab = NULL, Z.lab = NULL ) {
  
  # fix things so that t is the names, and not the indexes of data to extract
  
  # convert local bounds to extent object
  D.s = extent(D.s)
  
  # save time labels before they are transformed
  tLabs = names(Y)[t]
  
  
  # extract local data
  
  Y = extractRegion(Y, D.s, type.s, aggfact.s, mask.s)
  
  
  if(class(X)!='list')
    X = list(X)
  
  # extract regions and aggregate local covariates
  for(i in 1:length(X))
    X[[i]] = extractRegion(X[[i]], D.s, type.s, aggfact.s, mask.s)
  
  # build local design matrices for each timepoint
  X.mat = foreach(tt = t, .combine = 'abind3') %do% {
    
    # extract data from each predictor
    x = foreach(x = X, .combine='cbind') %do% { x@data@values[, tt] }
    
    # if requested, add intercept; return data
    if(intercept) {
      cbind(1, x)
    } else
      x
  }
  
  
  # extract remote data
  
  if(class(Z)!='list') {
    Z = list(Z)
    D.r = list(D.r)
    mask.r = list(mask.r)
  }
  
  # extract regions and aggregate remote covariates
  for(i in 1:length(Z))
    Z[[i]] = extractRegion(Z[[i]], D.r[[i]], type.r, aggfact.r, mask.r[[i]])
  
  # combine remote covariate data from each region
  Z.mat = foreach(z = Z, .combine = 'rbind') %do% { matrix(z@data@values[, t], 
                                                           ncol = length(t)) }

  
  # extract coordinates
  
  coords.s = coordinates(Y)
  coords.r = foreach(z = Z, .combine = 'rbind') %do% { coordinates(z) }
  
  # remove remote covariates that have NA data
  complete.data = complete.cases(Z.mat)
  Z.mat = matrix(Z.mat[complete.data,], ncol=length(t))
  coords.r = coords.r[complete.data,]
  
  
  # build return object
    
  res = list(
    tLabs = tLabs,
    coords.s = coords.s,
    coords.r = coords.r,
    X = X.mat,
    Y = Y@data@values[,t],
    Z = Z.mat
  )
  
  res$X.lab = X.lab
  res$Y.lab = Y.lab
  res$Z.lab = Z.lab
  
  class(res) = 'stData'
  
  res
}