#' Extract region from a SpatialGridDataFrame
#'
#' @return a modified SpatialGridDataFrame, sgdf, with the climatology for each
#'  location accessible via attr(sgdf@data@values, 'scaled:center') if anomalies
#'  were computed
#'  
#' @examples # Show an example where you plot the climatology grid too!
#'
#'
#' @export
#' 
#' @importFrom raster brick mask aggregate crop
#' @importFrom SDMTools aspect slope
#' 
#' @param type whether to return the raw data, anomalies (data minus temporal 
#' average at each location), or standardized anomalies (anomalies divided by
#' temporal standard deviation at each location)
#' @param extent raster::extent object featuring region to extract
#' @param aggfact if provided, will spatially average the data
#' @param mask if an sgdf is provided, the data will be masked before
#'  extraction, aggregation, and anomaly computation
#'  @param aspect TRUE to return the aspect of the surface at each location 
#'   instead of the value of the surface itself
#'  @param aspect.categories if aspect==TRUE, this specifies the number of 
#'   discrete categories to divide aspect numbers (0-360) into.  NULL if the
#'   original scale (0-360) should be kept. By design, the aspect categories
#'   will be centered on north in the first category.
#' 


extractRegion = function(sgdf, extent, 
                         type='response', 
                         aggfact=NULL, mask=NULL, 
                         aspect=F, aspect.categories=NULL,
                         slope=F) {
  
  # require(raster)
  # require(SDMTools)
  
  # convert sgdf to a raster object
  sgdf = brick(sgdf)
  
  # mask data
  if(!is.null(mask))
    sgdf = mask(sgdf, brick(mask))

  # extract and aggregate the data
  if(!is.null(aggfact) && aggfact > 1)
    sgdf = aggregate(sgdf, fact=aggfact)
  
  if(aspect) {
    for(i in 1:sgdf@data@nlayers) {
      # compute aspects
      sgdf[[i]] = aspect(sgdf[[i]], latlon = T)
      # classify aspects
      if(!is.null(aspect.categories)) {
        sgdf[[i]]@data@values = as.numeric(
            cut((sgdf[[i]]@data@values + 180/aspect.categories) %% 360, 
                seq(0, 360, length.out = aspect.categories+1))
          )
      }
    }
  } else if(slope) {
    for(i in 1:sgdf@data@nlayers) {
      # compute slopes
      sgdf[[i]] = slope(sgdf[[i]], latlon = T)
    }
  }
  
  # crop data
  sgdf = crop(sgdf, extent)
  
  # compute anomalies
  if(!is.na(pmatch(type, 'std.anomaly')))
    sgdf@data@values = t(scale(t(sgdf@data@values), center=T, scale=T))
  else if(!is.na(pmatch(type, 'anomaly')))
    sgdf@data@values = t(scale(t(sgdf@data@values), center=T, scale=F))
  
  sgdf
}