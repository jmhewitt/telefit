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
#' @importFrom raster crop brick aggregate scale mask
#' 
#' @param type whether to return the raw data, anomalies (data minus temporal 
#' average at each location), or standardized anomalies (anomalies divided by
#' temporal standard deviation at each location)
#' @param extent raster::extent object featuring region to extract
#' @param aggfact if provided, will spatially average the data
#' @param mask if an sgdf is provided, the data will be masked before
#'  extraction, aggregation, and anomaly computation
#' 


extractRegion = function(sgdf, extent, 
                         type='response', 
                         aggfact=NULL, mask=NULL) {
  
  # convert sgdf to a raster object
  sgdf = brick(sgdf)
  
  # mask data
  if(!is.null(mask))
    sgdf = mask(sgdf, brick(mask))
    
  # extract and aggregate the data
  sgdf = crop(sgdf, extent)
  if(!is.null(aggfact) && aggfact > 1)
    sgdf = aggregate(sgdf, fact=aggfact)
  
  # compute anomalies
  if(!is.na(pmatch(type, 'std.anomaly')))
    sgdf@data@values = t(scale(t(sgdf@data@values), center=T, scale=T))
  else if(!is.na(pmatch(type, 'anomaly')))
    sgdf@data@values = t(scale(t(sgdf@data@values), center=T, scale=F))
    
  
  sgdf
}