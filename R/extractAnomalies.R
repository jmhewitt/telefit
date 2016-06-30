#' Compute anomalies for a SpatialGridDataFrame and returns an extent
#'
#' @return a modified SpatialGridDataFrame, sgdf, with the climatology for each
#'  location accessible via attr(sgdf@data@values, 'scaled:center')
#'  
#' @examples # Show an example where you plot the climatology grid too!
#'
#'
#' @export
#' 
#' @importFrom raster crop brick aggregate scale mask
#' 
#' @param extent raster::extent object featuring region to extract
#' @param aggfact if provided, will spatially average the data
#' @param mask if an sgdf is provided, the data will be masked before
#'  extraction, aggregation, and anomaly computation
#' 


extractAnomalies = function(sgdf, extent, aggfact=NULL, mask=NULL) {
  
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
  sgdf@data@values = t(scale(t(sgdf@data@values), center=T, scale=F))
  
  sgdf
}