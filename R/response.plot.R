#' Plots teleconnection correlation maps
#'
#' This function provides basic plotting for analyses returned from cor.tel
#' 
#'
#' @export
#' 
#' @import ggplot2
#' @importFrom dplyr mutate filter "%>%"
#' 
#' @param Y a response field to plot; must be in long format with (lon.Y, lat.Y)
#'  coordinates
#' @param boxsize size of grid boxes plotted
#' 
#' @return a ggplot object with the response map
#'
#' 
#' 

response.plot = function(Y, boxsize=NULL ) {
  
  require(dplyr)
  
  # get country outlines ggplot format
  world = map_data("world")
  # duplicate countries for plotting with any map center
  world = rbind(world, world %>% mutate(long=long-360, group=group+max(group)+1))
  # add us state outlines to map data
  world = rbind(world, map_data('state') %>% mutate(group=group+max(world$group)+1))
  
  #
  # set commands to modify plotting options, if specified
  #
  
  if(is.null(boxsize))
    tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y)
  else
    tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y, width=boxsize, height=boxsize)
  
  
  # build base plot
  worldmap = ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_tile(tile.aes, data = Y  %>% 
                mutate(lon.Y = ifelse(lon.Y<=0, lon.Y, lon.Y-360)), 
              inherit.aes = F) +
    scale_fill_gradient2(low = "#a6611a", mid = '#f5f5f5', high = '#018571') +
    scale_x_continuous(trans = lon_trans()) +
    scale_y_continuous(trans = lat_trans()) +
    xlab('Longitude') +
    ylab('Latitude') + 
    geom_path() +
    theme_grey()
  
  # apply map projection and truncation
  lon.W = min(Y$lon.Y)
  lon.E = max(Y$lon.Y)
  lat.S = min(Y$lat.Y)
  lat.N = max(Y$lat.Y)
  worldmap + coord_fixed(xlim=c(lon.W, lon.E), ylim=c(lat.S, lat.N), ratio=1.3)
}