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
#' @param cor.tel an output object from cor.tel
#' @param coord.Y c(lon, lat) specifying a base point for which to plot remote 
#'  correlations
#' @param boxsize size of grid boxes plotted
#' @param lon.W westmost longitude to plot
#' @param lon.E eastmost longitude to plot
#' @param lat.S southmost latitude to plot
#' @param lat.N northmost latitude to plot
#' @param signif if TRUE, then cor.tel must have a column labeled 'signif' that
#'  indicates which correlations are significant.  These correlations will be
#'  printed in bold, and the rest will be printed more lightly
#' 
#' @return a ggplot object with the teleconnection map
#'
#' 
#' 

plot.cor.tel = function(cor.tel, coord.Y, lon.W, lon.E, lat.S, lat.N, 
                        signif=F, boxsize=NULL ) {
  
  require(dplyr)
  
  # get country outlines ggplot format
  world = map_data("world")
  # duplicate countries for plotting with any map center
  world = rbind(world, world %>% mutate(long=long-360, group=group+max(group)+1))
  # add us state outlines to map data
  world = rbind(world, map_data('state') %>% mutate(group=group+max(world$group)+1))
  
  lon_ = coord.Y[1]
  lat_ = coord.Y[2]
  
  #
  # set commands to modify plotting options, if specified
  #
  
  if(is.null(boxsize))
    tile.aes = aes(x=lon.Z, y=lat.Z, fill=cor)
  else
    tile.aes = aes(x=lon.Z, y=lat.Z, fill=cor, width=boxsize, height=boxsize)
  
  if(signif)
    alpha = .2
  else 
    alpha = 1
  
  # build base plot
  worldmap = ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_tile(tile.aes, data = cor.tel %>% filter(lon.Y==lon_, lat.Y==lat_) %>% 
                mutate(lon.Z = ifelse(lon.Z<=0, lon.Z, lon.Z-360)), 
              inherit.aes = F, alpha = alpha) +
    geom_point(aes(x=lon.Y, y=lat.Y),
               data = data.frame(lon.Y = lon_, lat.Y = lat_),
               inherit.aes = F, col=2) +
    scale_fill_gradient2(low = "#3679bd", high = '#bd3636') +
    xlab('Longitude') +
    ylab('Latitude') + 
    geom_path() +
    theme_grey()
  
  # add significant overlays, if applicable
  if(signif)
    if(nrow(cor.tel %>% filter(lon.Y==lon_, lat.Y==lat_) %>% 
            filter(signif==T)) > 0) {
      worldmap = worldmap + 
        geom_tile(tile.aes, data = cor.tel %>% 
                    filter(lon.Y==lon_, lat.Y==lat_) %>% filter(signif==T),
                  inherit.aes = F, color='black', lwd=1.75, alpha = 1)
    }
  
  # apply map projection and truncation
  worldmap + coord_fixed(xlim=c(lon.W, lon.E), ylim=c(lat.S, lat.N), ratio=1.3)
}