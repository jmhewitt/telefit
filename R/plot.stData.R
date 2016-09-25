#' Plot stData objects
#'
#' This function provides basic plotting for telefit package data.
#' 
#'
#' @export
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom stringr str_wrap
#' 
#' @param boxsize size of grid boxes plotted
#' @param t timepoint to plot.  Will automatically plot the first timepoint if
#'  t=NULL.
#' @param p column index of local covariate to plot if type='covariate'. Will 
#'  automatically assume the local covariate data includes an intercept and will
#'  plot the second column if p=NULL.
#' @param map name of map provided by the maps package. These include county, 
#'  france, italy, nz, state, usa, world, world2.  By default, all stData plots
#'  will include us state outlines.
#' @param region name of subregions to include. Defaults to . which includes 
#'  all subregions. See documentation for map for more details.
#' @param type Either 'response', 'cat.response', 'covariate', 'remote', 
#'  or 'teleconnection' 
#'  to specify which part of stData to plot.  Note that 'teleconnection' applies
#'  only if the stData object contains information about teleconnection effects,
#'  i.e., if it is a simulated dataset or otherwise modified to include 
#'  estimates of teleconnection effects.  Note that the value for type can
#'  be an abbreviation since partial matching is used during plotting.
#' @param stData Object of class stData to plot.
#' @param coord.s if plot type is 'teleconnection', specifies the longitude and 
#'  latitude of local coordinate for which to plot teleconnection effects. if 
#'  NULL, the middle local coordinate will be plotted.
#' @param zlim c(min, max) vector that specifies the colorscale limits
#' @param lab.teleconnection label used for fill scale in teleconnection plot
#' @param fill.lab.width line width for fill scale label
#' 
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stData = function( stData, type='response', t=NULL, boxsize=NULL, p=NULL,  
                        map='world', region='.', coord.s=NULL, zlim=NULL,
                        lab.teleconnection = expression(alpha),
                        fill.lab.width = 20 ) {

  if(is.null(t))
    t=stData$tLabs[1]
  
  if(is.null(p))
    p=2
  
  if(class(stData$Y)!='matrix')
    stData$Y = matrix(stData$Y, nrow = nrow(stData$coords.s))
  
  # extract dataset to plot
  match.opts = c('response', 'covariate', 'remote', 'teleconnection', 
                 'cat.response')
  type = match.opts[pmatch(type, match.opts)]
  if( type=='response' ) {
    Y = data.frame( Y = stData$Y[, match(t, stData$tLabs)],
                    lon.Y = stData$coords.s[,1], 
                    lat.Y = stData$coords.s[,2] )
    lab.col = stData$Y.lab
    scheme.col = list(low = "#a6611a", mid = '#f5f5f5', high = '#018571')
  } else if( type=='covariate' ) {
    Y = data.frame( Y = stData$X[, p, match(t, stData$tLabs)],
                    lon.Y = stData$coords.s[,1], 
                    lat.Y = stData$coords.s[,2] )
    lab.col = stData$X.lab
    scheme.col = list(low = "#008837", mid = '#f7f7f7', high = '#7b3294')
  } else if( type=='remote' ) {
    Y = data.frame( Y = stData$Z[, match(t, stData$tLabs)],
                    lon.Y = stData$coords.r[,1], 
                    lat.Y = stData$coords.r[,2] )
    lab.col = stData$Z.lab
    scheme.col = list(low = "#008837", mid = '#f7f7f7', high = '#7b3294')
  } else if( type=='teleconnection' ) {
    
    n = nrow(stData$coords.s)
    r = nrow(stData$coords.r)
    
    if(is.null(coord.s))
      coord.s = stData$coords.s[round(n/2),]
    
    Y = data.frame( Y = stData$alpha,
                    lon.Z = stData$coords.r[,1], 
                    lat.Z = stData$coords.r[,2],
                    lon.Y = rep(stData$coords.s[,1], rep(r,n)),
                    lat.Y = rep(stData$coords.s[,2], rep(r,n)) ) %>% 
      filter(lon.Y==coord.s[1], lat.Y==coord.s[2]) %>% 
      mutate(lon.Y=lon.Z, lat.Y=lat.Z )
    
    lab.col = lab.teleconnection
    scheme.col = list(low = "#0571b0", mid = '#f7f7f7', high = '#ca0020')
  } else if( type=='cat.response' ) {
    Y = data.frame( Y = stData$Y[, match(t, stData$tLabs)],
                    lon.Y = stData$coords.s[,1], 
                    lat.Y = stData$coords.s[,2] )
    lab.col = paste(stData$Y.lab, 'level')
    scheme.col = list(low = "#a6611a", mid = '#f5f5f5', high = '#018571')
  } 
  
  # compute truncations and apply wrapping
  if(type %in% c('remote', 'teleconnection')) {
    if(max(Y$lon.Y)>0) {
      if(min(Y$lon.Y)<0) {
        lon.E = max(Y %>% filter(lon.Y<=0) %>% select(lon.Y))
        lon.W = min(Y %>% filter(lon.Y>0) %>% select(lon.Y)) - 360
      } else {
        lon.E = max(Y$lon.Y) - 360
        lon.W = min(Y$lon.Y) - 360
      }
    } else {
      lon.E = max(Y$lon.Y)
      lon.W = min(Y$lon.Y)
    }
    lat.S = min(Y$lat.Y)
    lat.N = max(Y$lat.Y)
    
    Y = rbind(Y, Y %>% mutate(lon.Y=lon.Y-360))
  } else {
    lon.W = min(Y$lon.Y)
    lon.E = max(Y$lon.Y)
    lat.S = min(Y$lat.Y)
    lat.N = max(Y$lat.Y)
  }
  

  # get us state outlines ggplot format
  world = map_data('state', region=region)
  # get country outlines ggplot format
  if(map=='world') {
    # get raw outline data
    world.raw = map_data('world')
    # duplicate countries for plotting with any map center
    world.raw = rbind(world.raw, world.raw %>% 
                        mutate(long=long-360, group=group+max(group)+1))
    # add outline data to state outlines
    world = rbind(world, world.raw %>% mutate(group=group+max(world$group)+1))
  }
  
  
  #
  # set commands to modify plotting options, if specified
  #
  
  if(is.null(boxsize))
    tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y)
  else
    tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y, width=boxsize, height=boxsize)
  
  # wrap fill label
  lab.col = str_wrap(lab.col, width=fill.lab.width)
  
  if(type=='cat.response') {
    fillscale = scale_fill_brewer(lab.col, 
                                  type = 'div',
                                  palette = 'BrBG',
                                  direction = 1)
  } else if(is.null(zlim)) {
    fillscale = scale_fill_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high)
  } else  {
    fillscale = scale_fill_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high,
                                     limits = zlim)
  } 
  
  # build base plot
  worldmap = ggplot(world, aes(x=long, y=lat, group=group)) +
    geom_tile(tile.aes, data = Y  %>% 
                mutate(lon.Y = ifelse(lon.Y<=0, lon.Y, lon.Y-360)), 
              inherit.aes = F) +
    fillscale +
    scale_x_continuous(trans = lon_trans()) +
    scale_y_continuous(trans = lat_trans()) +
    xlab('Longitude') +
    ylab('Latitude') + 
    geom_path() +
    theme_grey() +
    ggtitle(t)
  
  # add coord.s to the plot and modify truncation
  if(!is.null(coord.s)) {
    worldmap = worldmap + geom_point(aes(x=lon.Y, y=lat.Y), 
                                     data = data.frame(lon.Y = coord.s[1],
                                                       lat.Y = coord.s[2]),
                                     col = 2, inherit.aes = F)
    
    lon.E = max(lon.E, coord.s[1])
    lon.W = min(lon.W, coord.s[1])
    lat.N = max(lat.N, coord.s[2])
    lat.S = min(lat.S, coord.s[2])
  }
  
  # apply map projection and truncation
  worldmap + coord_fixed(xlim=c(lon.W, lon.E), ylim=c(lat.S, lat.N), ratio=1.3)
}