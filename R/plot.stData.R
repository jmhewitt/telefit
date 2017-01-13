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
#'  'teleconnection', 'teleconnection_local', or 'teleconnection_knot_local'
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
#' @param category.breaks [ncoords x ncats] list of breakpoints used for binning
#'  responses into categories
#' @param coords.knots if plot type is 'remote', specifies the longitude and
#'  latitude of knot locations to overlay on the 'remote' plot
#' @param signif.telecon if TRUE, will highlight significant teleconnection
#'  effects when type=='teleconnection'
#' @param coord.r if plot type is 'teleconnection_local', specifes the longitude
#'  and latitude of remote coordinate for which to plot associated teleconnection
#'  effects.  if NULL, the middle remote coordinate will be plotted.
#'  
#' @return a ggplot object with the specified map
#'
#' 
#' 

plot.stData = function( stData, type='response', t=NULL, boxsize=NULL, p=NULL,  
                        map='world', region='.', coord.s=NULL, coord.r=NULL,
                        zlim=NULL,
                        lab.teleconnection = expression(alpha),
                        fill.lab.width = 20, category.breaks = NULL,
                        coords.knots = NULL, signif.telecon = F, dots=NULL, ...) {
  
  # merge unique list of dots
    dots = c(dots, list(...))
    dots = dots[!duplicated(dots)]
  # overwrite arguments to function if they exist in dots
    for(x in setdiff(names(formals(eval(match.call()[[1]]))), c('dots', '...'))) {
      if(x %in% names(dots)) {
        assign(eval(x), dots[[x]])
      }
    }
  
  if(!is.null(coord.s))
    coord.s=unlist(coord.s)
    
  if(is.null(t))
    t=stData$tLabs[1]
  
  if(is.null(p))
    p=2
  
  if(!is.null(stData$Y)) {
    if(class(stData$Y)!='matrix')
      stData$Y = matrix(stData$Y, nrow = nrow(stData$coords.s))  
  }
  
  if(!is.null(stData$Y.cat)) {
    if(class(stData$Y.cat)!='matrix')
      stData$Y.cat = matrix(stData$Y.cat, nrow = nrow(stData$coords.s))  
  }
  
  if(!is.null(coords.knots)) {
    coords.knots = data.frame(coords.knots)
    colnames(coords.knots) = c('lon', 'lat')
    coords.knots = rbind(coords.knots, coords.knots %>% mutate(lon=lon-360))
  }
  
  # extract dataset to plot
  match.opts = c('response', 'covariate', 'remote', 'teleconnection', 
                 'cat.response', 'teleconnection_knot', 'teleconnection_knot_local')
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
    # scheme.col = list(low = "#008837", mid = '#f7f7f7', high = '#7b3294')
    scheme.col = list(low = "#0571b0", mid = '#f7f7f7', high = '#ca0020')
  } else if( type=='teleconnection' ) {
    
    n = nrow(stData$coords.s)
    r = nrow(stData$coords.r)
    
    if(is.null(coord.s))
      coord.s = stData$coords.s[round(n/2),]
    
    Y = data.frame( Y = as.numeric(stData$alpha),
                    lon.Z = stData$coords.r[,1], 
                    lat.Z = stData$coords.r[,2],
                    lon.Y = rep(stData$coords.s[,1], rep(r,n)),
                    lat.Y = rep(stData$coords.s[,2], rep(r,n)) ) %>% 
      filter(lon.Y==coord.s[1], lat.Y==coord.s[2]) %>% 
      mutate(lon.Y=lon.Z, lat.Y=lat.Z )
    
    lab.col = lab.teleconnection
    scheme.col = list(low = "#0571b0", mid = '#f7f7f7', high = '#ca0020')
  } else if( type=='teleconnection_knot' ) {
    
    n = nrow(stData$coords.s)
    r_knots = nrow(stData$coords.knots)
    
    if(is.null(coord.s))
      coord.s = stData$coords.s[round(n/2),]
    
    Y = data.frame( Y = stData$alpha_knots,
                    signif = ifelse(stData$alpha_knots_signif, 3, 0),
                    lon.Z = stData$coords.knots[,1], 
                    lat.Z = stData$coords.knots[,2],
                    lon.Y = rep(stData$coords.s[,1], rep(r_knots,n)),
                    lat.Y = rep(stData$coords.s[,2], rep(r_knots,n)) ) %>% 
      filter(lon.Y==coord.s[1], lat.Y==coord.s[2]) %>% 
      mutate(lon.Y=lon.Z, lat.Y=lat.Z )
    
    lab.col = lab.teleconnection
    scheme.col = list(low = "#0571b0", mid = '#f7f7f7', high = '#ca0020')
  } else if( type=='teleconnection_knot_local' ) {
    
    n = nrow(stData$coords.s)
    r_knots = nrow(stData$coords.knots)
    
    if(is.null(coord.r))
      coord.r = stData$coords.s[round(r_knots/2),]
    
    Y = data.frame( Y = stData$alpha_knots,
                    lon.Z = stData$coords.knots[,1], 
                    lat.Z = stData$coords.knots[,2],
                    lon.Y = rep(stData$coords.s[,1], rep(r_knots,n)),
                    lat.Y = rep(stData$coords.s[,2], rep(r_knots,n)) ) %>% 
      filter(lon.Z==coord.r[1], lat.Z==coord.r[2])
    
    lab.col = lab.teleconnection
    scheme.col = list(low = "#0571b0", mid = '#f7f7f7', high = '#ca0020')
  } else if( type=='cat.response' ) {
    
    # categorize Y according to given breakpoints, if necessary
    if(is.null(stData$Y.cat)) {
      stData$Y.cat = stData$Y
      for(s in 1:nrow(stData$Y)) {
        stData$Y.cat[s,] = 1 + findInterval(stData$Y[s,], category.breaks[s,])
      }
    }
    
    # auto-determine labels based on the number of category breaks
    if(ncol(category.breaks)==1) {
      cat.labels = c('Below average', 'Above average')
      scheme.col = c('Below average'='#fc8d59', 'Above average'='#91bfdb')
    } else if(ncol(category.breaks)==2) {
      cat.labels = c('Below average', 'Near average', 'Above average')
      scheme.col = c('Below average'='#fc8d59', 'Near average'='#ffffbf', 
                     'Above average'='#91bfdb')
    }
    
    Y.cat = as.numeric(stData$Y.cat[, match(t, stData$tLabs)])
    Y.cat = factor(Y.cat, labels=cat.labels[sort(unique(Y.cat))])
    
    # build plotting frame
    Y = data.frame( Y = Y.cat,
                    lon.Y = stData$coords.s[,1], 
                    lat.Y = stData$coords.s[,2] )
    lab.col = paste(stData$Y.lab, 'level')
  } 
  
  # compute truncations and apply wrapping
  if(type %in% c('remote', 'teleconnection', 'teleconnection_knot')) {
    if(max(Y$lon.Y)>0) {
      if(min(Y$lon.Y)<0) {
        lon.E = max(Y %>% filter(lon.Y<=0) %>% dplyr::select(lon.Y))
        lon.W = min(Y %>% filter(lon.Y>0) %>% dplyr::select(lon.Y)) - 360
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
    world.raw = map_data('world') %>% filter(region!='USA')
    # duplicate countries for plotting with any map center
    world.raw = rbind(world.raw, world.raw %>% 
                        mutate(long=long-360, group=group+max(group)+1))
    # add outline data to state outlines
    world = rbind(world, world.raw %>% mutate(group=group+max(world$group)+1))
  }
  
  
  #
  # set commands to modify plotting options, if specified
  #
  
  if(type!='teleconnection_knot') {
    if(is.null(boxsize)) {
      tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y)
    }
    else {
      tile.aes = aes(x=lon.Y, y=lat.Y, fill=Y, width=boxsize, height=boxsize)
    }
  } else {
    if(signif.telecon) {
      point.aes = aes(x=lon.Y, y=lat.Y, fill=Y, stroke=signif)
    } else {
      point.aes = aes(x=lon.Y, y=lat.Y, fill=Y, stroke=0)
    }
  }
  
  # wrap fill label
  if(class(lab.col)!='expression') {
    lab.col = str_wrap(lab.col, width=fill.lab.width)
  }
  
  if(type=='cat.response') {
    fillscale = scale_fill_manual(lab.col, values = scheme.col)
  } else if(is.null(zlim)) {
    fillscale = scale_fill_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high)
    colscale = scale_color_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high)
  } else  {
    fillscale = scale_fill_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high,
                                     limits = zlim)
    colscale = scale_color_gradient2(lab.col,
                                     low = scheme.col$low, 
                                     mid = scheme.col$mid, 
                                     high = scheme.col$high,
                                     limits = zlim)
  } 
  
  # build base plot
  if(type!='teleconnection_knot') {
    
    worldmap = ggplot(world, aes(x=long, y=lat, group=group)) +
      geom_tile(tile.aes, data = Y  %>% 
                  mutate(lon.Y = ifelse(lon.Y<=0, lon.Y, lon.Y-360)), 
                inherit.aes = F) +
      fillscale +
      scale_x_continuous(trans = lon_trans()) +
      scale_y_continuous(trans = lat_trans()) +
      xlab('Longitude') +
      ylab('Latitude') 
    
    if(type %in% c('remote', 'teleconnection') ) {
      worldmap = worldmap + geom_polygon()
    } else {
      worldmap = worldmap + geom_path()
    }
    
    worldmap = worldmap + 
      theme_grey() +
      ggtitle(t)
    
    if(!is.null(coords.knots)) {
      worldmap = worldmap + geom_point(aes(x=lon, y=lat), data = coords.knots,
                                       col = 'black', fill = 'white', shape=21,
                                       inherit.aes = F)
    }
  } else {
    worldmap = ggplot(world, aes(x=long, y=lat, group=group)) +
      geom_point(point.aes, data = Y  %>% 
                  mutate(lon.Y = ifelse(lon.Y<=0, lon.Y, lon.Y-360)), 
                inherit.aes = F, size=4, shape=21) +
      fillscale +
      scale_color_manual('Significance', values=c('True'='black', 'False'='grey')) +
      scale_x_continuous(trans = lon_trans()) +
      scale_y_continuous(trans = lat_trans()) +
      xlab('Longitude') +
      ylab('Latitude') + 
      geom_path() +
      theme_grey() +
      ggtitle(t)
    
    if(!is.null(coords.knots)) {
      worldmap = worldmap + geom_point(aes(x=lon, y=lat), data = coords.knots,
                                       col = 'black', fill = 'white', shape=21,
                                       inherit.aes = F)
    }
  }
  
  
  # add coord.s to the plot and modify truncation
  if(!is.null(coord.s)) {
    worldmap = worldmap + geom_point(aes(x=lon.Y, y=lat.Y), 
                                     data = data.frame(lon.Y = coord.s[1],
                                                       lat.Y = coord.s[2]),
                                     col = 'black', fill = 'white', shape=21, 
                                     inherit.aes = F)
    
    lon.E = max(lon.E, coord.s[1])
    lon.W = min(lon.W, coord.s[1])
    lat.N = max(lat.N, coord.s[2])
    lat.S = min(lat.S, coord.s[2])
  }
  
  # apply map projection and truncation
  worldmap + coord_fixed(xlim=c(lon.W, lon.E), ylim=c(lat.S, lat.N), ratio=1.3)
}