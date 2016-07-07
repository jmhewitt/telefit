#' Formatting for longitude scales
#'
#' @importFrom scales trans_new
#' 

lon_trans = function() {
  
  trans_new( name = 'lon', transform = function(x){x}, 
             inverse = function(x){x}, 
             format = function(x) { 
               W = x < 0
               gsub('^0E$', '0', paste(gsub('-', '', as.character(x)), 
                                       ifelse(W, 'W', 'E'), sep=''))
             } )
  
}


#' Formatting for latidue scales
#'
#' @importFrom scales trans_new
#' 

lat_trans = function() {
  
  trans_new( name = 'lat', transform = function(x){x}, 
             inverse = function(x){x}, 
             format = function(x) { 
               S = x < 0
               gsub('^0N$', '0', paste(gsub('-', '', as.character(x)), 
                                       ifelse(S, 'S', 'N'), sep=''))
             } )
  
}