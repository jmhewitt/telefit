#' Standardized anomalies of CO Precipitation
#'
#' A dataset containing spatially-aggregated climate data from the ERA-Interim
#' and PRISM datasets.  The response comes from PRISM, average monthly 
#' precipitation in a DJF winter.  The covariates come from ERA-Interim, 
#' Colorado and Pacific Ocean (sea) surface temperatures.  All data has been
#' converted to standardized anomalies.
#'
#' @format A stData object with 33 years of observations:
#' \describe{
#'   \item{tLabs}{year labels for data columns}
#'   \item{coords.s}{centers of grid cells for Colorado data}
#'   \item{coords.r}{centers of grid cells for Pacific Ocean data}
#'   \item{X}{Array of design matrices for Colorado covariates}
#'   \item{Y}{Matrix of precipitation observations}
#'   \item{Z}{Matrix of Pacific Ocean data}
#'   \item{X.lab}{Label for covariate data, used by plotting functions}
#'   \item{Y.lab}{Label for response data, used by plotting functions}
#'   \item{Z.lab}{Label for covariate data, used by plotting functions}
#' }
#' 
#' @source \url{http://prism.oregonstate.edu}
#' @source \url{https://rda.ucar.edu/datasets/ds627.0/}
#' 
"coprecip"