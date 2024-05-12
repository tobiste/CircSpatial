#' Ocean Wind Data
#'
#' OceanWind provides a large dataset (495,688 observations) to further explore
#' the function CircDataimage. The OceanWind data was freely obtained from
#' ICOADS for the El Nino years 1972, 1976, 1982, 1987, 1991, 1994, and 1997,
#' January through April, in 1º increments for the area of longitude 0.5º E
#' to +359.5º E by latitude -59.5º N to +60.5º N.
#'
#' @docType data
#'
#' @usage data('OceanWind')
#' @aliases OceanWind
#'
#' @format `data.frame()`: \describe{
#' \item{Vector of time of observation = year + month/12, month in [1972.083, 1997.333]}
#' \item{x,y}{numeric vectors of longitude and latitude}
#' \item{u,v}{Vector of east and north component of wind (0.01 meters/second).}
#' }
#'
#' @source \url{http://dss.ucar.edu/datasets/ds540.1/data/msga.form.html}
#'
#' @keywords datasets
#'
#' @examples
#' data("OceanWind")
#' head(OceanWind)
"OceanWind"


#' World Land Mask
#'
#' WorldMask is used by [CircDataimage()] to restore land contours to the
#' circular dataimage of smoothed OceanWind. WorldMask was derived from the R package `fields` dataset world.dat.
#'
#' @docType data
#'
#' @usage data('WorldMask')
#' @aliases WorldMask
#'
#' @format `matrix()` of WorldMask is a matrix of 360 rows (0.5º to 359.5º
#' longitude) x 121 columns (-59.5º to +60.5º latitude) suited to OceanWind,
#' with elements NA where wind data is not missing and 1 where wind data is
#' missing.
#'
#' @keywords datasets
#'
#' @examples
#' data("WorldMask")
#' head(WorldMask)
"WorldMask"
