#' Test Pattern
#'
#' The function TestPattern makes an intuitive simple
#' test pattern to explore the function CircDataimage, which produces a GUI for
#' interactive imaging of circular-spatial data. TestPattern computes direction
#' such that the direction at any point is the angle between the line from
#' origin to point and the horizon.
#'
#' @return list with \describe{
#' \item{x,y}{Vectors of location coordinates}
#' \item{u}{Vector of horizontal component of cosines of direction.}
#' \item{v}{Vector of vertical component of sines of direction.}
#' }
#' @export
#'
#' @examples
#' TestPattern()
TestPattern <- function() {
    x <- seq(-1, 1, length = 251)
    y <- x
    x <- rep(x, 251)
    y <- rep(y, ea = 251)
    dir <- atan2(y, x)
    # dir[x == 0] <- NA
    # dir[y == 0] <- NA
    u <- cos(dir)
    v <- sin(dir)
    return(as.data.frame(list(x = x, y = y, u = u, v = v)))
  }
