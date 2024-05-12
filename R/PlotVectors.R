#' Traditional plots of vector-spatial data.
#'
#' @param x,y numeric. coordinates
#' @param h,v numeric vector of horizontal and vertical component.
#' Missing values are permitted.
#' @param UnitVector logical
#' @param TriIcon logical
#' @param AdjArrowLength Arrow length multiplier
#' @param AdjHeadLength Arrow head length multiplier
#' @param TriIconAdj Multiplies size of icons
#' @param TriRatio Length to width ratio of triangle icon.
#' @param JitterPlot logical. If `TRUE`, add jitter to location coordinates
#' @param Jitter Amount of jitter = `Jitter` * [stats::runif()] value.
#' @param ... Additional plot parameters
#'
#' @export
#'
#' @examples
#' data(OceanWind)
#' wind.1997.Jan <- OceanWind[OceanWind$year > 1997 & OceanWind$year < 1997.1, -1]
#'
#' ## Direction Only
#' PlotVectors(
#'   x = wind.1997.Jan$x, y = wind.1997.Jan$y, h = wind.1997.Jan$u, v = wind.1997.Jan$v,
#'   UnitVector = TRUE, AdjArrowLength = 0.75, AdjHeadLength = 0.75, xlim = c(320, 350), ylim = c(0, 30)
#' )
#'
#' ## Direction and Magnitude
#' # fields function arrows omits arrowheads with a warning on
#' # any arrow of length less than 0.001 inch.
#' PlotVectors(
#'   x = wind.1997.Jan$x, y = wind.1997.Jan$y, h = wind.1997.Jan$u, v = wind.1997.Jan$v,
#'   UnitVector = FALSE, TriIcon = FALSE, AdjArrowLength = 3, AdjHeadLength = 0.4, xlim = c(320, 350),
#'   ylim = c(0, 30)
#' )
#'
#' ## Triangle Icons
#' PlotVectors(
#'   x = wind.1997.Jan$x, y = wind.1997.Jan$y, h = wind.1997.Jan$u, v = wind.1997.Jan$v,
#'   UnitVector = FALSE, TriIcon = TRUE, TriIconAdj = 0.25, TriRatio = 4, xlim = c(320, 350),
#'   ylim = c(0, 30)
#' )
PlotVectors <- function(x, y, h, v, UnitVector = TRUE, TriIcon = FALSE, AdjArrowLength = 1, AdjHeadLength = 1, TriIconAdj = 1,
                        TriRatio = 4, JitterPlot = FALSE, Jitter = 1, ...) {
  if ((length(x) != length(y)) | (length(h) != length(v)) | (length(x) != length(h))) stop("lengths of vector inputs unequal")

  filter <- is.na(h) | is.na(v) | (h == 0 & v == 0)
  x <- x[!filter]
  y <- y[!filter]
  h <- h[!filter]
  v <- v[!filter]
  # fields function arrows omits arrowheads with a warning on any arrow of length less than 1/1000 inch.

  Dir <- atan2(v, h)
  Dir[Dir < 0] <- Dir[Dir < 0] + 2 * pi
  if (JitterPlot == TRUE) {
    x <- x + Jitter * runif(length(x))
    y <- y + Jitter * runif(length(y))
  }
  plot(x, y, ty = "n", asp = 1, ...)

  if (UnitVector) {
    arrow.plot(x, y, cos(Dir), sin(Dir), true.angle = TRUE, arrow.ex = AdjArrowLength * 0.05, length = AdjHeadLength * 0.125, angle = 20, xpd = FALSE)
  } else {
    if (TriIcon) {
      m <- sqrt(h^2 + v^2) # magnitude
      w <- sqrt(m / TriRatio)
      n <- length(x)
      xa <- x + TriIconAdj * w * cos(Dir + pi / 2)
      ya <- y + TriIconAdj * w * sin(Dir + pi / 2)
      xb <- x + TriIconAdj * TriRatio * w * cos(Dir)
      yb <- y + TriIconAdj * TriRatio * w * sin(Dir)
      xc <- x + TriIconAdj * w * cos(Dir - pi / 2)
      yc <- y + TriIconAdj * w * sin(Dir - pi / 2)

      for (i in 1:n) polygon(x = c(xa[i], xb[i], xc[i]), y = c(ya[i], yb[i], yc[i]), density = -1, col = 1)
    } else {
      arrow.plot(x, y, h, v, true.angle = TRUE, arrow.ex = AdjArrowLength * 0.05, length = AdjHeadLength * 0.125, angle = 20, xpd = FALSE)
    }
  }
}
