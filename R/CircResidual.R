#' #' Compute rotational residuals from the mean direction
#'
#' @param X,Y Vectors of horizontal and vertical coordinates of observation and
#' trend locations.
#' @param Raw Vector of horizontal coordinates of observation and trend locations.
#' @param Trend Vector of fitted direction in radians. NAs not allowed.
#' @param Plot If `FALSE` return output list. If `TRUE`, plot data (black),
#' model (tan), and residuals (dashed black) with asp=1.
#' @param AdjArrowLength Multiplies length of arrows in plots.
#' @param ... Additional plot parameters
#'
#' @return list with \describe{
#' \item{x,y}{Vectors of location coordinates}
#' \item{direction}{Vector of directions (in radians)}
#' }
#'
#' @importFrom fields arrow.plot
#' @export
#'
#' @examples
#' ## Construct Trend Model of 121 locations
#' xy <- expand.grid(1:11, 1:11) # grid
#' x1 <- xy[, 1]
#' y1 <- xy[, 2]
#' model.direction1 <- matrix(data = c(
#'   157, 141, 126, 113, 101, 90, 79, 67, 54, 40, 25, 152, 137, 123, 111, 100, 90, 80, 69, 57, 44, 30,
#'   147, 133, 120, 109, 99, 90, 81, 71, 60, 48, 35, 142, 129, 117, 107, 98, 90, 82, 73, 63, 52, 40,
#'   137, 125, 114, 105, 97, 90, 83, 75, 66, 56, 45, 132, 121, 111, 103, 96, 90, 84, 77, 69, 60, 50,
#'   127, 117, 108, 101, 95, 90, 85, 79, 72, 64, 55, 122, 113, 105, 99, 94, 90, 86, 81, 75, 68, 60,
#'   117, 109, 102, 97, 93, 90, 87, 83, 78, 72, 65, 112, 105, 99, 95, 92, 90, 88, 85, 81, 76, 70,
#'   107, 101, 96, 93, 91, 90, 89, 87, 84, 80, 75
#' ), ncol = 11, byrow = TRUE)
#' model.direction1 <- as.vector(model.direction1) * pi / 180
#'
#' ## Plot Trend Model, See Figure J-3 (a)
#' plot(x1, y1, type = "n", xlab = "", ylab = "", asp = 1)
#' fields::arrow.plot(x1, y1,
#'   u = cos(model.direction1), v = sin(model.direction1), arrow.ex = 0.1, xpd = TRUE,
#'   true.angle = TRUE, length = .1
#' )
#'
#' ## Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
#' ## from GRF (Range=4, spherical covariance).
#' set.seed(666)
#' crf1 <- SimulateCRF(
#'   CircDistr = "vM", Rho = sqrt(0.5), Range = 4, CovModel = "spherical",
#'   Grid = xy, OverFit = TRUE
#' )
#'
#' ## Plot CRF
#' par(mai = c(0.4, 0.35, .25, 0.25))
#' plot(crf1$x, crf1$y, type = "n", xlab = "", ylab = "", asp = 1)
#' fields::arrow.plot(
#'   a1 = crf1$x, a2 = crf1$y, u = cos(crf1$direction), v = sin(crf1$direction), arrow.ex = 0.1,
#'   xpd = TRUE, true.angle = TRUE, length = .1
#' )
#'
#' # Make sample
#' sample.direction1 <- model.direction1 + crf1$direction
#'
#' ## Plot Sample, See Figure J-3 (c)
#' sample.direction1 <- model.direction1 + crf1$direction
#' plot(x1, y1, type = "n", asp = 1)
#' fields::arrow.plot(
#'   a1 = x1, a2 = y1, u = cos(sample.direction1), v = sin(sample.direction1), arrow.ex = 0.125,
#'   xpd = TRUE, true.angle = TRUE, length = .1
#' )
#'
#' ## Fit An Appropriate Model
#' FitHoriz1 <- lm(cos(sample.direction1) ~ (x1 + y1))
#' FitVert1 <- lm(sin(sample.direction1) ~ (x1 + y1))
#' fitted.direction1 <- atan2(FitVert1$fitted.values, FitHoriz1$fitted.values)
#'
#' ## Plot Fitted Model
#' plot(x1, y1, type = "n", asp = 1, xlab = "", ylab = "")
#' fields::arrow.plot(x1, y1,
#'   u = cos(fitted.direction1), v = sin(fitted.direction1), arrow.ex = 0.1, xpd = TRUE,
#'   true.angle = TRUE, length = .1
#' )
#'
#' ## Compute Residuals
#' resids1 <- CircResidual(
#'   X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1,
#'   Plot = FALSE
#' )
#'
#' ## Plot Sample, Fitted Model, and Residual Rotations
#' CircResidual(
#'   X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1, Plot = TRUE,
#'   xlim = c(3, 7), ylim = c(3, 7)
#' )
CircResidual <- function(X, Y, Raw, Trend, Plot = FALSE, AdjArrowLength = 1, ...) {
  # 2008-11-10.2053
  # Assumptions: Raw may have NAs, trend has no NAs.  Trend locations and Raw locations are identical to compute residuals.
  # require(fields)
  if ((length(X) != length(Y)) | (length(X) != length(Raw)) | (length(X) != length(Trend)) | (length(Y) != length(Raw)) |
    (length(Y) != length(Trend)) | (length(Raw) != length(Trend))) {
    stop("lengths of vector inputs unequal")
  }
  if (AdjArrowLength <= 0) stop("AdjArrowLength invalid")
  if (sum(is.na(Trend)) > 0) stop("NAs not allowed in Trend")

  FilterNA <- is.na(Raw)
  x <- X[!FilterNA]
  y <- Y[!FilterNA]
  raw <- Raw[!FilterNA]
  trend <- Trend[!FilterNA]
  raw[raw < 0] <- raw[raw < 0] + 2 * pi # Like R1.Standardize in CircDataimage
  trend[trend < 0] <- trend[trend < 0] + 2 * pi
  circdist <- abs(raw - trend) # Linear distance in radians with NAs where raw has NAs
  circdist[circdist > pi] <- 2 * pi - circdist[circdist > pi] # Circular distance in radians
  resids <- circdist
  filter <- (trend > raw) & (trend - raw) < pi | (raw > trend) & (raw - trend) > pi
  resids[filter] <- -1 * circdist[filter]
  if (Plot == TRUE) {
    plot(X, Y, type = "n", xlab = "", ylab = "", asp = 1, ...)
    fields::arrow.plot(x, y, u = cos(raw), v = sin(raw), xpd = FALSE, true.angle = TRUE, arrow.ex = .15 * AdjArrowLength, length = .1, col = 1)
    fields::arrow.plot(X, Y, u = cos(Trend), v = sin(Trend), xpd = FALSE, true.angle = TRUE, arrow.ex = .15 * AdjArrowLength, length = .1, col = "tan", lwd = 3)
    fields::arrow.plot(x, y, u = cos(resids), v = sin(resids), xpd = FALSE, true.angle = TRUE, arrow.ex = .15 * AdjArrowLength, length = .1, col = 2, lty = 2)
  } else {
    return(list(x = x, y = y, direction = resids))
  }
}
