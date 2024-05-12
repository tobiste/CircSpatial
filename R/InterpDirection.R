#' Interpolate the model of mean direction at each kriging location.
#'
#' Interpolate models of direction cosines and sines, separately to avoid cross over.  Fit plane to triangular half of cell
#' (rectangular element of regular grid of measurement locations) in which interpolation location occurs.
#' Assumptions - Locations to interpolate are within range of (in.x, in.y), inputs have no missing.
#'
#' @param in.x,in.y vector of input horizontal and vertical coordinates
#' @param in.direction vector of input direction in radians
#' @param out.x,out.y vector of interpolation output horizontal and vertical coordinates
#'
#' @return list with \describe{
#' \item{out.x,out.y}{Vectors of location coordinates}
#' \item{direction}{Vector of directions (in radians)}
#' }
#' @export
#'
#' @examples
#' ## Construct Trend Model of 121 locations
#' x1 <- 1:11
#' y1 <- 1:11
#' y1 <- rep(y1, 11)
#' x1 <- rep(x1, each = 11)
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
#' ## Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
#' ## from GRF (Range=4, spherical covariance).
#' set.seed(666)
#' crf1 <- SimulateCRF(
#'   CircDistr = "vM", Rho = sqrt(0.5), Range = 4, CovModel = "spherical",
#'   Grid = cbind(x1, y1), OverFit = TRUE
#' )
#'
#' # Make sample
#' sample.direction1 <- model.direction1 + crf1$direction
#' ## Fit An Appropriate Model
#' ## Code for median polish is contained in Appendix K, Section K.12
#' FitHoriz1 <- lm(cos(sample.direction1) ~ (x1 + y1))
#' FitVert1 <- lm(sin(sample.direction1) ~ (x1 + y1))
#' fitted.direction1 <- atan2(FitVert1$fitted.values, FitHoriz1$fitted.values)
#' ## Compute Residuals
#' resids1 <- CircResidual(
#'   X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1,
#'   Plot = FALSE
#' )
#'
#' ## Fit cosine Models
#' CosinePlots(
#'   x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
#'   Plot = TRUE, Cloud = FALSE, Model = TRUE, nugget = 0, Range = 4.0, sill = 0.56, x.legend = 0.2,
#'   y.legend = 0.4, xlim = c(0, 8), ylim = c(0, 1)
#' )
#'
#' ## Krig to residuals using cosine Model
#' x2 <- seq(1, 11, by = 0.2)
#' n <- length(x2)
#' y2 <- x2
#' y2 <- rep(y2, n)
#' x2 <- rep(x2, each = n)
#' rm(n)
#' krig2 <- KrigCRF(
#'   krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
#'     resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = FALSE
#' )
#'
#' ## Interpolate Fitted Model
#' interp2 <- InterpDirection(
#'   in.x = x1, in.y = y1, in.direction = fitted.direction1, out.x = krig2$x,
#'   out.y = krig2$y
#' )
#'
#' ## Plot Interpolated Fitted Model and Overplot Fitted Model. See Figure J-17.
#' plot(interp2$x, interp2$y, type = "n", asp = 1, xlim = c(5, 8), ylim = c(5, 8), xlab = "", ylab = "")
#' fields::arrow.plot(interp2$x, interp2$y,
#'   u = cos(interp2$direction), v = sin(interp2$direction), arrow.ex = 0.09,
#'   xpd = FALSE, true.angle = TRUE, length = .1, col = "tan"
#' )
#' fields::arrow.plot(x1, y1,
#'   u = cos(fitted.direction1), v = sin(fitted.direction1), arrow.ex = 0.06, xpd = FALSE,
#'   true.angle = TRUE, length = .1, col = 1
#' )
#'
#' ## Plot Estimate Of Direction And Overplot Sample. See Figure J-18.
#' estimate2 <- interp2$direction + krig2$direction
#' plot(interp2$x, interp2$y, type = "n", xlab = "", ylab = "", asp = 1)
#' fields::arrow.plot(interp2$x, interp2$y,
#'   u = cos(estimate2), v = sin(estimate2), arrow.ex = 0.05, xpd = FALSE,
#'   true.angle = TRUE, length = .05, col = "tan"
#' )
#' fields::arrow.plot(x1, y1,
#'   u = cos(sample.direction1), v = sin(sample.direction1), arrow.ex = 0.05,
#'   xpd = FALSE, true.angle = TRUE, length = .05, col = 1
#' )
#'
#' ## Zoom
#' plot(interp2$x, interp2$y, type = "n", xlab = "", ylab = "", asp = 1, xlim = c(3, 6), ylim = c(3, 6))
#' fields::arrow.plot(interp2$x, interp2$y,
#'   u = cos(estimate2), v = sin(estimate2), arrow.ex = 0.075,
#'   xpd = FALSE, true.angle = TRUE, length = .05, col = "tan"
#' )
#' fields::arrow.plot(x1, y1,
#'   u = cos(sample.direction1), v = sin(sample.direction1), arrow.ex = 0.05,
#'   xpd = FALSE, true.angle = TRUE, length = .05, col = 1
#' )
InterpDirection <- function(in.x, in.y, in.direction, out.x, out.y) {
    minx.in <- min(in.x)
    maxx.in <- max(in.x)
    miny.in <- min(in.y)
    maxy.in <- max(in.y)
    minx.out <- min(out.x)
    maxx.out <- max(out.x)
    miny.out <- min(out.y)
    maxy.out <- max(out.y)
    if (minx.out < minx.in | maxx.out > maxx.in | miny.out < miny.in | maxy.out > maxy.in) stop("Interpolation range exceeds range of (in.x, in.y)")

    if ((length(in.x) != length(in.y)) | (length(in.x) != length(in.direction)) | (length(in.y) != length(in.direction))) {
      stop("lengths of vector inputs unequal")
    }

    if (length(out.x) != length(out.y)) stop("lengths of vector outputs unequal")


    # Organize model data
    X <- sort(unique(in.x), decreasing = FALSE) # Increases left to right
    m <- length(X)
    Y <- sort(unique(in.y), decreasing = TRUE) # Decreases top to bottom
    n <- length(Y)
    # Col of matrix of directions reflects the horiz or west to east component location
    xmin <- min(X)
    deltax <- X[2] - X[1]
    ymax <- max(Y)
    deltay <- Y[1] - Y[2]
    Col <- 1 + (in.x - xmin) / deltax
    Row <- 1 + (ymax - in.y) / deltay
    directions <- matrix(data = NA, nrow = n, ncol = m)
    directions[cbind(Row, Col)] <- in.direction # matrix of organized directions
    U <- cos(directions) # matrix of organized cosines of directions
    V <- sin(directions) # matrix of organized sines of directions

    n <- length(out.x)

    CosOut <- rep(NA, n) # for interpolated cosine
    SinOut <- CosOut # for interpolated sin

    p <- 1:length(X)
    q <- 1:length(Y)

    for (i in 1:n)
    {
      xx <- out.x[i]
      yy <- out.y[i]
      Vert <- FALSE
      Horiz <- FALSE
      if (sum(X == xx) == 1) Vert <- TRUE
      if (sum(Y == yy) == 1) Horiz <- TRUE

      if (Vert == FALSE & Horiz == FALSE) {
        west <- max(p[X <= xx])
        east <- west + 1
        south <- min(q[Y <= yy])
        north <- south - 1
        x.west <- X[west]
        x.east <- X[east]
        y.south <- Y[south]
        y.north <- Y[north]
        cos.nw <- U[north, west]
        cos.ne <- U[north, east]
        cos.sw <- U[south, west]
        cos.se <- U[south, east]
        sin.nw <- V[north, west]
        sin.ne <- V[north, east]
        sin.sw <- V[south, west]
        sin.se <- V[south, east]

        m <- (y.north - y.south) / (x.east - x.west) # 1 if vert res=horiz res
        b <- y.north - m * x.east
        ydiag <- m * xx + b
        if (yy <= ydiag) # On diagonal or in lower triangular
          {
            # Fit plane to lower triangular
            AB <- c(x.east - x.west, 0, cos.se - cos.sw)
            AC <- c(x.east - x.west, y.north - y.south, cos.ne - cos.sw)
            # Coefficients of cross product AB X AC
            a <- AB[2] * AC[3] - AB[3] * AC[2]
            b <- AB[3] * AC[1] - AB[1] * AC[3]
            c <- AB[1] * AC[2] - AB[2] * AC[1]
            CosOut[i] <- cos.sw + (a * (x.west - xx) + b * (y.south - yy)) / c
            # Fit plane to lower triangular
            AB <- c(x.east - x.west, 0, sin.se - sin.sw)
            AC <- c(x.east - x.west, y.north - y.south, sin.ne - sin.sw)
            # Coefficients of cross product AB X AC
            a <- AB[2] * AC[3] - AB[3] * AC[2]
            b <- AB[3] * AC[1] - AB[1] * AC[3]
            c <- AB[1] * AC[2] - AB[2] * AC[1]
            SinOut[i] <- sin.sw + (a * (x.west - xx) + b * (y.south - yy)) / c
          } else {
          # In upper triangular
          AC <- c(x.east - x.west, y.north - y.south, cos.ne - cos.sw)
          AD <- c(0, y.north - y.south, cos.nw - cos.sw)
          a <- AC[2] * AD[3] - AC[3] * AD[2]
          b <- AC[3] * AD[1] - AC[1] * AD[3]
          c <- AC[1] * AD[2] - AC[2] * AD[1]
          CosOut[i] <- cos.sw + (a * (x.west - xx) + b * (y.south - yy)) / c
          AC <- c(x.east - x.west, y.north - y.south, sin.ne - sin.sw)
          AD <- c(0, y.north - y.south, sin.nw - sin.sw)
          a <- AC[2] * AD[3] - AC[3] * AD[2]
          b <- AC[3] * AD[1] - AC[1] * AD[3]
          c <- AC[1] * AD[2] - AC[2] * AD[1]
          SinOut[i] <- sin.sw + (a * (x.west - xx) + b * (y.south - yy)) / c
        }
      } else if (Vert == TRUE & Horiz == FALSE) {
        p1 <- p[X == xx] # Column of vert grid line
        q1 <- min(q[Y < yy]) # Even spacing is not assumed
        q2 <- q1 - 1
        cos1 <- U[q1, p1]
        cos2 <- U[q2, p1]
        CosOut[i] <- cos1 + (cos2 - cos1) * (yy - Y[q1]) / (Y[q2] - Y[q1])
        sin1 <- V[q1, p1]
        sin2 <- V[q2, p1]
        SinOut[i] <- sin1 + (sin2 - sin1) * (yy - Y[q1]) / (Y[q2] - Y[q1])
      } else if (Vert == FALSE & Horiz == TRUE) {
        q1 <- q[Y == yy] # Row of horiz grid line
        p1 <- max(p[X < xx])
        p2 <- p1 + 1
        cos1 <- U[q1, p1]
        cos2 <- U[q1, p2]
        CosOut[i] <- cos1 + (cos2 - cos1) * (xx - X[p1]) / (X[p2] - X[p1])
        sin1 <- V[q1, p1]
        sin2 <- V[q1, p2]
        SinOut[i] <- sin1 + (sin2 - sin1) * (xx - X[p1]) / (X[p2] - X[p1])
      } else # Vert==TRUE & Horiz==TRUE
      {
        CosOut[i] <- U[q[Y == yy], p[X == xx]]
        SinOut[i] <- V[q[Y == yy], p[X == xx]]
      }
    }
    dir <- atan2(SinOut, CosOut)
    return(list(x = out.x, y = out.y, direction = dir))
  }
