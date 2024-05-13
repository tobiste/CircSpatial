#' Circular kriging
#'
#' Krige and smooth the circular-spatial residuals.
#'
#' @param krig.x,krig.y numeric vector of horizontal and vertical coordinates of kriging locations.
#' @param resid.x,resid.y numeric Vector of horizontal coordinates of rotational residuals or data. no NAs
#' @param resid.direction Vector of direction in radians of rotational residuals or data.
#' @param Model covariance model. Object of  class [RandomFields::RMmodel()],
#' e.g. [RandomFields::RMcauchy()], [RandomFields::RMexp()],
#' [RandomFields::RMgencauchy()], [RandomFields::RMgauss()],
#' [RandomFields::RMgneiting()],[RandomFields::RMmatern()],
#' [RandomFields::RMnugget()], [RandomFields::RMspheric()],
#' [RandomFields::RMstable()], [RandomFields::RMwhittle()]).
#' See [RandomFields::RFcov()] for more details.
#' @param Nugget nugget component of the variogram (this basically adds a nugget component to the model); if missing, nugget component is omitted
#' @param Range range parameter of the variogram model component; in case of anisotropy: major range
#' @param sill sill of the variogram model component, or model: Mean cosine at the Range.
#' @param params `list` that specifies free parameters in a formula description, see [RandomFields::RMformula()]
#' @param Smooth logical.
#' @param bandwidth Kernel smoothing bandwidth (>0)
#' @param Plot logical
#' @param Xlim,Ylim plot limits.
#' @param PlotVar logical. Plot circular kriging variance.
#' @param ... optional params to [RandomFields::RFcov()]
#'
#' @importFrom fields arrow.plot
#' @importFrom RandomFields RFcov RFoptions
#'
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
#' ## Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
#' ## from GRF (Range=4, spherical covariance).
#' set.seed(666)
#' crf1 <- SimulateCRF(
#'   CircDistr = "vM", Rho = sqrt(0.5), Range = 4, CovModel = "spherical",
#'   Grid = xy, OverFit = TRUE
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
#'
#' x2 <- seq(1, 11, by = 0.2)
#' y2 <- x2 ## Kriging locations
#'
#' m <- RandomFields::RMexp(var = 1, scale = 1)
#'
#' krig2 <- KrigCRF(
#'   krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y,
#'   resid.direction = resids1$direction, Model =m, Nugget = 0.0, Range = 4, sill = 0.56,
#'   Plot = FALSE
#' )
KrigCRF <- function(krig.x, krig.y, resid.x, resid.y, resid.direction, Model = RandomFields::RMexp(var = 1, scale = 1), Nugget = 0, Range, sill,
                    #params = list(M = 0,var = 1,nugg = 0,scale = 1),
                    Smooth = FALSE, bandwidth, Plot = FALSE, Xlim = NULL, Ylim = NULL, PlotVar = FALSE, ...) {
  # 2008-11-11.1213
  # select model from covariance models in R package Random Fields function CovarianceFct ## ??RFcov
  # resid.x and resid.y have no NAs

  if ((length(krig.x) != length(krig.y)) | (length(resid.x) != length(resid.y)) | (length(resid.x) != length(resid.direction)) |
    (length(resid.y) != length(resid.direction))) {
    stop("lengths of vector inputs unequal")
  }
  if ((Nugget < 0) | (Nugget > 1)) stop("Nugget invalid")

  # fix the order of the kriging coordinates
  xx <- sort(unique(krig.x))
  yy <- sort(unique(krig.y))
  nx <- length(xx)
  ny <- length(yy) # rectangular or square grid
  krig.y <- rep(yy, nx)
  krig.x <- rep(xx, each = ny)

  Distances <- as.matrix(dist(cbind(resid.x, resid.y)))
  Ncol <- ncol(Distances)
  K <- c()
  RandomFields::RFoptions("allow_duplicated_locations" = TRUE)

  for (i in 1:Ncol)
  {
    # K <- cbind(K, sill + (1 - Nugget - sill) * CovarianceFct(
    #   x = Distances[, i] / Range, model = Model,
    #   param = c(mean = 0, variance = 1, nugget = 0, scale = 1, ...), dim = 1, fctcall = "Cov"
    # ))
    K <- cbind(K, sill + (1 - Nugget - sill) * RandomFields::RFcov(
      x = Distances[, i] / Range, model = Model,
      #params = params,
      dim = 1
    ))
  }
  diag(K) <- 1 # TRUE even if nugget > 0 for any model
  Kinv <- solve(K)

  U <- t(cbind(cos(resid.direction), sin(resid.direction)))
  # V <- t(U) %*% U

  n <- length(krig.x) # krig.x=krig.y for square or rect grid
  krig.direction <- krig.variance <- vector(mode = "numeric", length = n)

  for (i in 1:n)
  {
    distances <- sqrt((krig.x[i] - resid.x)^2 + (krig.y[i] - resid.y)^2)
    # c <- sill + (1 - Nugget - sill) * CovarianceFct(
    #   x = distances / Range, model = Model,
    #   param = c(mean = 0, variance = 1, nugget = 0, scale = 1, ...), dim = 1, fctcall = "Covariance"
    # )
    # c <- sill + (1 - Nugget - sill) * CovarianceFct(
    #   x = distances / Range, model = Model,
    #   param = c(mean = 0, variance = 1, nugget = 0, scale = 1, ...), dim = 1, fctcall = "Cov"
    # )
    c <- sill + (1 - Nugget - sill) * RandomFields::RFcov(
      x = distances / Range, model = Model,
      #params = params,
      dim = 1
    )

    c[distances == 0] <- 1 # TRUE even if nugget > 0 for any model
    w <- Kinv %*% c

    # w <- (Kinv %*% c)/sqrt(as.numeric(t(c) %*% Kinv %*% V %*% Kinv %*% c)) # gives same directions
    u <- U %*% w
    krig.direction[i] <- atan2(u[2], u[1])
    krig.variance[i] <- 2 - 2 * sqrt(as.numeric(t(c) %*% Kinv %*% c))
  }

  if (Smooth) {
    xx.dx <- xx[2] - xx[1]
    yy.dy <- yy[2] - yy[1]
    # as.image loads the matrix by row
    ImageList.x <- fields::as.image(cos(krig.direction), x = data.frame(krig.x, krig.y), nrow = nx, ncol = ny, boundary.grid = FALSE)
    smooth.x <- fields::image.smooth(ImageList.x, theta = bandwidth)
    ImageList.y <- fields::as.image(sin(krig.direction), x = data.frame(krig.x, krig.y), nrow = nx, ncol = ny, boundary.grid = FALSE)
    smooth.y <- fields::image.smooth(ImageList.y, theta = bandwidth)
    krig.direction <- as.vector(t(atan2(smooth.y$z, smooth.x$z)))
  }

  if (Plot) {
    if (!PlotVar) {
      plot(krig.x, krig.y, ty = "n", xlab = "", ylab = "", asp = 1, xlim = Xlim, ylim = Ylim)
      fields::arrow.plot(
        a1 = krig.x, a2 = krig.y, u = cos(krig.direction), v = sin(krig.direction), arrow.ex = 0.06,
        xpd = FALSE, true.angle = TRUE, length = .05, col = 1
      )
    } else {
      krig.variance.matrix <- matrix(data = krig.variance, nrow = nx, ncol = ny, byrow = TRUE)
      if (!is.null(Xlim)) {
        filled.contour(
          x = xx, y = yy, z = krig.variance.matrix, xlim = Xlim, ylim = Ylim,
          color = terrain.colors, key.title = title(main = "Circ Krig \nVariance", cex.main = 0.8),
          asp = 1, plot.axes = {
            axis(1)
            axis(2)
            points(resid.x, resid.y, pch = 20, cex = .65)
          }
        )
      } else {
        filled.contour(
          x = xx, y = yy, z = krig.variance.matrix,
          color = terrain.colors, key.title = title(main = "Circ Krig \nVariance", cex.main = 0.8),
          asp = 1, plot.axes = {
            axis(1)
            axis(2)
            points(resid.x, resid.y, pch = 20, cex = .65)
          }
        )
      }
    }
  } else {
    return(list(x = krig.x, y = krig.y, direction = krig.direction, variance = krig.variance))
  }
}
