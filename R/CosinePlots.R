#' Cosine Plot
#'
#' Plot the empirical and fitted cosineograms of the spatial correlation of
#' circular-spatial data. Assumption: Isotropic circular random field
#'
#' @param x,y vectors of location coordinates, directions is a vector of
#' directions in radians.
#' @param directions Vector of direction of observations or residual rotations
#' in radians.
#' @param Lag Vector of ascending distances, beginning with zero, where mean
#' cosine is to be computed. `Lag.n.Adj > 0` multiplies the number of lag points.
#' @param Lag.n.Adj adjusts `nBins`. Multiplier (`> 0`) of the number of lag
#' points. Value `> 1` increases the number of points for more detail.
#' Value `< 1` decreases the number of points for less detail.
#' @param BinWAdj Multiplier (`>=1`) of bin width. Value `> 1` has a smoothing effect. Sturges rule determines `nBins`.
#' @param Plot logical. `TRUE` plot cosineocloud or cosineogram, else ouput list of points.
#' @param Cloud logical. `TRUE` plots cosineocloud, else cosineogram.
#' @param Model `TRUE` overplots exponential, gaussian, and spherical models with nugget, Range, and sill parameters.
#' @param nugget Model nugget or mean cosine near zero distance. `0 <= nugget <= 1`.
#' @param Range Model range.
#' @param sill Model sill
#' @param x.legend,y.legend adjust legend location
#' @param TrimMean Apply trimmed mean (`0.0` to `0.5`) in computing the mean cosine at a distance. See [base::mean()].
#' @param ... Additional plotting parameters
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
#'
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
#'
#' ## Fit An Appropriate Model
#' FitHoriz1 <- lm(cos(sample.direction1) ~ (x1 + y1))
#' FitVert1 <- lm(sin(sample.direction1) ~ (x1 + y1))
#' fitted.direction1 <- atan2(FitVert1$fitted.values, FitHoriz1$fitted.values)
#'
#' ## Compute Residuals
#' resids1 <- CircResidual(
#'   X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1,
#'   Plot = FALSE
#' )
#'
#' ## Output list of cosineogram coordinates for fitting analytically
#' cosineogram.out <- CosinePlots(
#'   x = resids1$x, y = resids1$y, directions = resids1$direction,
#'   Lag.n.Adj = 1, BinWAdj = 1, Plot = FALSE, Cloud = FALSE, Model = FALSE
#' )
#'
#' ## Cosineocloud
#' CosinePlots(
#'   x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
#'   Plot = TRUE, Cloud = TRUE
#' )
#'
#' ## Cosineogram
#' CosinePlots(
#'   x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
#'   Plot = TRUE, Cloud = FALSE, Model = FALSE
#' )
#'
#' abline(h = 0.56, col = 2)
#' abline(v = 4, col = 2)
#'
#' ## Fit cosine Models
#' CosinePlots(
#'   x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
#'   Plot = TRUE, Cloud = FALSE, Model = TRUE, nugget = 0, Range = 4.0, sill = 0.56, x.legend = .2,
#'   y.legend = 0.3, xlim = c(0, 8), ylim = c(0, 1)
#' )
CosinePlots <- function(x, y, directions, Lag = NULL, Lag.n.Adj = 1, BinWAdj = 1, Plot = TRUE,
           Cloud = FALSE, Model = FALSE, nugget = 0, Range = NULL, sill = NULL, x.legend = 0.6, y.legend = 1.0, TrimMean = 0.1, ...) {
  stopifnot(!((length(x) != length(y)) | (length(x) != length(directions)) | (length(y) != length(directions))))
  # if ((length(x) != length(y)) | (length(x) != length(directions)) | (length(y) != length(directions))) {
  #     stop("lengths of vector inputs unequal")
  #   }
    stopifnot(Lag.n.Adj > 0)
    #if (Lag.n.Adj <= 0) stop("Lag.n.Adj invalid")
    stopifnot(!((nugget < 0) | (nugget > 1)))
    #if ((nugget < 0) | (nugget > 1)) stop("nugget invalid")
    if (!is.null(Range)) {
      stopifnot(Range > 0)
      #if (Range <= 0) stop("Range negative")
    }
    if (!is.null(sill)) {
      stopifnot(!((sill < 0) | (sill >= 1)))
      #if ((sill < 0) | (sill >= 1)) stop("sill invalid")
      stopifnot(!(1 - nugget < sill))
      #if (1 - nugget < sill) stop("1-nugget < sill")
    }

    # Repair Input and Remove missings
    if (BinWAdj < 1) BinWAdj <- 1 # points will fall out of bins if adjust < 1
    filter <- !is.na(directions)
    x <- x[filter]
    y <- y[filter]
    directions <- directions[filter]

    # Pairwise cosines
    # Subroutine to compute circular distances in radians
    CircDist <- function(alpha, beta) {
      alpha[alpha < 0] <- 2 * pi + alpha[alpha < 0]
      beta[beta < 0] <- 2 * pi + beta[beta < 0]
      theta <- abs(alpha - beta)
      theta[theta > pi] <- 2 * pi - theta[theta > pi]
      return(theta)
    }
    Cosines <- cos(outer(directions, directions, FUN = "CircDist"))
    filter.tri <- upper.tri(Cosines)
    Cosines <- Cosines[filter.tri]

    # Pairwise distances
    Distances <- as.matrix(dist(cbind(x, y))) # Diagonal of zero distances
    X <- Distances[filter.tri] # vector of distances corresponding to vector of cosines Cosines

    if (Cloud) {
      Y <- Cosines
    } else {
      if (!is.null(Lag)) {
        # Assumes equally spaced lags, except for first lag point
        HalfBinWidth <- BinWAdj * 0.5 * (Lag[2] - Lag[1])
        Lag.n <- length(Lag)
      } else {
        nBins <- trunc(log2(sum(filter.tri)) + 1) # Sturges rule
        Lag.n <- trunc(Lag.n.Adj * (nBins + 1))
        nBins <- Lag.n - 1
        distance.max <- max(X)
        HalfBinWidth <- BinWAdj * 0.5 * distance.max / nBins
        Lag <- seq(0, distance.max, length.out = Lag.n)
      }

      Y <- vector(mode = "numeric", length = Lag.n)
      if (Lag[1] == 0) {
        Y[1] <- 1
        i1 <- 2
      } else {
        i1 <- 1
      }
      for (i in i1:Lag.n)
      {
        filter <- abs(X - Lag[i]) <= HalfBinWidth
        Y[i] <- mean(Cosines[filter], trim = TrimMean)
      }

      X <- Lag # Cloud=FALSE filtered distances replaced by lag vector
    }
    if (Plot) {
      plot(X, Y, xlab = "Distance", ylab = "Cosine", cex.main = 1, ...)
      if (Cloud == FALSE & Model == TRUE) {
        xx <- seq(min(X), max(X), length.out = 101)
        c.e <- 1 - nugget - (1 - nugget - sill) * (1 - exp(-3 * xx / Range))
        c.g <- 1 - nugget - (1 - nugget - sill) * (1 - exp(-3 * (xx / Range)^2))
        X1 <- xx[xx <= Range]
        c.s <- 1 - nugget - (1 - nugget - sill) * (1.5 * X1 / Range - 0.5 * (X1 / Range)^3)
        X2 <- xx[xx > Range]
        c.s <- c(c.s, rep(sill, length(X2)))
        lines(xx, c.e, col = 2, lty = 1, lwd = 1)
        lines(xx, c.g, col = "tan", lty = 1, lwd = 3)
        lines(xx, c.s, col = 4, lty = 2, lwd = 1)

        legend(
          x = x.legend * max(X), y = y.legend, c("Exponential", "Gaussian", "Spherical"),
          lty = c(1, 1, 2), col = c(2, "tan", 4), lwd = c(1, 3, 1), cex = 1.1
        )
      }
    } else {
      return(list(distance = X, cosine = Y))
    }
  }
