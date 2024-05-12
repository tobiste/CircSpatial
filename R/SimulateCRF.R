#' Simulate a circular random field
#'
#' @param N Number of spatial locations to simulate
#' @param CircDistr Circular distribution:
#' \describe{
#'  \item{"U"}{uniform}
#'  \item{"vM"}{von Mises}
#'  \item{"WrC"}{wrapped Cauchy distributions}
#'  \item{"Tri"}{triangular Cauchy distribution}
#'  \item{"Card}{cardioid Cauchy distributions.}
#' }
#' @param Rho Mean resultant length parameter
#' \describe{
#' \item{for triangular}{`0 < Rho <= 4/pi^2`}
#' \item{for cardioid}{`0 < Rho <= 0.5`}
#' \item{For vM and wrapped Cauchy}{`0 < Rho < 1`, `1 == `degenerate}
#' \item{For uniform}{`Rho = 0` (all directions have equal density)}
#' }
#' @param Mu Mean resultant direction of circular distribution (rad). \eqn{|\mu| \leq \pi}{mu <= pi}
#' @param Range Distance at which CRV are not correlated
#' @param Ext `Range*Ext` is horizontal and vertical length of sample space
#' @param CovModel Name of spatial correlation function, see [geoR::cov.spatial()]
#' @param Grid Regular or irregular `N x 2` matrix of simulation locations,
#' overrides `N` and `Ext`
#' @param Anisotropy Vector of geometric anisotropy (angle, ratio). Angle in
#' radians, ratio >= 1. See [geoR::grf()]
#' @param logical. \describe{
#' \item{`TRUE`}{standardization (centering and rescaling realization
#' of the GRV to `mean = 0` and `sd = 1`) results in closer fit for qualitative evaluation
#' of the CRV. Undesirable effects are loss of independence of the marginal
#' GRVs, biased GRF covariance, and biased testing.
#' Standardization is suitable for demonstration with closer fit,
#' visualization, and illustrations. Do not standardize for the purposes of
#' simulation and testing.}
#' \item{`FALSE`}{non-standardization (the default)
#' includes expected variation from transformation of variation in mean and
#' sd of sample of GRV.
#' }
#' }
#' @param Resolution For nonclosed form inverse CDF, circular quantiles are
#' interpolated at resolution Resolution. `0.001 <= Resolution <= 0.01`
#' recommended.
#'
#' @return list with \describe{
#' \item{x,y}{Vectors of location coordinates}
#' \item{direction}{Vector of directions (in radians)}
#' \item{Z}{Vector of simulated observations of the GRV of the GRF.}
#' }
#'
#' @importFrom geoR grf
#'
#' @export
#'
#' @examples
#' xy <- expand.grid(1:11, 1:11) # grid
#' SimulateCRF(CircDistr = "vM", Rho = sqrt(0.5), Range = 4, Grid = xy, OverFit = TRUE)
#'
#' SimulateCRF(N = 200, CircDistr = "Card", Rho = 0.4, Range = 5, Ext = 3, CovModel = "exponential")
#' SimulateCRF(CircDistr = "U", Range = 8, Ext = 3, CovModel = "gaussian")
#' SimulateCRF(CircDistr = "Tri", Rho = 0.5 * 4 / pi^2, Range = 8, Ext = 3, CovModel = "spherical")
#' SimulateCRF(CircDistr = "WrC", Rho = sqrt(0.8), Range = 8, Ext = 3, CovModel = "exponential")
#' SimulateCRF(N = 400, CircDistr = "WrC", Rho = sqrt(0.95), Range = 8, Ext = 3, CovModel = "spherical", Anisotropy = c(pi / 4, 3))
SimulateCRF <- function(N = 100,
                        CircDistr = c("U", "vM", "WrC", "Tri", "Card"),
                        Rho,
                        Mu = 0,
                        Range,
                        Ext = 1,
                        CovModel = c("circular", "matern", "exponential", "gaussian", "spherical", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget"),
                        Grid = NULL,
                        Anisotropy = NULL,
                        OverFit = FALSE,
                        Resolution = 0.01) {

  CircDistr <- match.arg(CircDistr)
  CovModel <- match.arg(CovModel)

  if (CircDistr == "U") Rho <- 0
  stopifnot(abs(Mu) <= pi)
  #if (abs(Mu) > pi) stop("abs(Mu) <= pi")

  if (!is.null(Grid)) {
    if (!inherits(Grid, "matrix")) Grid <- as.matrix(Grid)
    stopifnot(dim(Grid)[2] == 2)
    #if (dim(Grid)[2] != 2) stop("Grid not a N x 2 matrix")
    N <- dim(Grid)[1]
  }
  if (!is.null(Anisotropy)) {
    stopifnot(length(Anisotropy) == 2)
    #if (length(Anisotropy) != 2) stop("Anisotropy is not a 2 element vector. See geoR Help")
  }

  #if (N <= 0 | Rho < 0 | Range < 0 | Ext <= 0 | Resolution <= 0) stop("Improper numeric input")
  stopifnot(!(N <= 0 | Rho < 0 | Range < 0 | Ext <= 0 | Resolution <= 0))

  direction <- vector(mode = "numeric", length = N)

  # require(CircStats)
  # require(geoR)
  # Standard normal GRF, see Help geoR grf
  if (is.null(Grid)) {
    GRF <- grf(
      n = N,
      xlims = c(0, Range * Ext),
      ylims = c(0, Range * Ext),
      cov.model = CovModel,
      nugget = 0,
      cov.pars = c(1, Range),
      aniso.pars = Anisotropy,
      #RF = TRUE,
      messages = FALSE
    )
  } else {
    GRF <- grf(
      grid = Grid,
      cov.model = CovModel,
      nugget = 0,
      cov.pars = c(1, Range),
      aniso.pars = Anisotropy,
      # RF=TRUE,
      messages = FALSE
    )
  }

  Grid <- GRF$coords # N x 2 matrix
  x <- Grid[, 1]
  y <- Grid[, 2]
  Z <- GRF$data # Vector of GRV
  if (OverFit) {
    Z <- (Z - mean(Z)) / sd(Z)
    GRF$data <- Z
  }
  CumProbZ <- pnorm(Z, mean = 0, sd = 1, lower.tail = TRUE)

  if (CircDistr == "U") {
    direction <- -pi + 2 * pi * CumProbZ
  } else if (CircDistr == "Tri") {

    stopifnot(!(Rho == 0 | Rho > 4 / pi^2))
    #if (Rho == 0 | Rho > 4 / pi^2) stop("Tri: 0 < Rho <= 4/pi^2")
    filter <- CumProbZ < 0.5
    u1 <- CumProbZ[filter]
    a <- Rho / 8
    b <- (4 + pi^2 * Rho) / (8 * pi)
    c <- 0.5 - u1
    q <- -.5 * (b + sqrt(b^2 - 4 * a * c))
    direction[filter] <- c / q

    u2 <- CumProbZ[!filter]
    a <- -Rho / 8
    b <- (4 + pi^2 * Rho) / (8 * pi)
    c <- 0.5 - u2
    q <- -.5 * (b + sqrt(b^2 - 4 * a * c))
    direction[!filter] <- c / q
  } else {
    # For OTHER circular distributions compute table of circular CDF and interpolate
    CircScale <- seq(-pi, pi, length = 2 * pi / Resolution)
    # With resolution=.01, circular support from -pi to +pi has 629 elements, delta ~0.01000507, CircScale[315] = 0
    n <- length(CircScale)
    if (CircDistr == "vM") {

      stopifnot(!(Rho == 0 | Rho >= 1))
      #if (Rho == 0 | Rho >= 1) stop("vM: 0 < Rho < 1")
      CircProb <- rep(-1, n)
      Kappa <- CircStats::A1inv(Rho) # N. I Fisher, Statistical Analysis of Circular Data, 2000 p. 49
      # As direction increases from -pi, pvm increases from .5

      for (i in 1:length(CircScale)) CircProb[i] <- CircStats::pvm(CircScale[i], mu = 0, kappa = Kappa)
      filter <- CircScale < 0
      CircProb[filter] <- CircProb[filter] - 0.5
      CircProb[!filter] <- CircProb[!filter] + 0.5
    } else if (CircDistr == "Card") {

      stopifnot(!(Rho == 0 | Rho > 0.5))
      #if (Rho == 0 | Rho > 0.5) stop("Cardioid: 0 < Rho <= 0.5")
      CircProb <- (CircScale + pi + 2 * Rho * sin(CircScale)) / (2 * pi)
    } else if (CircDistr == "WrC") {

      stopifnot(!(Rho == 0 | Rho >= 1))
      #if (Rho == 0 | Rho >= 1) stop("Wrapped Cauchy: 0 < Rho < 1 ")
      Angles1 <- CircScale[CircScale < 0]
      Angles2 <- CircScale[CircScale >= 0]
      prob1 <- 0.5 - acos(((1 + Rho^2) * cos(Angles1) - 2 * Rho) / (1 + Rho^2 - 2 * Rho * cos(Angles1))) / (2 * pi)
      prob2 <- 0.5 + acos(((1 + Rho^2) * cos(Angles2) - 2 * Rho) / (1 + Rho^2 - 2 * Rho * cos(Angles2))) / (2 * pi)
      CircProb <- c(prob1, prob2)
    }
    CircProb[1] <- 0
    CircProb[n] <- 1

    # Interpolation
    DeltaTh <- CircScale[2] + pi
    for (i in 1:N)
    {
      p <- CumProbZ[i] # Cumulative prob of GRV
      a <- max((1:n)[CircProb <= p]) # Index
      if (a == n) {
        r <- 0
      } else {
        if (CircProb[a] == p) {
          r <- 0
        } else {
          r <- (p - CircProb[a]) / (CircProb[a + 1] - CircProb[a])
        }
      }
      direction[i] <- CircScale[a] + r * DeltaTh
    }
  }

  direction <- direction + Mu
  return(list(x = x, y = y, direction = direction, Z = Z))
}
