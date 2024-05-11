#' Simulate CRF ~ (Range, CircDistr, Rho, mu=0)
#'
#' @param N Number of spatial locations to simulate
#' @param CircDistr Circular distribution in {U, vM, WrC, Tri, Card},
#' @param Rho Mean resultant length parameter
#' \describe{
#' \item{for triangular}{`0 < Rho <= 4/pi^2`}
#' \item{for cardioid}{`0 < Rho <= 0.5`}
#' \item{For vM and wrapped Cauchy}{`0 < Rho < 1`, `1 == `degenerate}
#' \item{For uniform}{`Rho = 0`}
#' }
#' @param Mu argaethsdthsth
#' @param Range Distance at which CRV independent
#' @param Ext `Range*Ext` is horizontal and vertical length of sample space
#' @param CovModel Name of spatial correlation function, see package `geoR` Help `cov.spatial`
#' @param xy Regular or irregular `N x 2` matrix of simulation locations, overides `N` and `Ext`
#' @param Anisotropy asrgadh
#' @param logical. \describe{
#' \item{`TRUE`}{standardization (centering and rescaling realization
#' of the GRV to `mean = 0` and `sd = 1`) results in closer fit for qualitative evaluation
#' of the CRV. Undesirable effects are loss of independence of the marginal
#' GRVs, biased GRF covariance, and biased testing.
#' Standardization is suitable for demonstration with closer fit,
#' visualization, and illustrations. Do not standardize for the purposes of
#' simulation and testing.}
#'item{`FALSE`}{non-standardization (default)
#' includes expected variation from transformation of variation in mean and
#' sd of sample of GRV.
#' }
#' }
#' @param Resolution drh
#'
#' @return \describe{
#' \item{x,y}{Vectors of location coordinates}
#' \item{direction}{Vector of directions (in radians)}
#' \item{Z}{Vector of simulated observations of the associated GRV}
#' }
#'
#' @importFrom geoR grf
#'
#' @export
#'
#' @examples
#' xy <- expand.grid(1:11, 1:11) # grid
#'
#' model.direction1 <- c(
#'  157, 141, 126, 113, 101, 90, 79, 67, 54, 40, 25, 152, 137, 123, 111,
#'  100, 90, 80, 69, 57, 44, 30, 147, 133, 120, 109,  99, 90, 81, 71, 60,
#'  48, 35, 142, 129, 117, 107,  98, 90, 82, 73, 63, 52, 40, 137, 125, 114,
#'  105,  97, 90, 83, 75, 66, 56, 45, 132, 121, 111, 103,  96, 90, 84, 77,
#'  69, 60, 50, 127, 117, 108, 101,  95, 90, 85, 79, 72, 64, 55, 122, 113,
#'  105,  99,  94, 90, 86, 81, 75, 68, 60, 117, 109, 102,  97,  93, 90, 87,
#'  83, 78, 72, 65, 112, 105,  99,  95,  92, 90, 88, 85, 81, 76, 70, 107,
#'  101,  96,  93,  91, 90, 89, 87, 84, 80, 75)
#'  model.direction1 <- (model.direction1)*pi/180
#'
#'  # Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
#'  from GRF (Range=4, spherical covariance).
#'  set.seed(666)
#'  crf1 <- SimulateCRF(CircDistr = "vM", Rho = sqrt(0.5), Range = 4, xy = xy, OverFit = TRUE)
#'  crf1
SimulateCRF <- function(N=100,
                        CircDistr = c("U", "vM", "WrC", "Tri", "Card"),
                        Rho,
                        Mu=0,
                        Range,
                        Ext=1,
                        CovModel = c("circular", "matern", "exponential", "gaussian", "spherical",  "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget"),
                        xy=NULL,
                        Anisotropy=NULL,
                        OverFit=FALSE,
                        Resolution=0.01) {
  CircDistr <- match.arg(CircDistr)
  CovModel <- match.arg(CovModel)

	if(CircDistr == "U") Rho = 0
	if(abs(Mu) > pi) stop("abs(Mu) <= pi")

	if(!is.null(xy)) {
		if(!inherits(xy, "matrix")) xy <- as.matrix(xy)
		if(dim(xy)[2] !=2) stop("Grid not a N x 2 matrix")
		N <- dim(xy)[1]}
	if(!is.null(Anisotropy)) {
		if(length(Anisotropy) != 2) stop("Anisotropy is not a 2 element vector. See geoR Help")
	}

	if(N <=0 | Rho < 0 | Range < 0| Ext <=0 | Resolution <= 0) stop("Improper numeric input")

	direction <- vector(mode="numeric", length=N)

	# require(CircStats)
	# require(geoR)
	# Standard normal GRF, see Help geoR grf
	if (is.null(xy)) {
  GRF <- grf(
    n = N,
    xlims = c(0, Range * Ext),
    ylims = c(0, Range * Ext),
    cov.model = CovModel,
    nugget = 0,
    cov.pars = c(1, Range),
    aniso.pars = Anisotropy,
    RF = TRUE,
    messages = FALSE
  )
} else {
  GRF <- grf(
    grid = xy,
    cov.model = CovModel,
    nugget = 0,
    cov.pars = c(1, Range),
    aniso.pars = Anisotropy,
    # RF=TRUE,
    messages = FALSE
  )
}

	XY <- GRF$coords # N x 2 matrix
	x <- XY[,1]; y <- XY[,2]
	Z <- GRF$data # Vector of GRV
	if(OverFit) {Z <- (Z - mean(Z))/sd(Z); GRF$data <- Z}
	CumProbZ <- pnorm(Z, mean=0, sd=1, lower.tail = TRUE)

	if(CircDistr=="U") {direction <- -pi + 2*pi*CumProbZ} else
	if(CircDistr == "Tri")
	{
		if(Rho==0 | Rho > 4/pi^2) stop("Tri: 0 < Rho <= 4/pi^2")
		filter <- CumProbZ < 0.5
		u1 <- CumProbZ[filter]
		a <- Rho/8
		b <- (4+pi^2*Rho)/(8*pi)
		c <- 0.5 - u1
		q <- -.5*(b+sqrt(b^2-4*a*c))
		direction[filter] <- c/q

		u2 <- CumProbZ[!filter]
		a <- -Rho/8
		b <- (4+pi^2*Rho)/(8*pi)
		c <- 0.5 - u2
		q <- -.5*(b+sqrt(b^2-4*a*c))
		direction[!filter]<- c/q
	} else
	{
		# For OTHER circular distributions compute table of circular CDF and interpolate
		CircScale <- seq(-pi, pi, length=2*pi/Resolution)
		# With resolution=.01, circular support from -pi to +pi has 629 elements, delta ~0.01000507, CircScale[315] = 0
		n <- length(CircScale)
		if(CircDistr == "vM")
		{
			if(Rho==0 | Rho >= 1) stop("vM: 0 < Rho < 1")
			CircProb <- rep(-1, n)
			Kappa=A1inv(Rho) # N. I Fisher, Statistical Analysis of Circular Data, 2000 p. 49
			# As direction increases from -pi, pvm increases from .5
			for(i in 1:length(CircScale)) CircProb[i] <- pvm(CircScale[i], mu=0, kappa=Kappa)
			filter <- CircScale < 0
			CircProb[filter] <- CircProb[filter] - 0.5
			CircProb[!filter] <- CircProb[!filter] + 0.5
		} else
		if(CircDistr == "Card")
		{
			if(Rho==0 | Rho > 0.5) stop("Cardioid: 0 < Rho <= 0.5")
			CircProb <- (CircScale + pi + 2*Rho*sin(CircScale))/(2*pi)
		} else
		if(CircDistr == "WrC")
		{
			if(Rho==0 | Rho >= 1) stop("Wrapped Cauchy: 0 < Rho < 1 ")
			Angles1 <- CircScale[CircScale < 0]
			Angles2 <- CircScale[CircScale >= 0]
			prob1 <- 0.5 - acos(((1+Rho^2)*cos(Angles1) - 2*Rho)/(1 + Rho^2 - 2*Rho * cos(Angles1)))/(2*pi)
			prob2 <- 0.5 + acos(((1+Rho^2)*cos(Angles2) - 2*Rho)/(1 + Rho^2 - 2*Rho * cos(Angles2)))/(2*pi)
			CircProb <-c(prob1, prob2)
		}
		CircProb[1] <- 0; CircProb[n] <- 1

		# Interpolation
		DeltaTh <- CircScale[2] + pi
		for(i in 1:N)
		{
			p <- CumProbZ[i]  # Cumulative prob of GRV
			a <- max((1:n)[CircProb <= p]) # Index
			if(a==n) {r <- 0} else
			{
				if(CircProb[a]==p) {r <- 0} else {r <- (p -CircProb[a])/( CircProb[a+1] -CircProb[a])}
			}
			direction[i] <- CircScale[a] + r*DeltaTh
		}
	}

	direction <- direction + Mu
	return(list(x=x, y=y, direction=direction, Z=Z))
}
