## Construct Trend Model of 121 locations
xy <- expand.grid(1:11, 1:11) # grid
x1 <- xy[, 1]
y1 <- xy[, 2]
model.direction1 <- matrix(data = c(
  157, 141, 126, 113, 101, 90, 79, 67, 54, 40, 25, 152, 137, 123, 111, 100, 90, 80, 69, 57, 44, 30,
  147, 133, 120, 109, 99, 90, 81, 71, 60, 48, 35, 142, 129, 117, 107, 98, 90, 82, 73, 63, 52, 40,
  137, 125, 114, 105, 97, 90, 83, 75, 66, 56, 45, 132, 121, 111, 103, 96, 90, 84, 77, 69, 60, 50,
  127, 117, 108, 101, 95, 90, 85, 79, 72, 64, 55, 122, 113, 105, 99, 94, 90, 86, 81, 75, 68, 60,
  117, 109, 102, 97, 93, 90, 87, 83, 78, 72, 65, 112, 105, 99, 95, 92, 90, 88, 85, 81, 76, 70,
  107, 101, 96, 93, 91, 90, 89, 87, 84, 80, 75
), ncol = 11, byrow = TRUE)

model.direction1 <- as.vector(model.direction1) * pi / 180

## Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
## from GRF (Range=4, spherical covariance).
set.seed(666)
crf1 <- SimulateCRF(
  CircDistr = "vM", Rho = sqrt(0.5), Range = 4, CovModel = "spherical",
  Grid = xy, OverFit = TRUE
)

# Make sample
sample.direction1 <- model.direction1 + crf1$direction

## Fit An Appropriate Model
FitHoriz1 <- lm(cos(sample.direction1) ~ (x1 + y1))
FitVert1 <- lm(sin(sample.direction1) ~ (x1 + y1))
fitted.direction1 <- atan2(FitVert1$fitted.values, FitHoriz1$fitted.values)

## Compute Residuals
resids1 <- CircResidual(
  X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1,
  Plot = FALSE
)

## Output list of cosineogram coordinates for fitting analytically
cosineogram.out <- CosinePlots(
  x = resids1$x, y = resids1$y, directions = resids1$direction,
  Lag.n.Adj = 1, BinWAdj = 1, Plot = FALSE, Cloud = FALSE, Model = FALSE
)

## Cosineocloud
CosinePlots(
  x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
  Plot = TRUE, Cloud = TRUE
)

## Cosineogram
CosinePlots(
  x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
  Plot = TRUE, Cloud = FALSE, Model = FALSE
)

abline(h = 0.56, col = 2)
abline(v = 4, col = 2)

## Fit cosine Models
CosinePlots(
  x = resids1$x, y = resids1$y, directions = resids1$direction, Lag.n.Adj = 1, BinWAdj = 1,
  Plot = TRUE, Cloud = FALSE, Model = TRUE, nugget = 0, Range = 4.0, sill = 0.56, x.legend = .2,
  y.legend = 0.3, xlim = c(0, 8), ylim = c(0, 1)
)
