## Construct Trend Model of 121 locations
xy <- expand.grid(1:11, 1:11) # grid
x1 <- xy[, 1]
y1 <- xy[, 2]
model.direction1 <- c(
  157, 141, 126, 113, 101, 90, 79, 67, 54, 40, 25, 152, 137, 123, 111, 100, 90, 80, 69, 57, 44, 30,
  147, 133, 120, 109, 99, 90, 81, 71, 60, 48, 35, 142, 129, 117, 107, 98, 90, 82, 73, 63, 52, 40,
  137, 125, 114, 105, 97, 90, 83, 75, 66, 56, 45, 132, 121, 111, 103, 96, 90, 84, 77, 69, 60, 50,
  127, 117, 108, 101, 95, 90, 85, 79, 72, 64, 55, 122, 113, 105, 99, 94, 90, 86, 81, 75, 68, 60,
  117, 109, 102, 97, 93, 90, 87, 83, 78, 72, 65, 112, 105, 99, 95, 92, 90, 88, 85, 81, 76, 70,
  107, 101, 96, 93, 91, 90, 89, 87, 84, 80, 75
)
model.direction1 <- as.vector(model.direction1) * pi / 180

plot(x1, y1, ty = "n", xlab = "", ylab = "", asp = 1)
fields::arrow.plot(x1, y1, u = cos(model.direction1), v = sin(model.direction1), true.angle = TRUE, length = .05)


## Compute vM CRF of 121 observations, Rho=sqrt(0.5) so sill about 0.5,
## from GRF (Range=4, spherical covariance).
set.seed(666)
crf1 <- SimulateCRF(
  CircDistr = "vM", Rho = sqrt(0.5), Range = 4, CovModel = "spherical",
  Grid = xy, OverFit = TRUE
)
fields::arrow.plot(x1, y1, u = cos(crf1$direction), v = sin(crf1$direction), true.angle = TRUE, length = .05, col = "blue")


# Make sample
sample.direction1 <- model.direction1 + crf1$direction

## Fit An Appropriate Model
## Code for median polish is contained in Appendix K, Section K.12
FitHoriz1 <- lm(cos(sample.direction1) ~ (x1 + y1))
FitVert1 <- lm(sin(sample.direction1) ~ (x1 + y1))
fitted.direction1 <- atan2(FitVert1$fitted.values, FitHoriz1$fitted.values)

## Compute Residuals
resids1 <- CircResidual(
  X = x1, Y = y1, Raw = sample.direction1, Trend = fitted.direction1,
  Plot = FALSE
)


x2 <- seq(1, 11, by = 0.2)
y2 <- x2 ## Kriging locations
krig2 <- KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y,
  resid.direction = resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56,
  Plot = FALSE
)


## Plot Kriging, Residuals Overploted In Black
plot(krig2$x, krig2$y, ty = "n", xlab = "", ylab = "", xlim = c(5, 8), ylim = c(5, 8), asp = 1)
fields::arrow.plot(krig2$x, krig2$y,
  u = cos(krig2$direction), v = sin(krig2$direction), arrow.ex = 0.06,
  xpd = FALSE, true.angle = TRUE, length = .05, col = "tan"
)
fields::arrow.plot(resids1$x, resids1$y,
  u = cos(resids1$direction), v = sin(resids1$direction), arrow.ex =
    0.09, xpd = FALSE, true.angle = TRUE, length = .05, col = 1
)

KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10)
)


# Repeat with Nugget = 0.15 and 0.3
KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.44, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10)
)
fields::arrow.plot(resids1$x, resids1$y,
  u = cos(resids1$direction), v = sin(resids1$direction), arrow.ex =
    0.09, xpd = F, true.angle = T, length = .05, col = 2
)


KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10), Smooth = TRUE, bandwidth = 0.1
)

KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10), Smooth = TRUE, bandwidth = 2
)

KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10), Smooth = TRUE, bandwidth = 4
)

KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Xlim = c(7, 10), Ylim = c(7, 10), Smooth = TRUE, bandwidth = 10
)


## Plot kriging estimate variability at sample locations on a regular grid
KrigCRF(
  krig.x = x2, krig.y = y2, resid.x = resids1$x, resid.y = resids1$y, resid.direction =
    resids1$direction, Model = RandomFields::RMexp(), Nugget = 0.0, Range = 4, sill = 0.56, Plot = TRUE,
  Smooth = FALSE, PlotVar = TRUE
)

## Plot kriging estimate variability with random sample locations
set.seed(13)
crf6 <- SimulateCRF(
  N = 400, CircDistr = "Card", Rho = 0.4, Range = 4, Ext = 3,
  CovModel = "spherical"
)

## Best fit is spherical with range=2.85 and sill=0.15
CosinePlots(
  x = crf6$x, y = crf6$y, directions = crf6$direction, Lag.n.Adj = 1.5, BinWAdj = 1,
  Plot = TRUE, Cloud = FALSE, Model = TRUE, nugget = 0, Range = 2.85, sill = 0.15, x.legend = .14,
  y.legend = 0.75, xlim = c(0, 6), ylim = c(0, 1)
)

x6 <- y6 <- seq(4, 7, by = 0.02)
# This may take a long time
KrigCRF(
  krig.x = x6, krig.y = y6, resid.x = crf6$x, resid.y = crf6$y, resid.direction = crf6$direction,
  Model = RandomFields::RMspheric(), Nugget = 0.0, Range = 2.85, sill = 0.15, Plot = TRUE, PlotVar = TRUE
)

