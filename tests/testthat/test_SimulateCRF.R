xy <- expand.grid(1:11, 1:11) # grid
SimulateCRF(CircDistr = "vM", Rho = sqrt(0.5), Range = 4, Grid = xy, OverFit = TRUE)
SimulateCRF(N = 200, CircDistr = "Card", Rho = 0.4, Range = 5, Ext = 3, CovModel = "exponential")
SimulateCRF(CircDistr = "U", Range = 8, Ext = 3, CovModel = "gaussian")
SimulateCRF(CircDistr = "Tri", Rho = 0.5 * 4 / pi^2, Range = 8, Ext = 3, CovModel = "spherical")
SimulateCRF(CircDistr = "WrC", Rho = sqrt(0.8), Range = 8, Ext = 3, CovModel = "exponential")
SimulateCRF(N = 400, CircDistr = "WrC", Rho = sqrt(0.95), Range = 8, Ext = 3, CovModel = "spherical", Anisotropy = c(pi / 4, 3))
