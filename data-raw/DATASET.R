## code to prepare `DATASET` dataset goes here
load("data-raw/OceanWind.rda")
usethis::use_data(OceanWind, overwrite = TRUE)

load("data-raw/WorldMask.rda")
usethis::use_data(WorldMask, overwrite = TRUE)
