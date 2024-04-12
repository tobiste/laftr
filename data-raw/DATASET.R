## code to prepare `DATASET` dataset goes here
library(dplyr)

sample <- readxl::read_xlsx("inst/extdata/example.xlsx", sheet = "sample") |>
  rename(U = UCa, sU = sUCa)

# load the calibration measurements for zeta factor
standard <- readxl::read_xlsx("inst/extdata/example.xlsx", sheet = "standard") |>
  rename(U = UCa, sU = sUCa)


# load('data/standard.RData')
# load('data/sample.RData')


usethis::use_data(sample, overwrite = TRUE)
usethis::use_data(standard, overwrite = TRUE)
