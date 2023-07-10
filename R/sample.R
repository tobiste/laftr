#' Example dataset
#'
#' Taís's FCT measurements
#'
#' @docType data
#'
#' @usage data('sample')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("sample")
#' head("sample")
"sample"


#' Example zeta dataset
#'
#' Taís's DUR measurements
#'
#' @docType data
#'
#' @usage data('standard')
#'
#' @format An object of class \code{"tibble"}
#'
#' @keywords datasets
#' @examples
#' data("standard")
#' head("standard")
"standard"


# library(dplyr)
# sample <- readxl::read_xlsx('data/example.xlsx', sheet = 'sample')  %>%
#   rename(U = UCa, sU = sUCa)
#
# # load the calibration measurements for zeta factor
# standard <- readxl::read_xlsx('data/example.xlsx', sheet = 'standard') %>%
#   rename(U = UCa, sU = sUCa)
#
#
# save(sample, file = "data/sample.RData", compress = TRUE, ascii = TRUE)
# save(standard, file = "data/standard.RData", compress = TRUE, ascii = TRUE)
