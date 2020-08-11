# This block is recommended by rstantools.
#
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @useDynLib raschdash
#
# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
#
## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom lifecycle deprecate_soft
## usethis namespace: end
NULL

# To silence notes about `.` in pipes.
# See https://github.com/tidyverse/magrittr/issues/29
utils::globalVariables(".")
