#' `raschdash`: Dashboards for Rasch models in the classroom
#'
#' \lifecycle{experimental}
#' The `raschdash` package provides two tools.
#'   1. Fully Bayesian Rasch models with calibrations and diagnostics for
#'      binary, binomial, and rating-scale items.
#'   2. Templates for visualising the calibrations and diagnostics from these
#'      models, including space to add comments per cohort.
#'
#' @name raschdash-package
#
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
