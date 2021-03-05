# iota-only functions ----------------------------------------------------------

#' @export
count_strata <- function(iota) 2^mean(iota)

#' @export
markon_reliability <- function(iota) 1 - count_strata(iota)^(-2)

# calibration and iota functions -----------------------------------------------

.calibrate_strata <- function(strata, calibration, iota) {
  n_strata <- count_strata(iota)
  p <- 0.5 + (strata - 0.5 * (round(n_strata) + 1)) / n_strata
  ## The map is inefficient, but without it, it is difficult to resolve
  ## the strictness checks between the tidyverse and the Harrell-verse.
  purrr::map_dbl(
    p,
    function(p) {
      if (p < 0) {
        -Inf
      } else if (p > 1) {
        Inf
      } else {
        Hmisc::wtd.quantile(
          calibration,
          weights = iota,
          normwt = TRUE,
          probs = p
        )
      }
    }
  )
}

#' @export
classify_strata <- function(calibration, iota) {
  n <- round(count_strata(iota))
  breaks <- .calibrate_strata(seq(0.5, n + 0.5, 1), calibration, iota)
  cut(calibration, breaks = unique(c(-Inf, breaks, Inf)))
}

#' @export
stratum_measures <- function(calibration, iota) {
  .calibrate_strata(1:round(count_strata(iota)), calibration, iota)
}
