#' Constructor for \code{rdfit}
#'
#' Construct an \code{rdfit} object.
#'
#' @param data a data frame with columns for cohort, group, person, testlet,
#'             item, maximum score, and observed score
#' @param cohort_colname name of the cohort column
#' @param group_colname name of the group column
#' @param person_colname name of the person column
#' @param testlet_colname name of the testlet column
#' @param item_colname name of the item column
#' @param max_score_colname name of the maximum score column
#' @param obs_score_colname name of the observed score column
#'
#' @return An \code{rdfit} object.
#'
#' @export
new_rdfit <- function(data,
                      cohort_colname,
                      group_colname, person_colname,
                      testlet_colname, item_colname,
                      max_score_colname, obs_score_colname) {
        data <- tibble::as_tibble(data)
        stopifnot(
                tibble::has_name(
                        data,
                        c(
                                cohort_colname,
                                group_colname, person_colname,
                                testlet_colname, item_colname,
                                max_score_colname, obs_score_colname
                        )
                )
        )

}

rdfit <- function(data,
                  cohort_colname    = 'Cohort',
                  group_colname     = 'Group',
                  person_colname    = 'Student',
                  testlet_colname   = 'Testlet',
                  item_colname      = 'Item',
                  max_score_colname = 'Max',
                  obs_score_colname = 'Score')  {
}
