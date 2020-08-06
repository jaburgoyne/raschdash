#' @importFrom purrr is_character
#' @importFram tibble has_name
new_rdfit <- function(data,
                      cohort_colname,
                      group_colname, person_colname,
                      testlet_colname, item_colname,
                      max_score_colname, obs_score_colname) {
        data <- tibble::as_tibble(data)
        stopifnot(
                has_name(
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

new_rdfit <- function(data,
                      cohort_colname    = 'Cohort',
                      group_colname     = 'Group',
                      person_colname    = 'Student',
                      testlet_colname   = 'Testlet',
                      item_colname      = 'Item',
                      max_score_colname = 'Max',
                      obs_score_colname = 'Score')  {
