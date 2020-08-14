.heaviside_difference <- function(observed, expected) {
        0.5 * sign(observed - expected) + 0.5
}

.loo_statistics <- function(stanfit, grouped_data) {
        collapse
        loo <- loo::loo(stanfit, save_psis = TRUE)
        .loo_expectation <- function(x) {
                loo::E_loo(
                        x = x,
                        psis_object = loo$psis_object,
                        ## E_loo() issues an annoying warning in the absence
                        ## of log_ratios.
                        log_ratios = -as.matrix(stanfit, "log_lik")
                ) %>%
                        purrr::pluck("value")
        }
        dplyr::tibble(
                expected_score =
                        as.matrix(stanfit, "y_rep") %>%
                                .loo_expectation(),
                information_content =
                        purrr::pluck(loo, "pointwise") %>%
                                dplyr::as_tibble() %>%
                                purrr::pluck("elpd_loo") %>%
                                magrittr::divide_by(-log(2)),
                entropy =
                        as.matrix(stanfit, "log_lik_rep") %>%
                                .loo_expectation() %>%
                                magrittr::divide_by(-log(2)),
                p_score =
                        purrr::pluck(x, "data", "observed_score") %>%
                                matrix(
                                        nrow = nrow(loo), ncol = ncol(loo),
                                        byrow = TRUE
                                ) %>%
                                .heaviside_difference(
                                        as.matrix(stanfit, "y_rep")
                                ) %>%
                                .loo_expectation(),
                ## Subtract log_lik from log_lik_rep instead of the other way
                ## around because this statistic is about information content,
                ## not likelihood.
                p_information =
                        as.matrix(stanfit, "log_lik_rep") %>%
                                .heaviside_difference(
                                        as.matrix(stanfit, "log_lik")
                                ) %>%
                                .loo_expectation(),
                pareto_k = loo::pareto_k_values(loo)
        )
}

#' @importFrom rlang !!!
#' @export
summary.rdfit <- function(object, ...) {
        .add_group_observations <- function(df) {
                group_observations <-
                        dplyr::filter(
                                dplyr::ungroup(df),
                                vctrs::vec_equal_na(.data$person)
                        )
                person_observations <-
                        dplyr::filter(
                                dplyr::ungroup(df),
                                !vctrs::vec_equal_na(.data$person)
                        )
                group_observations %>%
                        dplyr::select(-.data$person) %>%
                        dplyr::inner_join(
                                dplyr::distinct(
                                        person_observations,
                                        .data$group,
                                        .data$person
                                )
                        ) %>%
                        dplyr::bind_rows(person_observations) %>%
                        dplyr::group_by(!!!dplyr::groups(df))
        }
        .regroup <- function(numbered_df) {
                if (dplyr::n_groups(numbered_df) <= 1) {
                        dplyr::select(numbered_df, -.data$row)
                } else {
                        group_vars <- dplyr::group_vars(numbered_df)
                        ## Add group observations to persons if both are
                        ## present in the grouping variables.
                        if (
                                all(
                                        vctrs::vec_in(
                                                needles = c("group", "person"),
                                                haystack = group_vars
                                        )
                                )
                        ) {
                                extended_df <-
                                        .add_group_observations(numbered_df)
                        } else {
                                extended_df <- numbered_df
                        }
                        extended_df %>%
                                dplyr::summarise(
                                        dplyr::across(
                                                c(
                                                        .data$max_score,
                                                        .data$observed_score
                                                ),
                                                sum
                                        ),
                                        dplyr::across(.data$row, list)
                                )
                }
        }
        ## TODO: Document the output.
        grouped_data <-
                purrr::pluck(object, "data") %>%
                dplyr::as_tibble() %>%
                dplyr::select(!dplyr::starts_with("stan")) %>%
                dplyr::mutate(row = dplyr::row_number()) %>%
                dplyr::group_by(...) %>%
                .regroup() %>%
                dplyr::ungroup()
        grouped_data
        # dplyr::bind_cols(
        #         grouped_data,
        #         .loo_statistics(purrr::pluck(object, "stanfit"), grouped_data))
}
