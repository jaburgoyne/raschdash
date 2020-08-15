.heaviside_difference <- function(observed, expected) {
        0.5 * sign(observed - expected) + 0.5
}

.loo_statistics <- function(x, indices = NULL, use_loo = FALSE) {
        stanfit <- purrr::pluck(x, "stanfit")
        .collapse <- function(par) {
                if (rlang::is_null(indices)) {
                        as.matrix(stanfit, par)
                } else {
                        sapply(
                                X = indices,
                                FUN = function(i) {
                                        apply(
                                                X = as.matrix(stanfit, par)[,i],
                                                MARGIN = 1,
                                                FUN = sum
                                        )
                                }
                        )
                }
        }
        log_lik <- .collapse("log_lik")
        observed_scores <-
                indices %>%
                purrr:::map_dbl(
                        function(i) {
                                x %>%
                                        purrr::pluck("data") %>%
                                        dplyr::slice(i) %>%
                                        dplyr::pull("observed_score") %>%
                                        sum()
                        }
                )
        y <-
                matrix(
                        data = observed_scores,
                        nrow = nrow(log_lik),
                        ncol = ncol(log_lik),
                        byrow = TRUE
                )
        log_lik_rep <- .collapse("log_lik_rep")
        y_rep <- .collapse("y_rep")
        ## Right now, using LOO is unreliable for grouped structures, except
        ## possibly for items. The issue is well known for hierarchical data
        ## (Vehtari et al., 2017). Using the simple mean is a solution for now,
        ## but in the future, full cross-validation or improved LOO may be
        ## worthwhile.
        if (use_loo) {
                loo <-
                        loo::loo(
                                log_lik,
                                r_eff =
                                        loo::relative_eff(
                                                x = exp(log_lik),
                                                chain_id =
                                                        rep(
                                                                1:ncol(stanfit),
                                                                nrow(stanfit)
                                                        )
                                        ),
                                save_psis = TRUE
                        )
                .loo_expectation <- function(x) {
                        loo::E_loo(
                                x = x,
                                psis_object = loo$psis_object,
                                ## E_loo() issues an annoying warning in the
                                ## absence of log_ratios.
                                log_ratios = -log_lik,
                        ) %>%
                                purrr::pluck("value")
                }
        } else {
                .loo_expectation <- function(x) apply(x, 2, mean)
        }
        dplyr::tibble(
                expected_score = y_rep %>% .loo_expectation(),
                information_content =
                        if (use_loo) {
                                loo %>%
                                        purrr::pluck("pointwise") %>%
                                        dplyr::as_tibble() %>%
                                        purrr::pluck("elpd_loo") %>%
                                        magrittr::divide_by(-log(2))
                        } else {
                                log_lik %>% .loo_expectation() %>%
                                magrittr::divide_by(-log(2))
                        },
                entropy =
                        log_lik_rep %>%
                        .loo_expectation() %>%
                        magrittr::divide_by(-log(2)),
                p_score =
                        .heaviside_difference(y, y_rep) %>%
                        .loo_expectation(),
                ## Subtract log_lik from log_lik_rep instead of the other way
                ## around because this statistic is about information content,
                ## not likelihood.
                p_information =
                        .heaviside_difference(log_lik_rep, log_lik) %>%
                        .loo_expectation()
        ) %>%
                dplyr::bind_cols(
                        if (use_loo) {
                                dplyr::tibble(
                                        pareto_k = loo::pareto_k_values(loo)
                                )
                        } else {
                                NULL
                        }
                )
}

#' @importFrom rlang !!!
#' @export
summary.rdfit <- function(object, ..., use_loo = NULL) {
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
                                ),
                                by = "group"
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
                        extended_df <-
                                if (vctrs::vec_in("group", group_vars)) {
                                        numbered_df %>%
                                                .add_group_observations()
                                } else if (
                                        vctrs::vec_in("person", group_vars)
                                ) {
                                        numbered_df %>%
                                                filter(
                                                        !vctrs::vec_equal_na(
                                                                .data$person
                                                        )
                                                )
                                } else {
                                        numbered_df
                                }
                        extended_df %>%
                                dplyr::summarise(
                                        n_obs = dplyr::n(),
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
        ## TODO: Document the output. Note that slop of infit vs qlogis(p_info)
        ##       is almost exactly 8, implying that plogis(4) and plogis(8)
        ##       should be infit cutoffs, and plogis(-4) at the bottom.
        grouped_data <-
                purrr::pluck(object, "data") %>%
                dplyr::as_tibble() %>%
                dplyr::select(!dplyr::starts_with("stan")) %>%
                dplyr::mutate(row = dplyr::row_number()) %>%
                dplyr::group_by(...) %>%
                .regroup() %>%
                dplyr::ungroup()
        dplyr::bind_cols(
                dplyr::select(grouped_data, !dplyr::any_of("row")),
                .loo_statistics(
                        x = object,
                        indices =
                                if (tibble::has_name(grouped_data, 'row')) {
                                        grouped_data$row
                                } else {
                                        NULL
                                },
                        ## By default, use LOO only for complete observations.
                        use_loo =
                                if (rlang::is_null(use_loo)) {
                                        !tibble::has_name(grouped_data, 'row')
                                } else {
                                        use_loo
                                }
                )
        )
}