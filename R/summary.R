.heaviside_difference <- function(observed, expected) {
        0.5 * sign(observed - expected) + 0.5
}

.loo_statistics <- function(x, indices = NULL, use_loo = FALSE) {
        stanfit <- purrr::pluck(x, "stanfit")
        .collapse <- function(par) {
                if (is_null(indices)) {
                        as.matrix(stanfit, par)
                } else {
                        sapply(
                                X = indices,
                                FUN = function(i) {
                                        apply(
                                                X = as.matrix(stanfit, par)[, i],
                                                MARGIN = 1,
                                                FUN = sum
                                        )
                                }
                        )
                }
        }
        log_lik <- .collapse("log_lik")
        observed_scores <-
                purrr:::map_dbl(
                        indices,
                        function(i) {
                                dplyr::summarise(
                                        dplyr::slice(x$data, i),
                                        dplyr::across(
                                                .data$observed_score,
                                                sum
                                        )
                                )
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
                                                exp(log_lik),
                                                chain_id =
                                                        rep(
                                                                1:ncol(stanfit),
                                                                nrow(stanfit)
                                                        )
                                        ),
                                save_psis = TRUE
                        )
                .expectation <- function(x) {
                        eloo(
                                loo::E_loo(
                                        x,
                                        psis_object = loo$psis_object,
                                        log_ratios = -log_lik,
                                )
                        )
                }
        } else {
                .expectation <- function(x) apply(x, 2, mean)
        }
        dplyr::tibble(
                expected_score = .expectation(y_rep),
                information_content =
                        if (use_loo) {
                                loo %>%
                                        purrr::pluck("pointwise") %>%
                                        dplyr::as_tibble() %>%
                                        purrr::pluck("elpd_loo") %>%
                                        magrittr::divide_by(-log(2))
                        } else {
                                .expectation(log_lik) / -log(2)
                        },
                entropy = .expectation(log_lik_rep) / -log(2),
                p_score = .expectation(.heaviside_difference(y, y_rep)),
                ## Subtract log_lik from log_lik_rep instead of the other way
                ## around because this statistic is about information content,
                ## not likelihood.
                p_information =
                        .expectation(
                                .heaviside_difference(log_lik_rep, log_lik)
                        ),
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

#' @export
summary.rdfit <- function(object, ..., use_loo = NULL) {
        .add_group_observations <- function(df) {
                group_observations <-
                        dplyr::filter(
                                dplyr::ungroup(df),
                                vec_equal_na(.data$person)
                        )
                person_observations <-
                        dplyr::filter(
                                dplyr::ungroup(df),
                                !vec_equal_na(.data$person)
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
                                if (vec_in("group", group_vars)) {
                                        .add_group_observations(numbered_df)
                                } else if (vec_in("person", group_vars)) {
                                        dplyr::filter(
                                                numbered_df,
                                                !vec_equal_na(.data$person)
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
                                if (tibble::has_name(grouped_data, "row")) {
                                        grouped_data$row
                                } else {
                                        NULL
                                },
                        ## By default, use LOO only for complete observations.
                        use_loo =
                                if (is_null(use_loo)) {
                                        !tibble::has_name(grouped_data, "row")
                                } else {
                                        use_loo
                                }
                )
        )
}
