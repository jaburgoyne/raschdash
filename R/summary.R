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
                                indices,
                                function(i) {
                                        apply(
                                                as.matrix(stanfit, par)[, i],
                                                1,
                                                sum
                                        )
                                }
                        )
                }
        }
        log_lik <- .collapse("log_lik")
        observed_scores <-
                if (is_null(indices)) {
                        x$data$observed_score
                } else {
                        purrr::map_dbl(
                                indices,
                                function(i) {
                                        sum(
                                                dplyr::pull(
                                                        dplyr::slice(x$data, i),
                                                        .data$observed_score
                                                )
                                        )
                                }
                        )
                }
        y <- matrix(observed_scores, nrow(log_lik), ncol(log_lik), byrow = TRUE)
        log_lik_rep <- .collapse("log_lik_rep")
        y_rep <- .collapse("y_rep")
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
                                new_eloo(
                                        loo$pointwise[, "elpd_loo"] / -log(2),
                                        loo::pareto_k_values(loo)
                                )
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
                        )
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
                                if (has_name(grouped_data, "row")) {
                                        grouped_data$row
                                } else {
                                        NULL
                                },
                        ## By default, use LOO only for complete observations.
                        use_loo =
                                if (is_null(use_loo)) {
                                        !has_name(grouped_data, "row")
                                } else {
                                        use_loo
                                }
                )
        )
}
