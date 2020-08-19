# Special functions ------------------------------------------------------------

.kl_entropy <- function(x) {
        n <- vec_size(x)
        y <- vec_sort(x)
        ## Return entropy in bits rather than nats.
        (mean(log(y[-1] - y[-n])) - digamma(1) + digamma(n)) / log(2)
}

.heaviside_difference <- function(observed, expected) {
        0.5 * sign(observed - expected) + 0.5
}

# Internal data manipulation  --------------------------------------------------

.add_calibrations <- function(df, stanfit) {
        par <-
                if (has_name(df, 'person') && !vec_duplicate_any(df$person)) {
                        'person'
                } else if (has_name(df, 'item')
                           && !vec_duplicate_any(df$item)) {
                        'item'
                } else if (has_name(df, 'group')
                           && !vec_duplicate_any(df$group)) {
                        'group'
                } else if (has_name(df, 'testlet')
                           && !vec_duplicate_any(df$testlet)) {
                        'testlet'
                }
        if (is_null(par)) {
                return(df)
        }
        stan_id <- stringr::str_c("stan_", par)
        stan_par <-
                if (vec_in(par, c("group", "person"))) {
                        stringr::str_c(par, "_ability")
                } else if (vec_in(par, c("testlet", "item"))) {
                        stringr::str_c(par, "_difficulty")
                }
        calibrations <-
                as.matrix(stanfit, stan_par)[, purrr::pluck(df, stan_id)]
        prior_calibrations <-
                as.matrix(stanfit, stringr::str_c("prior_", stan_par))
        dplyr::mutate(
                df,
                !!stan_par := apply(calibrations, 2, median),
                mad = apply(calibrations, 2, mad),
                relative_entropy =
                        magrittr::subtract(
                                apply(prior_calibrations, 2, .kl_entropy),
                                apply(calibrations, 2, .kl_entropy)
                        )
        )
}

.add_pp_statistics <- function(df, stanfit, use_loo = FALSE) {
        indices <- purrr::pluck(df, "row")
        .collapse <-
                if (is_null(indices)) {
                        function(par) as.matrix(stanfit, par)
                } else {
                        function(par) {
                                m <- as.matrix(stanfit, par)
                                sapply(
                                        indices,
                                        function(i) apply(m[, i], 1, sum)
                                )
                        }
                }
        log_lik <- .collapse("log_lik")
        y <-
                matrix(
                        df$observed_score,
                        nrow(log_lik), ncol(log_lik),
                        byrow = TRUE
                )
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
        dplyr::mutate(
                df,
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

.personify_group_observations <- function(df) {
        person_observations <-
                dplyr::ungroup(df) %>%
                dplyr::filter(!vec_equal_na(.data$person))
        group_observations <-
                dplyr::ungroup(df) %>%
                dplyr::filter(vec_equal_na(.data$person)) %>%
                dplyr::select(!dplyr::any_of(c("stan_person", "person"))) %>%
                dplyr::inner_join(
                        dplyr::distinct(
                                person_observations,
                                dplyr::across(
                                        c("group", "stan_person", "person")
                                )
                        ),
                        by = "group"
                )
        dplyr::group_by(
                dplyr::bind_rows(group_observations, person_observations),
                !!!dplyr::groups(df)
        )
}

.aggregate_scores <- function(df) {
        if (dplyr::n_groups(df) <= 1) {
                return(df)
        }
        group_vars <- dplyr::group_vars(df)
        .personify <-
                if (any(vec_in(c("stan_group", "group"), group_vars))) {
                        .personify_group_observations
                } else {
                        identity
                }
        .identify <-
                if (any(vec_in(c("stan_person", "person"), group_vars))) {
                        function(x) {
                                dplyr::filter(x, !vec_equal_na(.data$person))
                        }
                } else {
                        identity
                }

        ## In order to support information about persons with or without respect
        ## to their groups, replicate group observations per person whenever
        ## groups are present, and remove group observations entirely whenever
        ## persons are present without groups. The price of this behaviour is
        ## that it is not possible to look only at group-level observations.
        dplyr::summarise(
                df %>% .personify() %>% .identify(),
                n = dplyr::n(),
                dplyr::across(dplyr::ends_with("_score"), sum),
                dplyr::across("row", list),
                ## The stan_ variables are only useful for retrieving
                ## calibrations, in which case they will be unique. Choosing the
                ## first value is a cheap solution to do *something* for the
                ## other cases.
                dplyr::across(dplyr::starts_with("stan_"), dplyr::first),
        )
}

# Exports ----------------------------------------------------------------------

#' @export
summary.rdfit <- function(object, ..., use_loo = NULL) {
        ## TODO: Document the output. Note that slop of infit vs qlogis(p_info)
        ##       is almost exactly 8, implying that plogis(4) and plogis(8)
        ##       should be infit cutoffs, and plogis(-4) at the bottom.
        purrr::pluck(object, "data") %>%
                dplyr::mutate(row = dplyr::row_number()) %>%
                dplyr::group_by(...) %>%
                .aggregate_scores() %>%
                dplyr::ungroup() %>%
                .add_calibrations(object$stanfit) %>%
                .add_pp_statistics(
                        object$stanfit,
                        ## By default, use LOO only for complete observations.
                        use_loo =
                                if (is_null(use_loo)) {
                                        ...length() == 0
                                } else {
                                        use_loo
                                }
                ) %>%
                dplyr::select(
                        !dplyr::any_of("row") & !dplyr::starts_with("stan_")
                )
}
