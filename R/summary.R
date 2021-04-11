# Special functions ------------------------------------------------------------

.heaviside_difference <- function(observed, expected) {
  0.5 * sign(observed - expected) + 0.5
}

# Internal data manipulation  --------------------------------------------------

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
    return(
      dplyr::mutate(
        df,
        dplyr::across(dplyr::ends_with("_score"), ~{.x * weight})
      )
    )
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
    dplyr::across(
      dplyr::ends_with("_score"),
      ~sum(.x * .data$weight)
    ),
    dplyr::across(c("row", "weight"), list),
    ## The stan_ variables are only useful for retrieving
    ## calibrations, in which case they will be unique. Choosing the
    ## first value is a cheap solution to do *something* for the
    ## other cases.
    dplyr::across(dplyr::starts_with("stan_"), dplyr::first),
  )
}

.add_calibrations <- function(df, draws, mean, sd) {
  par <-
    if (has_name(df, "person") && !vec_duplicate_any(df$person)) {
      "person"
    } else if (has_name(df, "item") && !vec_duplicate_any(df$item)) {
      "item"
    } else if (has_name(df, "group") && !vec_duplicate_any(df$group)) {
      "group"
    } else if (has_name(df, "testlet") && !vec_duplicate_any(df$testlet)) {
      "testlet"
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
  calibrations <- posterior::subset_draws(draws, stan_par)
  prior_calibrations <-
    posterior::extract_variable(draws, stringr::str_c("prior_", stan_par))
  ## For a consistent estimator, k -> \infty as S -> \infty
  k <- round(sqrt(vec_size(prior_calibrations)))
  df %>%
    dplyr::inner_join(
      calibrations %>%
        summary(
          median = stats::median,
          mad = stats::mad,
          ~posterior::quantile2(.x, c(0.05, 0.25, 0.75, 0.95)),
          iota = ~{FNN::KL.divergence(c(.x), prior_calibrations, k)[k]}
        ) %>%
        tidyr::separate(
          "variable",
          c(NA, stan_id, NA),
          sep = "\\]|\\[",
          convert = TRUE
        ),
      by = stan_id
    ) %>%
    dplyr::mutate(
      dplyr::across(c("median", "q5", "q25", "q75", "q95"), ~{mean + sd * .x}),
      dplyr::across("mad", ~{sd * .x}),
      dplyr::across("iota", ~{.x / log(2)})
    ) %>%
    dplyr::rename(!!stan_par := "median")
}

.add_pp_statistics <- function(df, draws, use_loo, cores, mark) {
  indices <- purrr::pluck(df, "row")
  weights <- purrr::pluck(df, "weight")
  .collapse <-
    function(par) {
      m <- posterior::subset_draws(draws, par)
      mapply(
        function(i, w) apply(sweep(m[, , i], 3, w, `*`), 1:2, sum),
        indices,
        weights
      )
    }
  log_lik <- .collapse("log_lik")
  log_lik_rep <- .collapse("log_lik_rep")
  y_rep <- .collapse("y_rep")
  y <-
    matrix(
      purrr::pluck(df, "observed_score"),
      nrow(y_rep), ncol(y_rep),
      byrow = TRUE
    )
  if (use_loo) {
    loo <-
      loo::loo(
        log_lik,
        r_eff =
          loo::relative_eff(
            exp(log_lik),
            cores = cores,
            ## This ordering is dependent on how mapply() works earlier. :(
            ## But leaving draws in arrays does not work, because the PSIS
            ## object from loo currently reverts to matrix form regardless of
            ## the input.
            chain_id =
              rep(
                1:posterior::nchains(draws),
                each = posterior::niterations(draws)
              )
          ),
        cores = cores,
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
    expected_mark = mark(.data$expected_score / .data$max_score),
    information =
      if (use_loo) {
        new_eloo(
          loo$pointwise[, "elpd_loo"] / log(0.5),
          loo::pareto_k_values(loo)
        )
      } else {
        .expectation(log_lik) / log(0.5)
      },
    entropy = .expectation(log_lik_rep) / log(0.5),
    p_score = .expectation(.heaviside_difference(y, y_rep)),
    infit = .data$information / .data$entropy
  )
}

# Exports ----------------------------------------------------------------------

#' @export
print.rdfit <- function(x, ...) print(x$draws, ...)

#' @export
summary.rdfit <- function(object, ...,
                          use_loo = NULL, cores = getOption("mc.cores", 1),
                          mean = 0, sd = 1, mark = identity) {
  ## TODO: Document the output. Note that slope of infit vs qlogis(p_info)
  ##       is almost exactly 8, implying that plogis(4) and plogis(8)
  ##       should be infit cutoffs, and plogis(-4) at the bottom.
  purrr::pluck(object, "data") %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::group_by(...) %>%
    .aggregate_scores() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      observed_mark = mark(.data$observed_score / .data$max_score)
    ) %>%
    .add_calibrations(purrr::pluck(object, "draws"), mean, sd) %>%
    .add_pp_statistics(
      purrr::pluck(object, "draws"),
      ## By default, use LOO only for complete observations.
      use_loo =
        if (is_null(use_loo)) {
          ...length() == 0
        } else {
          use_loo
        },
      cores = cores,
      mark = mark
    ) %>%
    dplyr::select(
      !dplyr::any_of(c("row", "weight")) & !dplyr::starts_with("stan_")
    )
}
