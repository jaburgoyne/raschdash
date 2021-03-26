# Special functions ------------------------------------------------------------

.heaviside_difference <- function(observed, expected) {
  0.5 * sign(observed - expected) + 0.5
}

# Internal data manipulation  --------------------------------------------------

.add_entropies <- function(df, stanfit) {
  .row_entropies <- function(m) {
    apply(
      m,
      2,
      function(v)
        entropy::entropy(
          table(v),
          unit = "log2",
          verbose = FALSE
        )
    )
  }
  max_max_score <- max(df$max_score)
  n_groups <- abs(min(df$stan_person))
  n_persons <- max(df$stan_person)
  ## TODO: Can these be integrated into LOO somehow?
  posterior_entropy <- .row_entropies(as.matrix(stanfit, "y_rep"))
  prior_group_entropy <-
    .row_entropies(as.matrix(stanfit, "y_prior_group"))
  prior_person_entropy <-
    .row_entropies(as.matrix(stanfit, "y_prior_person"))
  prior_group_item_entropy <-
    .row_entropies(as.matrix(stanfit, "y_prior_group_item")) %>%
    matrix(max_max_score, n_groups)
  prior_item_entropy <-
    .row_entropies(as.matrix(stanfit, "y_prior_item")) %>%
    matrix(max_max_score, n_persons)
  df  %>%
    mutate(
      entropy = posterior_entropy,
      entropy_prior_person =
        case_when(
          stan_person < 0 ~ prior_group_entropy[stan_item],
          TRUE ~ prior_person_entropy[stan_item]
         ),
      entropy_prior_item =
        map2_dbl(
          max_score, stan_person,
          ~{
            if (.y < 0) {
              prior_group_item_entropy[.x, abs(.y)]
            } else {
              prior_item_entropy[.x, .y]
            }
          }
        ),
      iota_person = entropy_prior_person - entropy,
      iota_item = entropy_prior_item - entropy
    ) %>%
    select(-entropy_prior_person, -entropy_prior_item)
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
    return(
      dplyr::mutate(
        df,
        dplyr::across(
          c(dplyr::ends_with("_score"), entropy, dplyr::starts_with("iota_")),
          ~{.x * weight}
        )
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
      c(dplyr::ends_with("_score"), entropy, dplyr::starts_with("iota_")),
      ~sum(.x * weight)
    ),
    dplyr::across(c(row, weight), list),
    ## The stan_ variables are only useful for retrieving
    ## calibrations, in which case they will be unique. Choosing the
    ## first value is a cheap solution to do *something* for the
    ## other cases.
    dplyr::across(dplyr::starts_with("stan_"), dplyr::first),
  )
}

.add_calibrations <- function(df, object, mean, sd) {
  stanfit <- object$stanfit
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
  calibrations <- as.matrix(stanfit, stan_par)[, purrr::pluck(df, stan_id)]
  dplyr::mutate(
    df,
    !!stan_par := mean + sd * apply(calibrations, 2, stats::median),
    mad = sd * apply(calibrations, 2, stats::mad),
    `5%` = mean + sd * apply(calibrations, 2, stats::quantile, 0.05),
    `25%` = mean + sd * apply(calibrations, 2, stats::quantile, 0.25),
    `75%` = mean + sd * apply(calibrations, 2, stats::quantile, 0.75),
    `95%` = mean + sd * apply(calibrations, 2, stats::quantile, 0.95)
  )
}

.add_pp_statistics <- function(df, stanfit, use_loo, mark) {
  indices <- purrr::pluck(df, "row")
  weights <- purrr::pluck(df, "weight")
  .collapse <-
    function(par) {
      m <- as.matrix(stanfit, par)
      mapply(
        function(i, w) {
          if (vec_size(i) > 1) {
            apply(sweep(m[, i], 2, w, `*`), 1, sum)
          } else {
            w * m[, i]
          }
        },
        indices,
        weights
      )
    }
  log_lik <- .collapse("log_lik")
  y_rep <- .collapse("y_rep")
  y <-
    matrix(
      pluck(df, "observed_score"),
      nrow(log_lik), ncol(log_lik),
      byrow = TRUE
    )
  if (use_loo) {
    loo <-
      loo::loo(
        log_lik,
        r_eff =
          loo::relative_eff(
            exp(log_lik),
            chain_id =
              rep(
                1:ncol(object$stanfit),
                nrow(object$stanfit)
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
  p_score <- .expectation(.heaviside_difference(y, y_rep))
  dplyr::mutate(
    df,
    expected_score = .expectation(y_rep),
    expected_mark = mark(.data$expected_score / .data$max_score),
    information_content =
      if (use_loo) {
        new_eloo(
          loo$pointwise[, "elpd_loo"] / -log(2),
          loo::pareto_k_values(loo)
        )
      } else {
        .expectation(log_lik) / -log(2)
      },
    p_score =
      purrr::map2_dbl(indices, weights, ~{sum(p_score[.x] * .y) / sum(.y)}),
    infit = information_content / entropy
  )
}

# Exports ----------------------------------------------------------------------

#' @export
print.rdfit <- function(x, ...) print(x$stanfit, ...)

#' @export
summary.rdfit <- function(object, ..., use_loo = NULL,
                          mean = 0, sd = 1, mark = identity) {
  ## TODO: Document the output. Note that slope of infit vs qlogis(p_info)
  ##       is almost exactly 8, implying that plogis(4) and plogis(8)
  ##       should be infit cutoffs, and plogis(-4) at the bottom.
  purrr::pluck(object, "data") %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    .add_entropies(pluck(object, "stanfit")) %>%
    dplyr::group_by(...) %>%
    .aggregate_scores() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      observed_mark = mark(.data$observed_score / .data$max_score)
    ) %>%
    .add_calibrations(object, mean, sd) %>%
    .add_pp_statistics(
      pluck(object, "stanfit"),
      ## By default, use LOO only for complete observations.
      use_loo =
        if (is_null(use_loo)) {
          ...length() == 0
        } else {
          use_loo
        },
      mark = mark
    ) %>%
    dplyr::select(
      !dplyr::any_of(c("row", "weight")) & !dplyr::starts_with("stan_")
    )
}
