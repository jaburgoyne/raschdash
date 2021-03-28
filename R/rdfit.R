#' Fit a `raschdash` model
#'
#' Create `rdfit` objects from observed data.
#'
#' The `rdfit()` function is intended to be the primary method for users to create
#' `rdfit` objects. Most R users work with denormalised *dimensional data* in a
#' single data frame. In the case that the data are already stored in a
#' normalised structure, however, the raw `new_rdfit()` constructor will work.
#' The `new_rdfit()` constructor does not validate whether the data are properly
#' normalised, however, and Stan will throw sometimes cryptic error messages if
#' they are not.
#'
#' @param data a data frame with columns for cohort, group, person, testlet,
#'             item, maximum possible item score, observed score, and weight.
#'             For group-level observations, the person should be `NA`.
#' @param cohort <[tidyr::tidyr_tidy_select]> cohort column. Defaults to a
#'               generic *All* cohort.
#' @param group <[tidyr::tidyr_tidy_select]> group column. Defaults to a
#'              generic *All* group.
#' @param person <[tidyr::tidyr_tidy_select]> person column
#' @param testlet <[tidyr::tidyr_tidy_select]> testlet column. Defaults to a
#'                generic *All* testlet.
#' @param item <[tidyr::tidyr_tidy_select]> item column
#' @param max_score <[tidyr::tidyr_tidy_select]> maximum-score column
#' @param obs_score <[tidyr::tidyr_tidy_select]> observed score column
#' @param weight <[tidyr::tidyr_tidy_select]> weight column
#' @param k maximum score of rating-scale items. Defaults to no rating scales.
#' @param cohorts vector of unique cohort identifiers
#' @param groups vector of unique group identifiers
#' @param persons,person_groups vectors of unique person identifiers (`persons`)
#'                              and the identifiers of their corresponding
#'                              groups (`person_groups`)
#' @param testlets vector of unique testlet identifiers
#' @param items,item_testlets,item_max_scores vectors of unique item identifiers
#'                                            (`items`) and their corresponding
#'                                            testlets (`item_testlets`) and
#'                                            maximum scores(`max_scores`)
#' @param obs_group_cohorts,obs_groups,obs_group_items,obs_group_scores
#'            for the group-level observations, vectors of the corresponding
#'            cohorts (`obs_group_cohorts`), groups (`obs_groups`), items
#'            (`obs_group_items`), and observed scores (`obs_group_scores`)
#' @param group_obs_weights vector of non-negative regression weights per
#'            group-level observation
#' @param obs_person_cohorts,obs_persons,obs_person_items,obs_person_scores
#'            for the person-level observations, vectors of the corresponding
#'            cohorts (`obs_person_cohorts`), persons (`obs_persons`), items
#'            (`obs_person_items`), and observed scores (`obs_person_scores`)
#' @param person_obs_weights vector of non-negative regression weights per
#'            person-level observation
#' @param pars which parameters to save from the Stan model. Defaults to
#'             `standard`, which only saves parameters on the standard scale,
#'             but may also be `logit` to add logit-scale parameters,
#'             or `hyper` to add logit-scale parameters and hyper-parameters.
#'             For debugging, the `all` option adds other raw parameters
#'             specific to the Stan implementation.
#' @param ... additional parameters passed to [rstan::sampling()]
#'
#' @return An object of class `rdfit` with the following attributes.
#' \describe{
#'  \item{`data`}{The input data as passed to Stan.}
#'  \item{`stanfit`}{
#'   The [rstan::stanfit] object for the fitted model. It contains samples of
#'   the following parameters.
#'   * standard-scale parameters
#'       * `group_ability`: group abilities
#'       * `person_ability`: person abilities
#'       * `testlet_difficulty`: testlet difficulties
#'       * `item_difficulty`: item difficulties
#'       * `thresholds`: rating-scale threshold offsets from overall item
#'                       difficulty (if any rating-scale items are
#'                       present)
#'       * `sigma`: standard deviation of person abilities on the logit scale.
#'                  Used to convert between the logit and standard scales.
#'   * posterior predictive distribution
#'       * `y_rep`: predicted score for each observed combination of
#'                  group, person, and item
#'       * `prior_group_ability` (sampled under posterior hyper-parameters)
#'       * `prior_person_ability` (sampled under posterior hyper-parameters)
#'       * `prior_testlet_difficulty` (sampled under posterior hyper-parameters)
#'       * `prior_item_difficulty` (sampled under posterior hyper-parameters)
#'   * log likelihoods
#'       * `log_lik`: log likelihood of observed scores
#'       * `log_lik_rep`: log likelihood of predicted scores
#'       * `log_lik_prior_person`: log likelihood of predicted scores under the
#'                                 prior distribution of person and group
#'                                 parameters (sampled under posterior hyper-
#'                                 parameters)
#'       * `log_lik_prior_item`: log likelihood of predicted scores under the
#'                               prior distribution of item and testlet
#'                               parameters (sampled under posterior hyper-
#'                               parameters)
#'   * logit-scale parameters (when `pars = "logit"` or `pars = "hyper"`)
#'       * `xi`: group abilities
#'       * `eta`: person abilities
#'       * `epsilon`: testlet difficulties
#'       * `delta`: item difficulties
#'       * `tau`: rating-scale threshold offsets from overall item
#'                difficulty (if any rating-scale items are present)
#'    * hyper-parameters (when `pars = "hyper"`)
#'       * `nu`: mean difficulty
#'       * `psi`: scale of group abilities
#'       * `phi`: scale of person abilities (relative to group)
#'       * `theta_epsilon`: scale of testlet difficulties
#'       * `theta_upsilon`: scale of item/threshold difficulties (relative to testlet)
#'  }
#'
#'  \item{`loo`}{
#'   A [loo::loo] object for `stanfit`.
#'  }
#' }
#' @name rdfit-class
NULL

#' @describeIn rdfit-class Create an `rdfit` object from dimensional data.
#' @export
rdfit <- function(data,
                  cohort = "All",
                  group = "All", person,
                  testlet = "All", item, max_score,
                  obs_score,
                  weight,
                  k = 1,
                  ...) {
  obs <-
    ## rdfit objects return fit statistics as tibbles. If the data
    ## cannot be represented as a tibble, then the tidverse errors
    ## (or warnings) are all that is necessary.
    dplyr::as_tibble(data) %>%
    ## Skip column existence checks, because tibble will complain.
    dplyr::transmute(
      ## Fill in a generic All cohort if necessary.
      cohort = {{ cohort }},
      group = {{ group }},
      person = {{ person }},
      testlet = {{ testlet }},
      item = {{ item }},
      ## Rasch models require integral scores, and vctrs
      ## complains nicely.
      max_score =
        vec_cast(
          x = {{ max_score }},
          to = integer(),
          x_arg = "max_score"
        ),
      obs_score =
        vec_cast(
          x = {{ obs_score }},
          to = integer(),
          x_arg = "obs_score"
        ),
      weight =
        vec_cast(
          x = {{ weight }},
          to = numeric(),
          x_arg = "weight"
        )
    ) %>%
    ## Rasch models have no problem with missing observations, but
    ## Stan does not want to see them.
    dplyr::filter(!vec_equal_na(obs_score))
  ## Stan will complain about invalid scores and weights, but the messages will
  ## be opaque for R users.
  if (any(obs$max_score < 1)) {
    abort("Some maximum scores are non-positive.")
  }
  if (any(obs$obs_score > obs$max_score)) {
    abort("Some observed scores are greater than their maxima.")
  }
  if (any(obs$weight < 0)) {
    abort("Some weights are negative.")
  }
  ## Check for NAs and duplicates when preparing data for Stan. Stan's
  ## validation will catch these, but not with helpful messages.
  cohorts <- dplyr::distinct(obs, .data$cohort)
  if (any(vec_equal_na(cohorts$cohort))) {
    abort("Some cohorts are undefined.")
  }
  groups <- dplyr::distinct(obs, .data$group)
  if (any(vec_equal_na(groups$group))) {
    abort("Some groups are undefined.")
  }
  persons <-
    dplyr::filter(
      dplyr::distinct(obs, .data$person, .data$group),
      !vec_equal_na(.data$person)
    )
  ## Persons may be NA because at this point in the function, the group
  ## will not be (i.e., NA persons are guaranteed to refer to group
  ## observations).
  if (vec_duplicate_any(persons$person)) {
    abort("Some persons belong to multiple groups.")
  }
  testlets <- dplyr::distinct(obs, .data$testlet)
  if (any(vec_equal_na(testlets$testlet))) {
    abort("Some testlets are undefined.")
  }
  items <-
    dplyr::distinct(obs, .data$item, .data$testlet, .data$max_score)
  if (any(vec_equal_na(items$item))
  || any(vec_equal_na(items$max_score))) {
    abort("Some items are undefined.")
  }
  if (vec_duplicate_any(items$item)) {
    abort(
      stringr::str_c(
        "Some items belong to multiple testlets ",
        "or have inconsistent max scores."
      )
    )
  }
  ## A negative integer is no serious problem for K (it simply fails to
  ## select any items for rating scale), but Stan needs it to be
  ## non-negative in order to set the dimension of tau.
  k <- max(vec_cast(k, integer(), x_arg = "k"), 1)
  ## The constructor asks for group and individual observations to be
  ## separated.
  group_obs <- dplyr::filter(obs, vec_equal_na(.data$person))
  person_obs <- dplyr::filter(obs, !vec_equal_na(.data$person))
  new_rdfit(
    cohorts = cohorts$cohort,
    groups = groups$group,
    person_groups = persons$group, persons = persons$person,
    testlets = testlets$testlet,
    item_testlets = items$testlet, items = items$item,
    item_max_scores = items$max_score,
    obs_group_cohorts = group_obs$cohort,
    obs_groups = group_obs$group,
    obs_group_items = group_obs$item,
    obs_group_scores = group_obs$obs_score,
    group_obs_weights = group_obs$weight,
    obs_person_cohorts = person_obs$cohort,
    obs_persons = person_obs$person,
    obs_person_items = person_obs$item,
    obs_person_scores = person_obs$obs_score,
    person_obs_weights = person_obs$weight,
    k = k,
    ...
  )
}

#' @describeIn rdfit-class Create an `rdfit` object from normalised data.
#' @export
new_rdfit <- function(cohorts,
                      groups, persons, person_groups,
                      testlets, items, item_testlets, item_max_scores,
                      obs_group_cohorts, obs_groups, obs_group_items,
                      obs_group_scores, group_obs_weights,
                      obs_person_cohorts, obs_persons, obs_person_items,
                      obs_person_scores, person_obs_weights,
                      k = 1,
                      pars = c("standard", "logit", "hyper", "all"),
                      ...) {
  STANDARD_PARS <-
    c(
      "sigma",
      "group_ability", "person_ability",
      "testlet_difficulty", "item_difficulty",
      "thresholds",
      "prior_group_ability", "prior_person_ability",
      "prior_testlet_difficulty", "prior_item_difficulty",
      "y_rep",
      "log_lik", "log_lik_rep"
    )
  LOGIT_PARS <- c("xi", "eta", "epsilon", "delta", "tau")
  HYPER_PARS <-
    c("psi", "phi", "theta_epsilon", "theta_upsilon", "theta_tau")
  pars <- match.arg(pars)
  stan_pars <-
    switch(pars,
      "hyper" = c(STANDARD_PARS, LOGIT_PARS, HYPER_PARS),
      "logit" = c(STANDARD_PARS, LOGIT_PARS),
      "standard" = c(STANDARD_PARS),
      NA
    )
  if (length(testlets) == 1) {
    stan_pars <-
      setdiff(stan_pars, c("testlet_difficulty", "epsilon", "theta_epsilon"))
  }
  if (k == 1) {
    stan_pars <-
      setdiff(stan_pars, c("thresholds", "tau", "theta_tau"))
  }
  if (length(groups) == 1) {
    stan_pars <-
      setdiff(stan_pars, c("group_ability", "xi", "psi"))
  }
  ## Make tidy tibbles and sort for reliable indexing inside Stan.
  groups <-
    dplyr::tibble(
      stan_group = dplyr::dense_rank(groups),
      group = groups
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$stan_group)
  persons <-
    dplyr::tibble(
      stan_person = dplyr::dense_rank(persons),
      person = persons,
      group = person_groups
    ) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(groups, by = "group") %>%
    dplyr::arrange(.data$stan_person)
  testlets <-
    dplyr::tibble(
      stan_testlet = dplyr::dense_rank(testlets),
      testlet = testlets
    ) %>%
    dplyr::arrange(.data$stan_testlet)
  items <-
    dplyr::tibble(
      stan_item = dplyr::dense_rank(items),
      item = items,
      testlet = item_testlets,
      max_score = item_max_scores
    ) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(testlets, by = "testlet") %>%
    dplyr::arrange(.data$stan_item)
  ## Collate observations for Stan.
  group_observations <-
    dplyr::tibble(
      cohort = obs_group_cohorts,
      group = obs_groups,
      item = obs_group_items,
      observed_score = obs_group_scores,
      weight = group_obs_weights
    ) %>%
    dplyr::inner_join(
      dplyr::mutate(
        groups,
        ## The Stan model uses negative integers
        ##  to denote group observations.
        stan_person = -.data$stan_group,
        person = vec_cast(NA, to = obs_persons),
        .data$group
      ),
      by = "group"
    ) %>%
    dplyr::inner_join(items, by = "item")
  person_observations <-
    dplyr::tibble(
      cohort = obs_person_cohorts,
      person = obs_persons,
      item = obs_person_items,
      observed_score = obs_person_scores,
      weight = person_obs_weights
    ) %>%
    dplyr::inner_join(persons, by = "person") %>%
    dplyr::inner_join(items, by = "item")
  observations <-
    dplyr::bind_rows(group_observations, person_observations) %>%
    dplyr::select(
      .data$cohort,
      .data$stan_group, .data$group,
      .data$stan_person, .data$person,
      .data$stan_testlet, .data$testlet,
      .data$stan_item, .data$item,
      .data$max_score, .data$observed_score,
      .data$weight
    )
  ## Fit the model.
  stanfit <-
    rstan::sampling(
      stanmodels$raschdash,
      data =
        list(
          M = nrow(groups),
          N = nrow(persons),
          L = nrow(testlets),
          I = nrow(items),
          K = k,
          O = nrow(observations),
          mm = persons$stan_group,
          nn = observations$stan_person,
          ll = items$stan_testlet,
          ii = observations$stan_item,
          kk = items$max_score,
          y = observations$observed_score,
          w = observations$weight
        ),
      # Extra iterations would be necessary for reliable
      # 95% intervals (10K n_eff for η and δ).
      # iter = 10000,
      # The hierarchical models do not converge without
      # a high adapt_delta.
      control = list(adapt_delta = 0.99),
      pars = stan_pars,
      ...
    )
  structure(
    list(
      data = dplyr::as_tibble(observations),
      stanfit = stanfit
    ),
    class = "rdfit"
  )
}
