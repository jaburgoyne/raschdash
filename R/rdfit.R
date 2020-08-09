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
#'             item, maximum possible item score, and observed score. For
#'             group-level observations, the person should be `NA`.
#' @param cohort <[tidyr::tidyr_tidy_select]> cohort column. Defaults to a
#'               generic *All* cohort.
#' @param group <[tidyr::tidyr_tidy_select]> group column
#' @param person <[tidyr::tidyr_tidy_select]> person column
#' @param testlet <[tidyr::tidyr_tidy_select]> testlet column
#' @param item <[tidyr::tidyr_tidy_select]> item column
#' @param max_score <[tidyr::tidyr_tidy_select]> maximum-score column
#' @param obs_score <[tidyr::tidyr_tidy_select]> observed score column
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
#' @param obs_person_cohorts,obs_persons,obs_person_items,obs_person_scores
#'            for the person-level observations, vectors of the corresponding
#'            cohorts (`obs_person_cohorts`), persons (`obs_persons`), items
#'            (`obs_person_items`), and observed scores (`obs_person_scores`)
#' @param ... additional parameters passed to [rstan::sampling()]
#'
#' @return An object of class `rdfit` with the following attributes.
#' \describe{
#'  \item{`stanfit`}{
#'   The [rstan::stanfit] object for the fitted model. It contains samples of
#'   the following parameters.
#'   * logit-scale parameters
#'       * `xi`: group abilities
#'       * `eta`: person abilities
#'       * `epsilon`: testlet difficulties
#'       * `delta`: item difficulties
#'       * `tau`: rating-scale threshold offsets from overall item
#'                difficulty (if any rating-scale items are present)
#'       * `lambda`: standard deviation of person abilities. Used to
#'                   transform the logit scale to the standard scale.
#'   * standard-scale parameters
#'       * `group_ability`: group abilities
#'       * `person_ability`: person abilities
#'       * `testlet_difficulty`: testlet difficulties
#'       * `item_difficulty`: item difficulties
#'       * `thresholds`: rating-scale threshold offsets from overall item
#'                       difficulty (if any rating-scale items are
#'                       present)
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
#'  }
#'
#'  \item{`loo`}{
#'   A [loo::loo] object for `stanfit`.
#'  }
#' }
#' @name rdfit-class
NULL

#' @describeIn rdfit-class Create an `rdfit` object from dimensional data.
#' @importFrom dplyr as_tibble distinct filter if_else transmute
#' @importFrom rlang abort is_null
#' @importFrom vctrs vec_cast vec_duplicate_any vec_equal_na
#' @export
rdfit <- function(data,
                  cohort = "All",
                  group, person,
                  testlet, item, max_score,
                  obs_score,
                  k = 0,
                  ...) {
        obs <-
                ## rdfit objects return fit statistics as tibbles. If the data
                ## cannot be represented as a tibble, then the tidverse errors
                ## (or warnings) are all that is necessary.
                as_tibble(data) %>%
                ## Skip column existence checks, because tibble will complain.
                transmute(
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
                                )
                ) %>%
                ## Rasch models have no problem with missing observations, but
                ## Stan does not want to see them.
                filter(!vec_equal_na(obs_score))
        ## Stan will complain about invalid scores, but the messages will be
        ## opaque for R users.
        if (any(obs$max_score < 1)) {
                abort("Some maximum scores are non-positive.")
        }
        if (any(obs$obs_score > obs$max_score)) {
                abort("Some observed scores are greater than their maxima.")
        }
        ## Check for NAs and duplicates when preparing data for Stan. Stan's
        ## validation will catch these, but not with helpful messages.
        cohorts <- obs %>% distinct(.data$cohort)
        if (any(vec_equal_na(cohorts$cohort))) {
                abort("Some cohorts are undefined.")
        }
        groups <- obs %>% distinct(.data$group)
        if (any(vec_equal_na(groups$group))) {
                abort("Some groups are undefined.")
        }
        persons <-
                obs %>%
                distinct(.data$person, .data$group) %>%
                filter(!vec_equal_na(person))
        ## Persons may be NA because at this point in the function, the group
        ## will not be (i.e., NA persons are guaranteed to refer to group
        ## observations).
        if (vec_duplicate_any(persons$person)) {
                abort("Some persons belong to multiple groups.")
        }
        testlets <- obs %>% distinct(.data$testlet)
        if (any(vec_equal_na(testlets$testlet))) {
                abort("Some testlets are undefined.")
        }
        items <- obs %>% distinct(.data$item, .data$testlet, .data$max_score)
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
        k <- max(vec_cast(k, integer(), x_arg = "k"), 0)
        ## The constructor asks for group and individual observations to be
        ## separated.
        group_obs <- obs %>% filter(vec_equal_na(.data$person))
        person_obs <- obs %>% filter(!vec_equal_na(.data$person))
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
                obs_person_cohorts = person_obs$cohort,
                obs_persons = person_obs$person,
                obs_person_items = person_obs$item,
                obs_person_scores = person_obs$obs_score,
                k = k,
                ...
        )
}

#' @describeIn rdfit-class Create an `rdfit` object from normalised data.
#' @importFrom dplyr arrange bind_rows dense_rank distinct inner_join mutate
#' @importFrom rlang .data
#' @importFrom tibble new_tibble tibble
#' @export
new_rdfit <- function(cohorts,
                      groups, persons, person_groups,
                      testlets, items, item_testlets, item_max_scores,
                      obs_group_cohorts, obs_groups, obs_group_items,
                      obs_group_scores,
                      obs_person_cohorts, obs_persons, obs_person_items,
                      obs_person_scores,
                      k = 0,
                      ...) {
        ## Make tidy tibbles and sort for Stan.
        groups <-
                tibble(
                        stan_group = dense_rank(groups),
                        group = groups
                ) %>%
                distinct() %>%
                arrange(.data$stan_group)
        persons <-
                tibble(
                        stan_person = dense_rank(persons),
                        person = persons,
                        group = person_groups
                ) %>%
                distinct() %>%
                inner_join(groups, by = "group") %>%
                arrange(.data$stan_person)
        testlets <-
                tibble(
                        stan_testlet = dense_rank(testlets),
                        testlet = testlets
                ) %>%
                arrange(.data$stan_testlet)
        items <-
                tibble(
                        stan_item = dense_rank(items),
                        item = items,
                        testlet = item_testlets,
                        max_score = item_max_scores
                ) %>%
                distinct() %>%
                inner_join(testlets, by = "testlet") %>%
                arrange(.data$stan_item)
        # Collate observations for Stan.
        group_observations <-
                tibble(
                        cohort = obs_group_cohorts,
                        group = obs_groups,
                        item = obs_group_items,
                        obs_score = obs_group_scores,
                ) %>%
                inner_join(
                        mutate(
                                groups,
                                ## The Stan model uses negative integers
                                ##  to denote group observations.
                                stan_person = -.data$stan_group,
                                person = NA,
                                .data$group
                        ),
                        by = "group"
                ) %>%
                inner_join(items, by = "item")
        person_observations <-
                tibble(
                        cohort = obs_person_cohorts,
                        person = obs_persons,
                        item = obs_person_items,
                        obs_score = obs_person_scores,
                ) %>%
                inner_join(persons, by = "person") %>%
                inner_join(items, by = "item")
        observations <- bind_rows(group_observations, person_observations)
        # Fit the model.
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
                                        y = observations$obs_score
                                ),
                        # Extra iterations would be necessary for reliable
                        # 95% intervals (10K n_eff for η and δ).
                        # iter = 10000,
                        # The hierarchical model does not converge without
                        # a high adapt_delta.
                        control = list(adapt_delta = 0.99),
                        # The raw parameters are uninteresting.
                        pars = c(
                                "xi", "eta", "epsilon", "delta", "tau",
                                "lambda",
                                "group_ability", "person_ability",
                                "testlet_difficulty", "item_difficulty",
                                "thresholds",
                                "y_rep",
                                "prior_group_ability",
                                "prior_person_ability",
                                "prior_testlet_difficulty",
                                "prior_item_difficulty",
                                "log_lik", "log_lik_rep"
                        ),
                        ...
                )
        loo <- loo::loo(stanfit, save_psis = TRUE)
        ## TODO: Replace observations with posterior statistics.
        new_tibble(
                observations,
                stanfit = stanfit,
                loo = loo,
                class = "rdfit"
        )
}
