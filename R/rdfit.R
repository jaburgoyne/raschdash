#' Constructor for \code{rdfit}
#'
#' Construct an \code{rdfit} object.
#'
#' @param data a data frame with columns for cohort, group, person, testlet,
#'             item, maximum score, and observed score
#' @param cohort_colname name of the cohort column
#' @param group_colname name of the group column
#' @param person_colname name of the person column
#' @param testlet_colname name of the testlet column
#' @param item_colname name of the item column
#' @param max_score_colname name of the maximum score column
#' @param obs_score_colname name of the observed score column
#'
#' @return An \code{rdfit} object.
#'
#' @import dplyr
#' @importFrom rlang .data
#' @export
new_rdfit <- function(cohorts,
                      groups, person_groups, persons,
                      testlets, item_testlets, items, max_scores, K = 0,
                      obs_group_cohorts, obs_groups, obs_group_items,
                      obs_group_scores,
                      obs_person_cohorts, obs_persons, obs_person_items,
                      obs_person_scores,
                      ...) {
        ## Make tidy tibbles and sort for Stan..
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
                        max_score = max_scores
                ) %>%
                distinct() %>%
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
                                .data$group),
                        by = 'group'
                ) %>%
                inner_join(items, by = 'item')
        person_observations <-
                tibble(
                        cohort = obs_person_cohorts,
                        group = obs_persons,
                        item = obs_person_items,
                        obs_score = obs_person_scores,
                ) %>%
                inner_join(persons, by = 'person') %>%
                inner_join(items, by = 'item')
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
                                        K = K,
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
                        pars =
                                c(
                                        'epsilon_raw', 'upsilon_raw',
                                        'd_tau_raw', 'd_tau_err',
                                        'xi_raw', 'zeta_raw'
                                ),
                        include = FALSE,
                        ...
                )
        loo <- loo::loo(stanfit, save_psis = TRUE)
        structure(
                .Data = stanfit,
                loo = loo,
                class = c('rdfit', 'stanfit')
        )
}

#' @importFrom dplyr as_tibble distinct transmute
#' @importFrom rlang abort is_null
#' @importFrom vctrs vec_cast vec_duplicate_any vec_equal_na
rdfit <- function(data,
                  cohort = NULL,
                  group, person,
                  testlet, item, max_score, K = 0,
                  obs_score,
                  ...) {
        obs <-
                ## rdfit objects return fit statistics as tibbles. If the data
                ## cannot be represented as a tibble, then the tidverse errors
                ## (or warnings) are all that is necessary.
                as_tibble(data) %>%
                ## Skip column existence checks, because tibble will complain.
                transmute(
                        ## Fill in a generic All cohort if necessary.
                        cohort    =
                                ifelse(
                                        is_null(cohort),
                                        'All',
                                        {{cohort}}
                                ),
                        group     = {{group}},
                        person    = {{person}},
                        testlet   = {{testlet}},
                        item      = {{item}},
                        ## Rasch models require integral scores, and vctrs
                        ## complains nicely.
                        max_score =
                                vec_cast(
                                        x = {{max_score}},
                                        to = integer(),
                                        x_arg = 'max_score'
                                ),
                        obs_score =
                                vec_cast(
                                        x = {{obs_score}},
                                        to = integer(),
                                        x_arg = 'obs_score'
                                )
                ) %>%
                ## Rasch models have no problem with missing observations, but
                ## Stan does not want to see them.
                filter(!vec_equal_na(obs_score))
        ## Stan will complain about invalid scores, but the messages will be
        ## opaque for R users.
        if (any(obs$max_score < 1)) {
                abort('Some maximum scores are non-positive.')
        }
        if (any(obs$obs_score > obs$max_score)) {
                abort('Some observed scores are greater than their maxima.')
        }
        ## Check for NAs and duplicates when preparing data for Stan. Stan's
        ## validation will catch these, but not with helpful messages.
        cohorts <- obs %>% distinct(.data$cohort)
        if (any(vec_equal_na(cohorts$cohort))) {
                abort('Some cohorts are undefined.')
        }
        groups <- obs %>% distinct(.data$group)
        if (any(vec_equal_na(groups$group))) {
                abort('Some groups are undefined.')
        }
        persons <-
                obs %>%
                distinct(.data$person, .data$group) %>%
                filter(is.na(person))
        ## Persons may be NA because at this point in the function, the group
        ## will not be (i.e., NA persons are guaranteed to refer to group
        ## observations).
        if (vec_duplicate_any(persons$person)) {
                abort('Some persons belong to multiple groups.')
        }
        testlets <- obs %>% distinct(.data$testlet)
        if (any(vec_equal_na(testlets$testlet))) {
                abort('Some testlets are undefined.')
        }
        items <- obs %>% distinct(.data$item, .data$testlet, .data$max_score)
        if (any(vec_equal_na(items$item))
            || any(vec_equal_na(items$max_score))) {
                abort('Some items are undefined.')
        }
        if (vec_duplicate_any(items$item)) {
                abort(
                        stringr::str_c(
                                'Some items belong to multiple testlets ',
                                'or have inconsistent max scores.'
                                )
                )
        }
        ## A negative integer is no serious problem for K (it simply fails to
        ## select any items for rating scale), but Stan needs it to be non-
        ## negative in order to set the dimension of tau.
        K <- max(vec_cast(K, integer(), x_arg = 'K'), 0)
        ## The constructor asks for group and individual observations to be
        ## separated.
        group_obs <- obs %>% filter(vec_equal_na(.data$student))
        person_obs <- obs %>% filter(!vec_equal_na(.data$student))
        new_rdfit(
                cohorts = cohorts$cohort,
                groups = groups$group,
                person_groups = persons$group, persons = persons$person,
                testlets = testlets$testlet,
                item_testlets = items$testlet, items = items$item,
                max_scores = items$max_score, K = K,
                obs_group_cohorts = group_obs$cohort,
                obs_groups = group_obs$group,
                obs_group_items = group_obs$item,
                obs_group_scores = group_obs$obs_score,
                obs_person_cohorts = person_obs$cohort,
                obs_persons = person_obs$person,
                obs_person_items = person_obs$item,
                obs_person_scores = person_obs$score,
                ...
        )
}
