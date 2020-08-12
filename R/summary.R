.heaviside_difference <- function(observed, expected) {
        0.5 * sign(observed - expected) + 0.5
}

.loo_statistics <- function(x) {
        stanfit <- purrr::pluck(x, "stanfit")
        loo <- loo::loo(stanfit, save_psis = TRUE)
        .loo_expectation <- function(x) {
                loo::E_loo(
                        x = x,
                        psis_object = loo$psis_object,
                        ## Although the Pareto-k diagnostics are unused for now,
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
                        .heaviside_difference(as.matrix(stanfit, "y_rep")) %>%
                        .loo_expectation(),
                ## Subtract log_lik from log_lik_rep instead of the other way
                ## around because this statistic is about information content,
                ## not likelihood.
                p_information =
                        as.matrix(stanfit, "log_lik_rep") %>%
                        .heaviside_difference(as.matrix(stanfit, "log_lik")) %>%
                        .loo_expectation(),
                pareto_k = loo::pareto_k_values(loo)
        )
}

#' @export
summary.rdfit <- function(object, ...) {
        dplyr::bind_cols(
                purrr::pluck(object, "data"),
                .loo_statistics(object)
        ) %>%
                dplyr::as_tibble() %>%
                dplyr::select(!dplyr::starts_with("stan_"))
}
