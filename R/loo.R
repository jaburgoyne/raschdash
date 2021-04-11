#' @importFrom loo loo
#' @export
loo.rdfit <- function(x, ...) {
  log_lik <- purrr::pluck(x, "draws") %>% posterior::subset_draws("log_lik")
  r_eff <- loo::relative_eff(exp(log_lik), cores = getOption("mc.cores", 1))
  loo::loo(
    log_lik,
    r_eff = r_eff,
    cores = getOption("mc.cores", 1)
  )
}
