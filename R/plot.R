#' @export
plot.rdfit <- function(x, ..., pars = c("person_ability")) {
        rstan::plot(x$stanfit, ..., pars = pars)
}
