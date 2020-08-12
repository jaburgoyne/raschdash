#' @importFrom loo loo
#' @export
loo.rdfit <- function(x, ...) loo::loo(x$stanfit)
