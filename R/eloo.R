#' `eloo` vector
#'
#' Attach Pareto-*k* estimates to expectations from the `loo` package. Mostly
#' intended for internal development.
#'
#' The `eloo` class makes it easier to work with vectors that have Pareto *k*
#' values attached. For a conservative set of basic arithmetic operations,
#' `eloo` objects will retain their Pareto *k* values, but most operations will
#' fall gracefully back to standard numeric vector operations.
#'
#' @param value a numeric vector
#' @param pareto_k a numeric vector of Pareto-*k* shape parameters corresponding
#'                 to `x`
#' @param df a data frame with columns `value` and `pareto_k`
#'
#' @return An S3 vector of class `eloo`. The Pareto *k* parameters are
#'         accessible using [vctrs::field()].
#'
#' @name eloo-class
NULL

#' @describeIn eloo-class Create an `eloo` vector.
#' @export
new_eloo <- function(value = double(), pareto_k = double()) {
        vec_assert(value, ptype = double())
        vec_assert(pareto_k, ptype = double())
        new_rcrd(
                fields = list(value = value, pareto_k = pareto_k),
                class = "eloo"
        )
}

#' @describeIn eloo-class Create an `eloo` vector from the output of
#'                        [loo::E_loo()].
#' @export
eloo <- function(df) new_eloo(df$value, df$pareto_k)

#' @export
format.eloo <- function(x, ...) {
        dplyr::na_if(format(field(x, "value")), "NA")
}

#' @importFrom vctrs vec_ptype2
#' @export
vec_ptype2.eloo.eloo <- function(x, y, ...) new_eloo()

#' @export
vec_ptype2.double.eloo <- function(x, y, ...) double()

#' @export
vec_ptype2.eloo.double <- function(x, y, ...) double()

#' @importFrom vctrs vec_cast
#' @export
vec_cast.eloo.eloo <- function(x, to, ...) x

#' @export
vec_cast.double.eloo <- function(x, to, ...) field(x, "value")

#' @importFrom vctrs vec_math
#' @export
vec_math.eloo <- function(.fn, .x, ...) {
        vec_math_base(.fn, vec_cast(.x, double()))
}

# vec_arith has not yet received the update necessary to avoid the boilerplate
# functions and namespace tweaks (https://github.com/r-lib/vctrs/issues/1063).

#' @importFrom vctrs vec_arith
#' @method vec_arith eloo
#' @export
vec_arith.eloo <- function(op, x, y, ...) UseMethod("vec_arith.eloo", y)

#' @method vec_arith.eloo default
#' @export
vec_arith.eloo.default <- function(op, x, y, ...) stop_incompatible_op(op, x, y)

#' @method vec_arith.eloo eloo
#' @export
vec_arith.eloo.eloo <- function(op, x, y, ...) {
        vec_arith_base(op, vec_cast(x, double()), vec_cast(y, double()))
}

#' @method vec_arith.eloo numeric
#' @export
vec_arith.eloo.numeric <- function(op, x, y, ...) {
        if (vec_in(op, c("+", "-", "*", "/", "%%"))) {
                new_eloo(
                        value = vec_arith(op, field(x, "value"), y),
                        pareto_k = field(x, "pareto_k")
                )
        } else {
                vec_arith_base(op, vec_cast(x, double()), y)
        }
}

#' @method vec_arith.eloo MISSING
#' @export
vec_arith.eloo.MISSING <- function(op, x, y, ...) {
        if (op == "-") {
                -1 * x
        } else if (op == "+") {
                x
        } else {
                stop_incompatible_op(op, x, y)
        }
}

#' @export
vec_arith.numeric.eloo <- function(op, x, y, ...) {
        if (vec_in(op, c("+", "-", "*"))) {
                new_eloo(
                        value = vec_arith(op, x, field(y, "value")),
                        pareto_k = field(y, "pareto_k")
                )
        } else {
                vec_arith_base(op, x, vec_cast(y, double()))
        }
}
