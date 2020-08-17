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
        vctrs::vec_assert(value, ptype = double())
        vctrs::vec_assert(pareto_k, ptype = double())
        vctrs::new_rcrd(
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
        dplyr::na_if(format(vctrs::field(x, "value")), "NA")
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
vec_cast.double.eloo <- function(x, to, ...) vctrs::field(x, "value")

#' @importFrom vctrs vec_math
#' @export
vec_math.eloo <- function(.fn, .x, ...) {
        vctrs::vec_math_base(.fn, vctrs::vec_cast(.x, double()))
}

#' @importFrom vctrs vec_arith
#' @method vec_arith eloo
#' @export
vec_arith.eloo <- function(op, x, y, ...) {
        UseMethod("vec_arith.eloo", y)
}

#' @method vec_arith.eloo default
#' @export
vec_arith.eloo.default <- function(op, x, y, ...) {
        vctrs::stop_incompatible_op(op, x, y)
}

#' @method vec_arith.eloo eloo
#' @export
vec_arith.eloo.eloo <- function(op, x, y, ...) {
        vctrs::vec_math_base(
                op,
                vctrs::vec_cast(x, double()),
                vctrs::vec_cast(y, double())
        )
}

#' @method vec_arith.eloo numeric
#' @export
vec_arith.eloo.numeric <- function(op, x, y, ...) {
        switch(
                op,
                "+" = ,
                "-" = ,
                "*" = ,
                "/" = ,
                "%%" =
                        new_eloo(
                                value =
                                        vctrs::vec_arith(
                                                op = op,
                                                x = vctrs::field(x, "value"),
                                                y = y
                                        ),
                                pareto_k = vctrs::field(x, "pareto_k")
                        ),
                vctrs::vec_math_base(op, vctrs::vec_cast(x, double()), y)
        )
}

#' @method vec_arith.eloo MISSING
#' @export
vec_arith.eloo.MISSING <- function(op, x, y, ...) {
        switch(
                op,
                "+" = x,
                "-" = x * -1,
                vctrs::stop_incompatible_op(op, x, y)
        )
}

#' @export
vec_arith.numeric.eloo <- function(op, x, y, ...) {
        switch(
                op,
                "+" = ,
                "-" = ,
                "*" =
                        new_eloo(
                                value =
                                        vctrs::vec_arith(
                                                op = op,
                                                x = x,
                                                y = vctrs::field(y, "value")
                                        ),
                                pareto_k = vctrs::field(y, "pareto_k")
                        ),
                vctrs::vec_math_base(op, x, vctrs::vec_cast(y, double()))
        )
}
