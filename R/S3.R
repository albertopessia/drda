fn <- function(object, x, theta) {
  UseMethod("fn", object)
}

gradient_hessian <- function(object, theta) {
  UseMethod("gradient_hessian", object)
}

rss <- function(object) {
  UseMethod("rss", object)
}

rss_fixed <- function(object, known_param) {
  UseMethod("rss_fixed", object)
}

rss_gradient_hessian <- function(object) {
  UseMethod("rss_gradient_hessian", object)
}

rss_gradient_hessian_fixed <- function(object, known_param) {
  UseMethod("rss_gradient_hessian_fixed", object)
}

mle_asy <- function(object, theta) {
  UseMethod("mle_asy", object)
}

init <- function(object) {
  UseMethod("init", object)
}

fit <- function(object) {
  UseMethod("fit", object)
}

fit_constrained <- function(object) {
  UseMethod("fit_constrained", object)
}

fisher_info <- function(object, theta, sigma) {
  UseMethod("fisher_info", object)
}

curve_variance <- function(object, x) {
  UseMethod("curve_variance", object)
}

#' Area under the curve
#'
#' Evaluate the normalized area under the curve (NAUC).
#'
#' @details
#' The area under the curve (AUC) is the integral of the chosen model
#' `y(x; theta)` with respect to `x`.
#'
#' In real applications the response variable is usually contained within a
#' known interval. For example, if our response represents relative viability
#' against a control compound, the curve is expected to be between 0 and 1.
#' Let `ylim = c(yl, yu)` represent the admissible range of our function
#' `y(x; theta)`, that is `yl` is its lower bound and `yu` its upper bound.
#' Let `xlim = c(xl, xu)` represent the admissible range of the predictor
#' variable `x`. For example, when `x` represent the dose, the boundaries
#' are the minimum and maximum doses we can administer.
#'
#' To make the AUC value comparable between different compounds and/or studies,
#' this function sets a hard constraint on both the `x` variable and the
#' function `y`. The intervals can always be changed if needed.
#'
#' The integral calculated by this function is of the piece-wise function
#' `f(x; theta)` defined as
#'
#' `f(x; theta) = yl`, if `y(x; theta) < yl`
#'
#' `f(x; theta) = y(x; theta)`, if `yl <= y(x; theta) <= yu`
#'
#' `f(x; theta) = yu`, if `y(x; theta) > yu`
#'
#' The AUC is finally normalized by its maximum possible value, that is the
#' area of the rectangle with width `xu - xl` and height `yu - yl`.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param xlim numeric vector of length 2 with the lower and upped bound of the
#'  integration interval. Default is `c(-10, 10)` for the logistic function or
#'  `c(0, 1000)` for the log-logistic function.
#' @param ylim numeric vector of length 2 with the lower and upped bound of the
#'  allowed function values. Default is `c(0, 1)`.
#'
#' @return Numeric value representing the normalized area under the curve.
#'
#' @export
nauc <- function(object, xlim, ylim) {
  UseMethod("nauc", object)
}

#' Area above the curve
#'
#' Evaluate the normalized area above the curve (NAAC).
#'
#' @details
#' The area under the curve (AUC) is the integral of the chosen model
#' `y(x; theta)` with respect to `x`.
#'
#' In real applications the response variable is usually contained within a
#' known interval. For example, if our response represents relative viability
#' against a control compound, the curve is expected to be between 0 and 1.
#' Let `ylim = c(yl, yu)` represent the admissible range of our function
#' `y(x; theta)`, that is `yl` is its lower bound and `yu` its upper bound.
#' Let `xlim = c(xl, xu)` represent the admissible range of the predictor
#' variable `x`. For example, when `x` represent the dose, the boundaries
#' are the minimum and maximum doses we can administer.
#'
#' To make the AUC value comparable between different compounds and/or studies,
#' this function sets a hard constraint on both the `x` variable and the
#' function `y`. The intervals can always be changed if needed.
#'
#' The integral calculated by this function is of the piece-wise function
#' `f(x; theta)` defined as
#'
#' `f(x; theta) = yl`, if `y(x; theta) < yl`
#'
#' `f(x; theta) = y(x; theta)`, if `yl <= y(x; theta) <= yu`
#'
#' `f(x; theta) = yu`, if `y(x; theta) > yu`
#'
#' The AUC is finally normalized by its maximum possible value, that is the
#' area of the rectangle with width `xu - xl` and height `yu - yl`.
#'
#' The normalized area above the curve (NAAC) is simply `NAAC = 1 - NAUC`.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param xlim numeric vector of length 2 with the lower and upped bound of the
#'  integration interval. Default is `c(-10, 10)` for the logistic function or
#'  `c(0, 1000)` for the log-logistic function.
#' @param ylim numeric vector of length 2 with the lower and upped bound of the
#'  allowed function values. Default is `c(0, 1)`.
#'
#' @return Numeric value representing the normalized area above the curve.
#'
#' @export
naac <- function(object, xlim, ylim) {
  UseMethod("naac", object)
}
