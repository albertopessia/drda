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
#' Evaluate the normalized area under the curve (AUC).
#'
#' @details
#' The area under the curve (AUC) is simply the integral of the chosen model
#' `f(x; theta)` between `lower_bound` and `upper_bound` with respect to `x`.
#'
#' When the interval of integration is fixed, the curve `f(x; theta)` is
#' contained into the rectangle of height `omega` and width
#' `upper_bound - lower_bound`. The maximum area the curve can have is obviously
#' `(upper_bound - lower_bound) * omega`.
#'
#' We first shift the curve to set its minimum to 0. We then integrate the curve
#' and define the normalized AUC (NAUC) by dividing its value by the maximum
#' area.
#'
#' Default values of `lower_bound` and `upper_bound` were chosen based on common
#' dose ranges used in the literature. They are also symmetric around zero
#' so that `NAUC` is equal to `0.5` in the standard logistic model.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param lower_bound numeric value with the lower bound of the integration
#'   interval.
#' @param upper_bound numeric value with the upper bound of the integration
#'   interval.
#'
#' @return Numeric value representing the normalized area under the curve.
#'
#' @export
nauc <- function(object, lower_bound = -10, upper_bound = 10) {
  UseMethod("nauc", object)
}

#' Area above the curve
#'
#' Evaluate the normalized area above the curve (AAC).
#'
#' @details
#' The area under the curve (AUC) is simply the integral of the chosen model
#' `f(x; theta)` between `lower_bound` and `upper_bound` with respect to `x`.
#'
#' When the interval of integration is fixed, the curve `f(x; theta)` is
#' contained into the rectangle of height `omega` and width
#' `upper_bound - lower_bound`. The maximum area the curve can have is obviously
#' `(upper_bound - lower_bound) * omega`.
#'
#' We first shift the curve to set its minimum to 0. We then integrate the curve
#' and define the normalized AUC (NAUC) by dividing its value by the maximum
#' area. As a consequence, the normalized area above the curve is simply
#' `NAAC = 1 - NAUC`.
#'
#' Default values of `lower_bound` and `upper_bound` were chosen based on common
#' dose ranges used in the literature. They are also symmetric around zero
#' so that `NAAC` is equal to `0.5` in the standard logistic model.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param lower_bound numeric value with the lower bound of the integration
#'   interval.
#' @param upper_bound numeric value with the upper bound of the integration
#'   interval.
#'
#' @return Numeric value with the requested area.
#'
#' @export
naac <- function(object, lower_bound = -10, upper_bound = 10) {
  UseMethod("naac", object)
}
