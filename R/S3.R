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
#' The area under the curve (AUC) is simply the integral of the chosen model
#' `y(x; theta)` with respect to `x`.
#'
#' In real applications the response variable is usually contained within a
#' known interval. For example, if our response represents relative viability
#' against a control compound, the curve is then expected to be between 0 and 1.
#'
#' To make the AUC value comparable between different compounds and/or studies,
#' this function set a hard constraint on both the `x` variable and the function
#' `y`. The intervals can always be changed if needed.
#'
#' The function `f(x; theta)` that is integrated here is defined as
#'
#' `f(x; theta) = ylim[1], if y(x; theta) < ylim[1]`
#'
#' `f(x; theta) = y(x; theta), if ylim[1] <= y(x; theta) <= ylim[2]`
#'
#' `f(x; theta) = ylim[2], if y(x; theta) > ylim[2]`
#'
#' Finally, the AUC is normalized by its maximum possible value, that is the
#' area of the rectangle with width `xlim[2] - xlim[1]` and height
#' `ylim[2] - ylim[1]`.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param xlim numeric vector of length 2 with the lower and upped bound of the
#'  integration interval. Default is `c(-10, 10)`.
#' @param ylim numeric vector of length 2 with the lower and upped bound of the
#'  allowed function values. Default is `c(0, 1)`.
#'
#' @return Numeric value representing the normalized area under the curve.
#'
#' @export
nauc <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  UseMethod("nauc", object)
}

#' Area above the curve
#'
#' Evaluate the normalized area above the curve (NAAC).
#'
#' @details
#' The area under the curve (AUC) is simply the integral of the chosen model
#' `y(x; theta)` with respect to `x`.
#'
#' In real applications the response variable is usually contained within a
#' known interval. For example, if our response represents relative viability
#' against a control compound, the curve is then expected to be between 0 and 1.
#'
#' To make the AUC value comparable between different compounds and/or studies,
#' this function set a hard constraint on both the `x` variable and the function
#' `y`. The intervals can always be changed if needed.
#'
#' The function `f(x; theta)` that is integrated here is defined as
#'
#' `f(x; theta) = ylim[1], if y(x; theta) < ylim[1]`
#'
#' `f(x; theta) = y(x; theta), if ylim[1] <= y(x; theta) <= ylim[2]`
#'
#' `f(x; theta) = ylim[2], if y(x; theta) > ylim[2]`
#'
#' Finally, the AUC is normalized by its maximum possible value, that is the
#' area of the rectangle with width `xlim[2] - xlim[1]` and height
#' `ylim[2] - ylim[1]`.
#'
#' The normalized area above the curve is simply `NAAC = 1 - NAUC`.
#'
#' @param object fit object as returned by \code{\link[drda]{drda}}.
#' @param xlim numeric vector of length 2 with the lower and upped bound of the
#'  integration interval. Default is `c(-10, 10)`.
#' @param ylim numeric vector of length 2 with the lower and upped bound of the
#'  allowed function values. Default is `c(0, 1)`.
#'
#' @return Numeric value representing the normalized area above the curve.
#'
#' @export
naac <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  UseMethod("naac", object)
}
