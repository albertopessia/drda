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

find_optimum <- function(object) {
  UseMethod("find_optimum", object)
}

find_optimum_constrained <- function(object, constraint, known_param) {
  UseMethod("find_optimum_constrained", object)
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

nauc <- function(object, lower_bound, upper_bound) {
  UseMethod("nauc", object)
}

naac <- function(object, lower_bound, upper_bound) {
  UseMethod("naac", object)
}
