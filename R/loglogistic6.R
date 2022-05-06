# Fit dose-response data
#
# Use a Newton trust-region method to fit a log-logistic function to
# dose-response data.
#
# @param x numeric vector representing the fixed predictor variable.
# @param y numeric vector of observed values.
# @param w an optional vector of weights to be used in the fitting
#   process.
# @param start starting values for the parameters.
# @param max_iter maximum number of iterations in the optimization algorithm.
# @param lower_bound numeric vector with the minimum admissible values of the
#   parameters.
# @param upper_bound numeric vector with the maximum admissible values of the
#   parameters.
#
# @return An object of class `loglogistic*`.
loglogistic6_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 6) {
      stop("'start' must be of length 6", call. = FALSE)
    }

    if (start[3] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    if (start[4] <= 0) {
      stop("parameter 'phi' cannot be negative nor zero", call. = FALSE)
    }

    if (start[5] <= 0) {
      stop("parameter 'nu' cannot be negative nor zero", call. = FALSE)
    }

    if (start[6] <= 0) {
      stop("parameter 'xi' cannot be negative nor zero", call. = FALSE)
    }

    start[3:6] <- log(start[3:6])
  }

  object <- structure(
    list(
      x = x,
      y = y,
      w = w,
      n = length(y),
      stats = suff_stats(x, y, w),
      constrained = FALSE,
      start = start,
      max_iter = max_iter
    ),
    class = "loglogistic6"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      lower_bound <- rep(-Inf, 6)
    } else {
      if (length(lower_bound) != 6) {
        stop("'lower_bound' must be of length 6", call. = FALSE)
      }

      lower_bound[3] <- if (lower_bound[3] > 0) {
        log(lower_bound[3])
      } else {
        -Inf
      }

      lower_bound[4] <- if (lower_bound[4] > 0) {
        log(lower_bound[4])
      } else {
        -Inf
      }

      lower_bound[5] <- if (lower_bound[5] > 0) {
        log(lower_bound[5])
      } else {
        -Inf
      }

      lower_bound[6] <- if (lower_bound[6] > 0) {
        log(lower_bound[6])
      } else {
        -Inf
      }
    }

    if (is.null(upper_bound)) {
      upper_bound <- rep(Inf, 6)
    } else {
      if (length(upper_bound) != 6) {
        stop("'upper_bound' must be of length 6", call. = FALSE)
      }

      if (upper_bound[3] <= 0) {
        stop("'upper_bound[3]' cannot be negative nor zero.", call. = FALSE)
      }

      if (upper_bound[4] <= 0) {
        stop("'upper_bound[4]' cannot be negative nor zero.", call. = FALSE)
      }

      if (upper_bound[5] <= 0) {
        stop("'upper_bound[5]' cannot be negative nor zero.", call. = FALSE)
      }

      if (upper_bound[6] <= 0) {
        stop("'upper_bound[6]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[3:6] <- log(upper_bound[3:6])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 6-parameter log-logistic function
#'
#' Evaluate at a particular set of parameters the 6-parameter log-logistic
#' function.
#'
#' @details
#' The 6-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (xi * x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
#' `phi > 0`, `nu > 0`, and `xi > 0`. When `delta` is positive (negative) the
#' curve is monotonically increasing (decreasing).
#'
#' Parameter `alpha` is the value of the function when `x = 0`.
#' Parameter `delta` affects the value of the function when `x -> Inf`.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the mid-value of the function.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#' Parameter `xi` affects the value of the function when `x -> Inf`.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi, nu, xi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   log-logistic function.
#'
#' @export
loglogistic6_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  t1 <- x^eta
  t2 <- phi^eta

  alpha + delta * (t1 / (xi * t1 + nu * t2))^(1 / nu)
}

# @rdname loglogistic6_fn
fn.loglogistic6 <- function(object, x, theta) {
  loglogistic6_fn(x, theta)
}

# @rdname loglogistic6_fn
fn.loglogistic6_fit <- function(object, x, theta) {
  loglogistic6_fn(x, theta)
}

#' 6-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 6-parameter log-logistic function.
#'
#' @details
#' The 6-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (xi * x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
#' `phi > 0`, `nu > 0`, and `xi > 0`. When `delta` is positive (negative) the
#' curve is monotonically increasing (decreasing).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi, nu, xi)`.
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
loglogistic6_gradient <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xi * xe + nu * pe
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
  d <- delta * h / f

  a <- pe * lr
  q <- (eta * log(x) - log(f)) * f / nu
  s <- pe + q

  G <- matrix(1, nrow = k, ncol = 6)

  G[, 2] <- h
  G[, 3] <- a * d
  G[, 4] <- -eta * pe * d / phi
  G[, 5] <- -s * d / nu
  G[, 6] <- -xe * d / nu

  # gradient might not be defined when we plug x = 0 directly into the formula
  # however, the limits for x -> 0 are zero (not w.r.t. alpha)
  G[x_zero, -1] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  G
}

#' @rdname loglogistic6_gradient
loglogistic6_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xi * xe + nu * pe
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
  d <- delta * h / f

  a <- pe * lr
  p <- pe - xi * xe
  q <- (eta * log(x) - log(f)) * f / nu
  r <- d / f
  s <- pe + q

  H <- array(0, dim = c(k, 6, 6))

  H[, 3, 2] <- a * d / delta
  H[, 4, 2] <- -eta * pe * d / (delta * phi)
  H[, 5, 2] <- -s * d / (delta * nu)
  H[, 6, 2] <- -xe * d / (delta * nu)

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- lr * a * p * r
  H[, 4, 3] <- -(pe * f + eta * a * p) * r / phi
  H[, 5, 3] <- -(nu * pe + s) * a * r / nu
  H[, 6, 3] <- -(nu + 1) * xe * a * r / nu

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- eta * pe * (f + eta * p) * r / phi^2
  H[, 5, 4] <- eta * pe * (nu * pe + s) * r / (nu * phi)
  H[, 6, 4] <- (nu + 1) * eta * xe * pe * r / (nu * phi)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- (nu * pe^2 + (2 * f + s) * s) * r / nu^2
  H[, 6, 5] <- xe * (nu * pe + f + s) * r / nu^2

  H[, 2, 6] <- H[, 6, 2]
  H[, 3, 6] <- H[, 6, 3]
  H[, 4, 6] <- H[, 6, 4]
  H[, 5, 6] <- H[, 6, 5]
  H[, 6, 6] <- (nu + 1) * (xe / nu)^2 * r

  # Hessian might not be defined when we plug x = 0 directly into the formula
  # however, the limits for x -> 0 are zero
  H[x_zero, , ] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  H
}

#' @rdname loglogistic6_gradient
loglogistic6_gradient_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xi * xe + nu * pe
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
  d <- delta * h / f

  a <- pe * lr
  p <- pe - xi * xe
  q <- (eta * log(x) - log(f)) * f / nu
  r <- d / f
  s <- pe + q

  G <- matrix(1, nrow = k, ncol = 6)

  G[, 2] <- h
  G[, 3] <- a * d
  G[, 4] <- -eta * pe * d / phi
  G[, 5] <- -s * d / nu
  G[, 6] <- -xe * d / nu

  H <- array(0, dim = c(k, 6, 6))

  H[, 3, 2] <- a * d / delta
  H[, 4, 2] <- -eta * pe * d / (delta * phi)
  H[, 5, 2] <- -s * d / (delta * nu)
  H[, 6, 2] <- -xe * d / (delta * nu)

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- lr * a * p * r
  H[, 4, 3] <- -(pe * f + eta * a * p) * r / phi
  H[, 5, 3] <- -(nu * pe + s) * a * r / nu
  H[, 6, 3] <- -(nu + 1) * xe * a * r / nu

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- eta * pe * (f + eta * p) * r / phi^2
  H[, 5, 4] <- eta * pe * (nu * pe + s) * r / (nu * phi)
  H[, 6, 4] <- (nu + 1) * eta * xe * pe * r / (nu * phi)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- (nu * pe^2 + (2 * f + s) * s) * r / nu^2
  H[, 6, 5] <- xe * (nu * pe + f + s) * r / nu^2

  H[, 2, 6] <- H[, 6, 2]
  H[, 3, 6] <- H[, 6, 3]
  H[, 4, 6] <- H[, 6, 4]
  H[, 5, 6] <- H[, 6, 5]
  H[, 6, 6] <- (nu + 1) * (xe / nu)^2 * r

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero (not w.r.t. alpha)
  G[x_zero, -1] <- 0
  H[x_zero, , ] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  list(G = G, H = H)
}

#' 6-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 6-parameter log-logistic function.
#'
#' @details
#' The 6-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (xi * x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
#' `phi > 0`, `nu > 0`, and `xi > 0`.
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{loglogistic6_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)`, `phi2 = log(phi)`, `nu2 = log(nu)`, and
#' `xi2 = log(xi)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi, nu, xi)`.
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
loglogistic6_gradient_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  k1 <- eta / nu
  k2 <- xi / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- xi * c1 + nu * c2
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  c <- k2 * c1
  d <- g / f

  e <- log(x) - log(theta[4])

  p <- a * g
  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
  u <- q * s
  v <- q * t
  w <- q * c

  G <- matrix(1, nrow = k, ncol = 6)

  G[, 2] <- p
  G[, 3] <- delta * e * r
  G[, 4] <- -delta * r
  G[, 5] <- -delta * (v - u)
  G[, 6] <- -delta * w

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero (not w.r.t. alpha)
  G[x_zero, -1] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  G
}

#' @rdname loglogistic6_gradient_2
loglogistic6_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  k1 <- eta / nu
  k2 <- xi / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- xi * c1 + nu * c2
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  c <- k2 * c1
  d <- g / f

  e <- log(x) - log(theta[4])

  l <- (1 + nu) * c2 / f
  m <- (1 + nu) * xi * c1 / (nu * f)

  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
  u <- q * s
  v <- q * t
  w <- q * c
  y <- eta * log(x) - log(f)

  H <- array(0, dim = c(k, 6, 6))

  H[, 3, 2] <- e * r
  H[, 4, 2] <- -r
  H[, 5, 2] <- u - v
  H[, 6, 2] <- -w

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * e * (1 + eta * (l - 1) * e) * r
  H[, 4, 3] <- -delta * (1 + eta * (l - 1) * e) * r
  H[, 5, 3] <- -delta * e * (l + y / nu) * r
  H[, 6, 3] <- -delta * eta * e * l * w

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (l - 1) * r
  H[, 5, 4] <- delta * (l + y / nu) * r
  H[, 6, 4] <- delta * eta * l * w

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * ((l + y / nu) * c2 +
    (1 + y / nu) * (k1 * log(x) * f - s)) * q
  H[, 6, 5] <- delta * (1 + l + y / nu) * w

  H[, 2, 6] <- H[, 6, 2]
  H[, 3, 6] <- H[, 6, 3]
  H[, 4, 6] <- H[, 6, 4]
  H[, 5, 6] <- H[, 6, 5]
  H[, 6, 6] <- -delta * (1 - m) * w

  # Hessian might not be defined when we plug x = 0 directly into the formula
  # however, the limits for x -> 0 are zero
  H[x_zero, , ] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  H
}

#' @rdname loglogistic6_gradient_2
loglogistic6_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  k1 <- eta / nu
  k2 <- xi / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- xi * c1 + nu * c2
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  c <- k2 * c1
  d <- g / f

  e <- log(x) - log(theta[4])

  l <- (1 + nu) * c2 / f
  m <- (1 + nu) * xi * c1 / (nu * f)

  p <- a * g
  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
  u <- q * s
  v <- q * t
  w <- q * c
  y <- eta * log(x) - log(f)

  G <- matrix(1, nrow = k, ncol = 6)

  G[, 2] <- p
  G[, 3] <- delta * e * r
  G[, 4] <- -delta * r
  G[, 5] <- -delta * (v - u)
  G[, 6] <- -delta * w

  H <- array(0, dim = c(k, 6, 6))

  H[, 3, 2] <- e * r
  H[, 4, 2] <- -r
  H[, 5, 2] <- u - v
  H[, 6, 2] <- -w

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * e * (1 + eta * (l - 1) * e) * r
  H[, 4, 3] <- -delta * (1 + eta * (l - 1) * e) * r
  H[, 5, 3] <- -delta * e * (l + y / nu) * r
  H[, 6, 3] <- -delta * eta * e * l * w

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (l - 1) * r
  H[, 5, 4] <- delta * (l + y / nu) * r
  H[, 6, 4] <- delta * eta * l * w

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * ((l + y / nu) * c2 +
    (1 + y / nu) * (k1 * log(x) * f - s)) * q
  H[, 6, 5] <- delta * (1 + l + y / nu) * w

  H[, 2, 6] <- H[, 6, 2]
  H[, 3, 6] <- H[, 6, 3]
  H[, 4, 6] <- H[, 6, 4]
  H[, 5, 6] <- H[, 6, 5]
  H[, 6, 6] <- -delta * (1 - m) * w

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero (not w.r.t. alpha)
  G[x_zero, -1] <- 0
  H[x_zero, , ] <- 0

  # any NaN is because of corner cases where the derivatives are zero
  is_nan <- is.nan(G)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the gradient at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    G[is_nan] <- 0
  }

  is_nan <- is.nan(H)
  if (any(is_nan)) {
    warning(
      paste0(
        "issues while computing the Hessian at c(",
        paste(theta, collapse = ", "),
        ")"
      )
    )
    H[is_nan] <- 0
  }

  list(G = G, H = H)
}

# 6-parameter log-logistic function gradient and Hessian
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 6-parameter log-logistic function.
#
# @details
# The 6-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (xi * x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
# `phi > 0`, `nu > 0`, and `xi > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, `log(xi)`, and `log(nu)`.
#
# @param object object of class `loglogistic6`.
# @param theta numeric vector with the six parameters in the form
#   `c(alpha, delta, log(eta), log(phi), log(nu), log(xi))`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.loglogistic6 <- function(object, theta) {
  loglogistic6_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 6-parameter log-logistic model.
#
# @details
# The 6-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (xi * x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
# `phi > 0`, `nu > 0`, and `xi > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, `log(xi)`, and `log(nu)`.
#
# @param object object of class `loglogistic6`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.loglogistic6 <- function(object) {
  function(theta) {
    theta[3:6] <- exp(theta[3:6])
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.loglogistic6
rss_fixed.loglogistic6 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:6] <- exp(theta[3:6])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 6-parameter log-logistic model.
#
# @details
# The 6-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (xi * x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
# `phi > 0`, `nu > 0`, and `xi > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, `log(xi)`, and `log(nu)`.
#
# @param object object of class `loglogistic6`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.loglogistic6 <- function(object) {
  function(theta) {
    theta[3:6] <- exp(theta[3:6])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 6, 6))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- object$stats[, 2] * (r * H[, , 6] + G[, 6] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.loglogistic6
rss_gradient_hessian_fixed.loglogistic6 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:6] <- exp(theta[3:6])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 6, 6))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- object$stats[, 2] * (r * H[, , 6] + G[, 6] * G)

    list(
      G = apply(gradient[, idx, drop = FALSE], 2, sum),
      H = apply(hessian[, idx, idx, drop = FALSE], 2:3, sum)
    )
  }
}

# Maximum likelihood estimators
#
# Given a set of parameters, compute the maximum likelihood estimates of the
# linear parameters.
#
# @param object object of class `loglogistic6`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two parameters `alpha`
#   and `delta`.
mle_asy.loglogistic6 <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- exp(theta[4])
  nu <- exp(theta[5])
  xi <- exp(theta[6])

  g <- (x^eta / (xi * x^eta + nu * phi^eta))^(1 / nu)

  t1 <- 0
  t2 <- 0
  t3 <- 0
  t4 <- 0
  t5 <- 0

  for (i in seq_along(x)) {
    t1 <- t1 + w[i]
    t2 <- t2 + w[i] * g[i]
    t3 <- t3 + w[i] * g[i]^2
    t4 <- t4 + w[i] * y[i]
    t5 <- t5 + w[i] * g[i] * y[i]
  }

  denom <- t2^2 - t1 * t3

  if (denom != 0) {
    theta[1] <- (t2 * t5 - t3 * t4) / denom
    theta[2] <- (t2 * t4 - t1 * t5) / denom
  }

  theta
}

# Initialize vector of parameters
#
# Given the sufficient statistics, try to guess a good approximation to the
# Maximum Likelihood estimator of the six parameters of the log-logistic
# function.
#
# @param object object of class `loglogistic6`.
#
# @return Numeric vector of length 6 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.loglogistic6 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  theta <- if (is.null(object$start)) {
    # we initialize `nu = 1` and `xi = 1`, so that we start with a 4-parameter
    # log-logistic function
    #
    # define `beta = delta - alpha`, that is `beta` is the upper bound when
    # `alpha` is the lower bound.
    #
    # y = a + (b - a) * x^e / (x^e + p^e)
    # w = (y - a) / (b - a) = x^e / (x^e + p^e)
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(w / (1 - w)) = - e * log(p) + e * log(d)
    #
    # fit a linear model `z ~ u0 + u1 log(x)` and set `log_eta = log(u1)` and
    # `log_phi = -u0 / u1`
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- (stats[, 3] - min_value + 1.0e-8) / (max_value - min_value + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    lx <- log1p(stats[, 1])
    tmp <- lm(zv ~ lx)

    # the curve can either increase of decrease depending on the `alpha` and
    # `delta` parameter. However, we want `eta` to be positive. If `eta` is
    # negative we simply change its sign and switch curve direction.
    log_eta <- log(abs(tmp$coefficients[2]))
    log_phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    # find the maximum likelihood estimates of the linear parameters
    mle_asy(object, c(min_value, max_value, log_eta, log_phi, 0, 0))
  } else {
    mle_asy(object, object$start)
  }

  best_rss <- rss_fn(theta)

  # this is a latin hypercube design
  v <- 100
  param_set <- matrix(
    c(
      # log_eta
      -6.43, 2.39, -9.65, -0.13, -9.79, 1.69, 3.79, -7.97, -6.29, -4.47, 1.13,
      -9.37, -6.85, -9.23, 2.11, -5.87, -8.39, -0.55, -2.65, -2.23, -9.09,
      -4.89, 3.93, -2.37, -5.45, -7.41, -4.33, -8.25, 2.95, -2.09, 0.71, -2.93,
      -8.95, 3.09, -4.19, -9.51, -3.21, 0.15, -6.15, 1.97, -0.41, -1.53, -5.59,
      -5.17, -3.07, -1.25, 3.51, 0.29, -2.79, -3.91, 2.81, -7.13, 2.67, -4.05,
      -8.11, -8.53, -3.63, 1.27, 0.01, -6.71, 3.37, -6.57, -0.97, -1.11, 1.83,
      3.65, -1.67, -9.93, -4.75, -0.69, 0.43, -3.49, 1.55, -1.95, 0.85, -4.61,
      -8.81, -5.03, -6.99, -1.81, -2.51, -7.27, -1.39, 0.99, 1.41, 2.25, -7.69,
      0.57, -5.31, -7.55, 3.23, -5.73, -0.27, -3.35, 2.53, -8.67, -3.77, -7.83,
      -0.83, -6.01,
      # log_phi
      11.365, -16.105, 8.905, -6.675, -4.625, -16.515, 13.415, -2.165, -11.185,
      -14.055, 0.705, -8.725, 7.675, -17.745, -16.925, -13.645, 5.215, -7.085,
      -4.215, 6.445, 15.055, -12.825, 3.165, 7.265, -8.315, 18.335, -19.795,
      -3.395, -1.345, 10.545, 12.595, -5.035, 9.725, 18.745, 12.185, 3.575,
      -14.875, 17.105, -5.855, -17.335, 10.955, -2.985, 16.285, -10.365, 8.085,
      -14.465, 1.115, -15.285, -3.805, 6.855, 19.565, -19.385, 4.395, 9.315,
      1.935, -1.755, -9.955, 15.875, 10.135, -13.235, -15.695, 13.825, 17.925,
      -12.415, 8.495, -7.495, 14.235, -10.775, -9.545, -12.005, 15.465, 20.795,
      0.295, 6.035, 3.985, 2.755, 20.385, 16.695, 19.975, 4.805, -18.565,
      -2.575, 2.345, 5.625, -11.595, 14.645, -5.445, -9.135, 19.155, -7.905,
      -0.525, -0.935, -0.115, -18.155, 13.005, 1.525, 17.515, 11.775, -6.265,
      -18.975,
      # log_nu
      -0.74, -1.22, 1.5, -0.7, -0.18, 1.14, 0.5, -1.38, -1.66, -1.18, 0.18,
      -1.86, 0.54, -0.9, 1.1, -1.42, -0.1, 1.7, -1.26, 0.38, 1.66, -0.42, -1.98,
      -1.34, 1.62, 1.18, -0.22, -0.54, 0.74, -1.02, 0.66, 0.58, -0.86, 1.54,
      0.46, -0.58, 0.3, -1.62, -1.46, 1.82, 0.82, 0.34, -1.58, 1.22, -0.66,
      0.86, -0.5, -1.74, -0.94, 1.58, -1.82, 1.74, 0.9, -0.26, 1.02, 1.26, 0.22,
      0.94, -1.54, -0.78, 1.78, -0.02, -0.82, 1.86, 1.46, -1.1, -1.78, 0.62,
      1.06, 0.7, 0.06, -0.46, 1.34, 1.9, -1.5, 1.38, 1.98, -0.98, -1.9, -1.06,
      1.94, -1.3, 0.78, 1.3, -1.94, 0.1, -0.62, -0.14, -1.7, 0.26, -0.3, 0.02,
      -0.06, 0.14, -0.34, -0.38, 1.42, 0.42, 0.98, -1.14,
      # log_xi
      -1.94, -1.14, 1.5, -0.02, 1.14, -0.26, 1.82, -1.98, -1.34, 0.3, 0.46,
      -0.06, 1.78, -1.78, -0.66, 1.7, -1.1, 0.18, 0.9, 1.18, 0.22, -1.18, 0.1,
      -0.58, -0.18, -1.74, -0.82, 1.98, -1.9, 1.86, 0.06, -1.02, 0.34, 1.74,
      1.26, 1.46, -0.14, -0.62, 0.58, 0.86, -1.66, -1.54, 1.42, 1.02, 1.54,
      0.38, 1.62, -1.3, -1.22, 1.1, -0.86, 0.62, 1.34, -0.7, -0.3, -1.7, 1.66,
      1.58, -0.34, -0.98, -1.82, 1.3, -0.22, 0.5, 1.94, -1.06, 0.42, -0.9, 0.14,
      1.22, 0.54, -0.46, 0.26, 1.38, -0.54, -0.78, -1.5, -1.58, -1.38, -0.74,
      -1.86, 0.66, -0.1, 0.98, -1.46, 0.02, 0.7, -1.42, 0.78, -0.94, 0.94, 1.06,
      1.9, -1.26, 0.82, -0.42, -0.5, -0.38, 0.74, -1.62
    ),
    ncol = v, byrow = TRUE
  )

  theta_tmp <- matrix(nrow = 6, ncol = v)
  rss_tmp <- rep(10000, v)

  for (i in seq_len(v)) {
    current_par <- mle_asy(object, c(theta[1], theta[2], param_set[, i]))
    current_rss <- rss_fn(current_par)
    theta_tmp[, i] <- current_par
    rss_tmp[i] <- current_rss
  }

  # update the total iteration count
  # 1: initial crude estimation
  # 2: flat line approximation
  # v: total amount of grid points tested
  niter <- 2 + v

  ord <- order(rss_tmp)

  # select the best solution and other not so good solutions as starting points
  theta_1 <- theta_tmp[, ord[1]]
  theta_2 <- theta_tmp[, ord[5]]
  theta_3 <- theta_tmp[, ord[8]]

  if (object$constrained) {
    theta <- pmax(
      pmin(theta, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_1 <- pmax(
      pmin(theta_1, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_2 <- pmax(
      pmin(theta_2, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_3 <- pmax(
      pmin(theta_3, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
  }

  start <- cbind(theta, theta_1, theta_2, theta_3)

  tmp <- fit_nlminb(object, start, object$max_iter - niter)

  if (!is.infinite(tmp$rss) && (tmp$rss < best_rss)) {
    theta <- tmp$theta
    best_rss <- tmp$rss
    niter <- niter + tmp$niter
  }

  names(theta) <- NULL
  names(niter) <- NULL

  list(theta = theta, niter = niter)
}

# 6-parameter log-logistic fit
#
# Fit a 6-parameter log-logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 6-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (xi * x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
# `phi > 0`, `nu > 0`, and `xi > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, `log(xi)`, and `log(nu)`.
#
# @param object object of class `loglogistic6`.
#
# @return A list with the following components:
#   \describe{
#     \item{converged}{boolean. `TRUE` if the optimization algorithm converged,
#       `FALSE` otherwise.}
#     \item{iterations}{total number of iterations performed by the
#       optimization algorithm}
#     \item{constrained}{boolean. `TRUE` if optimization was constrained,
#       `FALSE` otherwise.}
#     \item{estimated}{boolean vector indicating which parameters were
#       estimated from the data.}
#     \item{coefficients}{maximum likelihood estimates of the model
#       parameters.}
#     \item{rss}{minimum value found of the residual sum of squares.}
#     \item{df.residual}{residual degrees of freedom.}
#     \item{fitted.values}{fitted mean values.}
#     \item{residuals}{residuals, that is response minus fitted values.}
#     \item{weights}{vector of weights used for the fit.}
#   }
fit.loglogistic6 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[3:6] <- exp(theta[3:6])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 6),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 6,
    fitted.values = loglogistic6_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic6_fit", "loglogistic")

  result
}

# @rdname fit.loglogistic6
fit_constrained.loglogistic6 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 6, 3)

  for (i in seq_len(6)) {
    lb_is_inf <- is.infinite(object$lower_bound[i])
    ub_is_inf <- is.infinite(object$upper_bound[i])

    if (object$lower_bound[i] == object$upper_bound[i]) {
      constraint[i, 2] <- TRUE
    } else if (!lb_is_inf || !ub_is_inf) {
      constraint[i, 3] <- TRUE
    } else {
      constraint[i, 1] <- TRUE
    }
  }

  known_param <- ifelse(constraint[, 2], object$lower_bound, NA_real_)

  solution <- find_optimum_constrained(object, constraint, known_param)

  # bring the parameters back to their natural scale
  theta <- object$lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta[3:6] <- exp(theta[3:6])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = loglogistic6_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic6_fit", "loglogistic")

  result
}

# 6-parameter log-logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 6-parameter log-logistic function. We assume that
# our observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 6-by-6 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function. We compute the
# observed Fisher information matrix because it has better finite sample
# properties.
#
# @param object object of class `loglogistic6`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.loglogistic6 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- loglogistic6_gradient_hessian(x, theta)

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 6)
  G[, 1] <- w * z * gh$G[, 1]
  G[, 2] <- w * z * gh$G[, 2]
  G[, 3] <- w * z * gh$G[, 3]
  G[, 4] <- w * z * gh$G[, 4]
  G[, 5] <- w * z * gh$G[, 5]
  G[, 6] <- w * z * gh$G[, 6]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 6, 6))

  H[, , 1] <- w * (z * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (z * gh$H[, , 2] + gh$G[, 2] * gh$G)
  H[, , 3] <- w * (z * gh$H[, , 3] + gh$G[, 3] * gh$G)
  H[, , 4] <- w * (z * gh$H[, , 4] + gh$G[, 4] * gh$G)
  H[, , 5] <- w * (z * gh$H[, , 5] + gh$G[, 5] * gh$G)
  H[, , 6] <- w * (z * gh$H[, , 6] + gh$G[, 6] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  v <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, v)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 6-parameter log-logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `loglogistic6_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.loglogistic6_fit <- function(object, x) {
  len <- length(x)

  V <- object$vcov[seq_len(6), seq_len(6)]

  if (any(is.na(V))) {
    return(rep(NA_real_, len))
  }

  G <- loglogistic6_gradient(x, object$coefficients)

  variance <- rep(NA_real_, len)

  for (i in seq_len(len)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 6-parameter log-logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 6-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (xi * x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu, xi)`, `eta > 0`,
# `phi > 0`, `nu > 0`, and `xi > 0`. The horizontal asymptote is
# `lambda = alpha + delta / xi^(1 / nu)`. When `xi = 1` it is
# obviously `lambda = alpha + delta`.
#
# The area under the curve (AUC) is the integral of `f(x; theta)` with respect
# to `x`.
#
#' @importFrom stats integrate
#'
#' @export
nauc.loglogistic6_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  if (length(xlim) != 2) {
    stop("'xlim' must be of length 2", call. = FALSE)
  }

  if (!is.numeric(xlim)) {
    stop("'xlim' must be a numeric vector of length 2", call. = FALSE)
  }

  # predictor cannot be negative, so we cannot integrate for x < 0
  if (xlim[1] < 0) {
    xlim[1] <- 0
  }

  if (xlim[1] >= xlim[2]) {
    stop("'xlim[1]' cannot be larger or equal to 'xlim[2]'", call. = FALSE)
  }

  if (length(ylim) != 2) {
    stop("'ylim' must be of length 2", call. = FALSE)
  }

  if (!is.numeric(ylim)) {
    stop("'ylim' must be a numeric vector of length 2", call. = FALSE)
  }

  if (ylim[1] >= ylim[2]) {
    stop("'ylim[1]' cannot be larger or equal to 'ylim[2]'", call. = FALSE)
  }

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]
  xi <- object$coefficients[6]

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- phi / (((delta / (ylim - alpha))^nu - xi) / nu)^(1 / eta)

  # value of the integral
  I <- 0

  # check if we really need to perform an integration
  flag <- TRUE

  # we might change the range of integration
  xlim_new <- xlim

  if (delta >= 0) {
    # curve is monotonically increasing
    lb <- alpha
    ub <- alpha + delta / xi^(1 / nu)

    if (lb < ylim[1]) {
      # the curve in `c(0, tmp[1])` is to be considered zero
      if (tmp[1] > xlim[2]) {
        # the integral is simply zero
        flag <- FALSE
      } else if (tmp[1] > xlim[1]) {
        xlim_new[1] <- tmp[1]
      }
    }

    if (ub > ylim[2]) {
      # the curve in `c(tmp[2], Inf)` is equal to `ylim[2]`
      if (tmp[2] < xlim[1]) {
        # not much to do in this case, the curve in the requested range is flat
        I <- I + (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

        # stop here
        flag <- FALSE
      } else if (tmp[2] < xlim[2]) {
        # standard before `tmp[2]` and flat in c(tmp[2], xlim[2])
        I <- I + (xlim[2] - tmp[2]) * (ylim[2] - ylim[1])

        # modify the range of integration
        xlim_new[2] <- tmp[2]
      }
    }
  } else {
    # curve is monotonically decreasing
    lb <- alpha + delta / xi^(1 / nu)
    ub <- alpha

    if (ub > ylim[2]) {
      # the first part of the curve in `c(0, tmp[2])` is equal to `ylim[2]`
      if (tmp[2] > xlim[2]) {
        # not much to do in this case, the curve in the requested range is flat
        I <- I + (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

        # stop here
        flag <- FALSE
      } else if (tmp[2] > xlim[1]) {
        # flat in c(xlim[1], tmp[2]) and standard after `tmp[2]`
        I <- I + (tmp[2] - xlim[1]) * (ylim[2] - ylim[1])

        # modify the range of integration
        xlim_new[1] <- tmp[2]
      }
    }

    if (lb < ylim[1]) {
      # the curve after `tmp[1]` is to be considered zero
      if (tmp[1] < xlim[1]) {
        # the integral is simply zero
        flag <- FALSE
      } else if (tmp[1] < xlim[2]) {
        xlim_new[2] <- tmp[1]
      }
    }
  }

  if (flag) {
    # we remove `ylim[1]` to shift the curve to zero
    f <- function(x) {
      fn(object, x, object$coefficients) - ylim[1]
    }

    I <- I + integrate(
      f, lower = xlim_new[1], upper = xlim_new[2],
      rel.tol = sqrt(.Machine$double.eps)
    )$value
  }

  nauc <- I / ((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))
  names(nauc) <- NULL

  nauc
}

#' @export
naac.loglogistic6_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  1 - nauc.loglogistic6_fit(object, xlim, ylim)
}

#' @export
effective_dose.loglogistic6_fit <- function(object, y, type = "relative") {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]
  xi <- object$coefficients[6]

  # value at -Inf is alpha
  # value at Inf is alpha + delta / xi^(1 / nu)
  fv <- if (type == "relative") {
    y[y <= 0 | y >= 1] <- NA_real_
    alpha + y * delta / xi^(1 / nu)
  } else if (type == "absolute") {
    y1 <- alpha
    y2 <- alpha + delta / xi^(1 / nu)

    if (delta > 0) {
      y[y < y1 | y > y2] <- NA_real_
    } else {
      y[y < y2 | y > y1] <- NA_real_
    }

    y
  } else {
    stop("invalid value for `type`", call. = FALSE)
  }

  z <- (fv - alpha) / delta
  x <- phi * (nu * z^nu / (1 - xi * z^nu))^(1 / eta)
  names(x) <- NULL

  x
}
