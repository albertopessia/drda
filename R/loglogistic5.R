# @rdname loglogistic6_new
loglogistic5_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 5) {
      stop("'start' must be of length 5", call. = FALSE)
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

    start[3:5] <- log(start[3:5])
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
    class = "loglogistic5"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      lower_bound <- rep(-Inf, 5)
    } else {
      if (length(lower_bound) != 5) {
        stop("'lower_bound' must be of length 5", call. = FALSE)
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
    }

    if (is.null(upper_bound)) {
      upper_bound <- rep(Inf, 5)
    } else {
      if (length(upper_bound) != 5) {
        stop("'upper_bound' must be of length 5", call. = FALSE)
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

      upper_bound[3:5] <- log(upper_bound[3:5])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 5-parameter log-logistic function
#'
#' Evaluate at a particular set of parameters the 5-parameter log-logistic
#' function.
#'
#' @details
#' The 5-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`,
#' `phi > 0`, and `nu > 0`.
#'
#' Parameter `alpha` is the value of the function when `x = 0`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the mid-value of the function.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   log-logistic function.
#'
#' @export
loglogistic5_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  t1 <- x^eta
  t2 <- phi^eta

  alpha + delta * (t1 / (t1 + nu * t2))^(1 / nu)
}

# @rdname loglogistic5_fn
fn.loglogistic5 <- function(object, x, theta) {
  loglogistic5_fn(x, theta)
}

# @rdname loglogistic5_fn
fn.loglogistic5_fit <- function(object, x, theta) {
  loglogistic5_fn(x, theta)
}

#' 5-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 5-parameter log-logistic function.
#'
#' @details
#' The 5-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`,
#' `phi > 0`, and `nu > 0`.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
loglogistic5_gradient <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + nu * pe
  h <- (x^eta / f)^(1 / nu)
  d <- delta * h / f

  a <- pe * lr
  q <- (eta * log(x) - log(f)) * f / nu
  s <- pe + q

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- h
  G[, 3] <- a * d
  G[, 4] <- -eta * pe * d / phi
  G[, 5] <- -s * d / nu

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

#' @rdname loglogistic5_gradient
loglogistic5_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + nu * pe
  h <- (x^eta / f)^(1 / nu)
  d <- delta * h / f

  a <- pe * lr
  p <- pe - xe
  q <- (eta * log(x) - log(f)) * f / nu
  r <- d / f
  s <- pe + q

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- a * d / delta
  H[, 4, 2] <- -eta * pe * d / (delta * phi)
  H[, 5, 2] <- -s * d / (delta * nu)

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- lr * a * p * r
  H[, 4, 3] <- -(pe * f + eta * a * p) * r / phi
  H[, 5, 3] <- -(nu * pe + s) * a * r / nu

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- eta * pe * (f + eta * p) * r / phi^2
  H[, 5, 4] <- eta * pe * (nu * pe + s) * r / (nu * phi)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- (nu * pe^2 + (2 * f + s) * s) * r / nu^2

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

#' @rdname loglogistic5_gradient
loglogistic5_gradient_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + nu * pe
  h <- (x^eta / f)^(1 / nu)
  d <- delta * h / f

  a <- pe * lr
  p <- pe - xe
  q <- (eta * log(x) - log(f)) * f / nu
  r <- d / f
  s <- pe + q

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- h
  G[, 3] <- a * d
  G[, 4] <- -eta * pe * d / phi
  G[, 5] <- -s * d / nu

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- a * d / delta
  H[, 4, 2] <- -eta * pe * d / (delta * phi)
  H[, 5, 2] <- -s * d / (delta * nu)

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- lr * a * p * r
  H[, 4, 3] <- -(pe * f + eta * a * p) * r / phi
  H[, 5, 3] <- -(nu * pe + s) * a * r / nu

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- eta * pe * (f + eta * p) * r / phi^2
  H[, 5, 4] <- eta * pe * (nu * pe + s) * r / (nu * phi)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- (nu * pe^2 + (2 * f + s) * s) * r / nu^2

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

#' 5-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 5-parameter log-logistic function.
#'
#' @details
#' The 5-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`,
#' `phi > 0`, and `nu > 0`.
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{loglogistic5_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)`, `phi2 = log(phi)`, and `nu2 = log(nu)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
loglogistic5_gradient_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  k1 <- eta / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + nu * c2
  e <- log(x) - log(phi)

  p <- (x^eta / f)^(1 / nu)
  q <- p / f
  r <- eta * c2 * q
  s <- f * log(f) / nu - c2
  t <- k1 * log(x) * f
  u <- q * s
  v <- q * t

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- p
  G[, 3] <- delta * e * r
  G[, 4] <- -delta * r
  G[, 5] <- -delta * (v - u)

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

#' @rdname loglogistic5_gradient_2
loglogistic5_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  k1 <- eta / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + nu * c2
  e <- log(x) - log(phi)

  l <- (1 + nu) * c2 / f

  p <- (x^eta / f)^(1 / nu)
  q <- p / f
  r <- eta * c2 * q
  s <- f * log(f) / nu - c2
  t <- k1 * log(x) * f
  u <- q * s
  v <- q * t
  y <- eta * log(x) - log(f)

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- e * r
  H[, 4, 2] <- -r
  H[, 5, 2] <- u - v

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * e * (1 + eta * (l - 1) * e) * r
  H[, 4, 3] <- -delta * (1 + eta * (l - 1) * e) * r
  H[, 5, 3] <- -delta * e * (l + y / nu) * r

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (l - 1) * r
  H[, 5, 4] <- delta * (l + y / nu) * r

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * ((l + y / nu) * c2 +
    (1 + y / nu) * (k1 * log(x) * f - s)) * q

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

#' @rdname loglogistic5_gradient_2
loglogistic5_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  k1 <- eta / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + nu * c2
  e <- log(x) - log(phi)

  l <- (1 + nu) * c2 / f

  p <- (x^eta / f)^(1 / nu)
  q <- p / f
  r <- eta * c2 * q
  s <- f * log(f) / nu - c2
  t <- k1 * log(x) * f
  u <- q * s
  v <- q * t
  y <- eta * log(x) - log(f)

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- p
  G[, 3] <- delta * e * r
  G[, 4] <- -delta * r
  G[, 5] <- -delta * (v - u)

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- e * r
  H[, 4, 2] <- -r
  H[, 5, 2] <- u - v

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * e * (1 + eta * (l - 1) * e) * r
  H[, 4, 3] <- -delta * (1 + eta * (l - 1) * e) * r
  H[, 5, 3] <- -delta * e * (l + y / nu) * r

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (l - 1) * r
  H[, 5, 4] <- delta * (l + y / nu) * r

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * ((l + y / nu) * c2 +
    (1 + y / nu) * (k1 * log(x) * f - s)) * q

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

# 5-parameter log-logistic function gradient and Hessian
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 5-parameter log-logistic function.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, `phi > 0`,
# and `nu > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, and `log(nu)`.
#
# @param object object of class `loglogistic5`.
# @param theta numeric vector with the five parameters in the form
#   `c(alpha, delta, log(eta), log(phi), log(nu))`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.loglogistic5 <- function(object, theta) {
  loglogistic5_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 5-parameter log-logistic model.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, `phi > 0`,
# and `nu > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, and `log(nu)`.
#
# @param object object of class `loglogistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.loglogistic5 <- function(object) {
  function(theta) {
    theta[3:5] <- exp(theta[3:5])
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.loglogistic5
rss_fixed.loglogistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:5] <- exp(theta[3:5])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 5-parameter log-logistic model.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, `phi > 0`,
# and `nu > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, and `log(nu)`.
#
# @param object object of class `loglogistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.loglogistic5 <- function(object) {
  function(theta) {
    theta[3:5] <- exp(theta[3:5])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 5, 5))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.loglogistic5
rss_gradient_hessian_fixed.loglogistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:5] <- exp(theta[3:5])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 5, 5))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)

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
# @param object object of class `loglogistic5`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two parameters `alpha`
#   and `delta`.
mle_asy.loglogistic5 <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- exp(theta[4])
  nu <- exp(theta[5])

  s1 <- x^eta
  s2 <- phi^eta

  g <- (s1 / (s1 + nu * s2))^(1 / nu)

  # when phi is extremely large or extremely small, the ratio can be problematic
  # in the case `x = 0` it must be set to zero
  g[x == 0] <- 0

  # when eta is extremely large the ratio might converge to Inf
  # this is the limit for s1 -> Inf
  g[is.infinite(s1)] <- 1

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
# @param object object of class `loglogistic5`.
#
# @return Numeric vector of length 5 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.loglogistic5 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  theta <- if (is.null(object$start)) {
    # we initialize `nu = 1` so that we start with a 4-parameter log-logistic
    # function
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
    mle_asy(object, c(min_value, max_value, log_eta, log_phi, 0))
  } else {
    mle_asy(object, object$start)
  }

  best_rss <- rss_fn(theta)

  # this is a space-filling design adapted from a max entropy grid
  v <- 200
  param_set <- matrix(
    c(
      # log_eta
      -9.9, -9.86, -9.85, -9.83, -9.79, -9.71, -9.69, -9.19, -9.11, -8.94,
      -8.69, -8.68, -8.16, -7.71, -7.47, -7.26, -7.13, -6.79, -6.78, -6.42,
      -6.26, -6.2, -5.94, -5.88, -5.41, -5.04, -4.74, -4.66, -4.21, -4.17,
      -3.63, -3.56, -3.52, -2.88, -2.83, -2.82, -2.36, -2.09, -1.88, -1.87,
      -1.7, -1.48, -1.36, -0.83, -0.71, -0.58, -0.57, -0.57, -0.56, -0.56,
      -0.56, -0.56, -0.55, -0.52, -0.52, -0.47, -0.36, -0.35, -0.32, -0.32,
      -0.29, -0.27, -0.22, -0.19, -0.18, -0.12, -0.11, -0.07, -0.07, -0.04, 0,
      0.03, 0.04, 0.09, 0.1, 0.16, 0.16, 0.21, 0.27, 0.27, 0.28, 0.28, 0.28,
      0.33, 0.34, 0.42, 0.44, 0.44, 0.5, 0.51, 0.51, 0.51, 0.57, 0.68, 0.7,
      0.72, 0.73, 0.76, 0.82, 0.86, 0.86, 0.9, 0.9, 0.94, 1.04, 1.04, 1.04,
      1.05, 1.09, 1.1, 1.13, 1.21, 1.21, 1.25, 1.27, 1.32, 1.32, 1.36, 1.37,
      1.43, 1.43, 1.44, 1.45, 1.47, 1.49, 1.55, 1.6, 1.63, 1.68, 1.72, 1.75,
      1.78, 1.83, 1.85, 1.87, 1.87, 1.88, 1.89, 1.95, 1.97, 2.01, 2.02, 2.09,
      2.1, 2.16, 2.19, 2.21, 2.24, 2.24, 2.26, 2.26, 2.27, 2.27, 2.27, 2.29,
      2.29, 2.29, 2.29, 2.3, 2.3, 2.33, 2.34, 2.35, 2.35, 2.35, 2.35, 2.36,
      2.37, 2.38, 2.64, 2.71, 2.71, 2.74, 2.75, 2.83, 2.92, 2.96, 3.02, 3.04,
      3.13, 3.14, 3.17, 3.19, 3.22, 3.28, 3.43, 3.46, 3.48, 3.5, 3.56, 3.62,
      3.71, 3.81, 3.87, 3.87, 3.89, 3.89, 3.9, 3.91, 3.91,
      # log_phi
      -14.04, 5.15, 4.8, 4.21, -7.23, 3.66, 4.57, -2.59, -18.41, 3.96, 5.29,
      -15.12, -10.19, -1.86, 5.04, -18.32, 4.7, 4.37, -10.62, 3.98, 3.02,
      -14.64, 5.29, -4.82, 4.96, -18.25, 4.67, -10.08, 4, -3.29, 4.35, -13.25,
      5.29, -1.32, 4.98, 3.85, 4.65, -6.74, -18.39, -12.4, 1.93, 3.96, 5.28,
      -18.37, -4.1, -0.17, -2.26, -1.48, -0.74, 0.44, 1.83, 3.69, 0.85, 2.55,
      4.24, 3.09, 1.45, 2.18, -1.59, -0.88, 4.53, 3.91, 0.26, -2.12, 3.41, 4.88,
      -0.41, 1.66, 2.66, 0.8, 3.73, -0.88, -1.7, 2.11, -11.26, 1.5, 3.4, -0.05,
      -2.27, 0.89, 2.43, 2.99, 3.77, 3.95, -0.66, 1.95, -1.32, 0.25, 3.9, -1.93,
      1.26, 5.29, 3.46, 1.92, -0.46, -1.32, 3.83, 2.83, 0.59, -2.28, 1.14,
      -1.06, 0.12, 2.33, -1.84, -0.57, 1.41, 2.88, 3.75, 0.9, -17.75, -7.57,
      0.32, 3.11, 1.93, -2.27, -0.95, 2.46, -0.53, -1.49, 0.21, 0.85, 4.37,
      3.89, 3.34, 1.4, -2.21, 2.15, -0.25, 0.48, -0.95, 3.48, 2.84, -1.79, 1.58,
      3.91, -2.28, 0.92, 2, -0.62, -1.27, 0.1, 3.53, 2.59, 4.66, -0.78, -5.01,
      -0.01, 1.44, -2.28, 5.02, 0.6, 4.12, 5.23, -13.52, -1.5, 3.13, 3.81, 1.91,
      2.4, 3.96, 4.72, -13.3, -7.79, -4.21, 0.91, -17.72, 3.85, 5.28, 3.17,
      -18.16, -12.58, -1.79, -7.38, 4.4, 3.47, -17.96, 3.95, 5.28, -8.67, -3.31,
      -13.02, 3.9, -18.39, 4.78, 0.26, -7.33, -18.03, -11.45, 3.88, -3.09,
      -14.78, 3.94, -11.04, 3.39, -17.21, -0.36, 5.3, -6.79, 4.65,
      # log_nu
      0.34, -1.5, 2.24, -0.62, -2.29, -1.52, -0.19, -1.71, -1, 1.21, 1.7, 0.21,
      -1.2, 2.19, -0.69, -1.07, -0.44, -0.62, 2.3, 1.28, -2.27, -2.24, -1.66,
      1.48, 1.39, -1.7, -1.15, 1.28, 1.59, -1.13, 0.13, -1.51, -2.29, -1.13,
      0.64, 1.91, 1.93, 1.95, 0.24, -0.52, -1.2, 0.36, -0.26, 0.21, 0.88, -0.46,
      1.96, 1.61, 0.46, 2.28, 2.24, -1.51, 1.06, -2.23, -1.82, -1.41, -2.07,
      -0.14, 1.36, 2.29, 1.06, -0.84, -1.31, -0.64, 2, 0.48, 1.28, 2.24, 2.19,
      1.69, -0.01, -0.35, -2.24, -1.92, 2.3, -1.88, -0.74, 0.07, 1.06, 1.74,
      1.8, -0.12, -0.08, 0.43, -2.15, -0.44, 1.96, -1.62, 1.57, 0.93, 0.34,
      -2.29, -1.18, 2.18, 2.14, 0.06, 1.05, -0.37, 0.73, 0.3, 1.39, 1.36,
      -0.84, 1.03, 0.44, 0.75, 1.62, -0.17, 0.74, -2.21, 0.25, 0, -1.96, 2.18,
      -2.03, -1.82, -1.89, 0.55, -0.4, 2.28, -1.78, 1.33, 0.07, -1.98, 1.95,
      0.7, 1.52, 0.79, -1.51, 0.77, 0.91, 0.39, -1.08, 1.53, 2.01, 1.65, 1.82,
      -2.02, 2.27, 2.24, -1.38, 0.85, 0.94, -0.85, 0.66, -2.26, 0.32, -0.25,
      2.26, -2.28, -2.29, 0.9, -0.73, 1.23, -0.51, -1.35, -0.1, 0.04, -1.38,
      0.74, 1.25, 0.26, -0.09, 2.15, 1.4, -1.58, 0.65, -1.94, -2.3, -1.18, 1.9,
      1.02, -2.08, 0.57, -0.2, -2.27, 1.53, -0.42, 1.75, -0.91, -0.79, -0.99,
      -0.9, 1.22, -1.87, -0.77, -1.23, -2.29, -1.33, -1.66, 1.78, -0.7, -0.23,
      0.52, -1.71, 1.83, -0.95, -1.3, -0.57, 1.13
    ),
    ncol = v, byrow = TRUE
  )

  theta_tmp <- matrix(nrow = 5, ncol = v)
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

# 5-parameter log-logistic fit
#
# Fit a 5-parameter log-logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, `phi > 0`,
# and `nu > 0`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)`,
# `log(phi)`, and `log(nu)`.
#
# @param object object of class `loglogistic5`.
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
fit.loglogistic5 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[3:5] <- exp(theta[3:5])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 5),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 5,
    fitted.values = loglogistic5_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic5_fit", "loglogistic")

  result
}

# @rdname fit.loglogistic5
fit_constrained.loglogistic5 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 5, 3)

  for (i in seq_len(5)) {
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
  theta[3:5] <- exp(theta[3:5])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = loglogistic5_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic5_fit", "loglogistic")

  result
}

# 5-parameter log-logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 5-parameter log-logistic function. We assume that
# our observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 5-by-5 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function. We compute the
# observed Fisher information matrix because it has better finite sample
# properties.
#
# @param object object of class `loglogistic5`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.loglogistic5 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- loglogistic5_gradient_hessian(x, theta)

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 5)
  G[, 1] <- w * z * gh$G[, 1]
  G[, 2] <- w * z * gh$G[, 2]
  G[, 3] <- w * z * gh$G[, 3]
  G[, 4] <- w * z * gh$G[, 4]
  G[, 5] <- w * z * gh$G[, 5]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 5, 5))

  H[, , 1] <- w * (z * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (z * gh$H[, , 2] + gh$G[, 2] * gh$G)
  H[, , 3] <- w * (z * gh$H[, , 3] + gh$G[, 3] * gh$G)
  H[, , 4] <- w * (z * gh$H[, , 4] + gh$G[, 4] * gh$G)
  H[, , 5] <- w * (z * gh$H[, , 5] + gh$G[, 5] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  v <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, v)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 5-parameter log-logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `loglogistic5_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.loglogistic5_fit <- function(object, x) {
  len <- length(x)

  V <- object$vcov[seq_len(5), seq_len(5)]

  if (any(is.na(V))) {
    return(rep(NA_real_, len))
  }

  G <- loglogistic5_gradient(x, object$coefficients)

  variance <- rep(NA_real_, len)

  for (i in seq_len(len)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 5-parameter log-logistic fit
#
# Find the dose that produced the observed response.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`,
# `phi > 0`, and `nu > 0`.
#
# This function evaluates the inverse function of `f(x; theta)`, that is
# if `y = fn(x; theta)` then `x = inverse_fn(y; theta)`.
inverse_fn.loglogistic5_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  phi / (((delta / (y - alpha))^nu - 1) / nu)^(1 / eta)
}

# 5-parameter log-logistic fit
#
# Evaluate at a particular point the gradient of the inverse log-logistic
# function.
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^(eta / nu) / (x^eta + nu * phi^eta)^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`,
# `phi > 0`, and `nu > 0`.
#
# This function evaluates the gradient of the inverse function.
inverse_fn_gradient.loglogistic5_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  h <- phi / eta
  z <- delta / (y - alpha)
  s <- z^nu
  u <- nu / (z^nu - 1)
  v <- u^(1 / eta)

  G <- matrix(0, nrow = length(y), ncol = 5)

  G[, 1] <- -h * z * s * u * v / delta
  G[, 2] <- -h * s * u * v / delta
  G[, 3] <- -h * log(u) * v / eta
  G[, 4] <- v
  G[, 5] <- -h * (log(z) * s * u - 1) * v / nu

  G
}

# 5-parameter log-logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 5-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, `phi > 0`,
# and `nu > 0`. The horizontal asymptote is `lambda = alpha + delta`.
#
# The area under the curve (AUC) is the integral of `f(x; theta)` with respect
# to `x`.
#
#' @importFrom stats integrate
#'
#' @export
nauc.loglogistic5_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
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

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- inverse_fn(object, ylim)

  # value of the integral
  I <- 0

  # check if we really need to perform an integration
  flag <- TRUE

  # we might change the range of integration
  xlim_new <- xlim

  if (delta >= 0) {
    # curve is monotonically increasing
    lb <- alpha
    ub <- alpha + delta

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
    lb <- alpha + delta
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
naac.loglogistic5_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  1 - nauc.loglogistic5_fit(object, xlim, ylim)
}

#' @export
effective_dose.loglogistic5_fit <- function(
  object, y, level = 0.95, type = "relative"
) {
  if (level <= 0 || level >= 1) {
    stop("Confidence level must be in the interval (0, 1)", call. = FALSE)
  }

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]

  # value at -Inf is alpha
  # value at Inf is alpha + delta / xi^(1 / nu)
  if (type == "relative") {
    y[y <= 0 | y >= 1] <- NA_real_
    y <- alpha + y * delta
  } else if (type == "absolute") {
    y1 <- alpha
    y2 <- alpha + delta

    if (delta > 0) {
      y[y < y1 | y > y2] <- NA_real_
    } else {
      y[y < y2 | y > y1] <- NA_real_
    }
  } else {
    stop("invalid value for `type`", call. = FALSE)
  }

  x <- inverse_fn(object, y)
  names(x) <- NULL

  V <- object$vcov[seq_len(5), seq_len(5)]
  G <- inverse_fn_gradient(object, y)

  std_err <- if (any(is.na(V))) {
    rep(NA_real_, length(y))
  } else{
    sqrt(diag(tcrossprod(crossprod(t(G), V), G)))
  }
  names(std_err) <- NULL

  q <- qnorm((1 - level) / 2)
  l <- round(level * 100)

  matrix(
    c(
      x,
      x + q * std_err,
      x - q * std_err
    ),
    nrow = length(y),
    ncol = 3,
    dimnames = list(
      round(y, digits = 2),
      c("Estimate", paste0(c("Lower .", "Upper ."), c(l, l)))
    )
  )
}
