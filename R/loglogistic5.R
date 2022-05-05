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
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
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
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
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
  g <- 1 / f^(1 / nu)
  h <- x^(eta / nu) * g
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
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  d <- g / f

  e <- log(x) - log(theta[4])

  p <- a * g
  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
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
  k2 <- 1 / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + nu * c2
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  c <- k2 * c1
  d <- g / f

  e <- log(x) - log(theta[4])

  l <- (1 + nu) * c2 / f

  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
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
  k2 <- 1 / nu

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + nu * c2
  g <- f^(-1 / nu)

  a <- x^k1
  b <- eta * c2
  c <- k2 * c1
  d <- g / f

  e <- log(x) - log(theta[4])

  l <- (1 + nu) * c2 / f

  p <- a * g
  q <- a * d
  r <- b * q
  s <- f * log(f) / nu - c2
  t <- log(a) * f
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
  # remove names in case they are set
  names(theta) <- NULL

  m <- object$m

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- exp(theta[4])
  nu <- exp(theta[5])

  g <- (x^eta / (x^eta + nu * phi^eta))^(1 / nu)

  t1 <- 0
  t2 <- 0
  t3 <- 0
  t4 <- 0
  t5 <- 0

  for (i in seq_len(m)) {
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

  # theta is so far a very crude estimation of our curve, but is the curve flat?
  y <- object$y
  x <- object$x
  w <- object$w

  idx <- !is.na(y) & !is.na(x) & !is.na(w) & !(w == 0)

  if (sum(idx) != length(object$y)) {
    y <- y[idx]
    x <- x[idx]
    w <- w[idx]
  }

  weighted_mean <- sum(w * y) / sum(w)
  linear_rss <- sum(w * (y - weighted_mean)^2)

  complete_rss <- sum(stats[, 2] * stats[, 4]) + best_rss

  loglik <- loglik_normal(c(linear_rss, complete_rss), m, sum(log(stats[, 2])))

  # variance + intercept -> 2 parameters estimated
  # variance + logistic4 -> 5 parameters estimated
  bic <- log(m) * c(2, 5) - 2 * loglik

  if (bic[1] <= bic[2]) {
    # we are in big problems as a flat horizontal line is likely the best model
    theta <- c(weighted_mean, 0, 0, 0, 0)
    best_rss <- rss_fn(theta)
  }

  # this is a maximum entropy design (generated by DiceDesign::dmaxDesign)
  # we only need a good design, no need to generate it every single time
  param_set <- matrix(
    c(
      -0.4521, -14.8242, 0.3443, -3.8785, 4.26, -0.0538, -3.1888, -8.4104,
      1.904, -5.9856, -16.6958, 0.9613, -8.3553, -1.4436, 1.5947, -6.1151,
      -4.0416, -0.4099, -1.3203, 16.5142, -1.9734, 1.6581, 1.2953, -0.9002,
      1.7268, 12.5047, 0.5134, -9.7705, 11.9881, -1.5569
    ),
    nrow = 3
  )

  v <- ncol(param_set)
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
  theta_2 <- theta_tmp[, ord[round(v / 3)]]
  theta_3 <- theta_tmp[, ord[round(2 * v / 3)]]

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

  # we remove `ylim[1]` to shift the curve to zero
  f <- function(x) {
    fn(object, x, object$coefficients) - ylim[1]
  }

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- phi / (((delta / (ylim - alpha))^nu - 1) / nu)^(1 / eta)

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
    I <- I + integrate(f, lower = xlim_new[1], upper = xlim_new[2])$value
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
effective_dose.loglogistic5_fit <- function(object, y, type = "relative") {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  # value at -Inf is alpha
  # value at Inf is alpha + delta / xi^(1 / nu)
  fv <- if (type == "relative") {
    y[y <= 0 | y >= 1] <- NA_real_
    alpha + y * delta
  } else if (type == "absolute") {
    y1 <- alpha
    y2 <- alpha + delta

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
  x <- phi * (nu * z^nu / (1 - z^nu))^(1 / eta)
  names(x) <- NULL

  x
}
