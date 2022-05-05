# @rdname logistic6_new
gompertz_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 4) {
      stop("'start' must be of length 4", call. = FALSE)
    }

    if (start[3] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    start[3] <- log(start[3])
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
    class = "gompertz"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      lower_bound <- rep(-Inf, 4)
    } else {
      if (length(lower_bound) != 4) {
        stop("'lower_bound' must be of length 4", call. = FALSE)
      }

      lower_bound[3] <- if (lower_bound[3] > 0) {
        log(lower_bound[3])
      } else {
        -Inf
      }
    }

    if (is.null(upper_bound)) {
      upper_bound <- rep(Inf, 4)
    } else {
      if (length(upper_bound) != 4) {
        stop("'upper_bound' must be of length 4", call. = FALSE)
      }

      if (upper_bound[3] <= 0) {
        stop("'upper_bound[3]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[3] <- log(upper_bound[3])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' Gompertz function
#'
#' Evaluate at a particular set of parameters the Gompertz function.
#'
#' @details
#' The Gompertz function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = exp(-exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)`, `alpha` is the value of the
#' function when `x -> -Inf`, `delta` is the (signed) height of the curve,
#' `eta > 0` is the steepness of the curve or growth rate, and `phi` is related
#' with the value of function at `x = 0`.
#'
#' When `delta < 0` the curve is monotonically decreasing while it is
#' monotonically increasing for `delta > 0`.
#'
#' @param x numeric vector at which the Gompertz function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   Gompertz function.
#'
#' @export
gompertz_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + delta * exp(-exp(-eta * (x - phi)))
}

# @rdname gompertz_fn
fn.gompertz <- function(object, x, theta) {
  gompertz_fn(x, theta)
}

# @rdname gompertz_fn
fn.gompertz_fit <- function(object, x, theta) {
  gompertz_fn(x, theta)
}

#' Gompertz function
#'
#' Evaluate at a particular set of parameters the Gompertz function.
#'
#' @details
#' The Gompertz function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = exp(-exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. When `delta` is
#' positive (negative) the curve is monotonically increasing (decreasing).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
gompertz_gradient <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  u <- b * g
  v <- expm1(-y)
  q <- 1 + y * v

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- delta * z * u
  G[, 4] <- -delta * eta * u

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

#' @rdname gompertz_gradient
gompertz_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  u <- b * g
  v <- expm1(-y)
  q <- 1 + y * v

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- z * u
  H[, 4, 2] <- -eta * u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * z^2 * v * u
  H[, 4, 3] <- -delta * q * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta^2 * v * u

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

#' @rdname gompertz_gradient
gompertz_gradient_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  u <- b * g
  v <- expm1(-y)
  q <- 1 + y * v

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- delta * z * u
  G[, 4] <- -delta * eta * u

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- z * u
  H[, 4, 2] <- -eta * u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * z^2 * v * u
  H[, 4, 3] <- -delta * q * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta^2 * v * u

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

#' Gompertz function
#'
#' Evaluate at a particular set of parameters the Gompertz function.
#'
#' @details
#' The Gompertz function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = exp(-exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. When `delta` is
#' positive (negative) the curve is monotonically increasing (decreasing).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{gompertz_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
gompertz_gradient_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  r <- -eta * b
  u <- r * g

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- -delta * z * u
  G[, 4] <- delta * u

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

#' @rdname gompertz_gradient_2
gompertz_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  r <- -eta * b
  u <- r * g
  v <- expm1(-y)
  q <- 1 + y * v

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- -z * u
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * z * q * u
  H[, 4, 3] <- delta * q * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- -delta * eta * u * v

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

#' @rdname gompertz_gradient_2
gompertz_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  z <- x - phi
  y <- eta * z

  b <- exp(-y)
  g <- exp(-b)

  r <- -eta * b
  u <- r * g
  v <- expm1(-y)
  q <- 1 + y * v

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- -delta * z * u
  G[, 4] <- delta * u

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- -z * u
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * z * q * u
  H[, 4, 3] <- delta * q * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- -delta * eta * u * v

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

# Gompertz function gradient and Hessian
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# Gompertz function.
#
# @details
# The Gompertz function `f(x; theta)` is defined here as
#
# `g(x; theta) = exp(-exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)`, `alpha` is the value of the
# function when `x -> -Inf`, `delta` is the (signed) height of the curve,
# `eta > 0` is the steepness of the curve or growth rate, and `phi` is related
# with the value of function at `x = 0`.
#
# @param object object of class `gompertz`.
# @param theta numeric vector with the four parameters in the form
#   `c(alpha, delta, eta, phi)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.gompertz <- function(object, theta) {
  gompertz_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# Gompertz model.
#
# @details
# The Gompertz function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) exp(-exp(-eta * (x - phi)))`
#
# where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
# asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
# steepness of the curve or growth rate, and `phi` related to the function
# value at `x = 0`.
#
# @param object object of class `gompertz`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.gompertz <- function(object) {
  function(theta) {
    theta[3] <- exp(theta[3])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.gompertz
rss_fixed.gompertz <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3] <- exp(theta[3])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a Gompertz model.
#
# @details
# The Gompertz function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) exp(-exp(-eta * (x - phi)))`
#
# where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
# asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
# steepness of the curve or growth rate, and `phi` related to the function
# value at `x = 0`.
#
# @param object object of class `gompertz`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.gompertz <- function(object) {
  function(theta) {
    theta[3] <- exp(theta[3])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 4, 4))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.gompertz
rss_gradient_hessian_fixed.gompertz <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3] <- exp(theta[3])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 4, 4))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)

    list(
      G = apply(gradient[, idx, drop = FALSE], 2, sum),
      H = apply(hessian[, idx, idx, drop = FALSE], 2:3, sum)
    )
  }
}

# Maximum likelihood estimators
#
# Given a set of parameters, compute the maximum likelihood estimates of the
# lower and upper horizontal asymptotes.
#
# @param object object of class `gompertz`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two asymptotes.
mle_asy.gompertz <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- theta[4]

  g <- exp(-exp(-eta * (x - phi)))

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
# Maximum Likelihood estimator of the four parameters of the Gompertz function.
#
# @param object object of class `gompertz`.
#
# @return Numeric vector of length 4 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.gompertz <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  theta <- if (is.null(object$start)) {
    # y = a + (b - a) * exp(-exp(-e * (x - p)))
    # w = (y - a) / (b - a) = exp(-exp(-e * (x - p)))
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(-log(w)) = -e * (x - p) = e * p - e * x
    #
    # fit `z ~ u0 + u1 x` and set `eta = -u1` and `phi = -u0 / u1`
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- (stats[, 3] - min_value + 1.0e-8) / (max_value - min_value + 2.0e-8)
    zv <- log(-log(zv))
    tmp <- lm(zv ~ stats[, 1])

    log_eta <- log(abs(tmp$coefficients[2]))
    phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    # find the maximum likelihood estimates of the linear parameters
    mle_asy(object, c(min_value, max_value, log_eta, phi))
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
    theta <- if (theta[2] >= 0) {
      c(
        0.9 * weighted_mean + 0.1 * min_value,
        1.0e-3,
        -5,
        object$stats[m, 1] + 100
      )
    } else {
      c(
        0.9 * weighted_mean + 0.1 * max_value,
        -1.0e-3,
        -5,
        object$stats[m, 1] + 100
      )
    }

    best_rss <- rss_fn(theta)
  }

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
      # phi
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
      -18.975
    ),
    ncol = v, byrow = TRUE
  )

  theta_tmp <- matrix(nrow = 4, ncol = v)
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
    # fix the candidates to be within the constraints
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
  }

  niter <- niter + tmp$niter

  names(theta) <- NULL
  names(niter) <- NULL

  list(theta = theta, niter = niter)
}

# Gompertz fit
#
# Fit a Gompertz function to observed data with a Maximum Likelihood approach.
#
# @details
# The Gompertz function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) exp(-exp(-eta * (x - phi)))`
#
# where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
# asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
# steepness of the curve or growth rate, and `phi` related to the function
# value at `x = 0`.
#
# @param object object of class `gompertz`.
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
fit.gompertz <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[3] <- exp(theta[3])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 4),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 4,
    fitted.values = gompertz_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("gompertz_fit", "logistic")

  result
}

# @rdname fit.gompertz
fit_constrained.gompertz <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 4, 3)

  for (i in seq_len(4)) {
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
  theta[3] <- exp(theta[3])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = gompertz_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("gompertz_fit", "logistic")

  result
}

# Gompertz fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the Gompertz function. We assume that our observations
# `y` are independent and such that `y = mu(x; theta) + sigma * epsilon`,
# where `epsilon` has a standard Normal distribution `N(0, 1)`.
#
# The 4-by-4 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function.
#
# @param object object of class `gompertz`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.gompertz <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- gompertz_gradient_hessian(x, theta)

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 4)
  G[, 1] <- w * z * gh$G[, 1]
  G[, 2] <- w * z * gh$G[, 2]
  G[, 3] <- w * z * gh$G[, 3]
  G[, 4] <- w * z * gh$G[, 4]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 4, 4))

  H[, , 1] <- w * (z * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (z * gh$H[, , 2] + gh$G[, 2] * gh$G)
  H[, , 3] <- w * (z * gh$H[, , 3] + gh$G[, 3] * gh$G)
  H[, , 4] <- w * (z * gh$H[, , 4] + gh$G[, 4] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# Gompertz fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `gompertz_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.gompertz_fit <- function(object, x) {
  m <- length(x)

  V <- object$vcov[1:4, 1:4]

  if (any(is.na(V))) {
    return(rep(NA_real_, m))
  }

  G <- gompertz_gradient(x, object$coefficients)

  variance <- rep(NA_real_, m)

  for (i in seq_len(m)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# Gompertz fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The Gompertz function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) exp(-exp(-eta * (x - phi)))`
#
# where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
# asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
# steepness of the curve or growth rate, and `phi` is related to the function
# value at `x = 0`.
#
# The area under the curve (AUC) is simply the integral of `f(x; theta)` with
# respect to `x`.
#
#' @importFrom stats integrate
#' @export
nauc.gompertz_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  if (length(xlim) != 2) {
    stop("'xlim' must be of length 2", call. = FALSE)
  }

  if (!is.numeric(xlim)) {
    stop("'xlim' must be a numeric vector of length 2", call. = FALSE)
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

  if (ylim[1] < 0) {
    stop("'ylim[1]' cannot be negative", call. = FALSE)
  }

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- delta / (ylim - alpha)
  tmp[tmp > 1] <- phi - log(log(tmp[tmp > 1])) / eta

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
      # the curve in `c(-Inf, tmp[1])` is to be considered zero
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
      # the first part of the curve in `c(-Inf, tmp[2])` is equal to `ylim[2]`
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
naac.gompertz_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  1 - nauc.gompertz_fit(object, xlim, ylim)
}

#' @export
effective_dose.gompertz_fit <- function(object, y, type = "relative") {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  # value at -Inf is alpha
  # value at Inf is alpha + delta
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

  x <- phi - log(log(delta / (fv - alpha))) / eta
  names(x) <- NULL

  x
}
