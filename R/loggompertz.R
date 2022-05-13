# @rdname loglogistic6_new
loggompertz_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 4) {
      stop("'start' must be of length 4", call. = FALSE)
    }

    if (start[3] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    if (start[4] <= 0) {
      stop("parameter 'phi' cannot be negative nor zero", call. = FALSE)
    }

    start[3:4] <- log(start[3:4])
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
    class = "loggompertz"
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

      lower_bound[4] <- if (lower_bound[4] > 0) {
        log(lower_bound[4])
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

      if (upper_bound[4] <= 0) {
        stop("'upper_bound[4]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[3:4] <- log(upper_bound[3:4])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' log-Gompertz function
#'
#' Evaluate at a particular set of parameters the log-Gompertz function.
#'
#' @details
#' The log-Gompertz function `f(x; theta)` is defined here as
#'
#' `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. By convention we set
#' `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   log-logistic function.
#'
#' @export
loggompertz_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- alpha + delta * exp(-(phi / x)^eta)
  f[x == 0] <- alpha

  f
}

# @rdname loggompertz_fn
fn.loggompertz <- function(object, x, theta) {
  loggompertz_fn(x, theta)
}

# @rdname loggompertz_fn
fn.loggompertz_fit <- function(object, x, theta) {
  loggompertz_fn(x, theta)
}

#' log-Gompertz function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' log-Gompertz function.
#'
#' @details
#' The log-Gompertz function `f(x; theta)` is defined here as
#'
#' `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. By convention we set
#' `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
loggompertz_gradient <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  d <- f * g
  q <- a * d

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- -delta * q
  G[, 4] <- -delta * eta * d / phi

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

#' @rdname loggompertz_gradient
loggompertz_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  b <- 1 - f
  d <- f * g
  l <- a * b
  q <- a * d

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- -q
  H[, 4, 2] <- -eta * d / phi

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * l * q
  H[, 4, 3] <- -delta * (1 + eta * l) * d / phi

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (1 - eta * b) * d / phi^2

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

#' @rdname loggompertz_gradient
loggompertz_gradient_hessian <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  b <- 1 - f
  d <- f * g
  l <- a * b
  q <- a * d

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- -delta * q
  G[, 4] <- -delta * eta * d / phi

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- -q
  H[, 4, 2] <- -eta * d / phi

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * l * q
  H[, 4, 3] <- -delta * (1 + eta * l) * d / phi

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * eta * (1 - eta * b) * d / phi^2

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

#' Log-Gompertz function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' log-Gompertz function.
#'
#' @details
#' The log-Gompertz function `f(x; theta)` is defined here as
#'
#' `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. By convention we set
#' `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{loggompertz_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)` and `phi2 = log(phi)`.
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
loggompertz_gradient_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  d <- eta * f * g
  q <- a * d

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- g
  G[, 3] <- -delta * q
  G[, 4] <- -delta * d

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

#' @rdname loggompertz_gradient_2
loggompertz_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  b <- 1 - f
  d <- eta * f * g
  l <- eta * a * b
  q <- a * d
  r <- b * d

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- -q
  H[, 4, 2] <- -d

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * (1 + l) * q
  H[, 4, 3] <- -delta * (1 + l) * d

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- -delta * eta * r

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

#' @rdname loggompertz_gradient_2
loggompertz_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  x_zero <- x == 0

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  f <- (phi / x)^eta
  g <- exp(-f)

  a <- log(phi / x)
  b <- 1 - f
  d <- eta * f * g
  l <- eta * a * b
  q <- a * d
  r <- b * d

  G <- matrix(1, nrow = k, ncol = 4)
  H <- array(0, dim = c(k, 4, 4))

  G[, 2] <- g
  G[, 3] <- -delta * q
  G[, 4] <- -delta * d

  H[, 3, 2] <- -q
  H[, 4, 2] <- -d

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * (1 + l) * q
  H[, 4, 3] <- -delta * (1 + l) * d

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- -delta * eta * r

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

# log-Gompertz function gradient and Hessian
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# log-Gompertz function.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loggompertz`.
# @param theta numeric vector with the four parameters in the form
#   `c(alpha, delta, log(eta), log(phi))`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.loggompertz <- function(object, theta) {
  loggompertz_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a log-Gompertz
# model.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loggompertz`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.loggompertz <- function(object) {
  function(theta) {
    theta[3:4] <- exp(theta[3:4])
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.loggompertz
rss_fixed.loggompertz <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:4] <- exp(theta[3:4])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 4-parameter log-logistic model.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loggompertz`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.loggompertz <- function(object) {
  function(theta) {
    theta[3:4] <- exp(theta[3:4])

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

# @rdname rss_gradient_hessian.loggompertz
rss_gradient_hessian_fixed.loggompertz <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[3:4] <- exp(theta[3:4])

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
# linear parameters.
#
# @param object object of class `loggompertz`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two parameters `alpha`
#   and `delta`.
mle_asy.loggompertz <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  g <- exp(-(exp(theta[4]) / x)^exp(theta[3]))

  # when theta[4] is extremely small and `x = 0` the ratio turns into 0 / 0
  # in such cases, `x = 0` has the priority and `g` must be set to zero
  g[x == 0] <- 0

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
# @param object object of class `loggompertz`.
#
# @return Numeric vector of length 4 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.loggompertz <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  # we now make a guess at the curve direction: is it increasing or decreasing?
  # if it is increasing we should observe a general positive trend in the finite
  # differences.
  pos_trnd <- mean(diff(object$stats[, 3]) >= 0)

  is_increasing <- if (pos_trnd > 0.6) {
    TRUE
  } else if (pos_trnd < 0.4) {
    FALSE
  } else {
    # it appears to be random noise so there is not much we can do about it
    object$stats[1, 3] < object$stats[object$m, 3]
  }

  theta <- if (is.null(object$start)) {
    # y = a + d * exp(-(p / x)^e)
    # w = (y - a) / d = exp(-(p / x)^e)
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(-log(w)) = e * log(p) - e * log(x)
    #
    # fit a linear model `z ~ u0 + u1 log(x)` and set `log_eta = log(-u1)` and
    # `log_phi = -u0 / u1`
    #
    # if the curve is increasing `a` is the minimum value, otherwise it is the
    # maximum.
    # `abs(d)` is the difference between the maximum and minimum value. `d` is
    # negative if the curve is decreasing.
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- if (is_increasing) {
      (stats[, 3] - min_value + 1.0e-8) / (max_value - min_value + 2.0e-8)
    } else {
      (stats[, 3] - max_value - 1.0e-8) / (min_value - max_value - 2.0e-8)
    }
    zv <- log(-log(zv))
    lx <- log1p(stats[, 1])
    tmp <- lm(zv ~ lx)

    # the curve can either increase of decrease depending on the `alpha` and
    # `delta` parameter. However, we want `eta` to be positive. If `eta` is
    # negative we simply change its sign and switch curve direction.
    log_eta <- log(abs(tmp$coefficients[2]))
    log_phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    # find the maximum likelihood estimates of the linear parameters
    mle_asy(object, c(min_value, max_value, log_eta, log_phi))
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
      -14.78, 3.94, -11.04, 3.39, -17.21, -0.36, 5.3, -6.79, 4.65
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

# log-Gompertz fit
#
# Fit a log-Gompertz function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loggompertz`.
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
fit.loggompertz <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[3:4] <- exp(theta[3:4])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 4),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 4,
    fitted.values = loggompertz_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loggompertz_fit", "loglogistic")

  result
}

# @rdname fit.loggompertz
fit_constrained.loggompertz <- function(object) {
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
  theta[3:4] <- exp(theta[3:4])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = loggompertz_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loggompertz_fit", "loglogistic")

  result
}

# log-Gompertz fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the log-Gompertz function. We assume that
# our observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 4-by-4 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function. We compute the
# observed Fisher information matrix because it has better finite sample
# properties.
#
# @param object object of class `loggompertz`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.loggompertz <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- loggompertz_gradient_hessian(x, theta)

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
  v <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, v)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# log-Gompertz fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `loggompertz_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.loggompertz_fit <- function(object, x) {
  len <- length(x)

  V <- object$vcov[seq_len(4), seq_len(4)]

  if (any(is.na(V))) {
    return(rep(NA_real_, len))
  }

  G <- loggompertz_gradient(x, object$coefficients)

  variance <- rep(NA_real_, len)

  for (i in seq_len(len)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# log-Gompertz fit
#
# Find the dose that produced the observed response.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# This function evaluates the inverse function of `f(x; theta)`, that is
# if `y = fn(x; theta)` then `x = inverse_fn(y; theta)`.
inverse_fn.loggompertz_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  phi / log(delta / (y - alpha))^(1 / eta)
}

# log-Gompertz fit
#
# Evaluate at a particular point the gradient of the inverse log-Gompertz
# function.
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# This function evaluates the gradient of the inverse function.
inverse_fn_gradient.loggompertz_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  h <- phi / eta
  z <- delta / (y - alpha)
  u <- 1 / log(z)
  v <- u^(1 / eta)

  G <- matrix(0, nrow = length(y), ncol = 5)

  G[, 1] <- -h * u * z * v / delta
  G[, 2] <- -h * u * v / delta
  G[, 3] <- -h * log(u) * v / eta
  G[, 4] <- v

  G
}

# log-Gompertz fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The log-Gompertz function `f(x; theta)` is defined here as
#
# `f(x; theta) = alpha + delta exp(-(phi / x)^eta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. By convention we set
# `f(0; theta) = lim_{x -> 0} f(x; theta) = alpha`.
#
# The area under the curve (AUC) is the integral of `f(x; theta)` with respect
# to `x`.
#
#' @importFrom stats integrate
#'
#' @export
nauc.loggompertz_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
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

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- c(NA_real_, NA_real_)
  z <- (ylim - alpha) / delta

  if (z[1] > 0 && z[1] < 1) {
    tmp[1] <- phi * (-log(z[1]))^(-1 / eta)
  }

  if (z[2] > 0 && z[2] < 1) {
    tmp[2] <- phi * (-log(z[2]))^(-1 / eta)
  }

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
naac.loggompertz_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  1 - nauc.loggompertz_fit(object, xlim, ylim)
}

#' @export
effective_dose.loggompertz_fit <- function(
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

  V <- object$vcov[seq_len(4), seq_len(4)]
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
