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
  # remove names in case they are set
  names(theta) <- NULL

  m <- object$m

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
    theta <- c(weighted_mean, 0, 0, 0)
    best_rss <- rss_fn(theta)
  }

  # this is a maximum entropy design (generated by DiceDesign::dmaxDesign)
  # we only need a good design, no need to generate it every single time
  param_set <- matrix(
    c(
      -0.4521, -14.8242, -3.8785, 4.26, -3.1888, -8.4104, -5.9856, -16.6958,
      -8.3553, -1.4436, -6.1151, -4.0416, -1.3203, 16.5142, 1.6581, 1.2953,
      1.7268, 12.5047, -9.7705, 11.9881
    ),
    nrow = 2
  )

  v <- ncol(param_set)
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

  # we remove `ylim[1]` to shift the curve to zero
  f <- function(x) {
    fn(object, x, object$coefficients) - ylim[1]
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
    I <- I + integrate(f, lower = xlim_new[1], upper = xlim_new[2])$value
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
effective_dose.loggompertz_fit <- function(object, y, type = "relative") {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

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
  x <- phi / (-log(z))^(1 / eta)
  names(x) <- NULL

  x
}
