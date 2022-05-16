# @rdname logistic6_new
logistic4_new <-  function(
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
    class = "logistic4"
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

#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)`, `alpha` is the value of the
#' function when `x -> -Inf`, `delta` is the (signed) height of the curve,
#' `eta > 0` is the steepness of the curve or growth rate (also known as the
#' Hill coefficient), and `phi` is the value of `x` at which the curve is equal
#' to its mid-point.
#'
#' When `delta < 0` the curve is monotonically decreasing while it is
#' monotonically increasing for `delta > 0`.
#'
#' The mid-point `f(phi; theta)` is equal to `alpha + delta / 2` while the value
#' of the function for `x -> Inf` is `alpha + delta`.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic4_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + delta / (1 + exp(-eta * (x - phi)))
}

# @rdname logistic4_fn
fn.logistic4 <- function(object, x, theta) {
  logistic4_fn(x, theta)
}

# @rdname logistic4_fn
fn.logistic4_fit <- function(object, x, theta) {
  logistic4_fn(x, theta)
}

#' 4-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
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
logistic4_gradient <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- 1 / f
  G[, 3] <- delta * t
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

#' @rdname logistic4_gradient
logistic4_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- t
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * q * t * (2 / f - 1 / b)
  H[, 4, 3] <- delta * (1 / eta + (2 - f / b) * t * f) * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * (2 / f - 1 / b) * r * u

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

#' @rdname logistic4_gradient
logistic4_gradient_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- 1 / f
  G[, 3] <- delta * t
  G[, 4] <- delta * u

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- t
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * q * t * (2 / f - 1 / b)
  H[, 4, 3] <- delta * (1 / eta + (2 - f / b) * t * f) * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * (2 / f - 1 / b) * r * u

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

#' 4-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. When `delta` is
#' positive (negative) the curve is monotonically increasing (decreasing).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{logistic4_gradient}}. To avoid the non-negative
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
logistic4_gradient_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- 1 / f
  G[, 3] <- delta * eta * t
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

#' @rdname logistic4_gradient_2
logistic4_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- eta * t
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 4, 3] <- delta * (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * (2 / f - 1 / b) * r * u

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

#' @rdname logistic4_gradient_2
logistic4_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  t <- (q / f) / f
  u <- (r / f) / f

  G <- matrix(1, nrow = k, ncol = 4)

  G[, 2] <- 1 / f
  G[, 3] <- delta * eta * t
  G[, 4] <- delta * u

  H <- array(0, dim = c(k, 4, 4))

  H[, 3, 2] <- eta * t
  H[, 4, 2] <- u

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 4, 3] <- delta * (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * (2 / f - 1 / b) * r * u

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

# 4-parameter logistic function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 4-parameter logistic function.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# @param object object of class `logistic4`.
# @param theta numeric vector with the four parameters in the form
#   `c(alpha, delta, eta, phi)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic4 <- function(object, theta) {
  logistic4_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 4-parameter logistic model.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# @param object object of class `logistic4`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.logistic4 <- function(object) {
  function(theta) {
    theta[3] <- exp(theta[3])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.logistic4
rss_fixed.logistic4 <- function(object, known_param) {
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
# against the mean of a 4-parameter logistic model.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# In our optimization algorithm, however, we consider the alternative
# parameterization and `z = log(eta)`.
#
# @param object object of class `logistic4`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic4 <- function(object) {
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

# @rdname rss_gradient_hessian.logistic4
rss_gradient_hessian_fixed.logistic4 <- function(object, known_param) {
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
# @param object object of class `logistic4`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two asymptotes.
mle_asy.logistic4 <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- theta[4]

  g <- 1 / (1 + exp(-eta * (x - phi)))

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
# Maximum Likelihood estimator of the four parameters of the logistic function.
#
# @param object object of class `logistic4`.
#
# @return Numeric vector of length 4 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.logistic4 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  theta <- if (is.null(object$start)) {
    # y = a + (b - a) / (1 + exp(-e * (x - p)))
    # w = (y - a) / (b - a) = 1 / (1 + exp(-e * (x - p)))
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(w / (1 - w)) = - e * p + e * x
    #
    # fit a linear model `z ~ u0 + u1 x` and set `eta = u1` and `phi = -u0 / u1`
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- (stats[, 3] - min_value + 1.0e-8) / (max_value - min_value + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    tmp <- lm(zv ~ stats[, 1])

    log_eta <- log(abs(tmp$coefficients[2]))
    phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    # find the maximum likelihood estimates of the linear parameters
    mle_asy(object, c(min_value, max_value, log_eta, phi))
  } else {
    mle_asy(object, object$start)
  }

  best_rss <- rss_fn(theta)

  # this is a space-filling design using a max entropy grid
  v <- 250
  param_set <- matrix(
    c(
      # log_eta
      -9.9, -9.85, -9.85, -9.74, -9.57, -9.4, -8.66, -8.48, -8.42, -8.37, -8.19,
      -7.96, -7.38, -7.31, -7.24, -7.1, -6.99, -6.88, -6.59, -6.47, -6.39,
      -5.62, -5.52, -5.5, -5.44, -5.4, -4.97, -4.5, -4.39, -4.31, -4.17, -4.08,
      -3.9, -3.49, -3.45, -3.23, -3.03, -2.93, -2.76, -2.72, -2.18, -1.98,
      -1.97, -1.83, -1.16, -0.92, -0.91, -0.85, -0.76, -0.71, -0.58, -0.57,
      -0.56, -0.55, -0.55, -0.54, -0.51, -0.48, -0.44, -0.41, -0.41, -0.37,
      -0.29, -0.29, -0.28, -0.26, -0.24, -0.24, -0.21, -0.18, -0.14, -0.13,
      -0.11, -0.08, -0.06, -0.05, -0.04, -0.03, -0.03, 0.02, 0.02, 0.02, 0.07,
      0.08, 0.1, 0.16, 0.16, 0.17, 0.18, 0.2, 0.22, 0.22, 0.24, 0.3, 0.31, 0.38,
      0.38, 0.4, 0.4, 0.4, 0.42, 0.43, 0.44, 0.45, 0.5, 0.5, 0.53, 0.57, 0.58,
      0.59, 0.62, 0.62, 0.66, 0.67, 0.69, 0.71, 0.71, 0.73, 0.73, 0.75, 0.78,
      0.8, 0.8, 0.85, 0.88, 0.9, 0.92, 0.93, 0.93, 0.93, 0.96, 0.98, 0.99, 1.02,
      1.03, 1.05, 1.06, 1.06, 1.07, 1.14, 1.16, 1.18, 1.23, 1.23, 1.24, 1.24,
      1.25, 1.27, 1.27, 1.28, 1.34, 1.34, 1.39, 1.4, 1.42, 1.42, 1.44, 1.45,
      1.46, 1.51, 1.51, 1.52, 1.55, 1.57, 1.59, 1.62, 1.63, 1.63, 1.67, 1.69,
      1.71, 1.72, 1.74, 1.81, 1.83, 1.84, 1.85, 1.87, 1.89, 1.92, 1.94, 1.95,
      1.97, 1.99, 2.02, 2.04, 2.08, 2.09, 2.1, 2.13, 2.16, 2.16, 2.17, 2.25,
      2.26, 2.26, 2.27, 2.28, 2.28, 2.3, 2.34, 2.35, 2.35, 2.36, 2.37, 2.43,
      2.5, 2.58, 2.58, 2.59, 2.65, 2.65, 2.71, 2.79, 2.84, 2.85, 2.88, 2.89,
      2.93, 2.93, 3.01, 3.02, 3.03, 3.11, 3.16, 3.2, 3.22, 3.25, 3.26, 3.27,
      3.39, 3.43, 3.43, 3.44, 3.45, 3.55, 3.6, 3.61, 3.61, 3.65, 3.65, 3.66,
      3.69, 3.85, 3.86, 3.87, 3.87, 3.89, 3.9, 3.9,
      # phi
      -12.25, 3.53, 8.28, -6.05, -17.21, -1.82, -0.02, -11.59, -16.72, -6.27,
      4.61, 9.07, -18.23, -9, 1.53, 5.53, -5.32, -1.69, 9.13, -12.86, -17.96,
      3.49, -9.37, 8.03, -5.24, 0.2, -13.39, -17.61, 0.88, -7.91, -4.03, 4.7, 9,
      -17.95, -13.12, -3.85, -8.09, 1.44, 9.45, 5.79, -11.62, -18.36, -3.68,
      0.78, 9.26, -7.58, -17.68, -2.31, 3.68, -12.74, -14.52, -4.81, 6.45,
      -12.14, 0.35, 2.89, -3.06, -17.99, -9.5, -5.94, -0.21, 4.97, 8.01, 1.61,
      -3.65, -7.53, -10.34, -13.82, -8.64, 0.18, 4.74, -2.09, -4.92, 7.67, 9.59,
      2.99, -0.07, -11.28, -6.83, -16.4, -13.69, 5.14, -9.66, -5.27, -18.29,
      8.41, -1.15, -2.76, 1.18, -9.91, -11.91, -6.4, 3.31, -4.13, -14.01, 6.97,
      -0.04, 1.79, 8.83, -8.63, -18.1, -2.51, -5.99, -14.78, 3.61, -10.88,
      -3.75, -0.87, -7.73, 2.14, 9.13, -13.17, -5.03, 6.38, -11.26, -2.16,
      -15.31, 1.11, -8.76, 8.2, 4.18, -17.97, -12.89, -6.3, -0.83, -4.8, 9.89,
      3.28, -9.76, -16.69, 1.34, 7.25, -14.92, -11.55, -3.61, -1.98, -17.51,
      5.39, -8.09, 0.95, -5.87, -15.77, 3, 7.98, -10.22, -0.05, -13.02, -3.4,
      -8.01, 5.44, -18.3, -6.65, -15.22, 4.64, -4.98, -11.01, 7.19, 1.11, -2.26,
      -12.01, -7.66, -16.52, -0.35, -4.32, -10.14, 3.35, 5.1, 9.65, -14.17,
      -6.28, -2.48, -17.94, 7.13, 0.83, -8.48, -10.85, -15.71, -1.41, 5.08,
      -3.85, 2.56, -13.28, 9.13, -6.3, -8.42, -2.57, 0.28, -15.34, 8.97, 2.19,
      5.12, -4.08, -8.62, -6.04, -14.71, -17.59, -12.24, -1.25, -10.56, 3.02,
      -5.34, -8.96, 3.05, 8.54, -13.61, -0.98, -17.5, -12.16, 8.26, 3.2, -0.31,
      -6.21, -16.95, -12.03, -8.32, 6.13, 2.53, -3.37, -17.61, 9.38, 5.03,
      -10.68, 0.05, -6.09, -17.42, -12.32, 9.65, 4.34, -2.26, -8.79, -17.26,
      2.25, -4.6, -10.98, 9.29, 6.59, -0.72, -18.35, -13.75, -7.78, 9.9, 3.45,
      -3.53, -0.87, -9.51, -17.21, 9.72, -6.03, 5.9, 2.5
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

# 4-parameter logistic fit
#
# Fit a 4-parameter logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# @param object object of class `logistic4`.
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
fit.logistic4 <- function(object) {
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
    fitted.values = logistic4_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic4_fit", "logistic")

  result
}

# @rdname fit.logistic4
fit_constrained.logistic4 <- function(object) {
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
    fitted.values = logistic4_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic4_fit", "logistic")

  result
}

# 4-parameter logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 4-parameter logistic function. We assume that our
# observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 4-by-4 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function.
#
# @param object object of class `logistic4`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic4 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- logistic4_gradient_hessian(x, theta)

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

# 4-parameter logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `logistic4_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.logistic4_fit <- function(object, x) {
  m <- length(x)

  V <- object$vcov[1:4, 1:4]

  if (any(is.na(V))) {
    return(rep(NA_real_, m))
  }

  G <- logistic4_gradient(x, object$coefficients)

  variance <- rep(NA_real_, m)

  for (i in seq_len(m)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 4-parameter logistic fit
#
# Find the dose that produced the observed response.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# This function evaluates the inverse function of `f(x; theta)`, that is
# if `y = fn(x; theta)` then `x = inverse_fn(y; theta)`.
inverse_fn.logistic4_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  x <- delta / (y - alpha)
  x[x > 1] <- phi - log(x[x > 1] - 1) / eta

  x
}

# 4-parameter logistic fit
#
# Evaluate at a particular point the gradient of the inverse logistic function.
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`.
#
# This function evaluates the gradient of the inverse function.
inverse_fn_gradient.logistic4_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]

  z <- delta / (y - alpha)
  u <- 1 / (z - 1)

  G <- matrix(1, nrow = length(y), ncol = 4)

  G[, 1] <- -z * z * u / (delta * eta)
  G[, 2] <- -z * u / (delta * eta)
  G[, 3] <- -log(u) / eta^2

  G
}

# 4-parameter logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 4-parameter logistic function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#
# where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
# asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
# steepness of the curve or growth rate (also known as the Hill coefficient),
# and `phi` is the value of `x` at which the curve is equal to its mid-point.
#
# The area under the curve (AUC) is simply the integral of `f(x; theta)` with
# respect to `x`.
#
#' @export
nauc.logistic4_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  if (length(xlim) != 2) {
    stop("'xlim' must be of length 2", call. = FALSE)
  }

  if (!is.numeric(xlim)) {
    stop("'xlim' must be a numeric vector", call. = FALSE)
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
    I <- I +
      (xlim_new[2] - xlim_new[1]) * (alpha + delta - ylim[1]) +
      delta * (
        log1p(exp(-eta * (xlim_new[2] - phi))) -
        log1p(exp(-eta * (xlim_new[1] - phi)))
      ) / eta
  }

  nauc <- I / ((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))
  names(nauc) <- NULL

  nauc
}

#' @export
naac.logistic4_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  1 - nauc.logistic4_fit(object, xlim, ylim)
}

#' @export
effective_dose.logistic4_fit <- function(
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
