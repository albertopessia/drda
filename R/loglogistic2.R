# @rdname loglogistic6_new
#
#' @importFrom stats lowess
loglogistic2_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  # 2-parameter log-logistic curve is tricky because according to our
  # parameterization we have two options:
  #
  # x^eta / (x^eta + phi^eta)
  #
  # or
  #
  # 1 - x^eta / (x^eta + phi^eta)
  #
  # this duality creates some problem while initializing the model because we
  # do not know yet if the curve is increasing or decreasing.
  #
  # We will try our best to infer it from the data.

  stats <- suff_stats(x, y, w)
  m <- nrow(stats)

  # depending on the quality of data, the guessing can be hard or easy
  #
  # we use `lowess` as a non-linear non-parametric smoother to get an idea of
  # the trend.
  xx <- sqrt(stats[, 2]) * stats[, 1]
  yy <- sqrt(stats[, 2]) * stats[, 3]
  fit <- lowess(xx, yy, f = 0.9)

  a <- 0
  d <- 1
  if (fit$y[1] > fit$y[m]) {
    a <- 1
    d <- -1
  }

  start <- if (!is.null(start)) {
    if (length(start) != 2) {
      stop("'start' must be of length 2", call. = FALSE)
    }

    if (start[1] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    if (start[2] <= 0) {
      stop("parameter 'phi' cannot be negative nor zero", call. = FALSE)
    }

    c(a, d, log(start))
  } else {
    c(a, d, NA_real_, NA_real_)
  }

  object <- structure(
    list(
      x = x,
      y = y,
      w = w,
      n = length(y),
      stats = stats,
      constrained = FALSE,
      start = start,
      max_iter = max_iter,
      m = m
    ),
    class = "loglogistic2"
  )

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      lower_bound <- rep(-Inf, 2)
    } else {
      if (length(lower_bound) != 2) {
        stop("'lower_bound' must be of length 2", call. = FALSE)
      }

      lower_bound[1] <- if (lower_bound[1] > 0) {
        log(lower_bound[1])
      } else {
        -Inf
      }

      lower_bound[2] <- if (lower_bound[2] > 0) {
        log(lower_bound[2])
      } else {
        -Inf
      }
    }

    if (is.null(upper_bound)) {
      upper_bound <- rep(Inf, 2)
    } else {
      if (length(upper_bound) != 2) {
        stop("'upper_bound' must be of length 2", call. = FALSE)
      }

      if (upper_bound[1] <= 0) {
        stop("'upper_bound[1]' cannot be negative nor zero.", call. = FALSE)
      }

      if (upper_bound[2] <= 0) {
        stop("'upper_bound[2]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound <- log(upper_bound)
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 2-parameter log-logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter log-logistic
#' function.
#'
#' @details
#' The 2-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = x^eta / (x^eta + phi^eta)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name) while
#' vector `c(alpha, delta)` is constrained to be either `c(0, 1)` (monotonically
#' increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#'
#' This function allows values other than {0, 1, -1} for `alpha` and `delta` but
#' will coerce them to their proper constraints.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`. `alpha` can only be equal to 0 or 1 while
#'   `delta` can only be equal to 1 or -1.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   log-logistic function.
#'
#' @export
loglogistic2_fn <- function(x, theta) {
  alpha <- 0
  delta <- 1
  if (theta[2] < 0) {
    alpha <- 1
    delta <- -1
  }

  eta <- theta[3]
  phi <- theta[4]

  t1 <- x^eta
  t2 <- phi^eta

  alpha + delta * t1 / (t1 + t2)
}

# @rdname loglogistic2_fn
fn.loglogistic2 <- function(object, x, theta) {
  loglogistic2_fn(x, c(object$start[1:2], theta))
}

# @rdname loglogistic2_fn
fn.loglogistic2_fit <- function(object, x, theta) {
  # within a fit, parameter theta is known exactly
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  t1 <- x^eta
  t2 <- phi^eta

  alpha + delta * t1 / (t1 + t2)
}

#' 2-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter log-logistic function.
#'
#' @details
#' The 2-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = x^eta / (x^eta + phi^eta)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
#' `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
#' increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
loglogistic2_gradient <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + pe
  h <- xe / f
  d <- h / f

  a <- pe * lr

  G <- matrix(0, nrow = k, ncol = 2)

  G[, 1] <- a * d
  G[, 2] <- -eta * pe * d / phi

  # gradient might not be defined when we plug x = 0 directly into the formula
  # however, the limits for x -> 0 are zero
  G[x_zero, ] <- 0

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

  sign(delta) * G
}

#' @rdname loglogistic2_gradient
loglogistic2_hessian <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + pe
  h <- xe / f
  d <- h / f

  a <- pe * lr
  p <- pe - xe
  r <- d / f

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- lr * a * p * r
  H[, 2, 1] <- -(pe * f + eta * a * p) * r / phi

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- eta * pe * (f + eta * p) * r / phi^2

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

  sign(delta) * H
}

#' @rdname loglogistic2_gradient
loglogistic2_gradient_hessian <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  pe <- phi^eta
  xe <- x^eta
  lr <- log(x / phi)

  f <- xe + pe
  h <- xe / f
  d <- h / f

  a <- pe * lr
  p <- pe - xe
  r <- d / f

  G <- matrix(0, nrow = k, ncol = 2)

  G[, 1] <- a * d
  G[, 2] <- -eta * pe * d / phi

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- lr * a * p * r
  H[, 2, 1] <- -(pe * f + eta * a * p) * r / phi

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- eta * pe * (f + eta * p) * r / phi^2

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero (not w.r.t. alpha)
  G[x_zero, ] <- 0
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

  list(G = sign(delta) * G, H = sign(delta) * H)
}

#' 2-parameter log-logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter log-logistic function.
#'
#' @details
#' The 2-parameter log-logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = x^eta / (x^eta + phi^eta)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
#' `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
#' `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
#' increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{loglogistic2_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)` and `phi2 = log(phi)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
loglogistic2_gradient_2 <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + c2
  e <- log(x) - log(phi)

  q <- (c1 / f) / f
  r <- eta * c2 * q

  G <- matrix(0, nrow = k, ncol = 2)

  G[, 1] <- e * r
  G[, 2] <- -r

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero
  G[x_zero, ] <- 0

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

  sign(delta) * G
}

#' @rdname loglogistic2_gradient_2
loglogistic2_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + c2
  e <- log(x) - log(phi)

  l <- 2 * c2 / f

  q <- (c1 / f) / f
  r <- eta * c2 * q

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- e * (1 + eta * (l - 1) * e) * r
  H[, 2, 1] <- -(1 + eta * (l - 1) * e) * r

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- eta * (l - 1) * r

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

  sign(delta) * H
}

#' @rdname loglogistic2_gradient_2
loglogistic2_gradient_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  x_zero <- x == 0

  eta <- theta[1]
  phi <- theta[2]

  xe <- x^eta
  pe <- phi^eta

  f <- xe + pe
  e <- log(x) - log(phi)

  l <- 2 * pe / f

  q <- (xe / f) / f
  r <- eta * pe * q

  G <- matrix(0, nrow = k, ncol = 2)

  G[, 1] <- e * r
  G[, 2] <- -r

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- e * (1 + eta * (l - 1) * e) * r
  H[, 2, 1] <- -(1 + eta * (l - 1) * e) * r

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- eta * (l - 1) * r

  # gradient and Hessian might not be defined when we plug x = 0 directly into
  # the formula
  # however, the limits for x -> 0 are zero
  G[x_zero, ] <- 0
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

  list(G = sign(delta) * G, H = sign(delta) * H)
}

# 2-parameter log-logistic function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 2-parameter log-logistic function.
#
# @details
# The 2-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^eta / (x^eta + phi^eta)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
# `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
# increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loglogistic2`.
# @param theta numeric vector with the two parameters in the form
#   `c(log(eta), log(phi))`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.loglogistic2 <- function(object, theta) {
  loglogistic2_gradient_hessian_2(object$stats[, 1], theta, object$start[2])
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 2-parameter log-logistic model.
#
# The 2-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^eta / (x^eta + phi^eta)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
# `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
# increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loglogistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.loglogistic2 <- function(object) {
  function(theta) {
    theta <- exp(theta)
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.loglogistic2
rss_fixed.loglogistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta <- exp(theta)

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 2-parameter log-logistic model.
#
# @details
# The 2-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^eta / (x^eta + phi^eta)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
# `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
# increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loglogistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.loglogistic2 <- function(object) {
  function(theta) {
    theta <- exp(theta)

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.loglogistic2
rss_gradient_hessian_fixed.loglogistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta <- exp(theta)

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

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
# @param object object of class `loglogistic2`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two parameters `alpha`
#   and `delta`.
mle_asy.loglogistic2 <- function(object, theta) {
  names(theta) <- NULL
  theta
}

# Initialize vector of parameters
#
# Given the sufficient statistics, try to guess a good approximation to the
# Maximum Likelihood estimator of the six parameters of the log-logistic
# function.
#
# @param object object of class `loglogistic2`.
#
# @return Numeric vector of length 2 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.loglogistic2 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  theta <- if (any(is.na(object$start))) {
    # data might not be compatible with a 2-parameter log-logistic function
    idx <- (stats[, 3] >= 0) & (stats[, 3] <= 1)
    xx <- stats[idx, 1]
    yy <- stats[idx, 3]

    # w = x^e / (x^e + p^e)
    #
    # by construction w is defined in (0, 1).
    #
    # z = log(w / (1 - w)) = - e * log(p) + e * log(d)
    #
    # fit a linear model `z ~ u0 + u1 log(x)` and set `log_eta = log(u1)` and
    # `log_phi = -u0 / u1`
    #
    # we add a very small number to avoid the logarithm of zero.
    zv <- (yy + 1.0e-8) / (1 + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    lx <- log1p(xx)
    tmp <- lm(zv ~ lx)

    # the curve can either increase of decrease depending on the `alpha` and
    # `delta` parameter. However, we want `eta` to be positive. If `eta` is
    # negative we simply change its sign and switch curve direction.
    log_eta <- log(abs(tmp$coefficients[2]))
    log_phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    c(log_eta, log_phi)
  } else {
    object$start[3:4]
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
      # log_phi
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

  theta_tmp <- matrix(nrow = 2, ncol = v)
  rss_tmp <- rep(10000, v)

  for (i in seq_len(v)) {
    theta_tmp[, i] <- param_set[, i]
    rss_tmp[i] <- rss_fn(param_set[, i])
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

# 2-parameter log-logistic fit
#
# Fit a 2-parameter log-logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 2-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^eta / (x^eta + phi^eta)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
# `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
# increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization `log(eta)` and
# `log(phi)`.
#
# @param object object of class `loglogistic2`.
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
fit.loglogistic2 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- c(object$start[1:2], exp(solution$optimum))

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(c(FALSE, TRUE), c(2, 2)),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 2,
    fitted.values = loglogistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic2_fit", "loglogistic")

  result
}

# @rdname fit.loglogistic2
fit_constrained.loglogistic2 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 2, 3)

  for (i in seq_len(2)) {
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
  theta <- c(object$start[1:2], exp(theta))

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = c(FALSE, FALSE, estimated),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = loglogistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("loglogistic2_fit", "loglogistic")

  result
}

# 2-parameter log-logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 2-parameter log-logistic function. We assume that
# our observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 2-by-2 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function. We compute the
# observed Fisher information matrix because it has better finite sample
# properties.
#
# @param object object of class `loglogistic2`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.loglogistic2 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta[3:4]) - y

  gh <- loglogistic2_gradient_hessian(x, theta[3:4], theta[2])

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 2)
  G[, 1] <- w * z * gh$G[, 1]
  G[, 2] <- w * z * gh$G[, 2]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 2, 2))

  H[, , 1] <- w * (z * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (z * gh$H[, , 2] + gh$G[, 2] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta[3:4])
  v <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, v)) / sigma^2

  lab <- c(names(theta)[-seq_len(2)], "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 2-parameter log-logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `loglogistic2_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.loglogistic2_fit <- function(object, x) {
  len <- length(x)

  V <- object$vcov[seq_len(2), seq_len(2)]

  if (any(is.na(V))) {
    return(rep(NA_real_, len))
  }

  theta <- object$coefficients
  G <- loglogistic2_gradient(x, theta[3:4], theta[2])

  variance <- rep(NA_real_, len)

  for (i in seq_len(len)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 2-parameter log-logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 2-parameter log-logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = x^eta / (x^eta + phi^eta)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `x >= 0`, `theta = c(alpha, delta, eta, phi)`, `eta > 0`, and
# `phi > 0`. Only `eta` and `phi` are free to vary (therefore the name), while
# `c(alpha, delta)` are constrained to be either `c(0, 1)` (monotonically
# increasing curve) or `c(1, -1)` (monotonically decreasing curve).
#
# The area under the curve (AUC) is the integral of `f(x; theta)` with respect
# to `x`.
#
#' @export
nauc.loglogistic2_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  nauc.loglogistic4_fit(object, xlim, ylim)
}

#' @export
naac.loglogistic2_fit <- function(object, xlim = c(0, 10), ylim = c(0, 1)) {
  1 - nauc.loglogistic4_fit(object, xlim, ylim)
}

#' @export
effective_dose.loglogistic2_fit <- function(
    object, y, level = 0.95, type = "relative"
) {
  effective_dose.loglogistic4_fit(object, y, level, type)
}
