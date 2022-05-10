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
  g <- 1 / f
  h <- xe * g
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
  g <- 1 / f
  h <- xe * g
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
  g <- 1 / f
  h <- xe * g
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
  g <- 1 / f

  d <- g / f
  e <- log(x) - log(phi)

  q <- c1 * d
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
  g <- 1 / f

  d <- g / f
  e <- log(x) - log(phi)

  l <- 2 * c2 / f

  q <- c1 * d
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

  c1 <- x^eta
  c2 <- phi^eta

  f <- c1 + c2
  g <- 1 / f

  d <- g / f
  e <- log(x) - log(phi)

  l <- 2 * c2 / f

  q <- c1 * d
  r <- eta * c2 * q

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
  z <- fn(object, x, theta) - y

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
effective_dose.loglogistic2_fit <- function(object, y, type = "relative") {
  effective_dose.loglogistic4_fit(object, y, type)
}
