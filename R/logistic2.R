# @rdname logistic6_new
logistic2_new <-  function(
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

    c(a, d, log(start[1]), start[2])
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
    class = "logistic2"
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

      upper_bound[1] <- log(upper_bound[1])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' This function allows values other than {0, 1, -1} for `alpha` and `delta` but
#' will coerce them to their proper constraints.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, delta, eta, phi)`. `alpha` can only be equal to 0 or 1 while
#'   `delta` can only be equal to 1 or -1.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic2_fn <- function(x, theta) {
  alpha <- 0
  delta <- 1
  if (theta[2] < 0) {
    alpha <- 1
    delta <- -1
  }

  eta <- theta[3]
  phi <- theta[4]

  alpha + delta / (1 + exp(-eta * (x - phi)))
}

# @rdname logistic2_fn
fn.logistic2 <- function(object, x, theta) {
  logistic2_fn(x, c(object$start[1:2], theta))
}

# @rdname logistic2_fn
fn.logistic2_fit <- function(object, x, theta) {
  # within a fit, parameter theta is known exactly
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + delta / (1 + exp(-eta * (x - phi)))
}

#' 2-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
logistic2_gradient <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b
  g <- 1 / f

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- t
  G[, 2] <- u

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

#' @rdname logistic2_gradient
logistic2_hessian <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b
  g <- 1 / f

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- q * t * (2 / f - 1 / b)
  H[, 2, 1] <- (1 / eta + (2 - f / b) * t / g) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

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

#' @rdname logistic2_gradient
logistic2_gradient_hessian <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b
  g <- 1 / f

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- t
  G[, 2] <- u

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- q * t * (2 / f - 1 / b)
  H[, 2, 1] <- (1 / eta + (2 - f / b) * t / g) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

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

#' 2-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
#' are free to vary (therefore the name) while vector `c(alpha, delta)` is
#' constrained to be either `c(0, 1)` (monotonically increasing curve) or
#' `c(1, -1)` (monotonically decreasing curve).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{logistic2_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, delta, eta, phi)`.
#' @param delta value of delta parameter (either 1 or -1).
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
logistic2_gradient_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  s <- 1 / f^2
  t <- q * s
  u <- r * s

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- eta * t
  G[, 2] <- u

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

#' @rdname logistic2_gradient_2
logistic2_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  s <- 1 / f^2
  u <- r * s

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- -y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 2, 1] <- (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

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

#' @rdname logistic2_gradient_2
logistic2_gradient_hessian_2 <- function(x, theta, delta) {
  k <- length(x)

  eta <- theta[1]
  phi <- theta[2]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + b

  q <- y * b
  r <- -eta * b

  s <- 1 / f^2
  t <- q * s
  u <- r * s

  G <- matrix(1, nrow = k, ncol = 2)

  G[, 1] <- eta * t
  G[, 2] <- u

  H <- array(0, dim = c(k, 2, 2))

  H[, 1, 1] <- -y * (1 + eta * (2 / f - 1 / b) * q) * u
  H[, 2, 1] <- (1 + eta * (2 / f - 1 / b) * q) * u

  H[, 1, 2] <- H[, 2, 1]
  H[, 2, 2] <- (2 / f - 1 / b) * r * u

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

# 2-parameter logistic function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 2-parameter logistic function.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name) while vector `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param theta numeric vector with the parameters in the form `c(eta, phi)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic2 <- function(object, theta) {
  logistic2_gradient_hessian_2(object$stats[, 1], theta, object$start[2])
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 2-parameter logistic model.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(p)` to evaluate the RSS associated to a particular
#   parameter choice `p`.
rss.logistic2 <- function(object) {
  function(theta) {
    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.logistic2
rss_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 2-parameter logistic model.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# To avoid issues with the non-negative constraints we consider in our
# optimization algorithm the alternative parameterization with `log(eta)`.
#
# @param object object of class `logistic2`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic2 <- function(object) {
  function(theta) {
    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

# @rdname rss_gradient_hessian.logistic2
rss_gradient_hessian_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[1] <- exp(theta[1])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
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
# lower and upper horizontal asymptotes.
#
# @param object object of class `logistic2`.
# @param theta vector of parameters.
#
# @return Numeric vector `theta`.
mle_asy.logistic2 <- function(object, theta) {
  theta
}

# Initialize vector of parameters
#
# Given the sufficient statistics, try to guess a good approximation to the
# Maximum Likelihood estimator of the four parameters of the logistic function.
#
# @param object object of class `logistic2`.
#
# @return Numeric vector of length 2 with a (hopefully) good starting point.
#
#' @importFrom stats lm
init.logistic2 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  theta <- if (any(is.na(object$start))) {
    # data might not be compatible with a 2-parameter log-logistic function
    idx <- (stats[, 3] >= 0) & (stats[, 3] <= 1)
    xx <- stats[idx, 1]
    yy <- stats[idx, 3]

    # we cannot guarantee that all values are between zero and one, so we assume
    # a 4-parameter logistic model as a starting point
    #
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
    zv <- (yy + 1.0e-8) / (1 + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    tmp <- lm(zv ~ xx)

    # the curve can either increase of decrease depending on the `alpha` and
    # `delta` parameter. However, we want `eta` to be positive. If `eta` is
    # negative we simply change its sign and switch curve direction.
    log_eta <- log(abs(tmp$coefficients[2]))
    phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    c(log_eta, phi)
  } else {
    object$start
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
    theta <- c(-5, object$stats[m, 1] + 100)
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
  theta_2 <- theta_tmp[, ord[round(v / 3)]]
  theta_3 <- theta_tmp[, ord[round(2 * v / 3)]]

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

# 2-parameter logistic fit
#
# Fit a 2-parameter logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# @param object object of class `logistic2`.
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
fit.logistic2 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- c(object$start[1:2], exp(solution$optimum[1]), solution$optimum[2])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(c(FALSE, TRUE), c(2, 2)),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - 2,
    fitted.values = logistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic2_fit", "logistic")

  result
}

# @rdname fit.logistic2
fit_constrained.logistic2 <- function(object) {
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
  theta <- c(object$start[1:2], exp(theta[1]), theta[2])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - sum(estimated),
    fitted.values = logistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic2_fit", "logistic")

  result
}

# 2-parameter logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 2-parameter logistic function. We assume that our
# observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 2-by-2 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function.
#
# @param object object of class `logistic2`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic2 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- logistic2_gradient_hessian(x, theta[3:4], theta[2])

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

  mu <- fn(object, object$x, theta)
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta)[-seq_len(2)], "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 2-parameter logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `logistic2_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.logistic2_fit <- function(object, x) {
  m <- length(x)

  V <- object$vcov[1:2, 1:2]

  if (any(is.na(V))) {
    return(rep(NA_real_, m))
  }

  theta <- object$coefficients
  G <- logistic2_gradient(x, theta[3:4], theta[2])

  variance <- rep(NA_real_, m)

  for (i in seq_len(m)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 2-parameter logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 2-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + exp(-eta * (x - phi)))`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi)` and `eta > 0`. Only `eta` and `phi`
# are free to vary (therefore the name), while `c(alpha, delta)` is
# constrained to be either `c(0, 1)` (monotonically increasing curve) or
# `c(1, -1)` (monotonically decreasing curve).
#
# The area under the curve (AUC) is simply the integral of `f(x; theta)` with
# respect to `x`.
#
#' @export
nauc.logistic2_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  nauc.logistic4_fit(object, xlim, ylim)
}

#' @export
naac.logistic2_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  1 - nauc.logistic4_fit(object, xlim, ylim)
}

#' @export
effective_dose.logistic2_fit <- function(object, y, type = "relative") {
  effective_dose.logistic4_fit(object, y, type)
}
