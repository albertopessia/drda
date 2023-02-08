# @rdname logistic6_new
logistic5_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 5) {
      stop("'start' must be of length 5", call. = FALSE)
    }

    if (start[3] <= 0) {
      stop("parameter 'eta' cannot be negative nor zero", call. = FALSE)
    }

    if (start[5] <= 0) {
      stop("parameter 'nu' cannot be negative nor zero", call. = FALSE)
    }

    start[c(3, 5)] <- log(start[c(3, 5)])
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
    class = "logistic5"
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

      if (upper_bound[5] <= 0) {
        stop("'upper_bound[5]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[c(3, 5)] <- log(upper_bound[c(3, 5)])
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 5-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#'
#' When `delta` is positive (negative) the curve is monotonically increasing
#' (decreasing). When `x -> -Inf` the value of the function is `alpha` while the
#' value of the function for `x -> Inf` is `alpha + delta `.
#'
#' Parameter `alpha` is the value of the function when `x -> -Inf`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the mid-value of the function.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic5_fn <- function(x, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  alpha + delta / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)
}

fn.logistic5 <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

fn.logistic5_fit <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

#' 5-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`. When
#' `delta` is positive (negative) the curve is monotonically increasing
#' (decreasing).
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Gradient or Hessian evaluated at the specified point.
#'
#' @export
logistic5_gradient <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- g
  G[, 3] <- delta * t
  G[, 4] <- delta * u
  G[, 5] <- delta * v / nu

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

gradient.logistic5_fit <- function(object, x) {
  logistic5_gradient(x, object$coefficients)
}

#' @export
#'
#' @rdname logistic5_gradient
logistic5_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- t
  H[, 4, 2] <- u
  H[, 5, 2] <- v / nu

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * q * t * ((1 + nu) / f - 1 / b)
  H[, 4, 3] <- delta * (1 / eta + (1 + nu - f / b) * t / g) * u
  H[, 5, 3] <- delta * (nu * u / eta + v) * t / (nu * g)

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * ((1 + nu) / f - 1 / b) * r * u
  H[, 5, 4] <- delta * (nu * u / eta + v) * u / (nu * g)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * (nu * (u / eta)^2 + v * (v - 2 * g)) / (nu^2 * g)

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

#' @export
#'
#' @rdname logistic5_gradient
logistic5_gradient_hessian <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- g
  G[, 3] <- delta * t
  G[, 4] <- delta * u
  G[, 5] <- delta * v / nu

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- t
  H[, 4, 2] <- u
  H[, 5, 2] <- v / nu

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- delta * q * t * ((1 + nu) / f - 1 / b)
  H[, 4, 3] <- delta * (1 / eta + (1 + nu - f / b) * t / g) * u
  H[, 5, 3] <- delta * (nu * u / eta + v) * t / (nu * g)

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * ((1 + nu) / f - 1 / b) * r * u
  H[, 5, 4] <- delta * (nu * u / eta + v) * u / (nu * g)

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * (nu * (u / eta)^2 + v * (v - 2 * g)) / (nu^2 * g)

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

#' 5-parameter logistic function gradient and Hessian
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined here as
#'
#' `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#' `f(x; theta) = alpha + delta g(x; theta)`
#'
#' where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`. When
#' `delta` is positive (negative) the curve is monotonically increasing
#' (decreasing).
#'
#' This set of functions use a different parameterization from
#' \code{link[drda]{logistic5_gradient}}. To avoid the non-negative
#' constraints of parameters, the gradient and Hessian computed here are for
#' the function with `eta2 = log(eta)` and `nu2 = log(nu)`.
#'
#' Note that argument `theta` is on the original scale and not on the log scale.
#'
#' @param x numeric vector at which the function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, delta, eta, phi, nu)`.
#'
#' @return Gradient or Hessian of the alternative parameterization evaluated at
#'   the specified point.
#'
#' @export
logistic5_gradient_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- y * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- g
  G[, 3] <- delta * eta * t
  G[, 4] <- delta * u
  G[, 5] <- delta * v

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

#' @export
#'
#' @rdname logistic5_gradient_2
logistic5_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- y * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- eta * t
  H[, 4, 2] <- u
  H[, 5, 2] <- v

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * y * (1 + eta * ((1 + nu) / f - 1 / b) * q) * u
  H[, 4, 3] <- delta * (1 + eta * ((1 + nu) / f - 1 / b) * q) * u
  H[, 5, 3] <- -delta * y * (nu * u / eta + v) * u / g

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * ((1 + nu) / f - 1 / b) * r * u
  H[, 5, 4] <- delta * (nu * u / eta + v) * u / g

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * (nu * (u / eta)^2 + v * (v - g)) / g

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

#' @export
#'
#' @rdname logistic5_gradient_2
logistic5_gradient_hessian_2 <- function(x, theta) {
  k <- length(x)

  delta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  y <- x - phi

  b <- exp(-eta * y)

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- y * b
  r <- -eta * b

  t <- (q * g) / f
  u <- (r * g) / f
  v <- u / eta + g * log(f) / nu

  G <- matrix(1, nrow = k, ncol = 5)

  G[, 2] <- g
  G[, 3] <- delta * eta * t
  G[, 4] <- delta * u
  G[, 5] <- delta * v

  H <- array(0, dim = c(k, 5, 5))

  H[, 3, 2] <- eta * t
  H[, 4, 2] <- u
  H[, 5, 2] <- v

  H[, 2, 3] <- H[, 3, 2]
  H[, 3, 3] <- -delta * y * (1 + eta * ((1 + nu) / f - 1 / b) * q) * u
  H[, 4, 3] <- delta * (1 + eta * ((1 + nu) / f - 1 / b) * q) * u
  H[, 5, 3] <- -delta * y * (nu * u / eta + v) * u / g

  H[, 2, 4] <- H[, 4, 2]
  H[, 3, 4] <- H[, 4, 3]
  H[, 4, 4] <- delta * ((1 + nu) / f - 1 / b) * r * u
  H[, 5, 4] <- delta * (nu * u / eta + v) * u / g

  H[, 2, 5] <- H[, 5, 2]
  H[, 3, 5] <- H[, 5, 3]
  H[, 4, 5] <- H[, 5, 4]
  H[, 5, 5] <- delta * (nu * (u / eta)^2 + v * (v - g)) / g

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

# 5-parameter logistic function gradient and Hessian
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 5-parameter logistic function.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# @param object object of class `logistic5`.
# @param theta numeric vector with the five parameters in the form
#   `c(alpha, delta, eta, phi, nu)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic5 <- function(object, theta) {
  logistic5_gradient_hessian_2(object$stats[, 1], theta)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 5-parameter logistic model.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider the alternative
# parameterization `v = log(nu)` and `z = log(eta)`.
#
# @param object object of class `logistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.logistic5 <- function(object) {
  function(theta) {
    theta[c(3, 5)] <- exp(theta[c(3, 5)])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.logistic5
rss_fixed.logistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[c(3, 5)] <- exp(theta[c(3, 5)])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# Residual sum of squares
#
# Evaluate the gradient and Hessian of the residual sum of squares (RSS)
# against the mean of a 5-parameter logistic model.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider the alternative
# parameterization `v = log(nu)` and `z = log(eta)`.
#
# @param object object of class `logistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic5 <- function(object) {
  function(theta) {
    theta[c(3, 5)] <- exp(theta[c(3, 5)])

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

# @rdname rss_gradient_hessian.logistic5
rss_gradient_hessian_fixed.logistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[c(3, 5)] <- exp(theta[c(3, 5)])

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
# lower and upper horizontal asymptotes.
#
# @param object object of class `logistic5`.
# @param theta vector of parameters.
#
# @return Numeric vector of length 2 with the MLE of the two asymptotes.
mle_asy.logistic5 <- function(object, theta) {
  names(theta) <- NULL

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- exp(theta[3])
  phi <- theta[4]
  nu <- exp(theta[5])

  g <- (1 + nu * exp(-eta * (x - phi)))^(-1 / nu)

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
# Maximum Likelihood estimator of the five parameters of the logistic function.
#
# @param object object of class `logistic5`.
#
# @return Numeric vector of length 5 with a (hopefully) good starting point.
#
#' @importFrom stats lm
#'
#' @noRd
init.logistic5 <- function(object) {
  stats <- object$stats
  rss_fn <- rss(object)

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  theta <- if (is.null(object$start)) {
    # we initialize `xi = 1` so that we start with a 4-parameter logistic
    # function
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
    zv <- (stats[, 3] - min_value + 1.0e-8) / (max_value - min_value + 2.0e-8)
    zv <- log(zv) - log1p(-zv)
    tmp <- lm(zv ~ stats[, 1])

    log_eta <- log(abs(tmp$coefficients[2]))
    phi <- -tmp$coefficients[1] / tmp$coefficients[2]

    # find the maximum likelihood estimates of the linear parameters
    mle_asy(object, c(min_value, max_value, log_eta, phi, 0))
  } else {
    mle_asy(object, object$start)
  }

  best_rss <- rss_fn(theta)

  # this is a space-filling design using a max entropy grid
  v <- 250
  param_set <- matrix(
    c(
      # log_eta
      -9.82, -9.78, -9.74, -9.71, -9.64, -9.55, -9.55, -9.51, -9.48, -9.46,
      -9.25, -8.95, -8.83, -8.75, -8.73, -8.7, -8.6, -8.58, -8.57, -8.56, -8.55,
      -8.34, -8.28, -8.19, -8.18, -8.08, -8.06, -8.01, -7.95, -7.8, -7.74,
      -7.73, -7.67, -7.56, -7.5, -7.46, -7.4, -7.25, -7.18, -7.17, -7.17, -7.05,
      -7.05, -7.04, -6.95, -6.92, -6.92, -6.92, -6.89, -6.66, -6.63, -6.56,
      -6.54, -6.38, -6.32, -6.26, -6.26, -6.19, -6.14, -6.14, -6.08, -5.99,
      -5.9, -5.89, -5.88, -5.88, -5.77, -5.75, -5.72, -5.65, -5.54, -5.51,
      -5.36, -5.35, -5.32, -5.28, -5.1, -5.08, -5.03, -5, -4.98, -4.94, -4.93,
      -4.92, -4.9, -4.9, -4.82, -4.73, -4.67, -4.48, -4.46, -4.43, -4.37, -4.32,
      -4.29, -4.27, -4.23, -4.17, -4.17, -4.06, -4.02, -3.94, -3.94, -3.92,
      -3.92, -3.87, -3.77, -3.7, -3.65, -3.59, -3.57, -3.57, -3.55, -3.3, -3.3,
      -3.28, -3.27, -3.18, -3.15, -3.09, -3.04, -2.97, -2.93, -2.91, -2.88,
      -2.88, -2.85, -2.82, -2.77, -2.72, -2.64, -2.64, -2.61, -2.61, -2.54,
      -2.3, -2.25, -2.21, -2.18, -2.17, -2.15, -2.09, -2, -1.97, -1.93, -1.88,
      -1.84, -1.83, -1.82, -1.79, -1.6, -1.6, -1.57, -1.56, -1.55, -1.52, -1.36,
      -1.35, -1.23, -1.23, -1.14, -1.13, -1.11, -1.02, -1, -0.99, -0.95, -0.89,
      -0.78, -0.73, -0.62, -0.61, -0.58, -0.5, -0.46, -0.41, -0.39, -0.23, -0.2,
      -0.17, -0.15, -0.14, -0.1, -0.05, 0.03, 0.04, 0.1, 0.1, 0.16, 0.19, 0.31,
      0.36, 0.36, 0.37, 0.49, 0.49, 0.51, 0.61, 0.63, 0.75, 0.81, 0.81, 0.84,
      0.87, 0.93, 1.12, 1.24, 1.27, 1.32, 1.49, 1.52, 1.56, 1.58, 1.68, 1.72,
      1.75, 1.92, 2, 2.09, 2.1, 2.1, 2.18, 2.19, 2.23, 2.23, 2.46, 2.53, 2.55,
      2.67, 2.71, 2.72, 2.79, 2.86, 2.88, 2.94, 3.03, 3.04, 3.15, 3.21, 3.52,
      3.54, 3.63, 3.66, 3.72, 3.74, 3.74, 3.75, 3.88, 3.9, 3.91,
      # phi
      -10.72, 7.43, -4.6, -6.85, -0.37, -2.6, -8.5, 3.38, 9.46, -13.59, -16.39,
      -3.49, -14.78, 0.74, -10.22, -6.25, -1.7, -12.81, 4.79, -8.57, 3.19, 1.81,
      -7.21, -4.09, -11.67, 0.15, -16.69, -2.81, -9.15, 6.01, -1.19, -5.17,
      -14.25, -7.96, 2.49, -4.26, -6.55, -11.14, 1.26, 6.58, -12.95, -2.64,
      4.84, -9.6, 8.9, 3.17, -7.53, -18.06, 0.26, -4.35, -10.42, -6.47, -13.03,
      -1.33, -14.2, 2.12, 6.76, -8.75, 0.46, 4.01, -2.73, -11.54, -15.22, -4.41,
      -17.91, -1.17, 9, -10.2, -6.52, 1.4, -7.86, -12.49, -2.92, -13.63, 5.29,
      -0.02, -6.47, -15.27, -10.03, -8.35, 6.57, -1.93, 8.68, -4.72, 3.96, 2.14,
      -17.65, -13.44, -11.54, -0.08, 5.82, -9.7, 8.17, -15.14, -3.17, -17.6,
      1.89, -7.01, -3.98, -1.19, 3.47, 9.37, -5.03, -11.52, -14.03, 0.64, 6.51,
      -8.92, 4.92, 7.72, -2.57, -18.31, -0.95, -11.32, -16.32, 1.84, -9.45,
      -7.06, 2.95, 4.18, -4.35, 7.3, 9.69, -14.05, 8.25, -12.54, -8.7, -2.83,
      -5.73, -0.28, 1.32, 6.36, -1.26, -11.39, -4.62, -9.87, -14.09, -18.06,
      7.83, 4.47, -15.74, -8.37, -12.23, -2.55, 2.96, -5.32, 0.52, -4.18, -6.71,
      5.89, -13.33, -1.48, -9.77, -18.29, 7.12, 9.22, -16.75, 2.13, -5.8, -7.78,
      -14.42, 3.55, -2.58, -4.49, -11.28, -0.51, -9.87, 5.09, 7.42, -18.31,
      -1.37, -12.31, 3.21, 0.9, -7.48, -2.86, -16.28, -11.05, 6.88, -5.19,
      -14.1, -3.68, -9.66, 1.38, -15.52, 9.36, -1.68, 3.95, -6.08, -17.54, 5.98,
      2.22, -0.54, -13.56, -7.35, -11.3, 0.78, -9.39, -3.76, -1.75, 9.39, -5.18,
      8.01, 4.42, -15.16, -8.3, -6.64, -10.14, 2.86, -3.93, -12.67, 0.16,
      -11.42, -14.22, -2.44, 6.6, 3.69, -6.1, 9.84, -17.15, -0.97, 1.77, -10.32,
      -3.67, 8.7, 5.01, -4.73, -7.33, -8.86, -2.51, -13.49, 4.2, -1.12, -6.06,
      -12.29, -15.43, -10.48, -16.71, 0.97, 2.51, 7.63, -3.8, -5.32, 5.64,
      -7.34, -1.17, -10.01, -14.72, -13.37, -8.87,
      # log_nu
      0.75, -0.23, 0.38, 0.15, -0.73, -0.21, 0.25, 1.75, -2.14, 0.14, -0.76,
      -0.05, -0.43, -1.42, 1.09, 0.08, 0.63, 1.5, 1.64, -1.76, -0.14, -1.82,
      -2.08, 0.25, 2.26, 1.94, -1.91, -2.2, -0.06, -0.46, 1.77, -0.61, 0.5,
      -1.14, -1.34, -1.87, -0.25, 0.82, -0.19, 1.44, 1.78, 1.18, 1.79, -1.18,
      -0.6, 1.07, 1.46, 1.21, 0.83, -1.9, 1.12, -1.49, 0.06, -0.69, 0.48, 0.71,
      -1.1, -1.97, 1, 1.95, -2.25, 0.6, -1.5, 0.39, 1.04, -0.33, 0.75, 1.98,
      0.7, 1.12, -1.66, -0.28, -0.35, -0.52, -0.25, -1.09, 1.54, 2.26, 0, -0.38,
      1.72, 1.77, -0.76, 1.52, 1.24, 1.74, 1.87, 1.47, 0.15, -0.72, 0.54, -0.86,
      0.12, -0.97, 1.55, 0.21, -0.11, 0.85, -1.93, 1.18, 0.42, -1.81, -0.61,
      2.3, 1.47, 1.07, 0.04, -0.18, -1.8, -1.59, 0.28, -0.44, -2.22, 1.27, 1.05,
      0.63, -1.22, -1.11, 0.78, 1.26, 1.25, 2.11, -1.42, -1.01, 0.09, -0.3,
      2.27, 0.95, -1.5, 0.41, 0.05, 0.75, 0.21, -2.29, -0.88, 1.34, 0.84, -1.09,
      -0.96, 0.63, -2.17, -0.57, 0.59, -1.23, 0.07, -0.07, -0.65, 1.1, -1.09,
      -1.05, -1.39, 1.84, 0.31, 1.42, 0.97, -0.97, 1.12, -1.39, -0.52, -1.94,
      0.35, -0.93, -1.75, 2.11, -0.44, -1.4, -2.29, 0.89, -0.47, -1.99, 0.72,
      -1.96, -1.52, 0.1, 0.35, -0.78, 2.11, 0.48, -1.77, -2.16, -0.97, -1.61,
      0.57, 1.75, 2.21, 0.77, 0.66, -0.82, -2.29, 1.87, 1.32, 0.83, -0.5, -0.51,
      -0.2, 0.34, 0.97, 1.72, -0.09, -0.31, 1.32, -0.92, 0.03, 0.12, 1.26, 0.25,
      -1.29, 0.84, -0.68, -0.23, -0.04, 0.48, 0.4, -1.77, 1.34, 1.54, -1.49,
      -1.61, -0.49, -0.32, 0.87, -0.8, -1.3, 1.32, -1.73, -1.59, 1.94, -2.24,
      -0.48, 0.49, 2.19, 1.62, -1.37, -0.71, -1.54, -0.89, 2.05, -1.53, -1.96,
      2.16, -0.79, 0.03, -1.31, -1.27, 2.05, -1.15, 1.68, -0.33, 0.98, 1.05
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

  # many points might have the same approximate RSS therefore we must choose the
  # most likely ones
  ord <- order(
    round(rss_tmp, digits = 4),
    apply(theta_tmp, 2, function(x) sqrt(crossprod(x)))
  )

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

# 5-parameter logistic fit
#
# Fit a 5-parameter logistic function to observed data with a Maximum
# Likelihood approach.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider the alternative
# parameterization `v = log(nu)` and `z = log(eta)`.
#
# @param object object of class `logistic5`.
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
fit.logistic5 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[c(3, 5)] <- exp(theta[c(3, 5)])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 5),
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 5,
    fitted.values = logistic5_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic5_fit", "logistic")

  result
}

# @rdname fit.logistic5
fit_constrained.logistic5 <- function(object) {
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
  theta[c(3, 5)] <- exp(theta[c(3, 5)])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = logistic5_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "delta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- c("logistic5_fit", "logistic")

  result
}

# 5-parameter logistic fit
#
# Evaluate the Fisher information matrix at the maximum likelihood estimate.
#
# @details
# Let `mu(x; theta)` be the 5-parameter logistic function. We assume that our
# observations `y` are independent and such that
# `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
# distribution `N(0, 1)`.
#
# The 5-by-5 (symmetric) Fisher information matrix is the expected value of
# the negative Hessian matrix of the log-likelihood function. We compute the
# observed Fisher information matrix because it has better finite sample
# properties.
#
# @param object object of class `logistic5`.
# @param theta numeric vector with the model parameters.
# @param sigma estimate of the standard deviation.
#
# @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic5 <- function(object, theta, sigma) {
  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]
  z <- fn(object, x, theta) - y

  gh <- logistic5_gradient_hessian(x, theta)

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
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - sum(object$w > 0)

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 5-parameter logistic fit
#
# Find the dose that produced the observed response.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# This function evaluates the inverse function of `f(x; theta)`, that is
# if `y = fn(x; theta)` then `x = inverse_fn(y; theta)`.
inverse_fn.logistic5_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  x <- (delta / (y - alpha))^nu
  x[!is.na(x) & (x > 1)] <- phi - log((x[!is.na(x) & (x > 1)] - 1) / nu) / eta

  x
}

# 5-parameter logistic fit
#
# Evaluate at a particular point the gradient of the inverse logistic function.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `g(x; theta) = 1 / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
# `f(x; theta) = alpha + delta g(x; theta)`
#
# where `theta = c(alpha, delta, eta, phi, nu)`, `eta > 0`, and `nu > 0`.
#
# This function evaluates the gradient of the inverse function.
inverse_fn_gradient.logistic5_fit <- function(object, y) {
  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  eta <- object$coefficients[3]
  nu <- object$coefficients[5]

  z <- delta / (y - alpha)
  s <- z^nu
  u <- nu / (s - 1)

  G <- matrix(1, nrow = length(y), ncol = 5)

  G[, 1] <- -z * s * u / (delta * eta)
  G[, 2] <- -s * u / (delta * eta)
  G[, 3] <- -log(u) / eta^2
  G[, 5] <- (1 - s * log(z) * u) / (eta * nu)

  G
}
