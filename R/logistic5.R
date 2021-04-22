# @rdname logistic6_new
logistic5_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 5) {
      stop("'start' must be of length 5", call. = FALSE)
    }

    if (start[2] <= start[1]) {
      stop("parameter 'beta' cannot be smaller than 'alpha'", call. = FALSE)
    }

    if (start[3] == 0) {
      stop("parameter 'eta' cannot be initialized to zero", call. = FALSE)
    }

    if (start[5] <= 0) {
      stop("parameter 'nu' cannot be negative nor zero", call. = FALSE)
    }

    start[5] <- log(start[5])
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
      rep(-Inf, 5)
    } else {
      if (length(lower_bound) != 5) {
        stop("'lower_bound' must be of length 5", call. = FALSE)
      }

      if (!is.infinite(lower_bound[2]) && (lower_bound[2] < lower_bound[1])) {
        # lower bound on alpha is a stronger constraint because beta > alpha
        lower_bound[2] <- lower_bound[1]
      }

      lower_bound[5] <- if (lower_bound[5] > 0) {
        log(lower_bound[5])
      } else {
        -Inf
      }
    }

    if (is.null(upper_bound)) {
      rep(Inf, 5)
    } else {
      if (length(upper_bound) != 5) {
        stop("'upper_bound' must be of length 5", call. = FALSE)
      }

      if (!is.infinite(upper_bound[1]) && (upper_bound[1] > upper_bound[2])) {
        # upper bound on beta is a stronger constraint because alpha < beta
        upper_bound[1] <- upper_bound[2]
      }

      if (upper_bound[5] <= 0) {
        stop("'upper_bound[5]' cannot be negative nor zero.", call. = FALSE)
      }

      upper_bound[5] <- log(upper_bound[5])
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
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, beta, eta, phi, nu)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic5_fn <- function(x, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)
}

# 5-parameter logistic function
#
# Evaluate at a particular set of parameters the 5-parameter logistic function.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# @param object object of class `logistic5`.
# @param x numeric vector at which the logistic function is to be evaluated.
# @param theta numeric vector with the five parameters in the form
#   `c(alpha, beta, eta, phi, nu)`.
#
# @return Numeric vector with the values of the logistic function.
fn.logistic5 <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

# @rdname fn.logistic5
fn.logistic5_fit <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

# 5-parameter logistic function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# 5-parameter logistic function.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# @param object object of class `logistic5`.
# @param theta numeric vector with the five parameters in the form
#   `c(alpha, beta, eta, phi, nu)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic5 <- function(object, theta) {
  x <- object$stats[, 1]

  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  omega <- beta - alpha

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s
  v <- u / eta + g * log(f) / nu

  gradient <- matrix(0, nrow = length(x), ncol = 5)
  hessian <- array(0, dim = c(length(x), 5, 5))

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- omega * t
  gradient[, 4] <- omega * u
  gradient[, 5] <- omega * v

  hessian[, 1, 1] <- 0
  hessian[, 2, 1] <- 0
  hessian[, 3, 1] <- -t
  hessian[, 4, 1] <- -u
  hessian[, 5, 1] <- -v

  hessian[, 1, 2] <- 0
  hessian[, 2, 2] <- 0
  hessian[, 3, 2] <- t
  hessian[, 4, 2] <- u
  hessian[, 5, 2] <- v

  hessian[, 1, 3] <- hessian[, 3, 1]
  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- omega * q * t * ((1 + nu) / f - 1 / b)
  hessian[, 4, 3] <- omega * (1 / eta + (1 + nu - f / b) * t / g) * u
  hessian[, 5, 3] <- omega * (nu * u / eta + v) * t / g

  hessian[, 1, 4] <- hessian[, 4, 1]
  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- omega * ((1 + nu) / f - 1 / b) * r * u
  hessian[, 5, 4] <- omega * (nu * u / eta + v) * u / g

  hessian[, 1, 5] <- hessian[, 5, 1]
  hessian[, 2, 5] <- hessian[, 5, 2]
  hessian[, 3, 5] <- hessian[, 5, 3]
  hessian[, 4, 5] <- hessian[, 5, 4]
  hessian[, 5, 5] <- omega * (nu * (u / eta)^2 + v * (v - g)) / g

  # When `b` is infinite, gradient and Hessian show NaNs
  # these are the limits for b -> Inf
  if (any(is.nan(gradient))) {
    gradient[, 1][is.nan(gradient[, 1])] <- 1
    gradient[, 2:5][is.nan(gradient[, 2:5])] <- 0
  }

  if (any(is.nan(hessian))) {
    hessian[is.nan(hessian)] <- 0
  }

  list(G = gradient, H = hessian)
}

# Residual sum of squares
#
# Evaluate the residual sum of squares (RSS) against the mean of a
# 5-parameter logistic model.
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider instead the equivalent
# function `f(x; theta')`
#
# `alpha + (beta - alpha) / (1 + exp(v - eta * (x - phi)))^(-exp(-v))`
#
# @param object object of class `logistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the RSS associated to a
#   particular parameter choice `theta`.
rss.logistic5 <- function(object) {
  function(theta) {
    theta[5] <- exp(theta[5])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.logistic5
rss_fixed.logistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[5] <- exp(theta[5])

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
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider instead the equivalent
# function `f(x; theta')`
#
# `alpha + (beta - alpha) / (1 + exp(v - eta * (x - phi)))^(-exp(-v))`
#
# @param object object of class `logistic5`.
# @param known_param numeric vector with the known fixed values of the model
#   parameters, if any.
#
# @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic5 <- function(object) {
  function(theta) {
    theta[5] <- exp(theta[5])

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
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[5] <- exp(theta[5])

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
  m <- object$m

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- theta[3]
  phi <- theta[4]
  log_nu <- theta[5]

  g <- 1 / (1 + exp(log_nu - eta * (x - phi)))^exp(-log_nu)

  t1 <- 0
  t2 <- 0
  t3 <- 0
  t4 <- 0
  t5 <- 0

  for (i in seq_len(m)) {
    t1 <- t1 + w[i] * g[i] * (g[i] - 1)
    t2 <- t2 + w[i] * (g[i] - 1)^2
    t3 <- t3 + w[i] * g[i]^2
    t4 <- t4 + w[i] * y[i] * (g[i] - 1)
    t5 <- t5 + w[i] * y[i] * g[i]
  }

  denom <- t1^2 - t2 * t3

  if (denom != 0) {
    alpha <- -(t1 * t5 - t3 * t4) / denom
    beta <- (t1 * t4 - t2 * t5) / denom

    if (beta > alpha) {
      theta[1] <- alpha
      theta[2] <- beta
    }
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
init.logistic5 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  linear_fit <- summary(lm(stats[, 3] ~ stats[, 1], weights = stats[, 2]))
  linear_coef <- linear_fit$coefficients

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  # we initialize nu = 1, so that we start with a 4-parameter logistic function
  #
  # y = a + (b - a) / (1 + exp(-e * (x - p)))
  # w = (y - a) / (b - a)
  # z = log(w / (1 - w)) = e * (x - p) = - e * p + e * x
  #
  # fit a linear model `z ~ u0 + u1 x` and set `eta = u1` and `phi = -u0 / u1`
  # we add a small number to avoid the logarithm of zero
  zv <- (stats[, 3] - min_value + 1.0e-5) / (max_value - min_value + 2.0e-5)
  zv <- log(zv) - log1p(-zv)
  tmp <- lm(zv ~ stats[, 1])

  eta <- tmp$coefficients[2]
  phi <- -tmp$coefficients[1] / eta

  theta <- mle_asy(object, c(min_value, max_value, eta, phi, 0))

  best_rss <- rss_fn(theta)
  niter <- 1

  if (linear_coef[2, 4] > 0.2) {
    # we are in big problems as a flat horizontal line is likely the best model
    weighted_mean <- sum(object$w * object$y) / sum(object$w)

    theta <- c(
      0.9 * weighted_mean + 0.1 * min_value,
      0.9 * weighted_mean + 0.1 * max_value,
      if (linear_coef[2, 1] <= 0) -1.0e-3 else 1.0e-3,
      object$stats[m, 1] + 100,
      0
    )

    best_rss <- rss_fn(theta)
  }

  delta <- mean(diff(stats[, 1]))

  v1 <- 20L
  v2 <- 20L
  v3 <- 5L
  v <- v1 * v2 * v3
  eta_set <- seq(-5, 5, length.out = v1)
  phi_set <- seq(
    stats[1, 1] - delta, stats[m, 1] + delta, length.out = v2
  )
  log_nu_set <- seq(-1, 0.5, length.out = v3)

  theta_tmp <- matrix(nrow = 5, ncol = v)
  rss_tmp <- rep(10000, v)

  # we check extreme values in case of problematic data
  theta_eta_1 <- matrix(nrow = 5, ncol = v2 * v3)
  rss_eta_1 <- rep(10000, v2 * v3)

  theta_eta_2 <- matrix(nrow = 5, ncol = v2 * v3)
  rss_eta_2 <- rep(10000, v2 * v3)

  i <- 0
  j <- 0
  for (phi in phi_set) {
    for (log_nu in log_nu_set) {
      j <- j + 1

      for (eta in eta_set) {
        i <- i + 1

        current_par <- mle_asy(object, c(theta[1], theta[2], eta, phi, log_nu))
        current_rss <- rss_fn(current_par)
        theta_tmp[, i] <- current_par
        rss_tmp[i] <- current_rss
      }

      current_par <- mle_asy(object, c(theta[1], theta[2], -20, phi, log_nu))
      current_rss <- rss_fn(current_par)
      theta_eta_1[, j] <- current_par
      rss_eta_1[j] <- current_rss

      current_par <- mle_asy(object, c(theta[1], theta[2], 20, phi, log_nu))
      current_rss <- rss_fn(current_par)
      theta_eta_2[, j] <- current_par
      rss_eta_2[j] <- current_rss
    }
  }

  ord <- order(rss_tmp)

  theta_1 <- theta_tmp[, ord[1]]
  theta_2 <- theta_tmp[, ord[round(v / 3)]]
  theta_3 <- theta_tmp[, ord[round(2 * v / 3)]]
  theta_4 <- theta_eta_1[, order(rss_eta_1)[1]]
  theta_5 <- theta_eta_2[, order(rss_eta_2)[1]]

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
    theta_4 <- pmax(
      pmin(theta_4, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
    theta_5 <- pmax(
      pmin(theta_5, object$upper_bound, na.rm = TRUE),
      object$lower_bound, na.rm = TRUE
    )
  }

  start <- cbind(theta, theta_1, theta_2, theta_3, theta_4, theta_5)

  tmp <- fit_nlminb(object, rss_fn, start)

  if (!is.infinite(tmp$rss) && (tmp$rss < best_rss)) {
    theta <- tmp$theta
    best_rss <- tmp$rss
    niter <- niter + tmp$niter
  }

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
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# In our optimization algorithm, however, we consider instead the equivalent
# function `f(x; theta')`
#
# `alpha + (beta - alpha) / (1 + exp(v - eta * (x - phi)))^(-exp(-v))`
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
  theta[5] <- exp(theta[5])

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

  param_names <- c("alpha", "beta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic5_fit"

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
  theta[5] <- exp(theta[5])

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

  param_names <- c("alpha", "beta", "eta", "phi", "nu")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic5_fit"

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

  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  omega <- beta - alpha

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s
  v <- u / eta + g * log(f) / nu

  d <- fn(object, x, theta) - y

  gradient <- matrix(0, nrow = object$m, ncol = 5)

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- omega * t
  gradient[, 4] <- omega * u
  gradient[, 5] <- omega * v / nu

  # When `b` is infinite, gradient shows NaNs
  if (any(is.nan(gradient))) {
    # these are the limits for b -> Inf
    gradient[, 1][is.nan(gradient[, 1])] <- 1
    gradient[, 2:5][is.nan(gradient[, 2:5])] <- 0
  }

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 5)
  G[, 1] <- w * d * gradient[, 1]
  G[, 2] <- w * d * gradient[, 2]
  G[, 3] <- w * d * gradient[, 3]
  G[, 4] <- w * d * gradient[, 4]
  G[, 5] <- w * d * gradient[, 5]

  G <- apply(G, 2, sum)

  hessian <- array(0, dim = c(object$m, 5, 5))

  hessian[, 1, 1] <- 0
  hessian[, 2, 1] <- 0
  hessian[, 3, 1] <- -t
  hessian[, 4, 1] <- -u
  hessian[, 5, 1] <- -v / nu

  hessian[, 1, 2] <- 0
  hessian[, 2, 2] <- 0
  hessian[, 3, 2] <- t
  hessian[, 4, 2] <- u
  hessian[, 5, 2] <- v / nu

  hessian[, 1, 3] <- hessian[, 3, 1]
  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- omega * q * t * ((1 + nu) / f - 1 / b)
  hessian[, 4, 3] <- omega * (1 / eta + (1 + nu - f / b) * t / g) * u
  hessian[, 5, 3] <- omega * (nu * u / eta + v) * t / (nu * g)

  hessian[, 1, 4] <- hessian[, 4, 1]
  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- omega * ((1 + nu) / f - 1 / b) * r * u
  hessian[, 5, 4] <- omega * (nu * u / eta + v) * u / (nu * g)

  hessian[, 1, 5] <- hessian[, 5, 1]
  hessian[, 2, 5] <- hessian[, 5, 2]
  hessian[, 3, 5] <- hessian[, 5, 3]
  hessian[, 4, 5] <- hessian[, 5, 4]
  hessian[, 5, 5] <- omega * (nu * (u / eta)^2 + v * (v - 2 * g)) / (nu^2 * g)

  # When `b` is infinite, gradient shows NaNs
  if (any(is.nan(hessian))) {
    # these are the limits for b -> Inf
    hessian[is.nan(hessian)] <- 0
  }

  H <- array(0, dim = c(object$m, 5, 5))

  H[, , 1] <- w * (d * hessian[, , 1] + gradient[, 1] * gradient)
  H[, , 2] <- w * (d * hessian[, , 2] + gradient[, 2] * gradient)
  H[, , 3] <- w * (d * hessian[, , 3] + gradient[, 3] * gradient)
  H[, , 4] <- w * (d * hessian[, , 4] + gradient[, 4] * gradient)
  H[, , 5] <- w * (d * hessian[, , 5] + gradient[, 5] * gradient)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - object$n

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

# 5-parameter logistic fit
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object of class `logistic5_fit`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.logistic5_fit <- function(object, x) {
  m <- length(x)

  V <- object$vcov[1:5, 1:5]

  if (any(is.na(V))) {
    return(rep(NA_real_, m))
  }

  alpha <- object$coefficients[1]
  beta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  omega <- beta - alpha

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)

  q <- (x - phi) * b
  r <- -eta * b

  s <- g / f
  t <- q * s
  u <- r * s
  v <- u / eta + g * log(f) / nu

  G <- matrix(0, nrow = m, ncol = 5)
  G[, 1] <- 1 - g
  G[, 2] <- g
  G[, 3] <- omega * t
  G[, 4] <- omega * u
  G[, 5] <- omega * v / nu

  # When `b` is infinite, gradient shows NaNs
  if (any(is.nan(G))) {
    # these are the limits for b -> Inf
    G[, 1][is.nan(G[, 1])] <- 1
    G[, 2:5][is.nan(G[, 2:5])] <- 0
  }

  variance <- rep(NA_real_, m)

  for (i in seq_len(m)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

# 5-parameter logistic fit
#
# Evaluate the normalized area under the curve (AUC) and area above the curve
# (AAC).
#
# @details
# The 5-parameter logistic function `f(x; theta)` is defined here as
#
# `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#
# where `theta = c(alpha, beta, eta, phi, nu)`, `beta > alpha`, and `nu > 0`.
#
# The area under the curve (AUC) is simply the integral of `f(x; theta)` with
# respect to `x`.
#
#' @importFrom stats integrate
#' @export
nauc.logistic5_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
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
  beta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  I <- 0
  xlim_new <- xlim

  if (alpha < ylim[1]) {
    tmp <- phi - log((((beta - alpha) / (ylim[1] - alpha))^nu - 1) / nu) / eta

    # if the curve is decreasing we change the upper bound of integration,
    # otherwise the lower bound
    if (eta < 0) {
      if (tmp < xlim[2]) {
        xlim_new[2] <- tmp
      }
    } else {
      if (tmp > xlim[1]) {
        xlim_new[1] <- tmp
      }
    }
  }

  if (beta > ylim[2]) {
    tmp <- phi - log((((beta - alpha) / (ylim[2] - alpha))^nu - 1) / nu) / eta

    # if the curve is decreasing we change the lower bound of integration,
    # otherwise the upper bound
    # in any case, we must now consider the area of the rectangle
    if (eta < 0) {
      if (tmp > xlim[1]) {
        I <- I + (tmp - xlim[1]) * (ylim[2] - ylim[1])
        xlim_new[1] <- tmp
      }
    } else {
      if (tmp < xlim[2]) {
        I <- I + (xlim[2] - tmp) * (ylim[2] - ylim[1])
        xlim_new[2] <- tmp
      }
    }
  }

  f <- if (ylim[1] == 0) {
    function(x) {
      fn(object, x, object$coefficients)
    }
  } else {
    function(x) {
      fn(object, x, object$coefficients) - ylim[1]
    }
  }

  I <- I + integrate(f, lower = xlim_new[1], upper = xlim_new[2])$value

  nauc <- I / ((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))
  names(nauc) <- NULL

  nauc
}

#' @export
naac.logistic5_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  1 - nauc(object, xlim, ylim)
}
