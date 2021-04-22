# @rdname logistic6_new
gompertz_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 4) {
      stop("'start' must be of length 4", call. = FALSE)
    }

    if (start[2] <= start[1]) {
      stop("parameter 'beta' cannot be smaller than 'alpha'", call. = FALSE)
    }

    if (start[3] == 0) {
      stop("parameter 'eta' cannot be initialized to zero", call. = FALSE)
    }
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
      rep(-Inf, 4)
    } else {
      if (length(lower_bound) != 4) {
        stop("'lower_bound' must be of length 4", call. = FALSE)
      }

      if (lower_bound[2] < lower_bound[1]) {
        # lower bound on alpha is a stronger constraint because beta > alpha
        lower_bound[2] <- lower_bound[1]
      }
    }

    if (is.null(upper_bound)) {
      rep(Inf, 4)
    } else {
      if (length(upper_bound) != 4) {
        stop("'upper_bound' must be of length 4", call. = FALSE)
      }

      if (upper_bound[1] > upper_bound[2]) {
        # upper bound on beta is a stronger constraint because alpha < beta
        upper_bound[1] <- upper_bound[2]
      }
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
#' `alpha + (beta - alpha) exp(-exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta > alpha` is the upper horizontal asymptote, `eta` is the
#' steepness of the curve or growth rate, and `phi` related to the function
#' value at `x = 0`.
#'
#' @param x numeric vector at which the Gompertz function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   Gompertz function.
#'
#' @export
gompertz_fn <- function(x, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + (beta - alpha) * exp(-exp(-eta * (x - phi)))
}

# Gompertz function
#
# Evaluate at a particular set of parameters the Gompertz function.
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
# @param x numeric vector at which the Gompertz function is to be evaluated.
# @param theta numeric vector with the four parameters in the form
#   `c(alpha, beta, eta, phi)`.
#
# @return Numeric vector of the same length of `x` with the values of the
#   Gompertz function.
fn.gompertz <- function(object, x, theta) {
  gompertz_fn(x, theta)
}

# @rdname fn.gompertz
fn.gompertz_fit <- function(object, x, theta) {
  gompertz_fn(x, theta)
}

# Gompertz function
#
# Evaluate at a particular set of parameters the gradient and Hessian of the
# Gompertz function.
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
# @param theta numeric vector with the four parameters in the form
#   `c(alpha, beta, eta, phi)`.
#
# @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.gompertz <- function(object, theta) {
  x <- object$stats[, 1]

  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  omega <- beta - alpha

  b <- exp(-eta * (x - phi))

  g <- exp(-b)

  q <- (x - phi) * b
  r <- -eta * b

  t <- q * g
  u <- r * g

  v <- expm1(-eta * (x - phi))

  gradient <- matrix(0, nrow = object$m, ncol = 4)
  hessian <- array(0, dim = c(object$m, 4, 4))

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- omega * t
  gradient[, 4] <- omega * u

  hessian[, 1, 1] <- 0
  hessian[, 2, 1] <- 0
  hessian[, 3, 1] <- -t
  hessian[, 4, 1] <- -u

  hessian[, 1, 2] <- 0
  hessian[, 2, 2] <- 0
  hessian[, 3, 2] <- t
  hessian[, 4, 2] <- u

  hessian[, 1, 3] <- hessian[, 3, 1]
  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- omega * (x - phi) * t * v
  hessian[, 4, 3] <- omega * (u / eta - eta * t * v)

  hessian[, 1, 4] <- hessian[, 4, 1]
  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- -omega * eta * u * v

  # When `b` is infinite, gradient and Hessian show NaNs
  if (any(is.nan(gradient))) {
    # these are the limits for b -> Inf
    gradient[, 1][is.nan(gradient[, 1])] <- 1
    gradient[, 2:4][is.nan(gradient[, 2:4])] <- 0
  }

  if (any(is.nan(hessian))) {
    hessian[is.nan(hessian)] <- 0
  }

  list(G = gradient, H = hessian)
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
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

# @rdname rss.gompertz
rss_fixed.gompertz <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

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
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

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
  m <- object$m

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- theta[3]
  phi <- theta[4]

  g <- exp(-exp(-eta * (x - phi)))

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

  linear_fit <- summary(lm(stats[, 3] ~ stats[, 1], weights = stats[, 2]))
  linear_coef <- linear_fit$coefficients

  min_value <- min(stats[, 3])
  max_value <- max(stats[, 3])

  # y = a + (b - a) * exp(-exp(-e * (x - p)))
  # w = (y - a) / (b - a)
  # z = log(-log(w)) = -e * (x - p) = e * p - e * x
  #
  # fit a linear model `z ~ u0 + u1 x` and set `eta = -u1` and `phi = -u0 / u1`
  # we add a small number to avoid the logarithm of zero
  zv <- (stats[, 3] - min_value + 1.0e-5) / (max_value - min_value + 2.0e-5)
  zv <- log(-log(zv))
  tmp <- lm(zv ~ stats[, 1])

  eta <- -tmp$coefficients[2]
  phi <- -tmp$coefficients[1] / eta

  theta <- mle_asy(object, c(min_value, max_value, eta, phi))

  best_rss <- rss_fn(theta)
  niter <- 1

  if (linear_coef[2, 4] > 0.2) {
    # we are in big problems as a flat horizontal line is likely the best model
    weighted_mean <- sum(object$w * object$y) / sum(object$w)

    theta <- c(
      0.9 * weighted_mean + 0.1 * min_value,
      0.9 * weighted_mean + 0.1 * max_value,
      if (linear_coef[2, 1] <= 0) -1.0e-3 else 1.0e-3,
      object$stats[m, 1] + 100
    )

    best_rss <- rss_fn(theta)
  }

  delta <- mean(diff(stats[, 1]))

  v1 <- 20L
  v2 <- 20L
  v <- v1 * v2
  eta_set <- seq(-5, 5, length.out = v1)
  phi_set <- seq(
    stats[1, 1] - delta, stats[m, 1] + delta, length.out = v2
  )

  theta_tmp <- matrix(nrow = 4, ncol = v)
  rss_tmp <- rep(10000, v)

  # we check extreme values in case of problematic data
  theta_eta_1 <- matrix(nrow = 4, ncol = v2)
  rss_eta_1 <- rep(10000, v2)

  theta_eta_2 <- matrix(nrow = 4, ncol = v2)
  rss_eta_2 <- rep(10000, v2)

  i <- 0
  j <- 0
  for (phi in phi_set) {
    j <- j + 1

    for (eta in eta_set) {
      i <- i + 1

      current_par <- mle_asy(object, c(theta[1], theta[2], eta, phi))
      current_rss <- rss_fn(current_par)
      theta_tmp[, i] <- current_par
      rss_tmp[i] <- current_rss
    }

    current_par <- mle_asy(object, c(theta[1], theta[2], -20, phi))
    current_rss <- rss_fn(current_par)
    theta_eta_1[, j] <- current_par
    rss_eta_1[j] <- current_rss

    current_par <- mle_asy(object, c(theta[1], theta[2], 20, phi))
    current_rss <- rss_fn(current_par)
    theta_eta_2[, j] <- current_par
    rss_eta_2[j] <- current_rss
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

  if (theta[2] < theta[1]) {
    # this is the dual solution, not the one we want
    # there is no easy formula to translate between the two, so we set the
    # solution to a sub-optimal one and let `ntrm` do the rest
    tmp <- theta[1]
    theta[1] <- theta[2]
    theta[2] <- tmp
    theta[3] <- -theta[3]
    theta[4] <- -theta[4]
  }

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

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "gompertz_fit"

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

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "gompertz_fit"

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
  w <- object$stats[, 2]
  d <- fn(object, object$stats[, 1], theta) - object$stats[, 3]

  gh <- gradient_hessian(object, theta)

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 4)
  G[, 1] <- w * d * gh$G[, 1]
  G[, 2] <- w * d * gh$G[, 2]
  G[, 3] <- w * d * gh$G[, 3]
  G[, 4] <- w * d * gh$G[, 4]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 4, 4))

  H[, , 1] <- w * (d * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (d * gh$H[, , 2] + gh$G[, 2] * gh$G)
  H[, , 3] <- w * (d * gh$H[, , 3] + gh$G[, 3] * gh$G)
  H[, , 4] <- w * (d * gh$H[, , 4] + gh$G[, 4] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - object$n

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

  alpha <- object$coefficients[1]
  beta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  omega <- beta - alpha

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  s <- 1 / f^2
  t <- q * s
  u <- r * s

  G <- matrix(0, nrow = m, ncol = 4)

  G[, 1] <- b / f
  G[, 2] <- 1 / f
  G[, 3] <- omega * t
  G[, 4] <- omega * u

  # When `b` is infinite, gradient shows NaNs
  if (any(is.nan(G))) {
    # these are the limits for b -> Inf
    G[, 1][is.nan(G[, 1])] <- 1
    G[, 2:4][is.nan(G[, 2:4])] <- 0
  }

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
  beta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]

  I <- 0
  xlim_new <- xlim

  if (alpha < ylim[1]) {
    tmp <- phi - log(log((beta - alpha) / (ylim[1] - alpha))) / eta

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
    tmp <- phi - log(log((beta - alpha) / (ylim[2] - alpha))) / eta

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
naac.gompertz_fit <- function(object, xlim = c(-10, 10), ylim = c(0, 1)) {
  1 - nauc(object, xlim, ylim)
}
