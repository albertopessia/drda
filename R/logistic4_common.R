#' Sufficient statistics
#'
#' Compute sufficient statistics for each unique value of the predictor
#' variable.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#'
#' @return Numeric matrix where the first column `x` are the sorted unique
#'   values of input vector `x`, second column `n` are the corresponding number
#'   of observations, third column `m` are the sample means, and fourth column
#'   `v` are the (uncorrected) sample variances.
logistic4_suff_stat <- function(
  x,
  y
) {
  unique_x <- sort(unique(x))
  k <- length(unique_x)

  s <- by(y, x, function(z) {
    n <- length(z)
    m <- mean(z)
    v <- sum((z - m)^2) / n
    c(n, m, v)
  })

  result <- cbind(
    unique_x,
    matrix(unlist(s), nrow = k, ncol = 3, byrow = TRUE)
  )
  colnames(result) <- c("x", "n", "m", "v")

  result
}

#' Initialize vector of parameters
#'
#' Given the sufficient statistics, try to guess a good approximation to the
#' Maximum Likelihood estimator of the four parameters of the logistic function.
#'
#' @param stats numeric matrix of sufficient statistics.
#'
#' @return Numeric vector of length 4 with a (hopefully) good starting point.
#'
#' @importFrom stats cov var
logistic4_init <- function(
  stats
) {
  alpha <- min(stats[, 3])
  beta <- max(stats[, 3])

  # perform an approximate linearisation and solve for the remaining two
  # parameters
  a <- alpha - 1.0e-10
  b <- beta + 1.0e-10

  z <- (stats[, 3] - a) / (b - a)
  w <- log(z / (1 - z))

  # we now assume that `w` is approximately eta * (x - phi)
  # the following two formulas are the least square estimators, i.e.
  # the variables that minimize sum((w - eta * (x - phi))^2)
  eta <- cov(stats[, 1], w) / var(stats[, 1])
  phi <- mean(stats[, 1]) - mean(w) / eta

  c(alpha, beta, eta, phi)
}

#' Variance estimator
#'
#' Compute the corrected variance estimate associated to a 4-parameter logistic
#' function.
#'
#' @param rss value of the residual sum of squares, that is its minimum possible
#'   value because we are using maximum likelihood estimation.
#' @param n total sample size.
#' @param k number of estimated parameters, with a default of 4 in case of
#'   unconstrained optimization.
#'
#' @return Corrected maximum likelihood estimate of the variance.
logistic4_variance <- function(
  rss,
  n,
  k = 4
) {
  df <- n - k
  rss / df
}

#' Log-likelihood value
#'
#' Evaluate the maximum of the log-likelihood function.
#'
#' @param variance corrected maximum likelihood estimate of the variance.
#' @param n total sample size.
#' @param k number of estimated parameters, with a default of 4 in case of
#'   unconstrained optimization.
#' @param log_w sum of log-weights, with a default value of 0 in case of weights
#'   all equal to 1.
#'
#' @return Corrected maximum likelihood estimator of the variance.
logistic4_loglik <- function(
  variance,
  n,
  k = 4,
  log_w = 0
) {
  - (n * (1 + log(2 * pi * variance) + log((n - k) / n)) - log_w) / 2
}

#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function is defined in this package as
#'
#' `f(x; theta) = alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param param numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic4_function <- function(
  x,
  param
) {
  param[1] + (param[2] - param[1]) / (1 + exp(-param[3] * (x - param[4])))
}

#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient of the 4-parameter
#' logistic function.
#'
#' @details
#' The 4-parameter logistic function is defined in this package as
#'
#' `f(x; theta) = alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param param numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return Numeric matrix of dimension length(x)-by-4, where each row is the
#'   gradient of the logistic function at the corresponding element of `x`.
logistic4_gradient <- function(
  x,
  param
) {
  d <- param[2] - param[1]
  y <- x - param[4]
  w <- -param[3] * y

  log_denom <- log1p(exp(w))

  t1 <- exp(-log_denom)
  t2 <- exp(w - log_denom)

  gradient <- matrix(1, nrow = length(x), ncol = 4)
  gradient[, 1] <- t2
  gradient[, 2] <- t1
  gradient[, 3] <- y * d * t1 * t2
  gradient[, 4] <- -param[3] * d * t1 * t2

  gradient
}

#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the Hessian matrix of the
#' 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function is defined in this package as
#'
#' `f(x; theta) = alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param param numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return Array `H` of dimension length(x)-by-4-by-4, where `H[i, , ]` is the
#'   4-by-4 Hessian matrix at `x[i]`.
logistic4_hessian <- function(
  x,
  param
) {
  d <- param[2] - param[1]
  y <- x - param[4]
  w <- -param[3] * y

  log_denom <- log1p(exp(w))

  t1 <- exp(-log_denom)
  t2 <- exp(w - log_denom)
  t3 <- t1 * t2
  t4 <- 2 * t2 - 1

  hessian <- array(0, dim = c(length(x), 4, 4))

  hessian[, 3, 1] <- -y * t3
  hessian[, 4, 1] <- param[3] * t3

  hessian[, 3, 2] <- y * t3
  hessian[, 4, 2] <- -param[3] * t3

  hessian[, 1, 3] <- hessian[, 3, 1]
  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- d * y^2 * t3 * t4
  hessian[, 4, 3] <- -d * t3 * (1 - w * t4)

  hessian[, 1, 4] <- hessian[, 4, 1]
  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- d * param[3]^2 * t3 * t4

  hessian
}
