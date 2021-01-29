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
#' @importFrom stats cov median var
logistic4_init <- function(
  stats
) {
  k <- nrow(stats)
  delta <- mean(diff(stats[, 1]))

  rss <- logistic4_rss(stats)

  n <- stats[, 2]

  growth_rate <- seq(-2, -0.01, length.out = 15)

  logx_midpoint <- seq(
    stats[1, 1] - 0.5 * delta, stats[k, 1] + 0.5 * delta, length.out = 15
  )

  theta <- c(
    min(stats[, 3]),
    max(stats[, 3]),
    -1,
    stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1]
  )

  if (stats[k, 3] > stats[k, 1]) {
    theta[3] <- 1
  }

  best_rss <- rss(theta)

  for (xm in logx_midpoint) {
    for (gr in growth_rate) {
      f <- exp(-log1p(exp(-gr * (stats[, 1] - xm))))

      minimum_numer <- 0
      maximum_numer <- 0
      denom <- 0

      for (i1 in seq_len(k)) {
        x1 <- n[i1] * f[i1]
        x2 <- n[i1] * (f[i1] - 1)

        tmp_numer <- 0
        tmp_denom <- 0
        for (i2 in seq_len(k)) {
          x3 <- n[i2] * stats[i2, 3]
          x4 <- f[i1] - f[i2]

          tmp_numer <- tmp_numer + x3 * x4
          tmp_denom <- tmp_denom + n[i2] * x4
        }

        minimum_numer <- minimum_numer + x1 * tmp_numer
        maximum_numer <- maximum_numer + x2 * tmp_numer
        denom <- denom + x1 * tmp_denom
      }

      minimum <- minimum_numer / denom
      maximum <- maximum_numer / denom

      current_par <- c(minimum, maximum, gr, xm)
      current_rss <- rss(current_par)

      if (!is.nan(current_rss) && (current_rss < best_rss)) {
        theta <- current_par
        best_rss <- current_rss
      }
    }
  }

  theta
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
