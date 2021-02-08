#' Sufficient statistics
#'
#' Compute sufficient statistics for each unique value of the predictor
#' variable.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of weights.
#'
#' @return Numeric matrix where the first column are the sorted unique values of
#'   input vector `x`, second column are the corresponding number of
#'   observations (or total weights), third column are the (weighted) sample
#'   means, and fourth column are the uncorrected (weighted) sample variances.
suff_stats <- function(x, y) {
  unique_x <- sort(unique(x))
  k <- length(unique_x)

  stats <- matrix(NA_real_, nrow = k, ncol = 4)
  colnames(stats) <- c("x", "n", "m", "v")

  for (i in seq_len(k)) {
    idx <- x == unique_x[i]

    z <- y[idx]

    n <- sum(idx)
    m <- mean(z)
    v <- sum((z - m)^2) / n

    stats[i, 1] <- unique_x[i]
    stats[i, 2] <- n
    stats[i, 3] <- m
    stats[i, 4] <- v
  }

  stats
}

#' @rdname suff_stats
suff_stats_weighted <- function(x, y, w) {
  unique_x <- sort(unique(x))
  k <- length(unique_x)

  stats <- matrix(NA_real_, nrow = k, ncol = 4)
  colnames(stats) <- c("x", "w", "m", "v")

  for (i in seq_len(k)) {
    idx <- x == unique_x[i]

    z <- y[idx]
    q <- w[idx]

    t <- sum(q)
    m <- sum(q * z) / t
    v <- sum(q * (z - m)^2) / t

    stats[i, 1] <- unique_x[i]
    stats[i, 2] <- t
    stats[i, 3] <- m
    stats[i, 4] <- v
  }

  stats
}

#' Variance estimator
#'
#' Compute the corrected variance estimate associated with a Normal
#' distribution.
#'
#' @param rss value of the residual sum of squares.
#' @param df residual degrees of freedom.
#'
#' @return Corrected maximum likelihood estimate of the variance.
variance_normal <- function(rss, df) {
  rss / df
}

#' Log-likelihood value
#'
#' Evaluate the maximum of the log-likelihood function associated with a Normal
#' distribution.
#'
#' @param deviance deviance value, i.e. residual sum of squares.
#' @param n total sample size.
#' @param log_w sum of log-weights.
#'
#' @return Value of the log-likelihood function at the MLE.
loglik_normal <- function(deviance, n, log_w = 0) {
  - (n * (1 + log(2 * pi * deviance / n)) - log_w) / 2
}

#' Approximate variance-covariance matrix
#'
#' Using the Fisher information matrix, compute an approximate
#' variance-covariance matrix of the estimated model parameters.
#'
#' @param fim Fisher information matrix.
#'
#' @return Approximate variance-covariance matrix.
approx_vcov <- function(fim) {
  vcov <- tryCatch({
      chol2inv(chol(fim))
    },
    error = function(e) {
      matrix(NA, nrow = nrow(fim), ncol = nrow(fim))
    }
  )

  lab <- rownames(fim)
  rownames(vcov) <- lab
  colnames(vcov) <- lab

  vcov
}
