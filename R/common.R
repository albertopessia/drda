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
#'   input vector `x`, second column are the total weights (sample size if no
#'   weights), third column are the (weighted) sample means, and fourth column
#'   are the uncorrected (weighted) sample variances.
suff_stats <- function(x, y, w) {
  unique_x <- sort(unique(x))
  k <- length(unique_x)

  stats <- matrix(NA_real_, nrow = k, ncol = 4)
  colnames(stats) <- c("x", "n", "m", "v")

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
      solve(fim)
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

#' Fit a function to observed data
#'
#' Use a Newton trust-region method to fit a function to observed data.
#'
#' @param object object of some model class.
#' @param constraint boolean matrix with constraint specifications.
#' @param known_param numeric vector of fixed known parameters.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{optimum}{maximum likelihood estimates of the model parameters.}
#'     \item{minimum}{(local) minimum of the residual sum of squares around the
#'      means.}
#'     \item{converged}{boolean. `TRUE` if the optimization algorithm converged,
#'       `FALSE` otherwise.}
#'     \item{iterations}{total number of iterations performed by the
#'       optimization algorithm}
#'   }
find_optimum <- function(object) {
  start <- if (!is.null(object$start)) {
    object$start
  } else {
    init(object)
  }

  rss_fn <- rss(object)
  rss_gh <- rss_gradient_hessian(object)

  update_fn <- function(theta) {
    mle_asy(object, theta)
  }

  ntrm(rss_fn, rss_gh, start, object$max_iter, update_fn)
}

#' @rdname find_optimum
find_optimum_constrained <- function(object, constraint, known_param) {
  start <- if (!is.null(object$start)) {
    # equality constraints have the priority over the provided starting values
    ifelse(is.na(known_param), object$start, known_param)
  } else {
    ifelse(is.na(known_param), init(object), known_param)
  }

  if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and Hessian
    idx <- which(!constraint[, 2])

    rss_fn <- rss_fixed(object, known_param)
    rss_gh <- rss_gradient_hessian_fixed(object, known_param)

    if (all(constraint[idx, 1])) {
      # we only have equality constraints, so after fixing the parameters what
      # remains is an unconstrained optimization
      ntrm(rss_fn, rss_gh, start[idx], object$max_iter)
    } else {
      ntrm_constrained(
        rss_fn, rss_gh, start[idx], object$max_iter, object$lower_bound[idx],
        object$upper_bound[idx]
      )
    }
  } else {
    rss_fn <- rss(object)
    rss_gh <- rss_gradient_hessian(object)

    ntrm_constrained(
      rss_fn, rss_gh, start, object$max_iter, object$lower_bound,
      object$upper_bound
    )
  }
}
