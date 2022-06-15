# Sufficient statistics
#
# Compute sufficient statistics for each unique value of the predictor
# variable.
#
# @param x numeric vector representing the fixed predictor variable.
# @param y numeric vector of observed values.
# @param w numeric vector of weights.
#
# @return Numeric matrix where the first column are the sorted unique values of
#   input vector `x`, second column are the total weights (sample size if all
#   weights are 1), third column are the (weighted) sample means, and fourth
#   column are the uncorrected (weighted) sample variances.
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

# Variance estimator
#
# Compute the corrected variance estimate associated with a Normal distribution.
#
# @param rss value of the residual sum of squares.
# @param df residual degrees of freedom.
#
# @return Corrected maximum likelihood estimate of the variance.
variance_normal <- function(rss, df) {
  rss / df
}

# Log-likelihood value
#
# Evaluate the maximum of the log-likelihood function associated with a Normal
# distribution.
#
# @param deviance deviance value, i.e. residual sum of squares.
# @param n total sample size.
# @param log_w sum of log-weights.
#
# @return Value of the log-likelihood function at the MLE.
loglik_normal <- function(deviance, n, log_w = 0) {
  -(n * (1 + log(2 * pi * deviance / n)) - log_w) / 2
}

# Approximate variance-covariance matrix
#
# Using the Fisher information matrix, compute an approximate
# variance-covariance matrix of the estimated model parameters.
#
# @param fim Fisher information matrix.
#
# @return Approximate variance-covariance matrix.
approx_vcov <- function(fim) {
  # total number of parameters (last one is always the standard deviation)
  p <- nrow(fim)

  vcov <- tryCatch(chol2inv(chol(fim)), error = function(e) NULL)

  if (is.null(vcov)) {
    # row/column associated to the standard deviation can create problems
    vcov <- tryCatch(chol2inv(chol(fim[-p, -p])), error = function(e) NULL)

    if (!is.null(vcov)) {
      # we succeeded
      vcov <- cbind(rbind(vcov, rep(NA_real_, p - 1)), rep(NA_real_, p))
    }
  }

  if (is.null(vcov)) {
    # sometimes `solve` succeeds where the previous one instead fails.
    # However, if the Cholesky decomposition failed, it is most likely a
    # sign that the solution is not admissible.
    # We will check later if all variances are non-negative.
    vcov <- tryCatch(solve(fim), error = function(e) NULL)

    if (is.null(vcov)) {
      vcov <- tryCatch(solve(fim[-p, -p]), error = function(e) NULL)

      if (!is.null(vcov)) {
        # we succeeded
        vcov <- cbind(rbind(vcov, rep(NA_real_, p - 1)), rep(NA_real_, p))
      }
    }
  }

  if (!is.null(vcov)) {
    d <- diag(vcov)

    # variances cannot be negative, but small values within tolerable numerical
    # error can be considered zero
    eps <- sqrt(.Machine$double.eps)
    d[d > -eps & d < 0] <- 0

    idx <- which(d < 0)

    diag(vcov) <- d

    if (length(idx) > 0) {
      # is the negative variance only the one associated with sigma?
      # In this case, we can be a bit more flexible in accepting the solution
      # because the variance-covariance matrix of the core parameters is ok
      if (length(idx) == 1 && idx[1] == p) {
        vcov[p, ] <- NA_real_
        vcov[, p] <- NA_real_
      } else {
        # some variances are negative, we cannot trust this approximation
        vcov <- matrix(NA_real_, nrow = p, ncol = p)
      }
    }
  } else {
    vcov <- matrix(NA_real_, nrow = p, ncol = p)
  }

  lab <- rownames(fim)
  rownames(vcov) <- lab
  colnames(vcov) <- lab

  vcov
}

# Algorithm initialization
#
# Use `nlminb` to find a good initial value.
#
# @param object object of some model class.
# @param start matrix of candidate starting points.
# @param niter maximum number of iterations
#
#' @importFrom stats nlminb
#'
#' @noRd
fit_nlminb <- function(object, start, max_iter) {
  # define the objective function
  f <- rss(object)

  # we give `nlminb` also the gradient and Hessian functions, however it wants
  # them as separate functions
  gh <- rss_gradient_hessian(object)

  g <- function(x) {
    gh(x)$G
  }

  h <- function(x) {
    gh(x)$H
  }

  fit_fn <- if (!object$constrained) {
    function(x, k) {
      y <- nlminb(
        start = x, objective = f, gradient = g, hessian = h,
        control = list(eval.max = k, iter.max = k)
      )

      list(
        par = mle_asy(object, y$par), niter = y$iterations
      )
    }
  } else {
    function(x, k) {
      y <- nlminb(
        start = x, objective = f, gradient = g, hessian = h,
        control = list(eval.max = k, iter.max = k),
        lower = object$lower_bound, upper = object$upper_bound
      )

      list(par = y$par, niter = y$iterations)
    }
  }

  best_par <- rep(NA_real_, nrow(start))
  best_rss <- Inf

  niter <- 0
  for (i in seq_len(ncol(start))) {
    tmp <- tryCatch(
      suppressWarnings(fit_fn(start[, i], max_iter)),
      error = function(e) NULL
    )

    if (!is.null(tmp)) {
      current_rss <- f(tmp$par)
      max_iter <- max(0, max_iter - tmp$niter)
      niter <- niter + tmp$niter

      if (!is.nan(current_rss) && (current_rss < best_rss)) {
        best_par <- tmp$par
        best_rss <- current_rss
      }
    }

    if (max_iter == 0) {
      break
    }
  }

  list(theta = best_par, rss = best_rss, niter = niter)
}

# Fit a function to observed data
#
# Use a Newton trust-region method to fit a function to observed data.
#
# @param object object of some model class.
# @param constraint boolean matrix with constraint specifications.
# @param known_param numeric vector of fixed known parameters.
#
# @return A list with the following components:
#   \describe{
#     \item{optimum}{maximum likelihood estimates of the model parameters.}
#     \item{minimum}{(local) minimum of the residual sum of squares around the
#      means.}
#     \item{converged}{boolean. `TRUE` if the optimization algorithm converged,
#       `FALSE` otherwise.}
#     \item{iterations}{total number of iterations performed by the
#       optimization algorithm}
#   }
find_optimum <- function(object) {
  start <- init(object)

  rss_fn <- rss(object)
  rss_gh <- rss_gradient_hessian(object)

  update_fn <- function(theta) {
    mle_asy(object, theta)
  }

  max_iter <- max(0, object$max_iter - start$niter)

  solution <- ntrm(rss_fn, rss_gh, start$theta, max_iter, update_fn)
  solution$iterations <- solution$iterations + start$niter

  solution
}

# @rdname find_optimum
find_optimum_constrained <- function(object, constraint, known_param) {
  start <- init(object)

  max_iter <- max(0, object$max_iter - start$niter)

  # equality constraints have the priority over the provided starting values
  theta <- ifelse(is.na(known_param), start$theta, known_param)

  solution <- if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and Hessian
    idx <- which(!constraint[, 2])

    rss_fn <- rss_fixed(object, known_param)
    rss_gh <- rss_gradient_hessian_fixed(object, known_param)

    if (all(constraint[idx, 1])) {
      # we only have equality constraints, so after fixing the parameters what
      # remains is an unconstrained optimization
      ntrm(rss_fn, rss_gh, theta[idx], max_iter)
    } else {
      ntrm_constrained(
        rss_fn, rss_gh, theta[idx], max_iter, object$lower_bound[idx],
        object$upper_bound[idx]
      )
    }
  } else {
    rss_fn <- rss(object)
    rss_gh <- rss_gradient_hessian(object)

    ntrm_constrained(
      rss_fn, rss_gh, theta, max_iter, object$lower_bound, object$upper_bound
    )
  }

  solution$iterations <- solution$iterations + start$niter

  solution
}
