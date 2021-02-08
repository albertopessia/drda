#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic4_function <- function(x, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))
}

#' 4-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 4-parameter logistic function.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(-eta * (x - phi)))`
#'
#' where `omega` is the log-width of the curve. Note that it is simply
#' `beta = alpha + exp(omega)`.
#'
#' NOTE: This function expects the input vector
#' `theta = c(alpha, beta, eta, phi)` but computes the gradient and Hessian of
#' the alternative parametrization `theta' = (alpha, omega, eta, phi)`.
#' This duality is unfortunately a trade-off between what users expect (first
#' version is easier to interpret) and the optimization stability of the second
#' version.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(alpha, beta, eta, phi)`.
#'
#' @return List of two elements. Element `G` is a numeric matrix of dimension
#'   length(x)-by-4, where each row is the gradient of the logistic function at
#'   the corresponding element of `x`. Element `H` is an array of dimension
#'   length(x)-by-4-by-4, where `H[k, , ]` is the 4-by-4 Hessian matrix at
#'   `x[k]`.
logistic4_gradient_hessian <- function(x, theta) {
  eta <- theta[3]
  phi <- theta[4]

  b <- exp(-eta * (x - phi))

  f <- 1 + b
  h <- (theta[2] - theta[1]) / f

  q <- (x - phi) * b
  r <- -eta * b

  s <- h / f
  t <- q / f
  u <- r / f

  gradient <- matrix(1, nrow = length(x), ncol = 4)
  hessian <- array(0, dim = c(length(x), 4, 4))

  gradient[, 2] <- h
  gradient[, 3] <- q * s
  gradient[, 4] <- r * s

  hessian[, 2, 2] <- gradient[, 2]
  hessian[, 3, 2] <- gradient[, 3]
  hessian[, 4, 2] <- gradient[, 4]

  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- (2 * t - (x - phi)) * gradient[, 3]
  hessian[, 4, 3] <- (2 * t - (x - phi) + 1 / eta) * gradient[, 4]

  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- (2 * u + eta) * gradient[, 4]

  list(G = gradient, H = hessian)
}

#' Residual sum of squares
#'
#' Evaluate the residual sum of squares (RSS) against the mean of a
#' 4-parameter logistic model.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(-eta * (x - phi)))`
#'
#' where `omega` is the log-width of the curve. Note that it is simply
#' `beta = alpha + exp(omega)`.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
logistic4_rss <- function(stats) {
  function(z) {
    theta <- z
    theta[2] <- theta[1] + exp(theta[2])

    mu <- logistic4_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' @rdname logistic4_rss
logistic4_rss_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])

    mu <- logistic4_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient and Hessian of the residual sum of squares (RSS)
#' against the mean of a 4-parameter logistic model.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(-eta * (x - phi)))`
#'
#' where `omega` is the log-width of the curve. Note that it is simply
#' `beta = alpha + exp(omega)`.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
logistic4_rss_gradient_hessian <- function(stats) {
  function(z) {
    theta <- z
    theta[2] <- theta[1] + exp(theta[2])

    mu <- logistic4_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic4_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 4, 4))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- stats[, 2] * (r * H[, , 4] + G[, 4] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

#' @rdname logistic4_rss_gradient_hessian
logistic4_rss_gradient_hessian_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])

    mu <- logistic4_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic4_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 4, 4))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- stats[, 2] * (r * H[, , 4] + G[, 4] * G)

    list(
      G = apply(gradient[, idx, drop = FALSE], 2, sum),
      H = apply(hessian[, idx, idx, drop = FALSE], 2:3, sum)
    )
  }
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
#' @importFrom stats median
logistic4_init <- function(stats) {
  k <- nrow(stats)
  delta <- mean(diff(stats[, 1]))

  rss <- logistic4_rss(stats)

  n <- stats[, 2]

  eta_set <- seq(-2, -0.01, length.out = 15)
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[k, 1] + 0.5 * delta, length.out = 15
  )

  theta <- c(
    min(stats[, 3]),
    max(stats[, 3]),
    -1,
    stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1]
  )

  if (stats[k, 3] > stats[1, 3]) {
    theta[3] <- 1
  }

  best_rss <- rss(theta)

  for (phi in phi_set) {
    for (eta in eta_set) {
      f <- exp(-log1p(exp(-eta * (stats[, 1] - phi))))

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

      current_par <- c(minimum, log(maximum - minimum), eta, phi)
      current_rss <- rss(current_par)

      if (!is.nan(current_rss) && (current_rss < best_rss)) {
        theta <- current_par
        best_rss <- current_rss
      }
    }
  }

  names(theta) <- NULL

  theta
}

#' 4-parameter logistic fit
#'
#' Evaluate the Fisher information matrix at the maximum likelihood estimate.
#'
#' @details
#' Let `mu(x; theta)` be the 4-parameter logistic function. We assume that our
#' observations `y` are independent and such that
#' `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
#' distribution `N(0, 1)`.
#'
#' The 4-by-4 (symmetric) Fisher information matrix is the expected value of
#' the negative Hessian matrix of the log-likelihood function.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param n effective sample size.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
logistic4_fisher_info_normal <- function(stats, n, theta, sigma) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]

  b <- exp(-eta * (stats[, 1] - phi))

  f <- 1 + b
  g <- 1 / f
  h <- (beta - alpha) * g

  q <- (stats[, 1] - phi) * b
  r <- -eta * b

  t <- q / f
  u <- r / f

  gradient <- matrix(1, nrow = nrow(stats), ncol = 4)

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- t * h
  gradient[, 4] <- u * h

  tmp <- array(0, dim = c(nrow(stats), 4, 4))
  tmp[, , 1] <- stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- stats[, 2] * gradient[, 2] * gradient
  tmp[, , 3] <- stats[, 2] * gradient[, 3] * gradient
  tmp[, , 4] <- stats[, 2] * gradient[, 4] * gradient

  fim <- matrix(0, nrow = 5, ncol = 5)
  fim[1:4, 1:4] <- apply(tmp, 2:3, sum)
  fim[5, 5] <- n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' Fit a 4-parameter logistic function
#'
#' Use a Newton trust-region method to fit a 4-parameter logistic function to
#' observed data.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param constraint boolean matrix with a representation of the constraints.
#' @param lower_bound numeric vector of length 4 with the minimum admissible
#'   values.
#' @param upper_bound numeric vector of length 4 with the maximum admissible
#'   values.
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
logistic4_ntrm <- function(stats, start, max_iter) {
  init <- if (!is.null(start)) {
    start
  } else {
    logistic4_init(stats)
  }

  fn <- logistic4_rss(stats)
  gh <- logistic4_rss_gradient_hessian(stats)

  ntrm(fn, gh, init, max_iter)
}

#' @rdname logistic4_ntrm
logistic4_ntrm_constrained <- function(
  stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
) {
  init <- if (!is.null(start)) {
    # equality constraints have the priority over the provided starting values
    ifelse(is.na(known_param), start, known_param)
  } else {
    ifelse(is.na(known_param), logistic4_init(stats), known_param)
  }

  if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and Hessian
    idx <- which(!constraint[, 2])

    fn <- logistic4_rss_fixed(stats, known_param)
    gh <- logistic4_rss_gradient_hessian_fixed(stats, known_param)

    if (all(constraint[idx, 1])) {
      # we only have equality constraints, so after fixing the parameters what
      # remains is an unconstrained optimization
      ntrm(fn, gh, init[idx], max_iter)
    } else {
      ntrm_constrained(
        fn, gh, init[idx], max_iter, lower_bound[idx], upper_bound[idx]
      )
    }
  } else {
    fn <- logistic4_rss(stats)
    gh <- logistic4_rss_gradient_hessian(stats)

    ntrm_constrained(
      fn, gh, init, max_iter, lower_bound, upper_bound
    )
  }
}

#' 4-parameter logistic fit
#'
#' Fit a 4-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 4-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of optional weights.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param lower_bound numeric vector of length 4 with the minimum admissible
#'   value of `alpha`, `beta`, `eta`, and `phi` respectively. Values can be
#'   equal to `-Inf`.
#' @param upper_bound numeric vector of length 4 with the maximum admissible
#'   value of `alpha`, `beta`, `eta`, and `phi` respectively. Values can be
#'   equal to `Inf`.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{converged}{boolean. `TRUE` if the optimization algorithm converged,
#'       `FALSE` otherwise.}
#'     \item{iterations}{total number of iterations performed by the
#'       optimization algorithm}
#'     \item{constrained}{boolean. `TRUE` if optimization was constrained,
#'       `FALSE` otherwise.}
#'     \item{estimated}{boolean vector indicating which parameters were
#'       estimated from the data.}
#'     \item{coefficients}{maximum likelihood estimates of the model
#'       parameters.}
#'     \item{rss}{minimum value found of the residual sum of squares.}
#'     \item{df.residual}{residual degrees of freedom.}
#'     \item{fitted.values}{fitted mean values.}
#'     \item{residuals}{residuals, that is response minus fitted values.}
#'     \item{weights}{vector of weights used for the fit.}
#'   }
logistic4_fit <- function(x, y, start, max_iter) {
  stats <- suff_stats(x, y)

  solution <- logistic4_ntrm(stats, start, max_iter)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 4),
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 4,
    fitted.values = logistic4_function(x, theta)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic4_fit
logistic4_fit_constrained <- function(
  x, y, start, max_iter, lower_bound, upper_bound
) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 4, 3)

  for (i in seq_len(4)) {
    lb_is_inf <- is.infinite(lower_bound[i])
    ub_is_inf <- is.infinite(upper_bound[i])

    if (lower_bound[i] == upper_bound[i]) {
      constraint[i, 2] <- TRUE
    } else if (!lb_is_inf || !ub_is_inf) {
      constraint[i, 3] <- TRUE
    } else {
      constraint[i, 1] <- TRUE
    }
  }

  known_param <- ifelse(constraint[, 2], lower_bound, NA_real_)

  stats <- suff_stats(x, y)

  solution <- logistic4_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  # bring the parameters back to their natural scale
  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic4_function(x, theta)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic4_fit
logistic4_fit_weighted <- function(x, y, w, start, max_iter) {
  w_zero <- w == 0
  w_positive <- !w_zero

  if (any(w_zero)) {
    w <- w[w_positive]
    x <- x[w_positive]
    y <- y[w_positive]

    if (length(y) == 0) {
      # all cases have weight zero
      return(
        list(
          converged = FALSE,
          iterations = 0L,
          constrained = FALSE,
          estimated = rep(NA_real_, 4),
          coefficients = rep(NA_real_, 4),
          rss = NA_real_,
          df.residual = 0L,
          fitted.values = numeric(0),
          residuals = numeric(0),
          weights = numeric(0)
        )
      )
    }
  }

  stats <- suff_stats_weighted(x, y, w)

  solution <- logistic4_ntrm(stats, start, max_iter)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 4),
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 4,
    fitted.values = logistic4_function(x, theta),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic4_fit
logistic4_fit_weighted_constrained <- function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  w_zero <- w == 0
  w_positive <- !w_zero

  if (any(w_zero)) {
    w <- w[w_positive]
    x <- x[w_positive]
    y <- y[w_positive]

    if (length(y) == 0) {
      # all cases have weight zero
      return(
        list(
          converged = FALSE,
          iterations = 0L,
          constrained = TRUE,
          estimated = rep(NA_real_, 4),
          coefficients = rep(NA_real_, 4),
          rss = NA_real_,
          df.residual = 0L,
          fitted.values = numeric(0),
          residuals = numeric(0),
          weights = numeric(0)
        )
      )
    }
  }

  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 4, 3)

  for (i in seq_len(4)) {
    lb_is_inf <- is.infinite(lower_bound[i])
    ub_is_inf <- is.infinite(upper_bound[i])

    if (lower_bound[i] == upper_bound[i]) {
      constraint[i, 2] <- TRUE
    } else if (!lb_is_inf || !ub_is_inf) {
      constraint[i, 3] <- TRUE
    } else {
      constraint[i, 1] <- TRUE
    }
  }

  known_param <- ifelse(constraint[, 2], lower_bound, NA_real_)

  stats <- suff_stats_weighted(x, y, w)

  solution <- logistic4_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  # bring the parameters back to their natural scale
  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic4_function(x, theta),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}
