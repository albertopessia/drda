#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic2_function <- function(x, theta) {
  eta <- theta[1]
  phi <- theta[2]

  1 / (1 + exp(-eta * (x - phi)))
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the four parameters in the form
#'   `c(eta, phi)`.
#'
#' @return List of two elements. Element `G` is a numeric matrix of dimension
#'   length(x)-by-2, where each row is the gradient of the logistic function at
#'   the corresponding element of `x`. Element `H` is an array of dimension
#'   length(x)-by-2-by-2, where `H[k, , ]` is the 2-by-2 Hessian matrix at
#'   `x[k]`.
logistic2_gradient_hessian <- function(x, theta) {
  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  t <- q / f
  u <- r / f

  gradient <- matrix(1, nrow = length(x), ncol = 2)
  hessian <- array(0, dim = c(length(x), 2, 2))

  gradient[, 1] <- q / f^2
  gradient[, 2] <- r / f^2

  hessian[, 1, 1] <- (2 * t - (x - phi)) * gradient[, 1]
  hessian[, 2, 1] <- (2 * t - (x - phi) + 1 / eta) * gradient[, 2]

  hessian[, 1, 2] <- hessian[, 2, 1]
  hessian[, 2, 2] <- (2 * u + eta) * gradient[, 2]

  list(G = gradient, H = hessian)
}

#' Residual sum of squares
#'
#' Evaluate the residual sum of squares (RSS) against the mean of a
#' 2-parameter logistic model.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
logistic2_rss <- function(stats) {
  function(theta) {
    mu <- logistic2_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' @rdname logistic2_rss
logistic2_rss_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    mu <- logistic2_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient and Hessian of the residual sum of squares (RSS)
#' against the mean of a 2-parameter logistic model.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
logistic2_rss_gradient_hessian <- function(stats) {
  function(theta) {
    mu <- logistic2_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic2_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 2, 2))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

#' @rdname logistic2_rss_gradient_hessian
logistic2_rss_gradient_hessian_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    mu <- logistic2_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic2_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 2, 2))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)

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
#' @return Numeric vector of length 2 with a (hopefully) good starting point.
#'
#' @importFrom stats median
logistic2_init <- function(stats) {
  k <- nrow(stats)
  delta <- mean(diff(stats[, 1]))

  rss <- logistic2_rss(stats)

  eta_set <- seq(-2, -0.01, length.out = 15)
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[k, 1] + 0.5 * delta, length.out = 15
  )

  theta <- c(-1, stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1])

  if (stats[k, 3] > stats[1, 3]) {
    theta[1] <- 1
  }

  best_rss <- rss(theta)

  for (phi in phi_set) {
    for (eta in eta_set) {
      current_par <- c(eta, phi)
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

#' 2-parameter logistic fit
#'
#' Evaluate the Fisher information matrix at the maximum likelihood estimate.
#'
#' @details
#' Let `mu(x; theta)` be the 2-parameter logistic function. We assume that our
#' observations `y` are independent and such that
#' `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
#' distribution `N(0, 1)`.
#'
#' The 2-by-2 (symmetric) Fisher information matrix is the expected value of
#' the negative Hessian matrix of the log-likelihood function.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param n effective sample size.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
logistic2_fisher_info_normal <- function(stats, n, theta, sigma) {
  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (stats[, 1] - phi))

  f <- 1 + b

  q <- (stats[, 1] - phi) * b
  r <- -eta * b

  gradient <- matrix(1, nrow = nrow(stats), ncol = 2)

  gradient[, 1] <- q / f^2
  gradient[, 2] <- r / f^2

  tmp <- array(0, dim = c(nrow(stats), 2, 2))
  tmp[, , 1] <- stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- stats[, 2] * gradient[, 2] * gradient

  fim <- matrix(0, nrow = 3, ncol = 3)
  fim[1:2, 1:2] <- apply(tmp, 2:3, sum)
  fim[3, 3] <- n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' Fit a 2-parameter logistic function
#'
#' Use a Newton trust-region method to fit a 2-parameter logistic function to
#' observed data.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param constraint boolean matrix with a representation of the constraints.
#' @param lower_bound numeric vector of length 2 with the minimum admissible
#'   values.
#' @param upper_bound numeric vector of length 2 with the maximum admissible
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
logistic2_ntrm <- function(stats, start, max_iter) {
  init <- if (!is.null(start)) {
    start
  } else {
    logistic2_init(stats)
  }

  fn <- logistic2_rss(stats)
  gh <- logistic2_rss_gradient_hessian(stats)

  ntrm(fn, gh, init, max_iter)
}

#' @rdname logistic2_ntrm
logistic2_ntrm_constrained <- function(
  stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
) {
  init <- if (!is.null(start)) {
    # equality constraints have the priority over the provided starting values
    ifelse(is.na(known_param), start, known_param)
  } else {
    ifelse(is.na(known_param), logistic2_init(stats), known_param)
  }

  if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and Hessian
    idx <- which(!constraint[, 2])

    fn <- logistic2_rss_fixed(stats, known_param)
    gh <- logistic2_rss_gradient_hessian_fixed(stats, known_param)

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
    fn <- logistic2_rss(stats)
    gh <- logistic2_rss_gradient_hessian(stats)

    ntrm_constrained(
      fn, gh, init, max_iter, lower_bound, upper_bound
    )
  }
}

#' 2-parameter logistic fit
#'
#' Fit a 2-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of optional weights.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param lower_bound numeric vector of length 2 with the minimum admissible
#'   value of `eta` and `phi` respectively. Values can be equal to `-Inf`.
#' @param upper_bound numeric vector of length 4 with the maximum admissible
#'   value of `eta` and `phi` respectively. Values can be equal to `Inf`.
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
logistic2_fit <- function(x, y, start, max_iter) {
  stats <- suff_stats(x, y)

  solution <- logistic2_ntrm(stats, start, max_iter)

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 2),
    coefficients = solution$optimum,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 2,
    fitted.values = logistic2_function(x, solution$optimum)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic2_fit
logistic2_fit_constrained <- function(
  x, y, start, max_iter, lower_bound, upper_bound
) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 2, 3)

  for (i in seq_len(2)) {
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

  solution <- logistic2_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic2_function(x, theta)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic2_fit
logistic2_fit_weighted <- function(x, y, w, start, max_iter) {
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
          estimated = rep(NA_real_, 2),
          coefficients = rep(NA_real_, 2),
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

  solution <- logistic2_ntrm(stats, start, max_iter)

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 2),
    coefficients = solution$optimum,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 2,
    fitted.values = logistic2_function(x, solution$optimum),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic2_fit
logistic2_fit_weighted_constrained <- function(
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
          estimated = rep(NA_real_, 2),
          coefficients = rep(NA_real_, 2),
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
  constraint <- matrix(FALSE, 2, 3)

  for (i in seq_len(2)) {
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

  solution <- logistic2_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic2_function(x, theta),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}
