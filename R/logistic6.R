#' 6-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 6-parameter logistic function.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic6_function <- function(x, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)
}

#' 6-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 6-parameter logistic function.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (exp(xi') + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)`, `nu = exp(nu')`,
#' `xi = exp(xi')`.
#'
#' NOTE: This function expects the input vector
#' `theta = c(alpha, beta, eta, phi, nu, xi)` but computes the gradient and
#' Hessian of the alternative parametrization
#' `theta' = (alpha, omega, eta, phi, nu', xi')`.
#' This duality is unfortunately a trade-off between what users expect (first
#' version is easier to interpret) and the optimization stability of the second
#' version.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @return List of two elements. Element `G` is a numeric matrix of dimension
#'   length(x)-by-6, where each row is the gradient of the logistic function at
#'   the corresponding element of `x`. Element `H` is an array of dimension
#'   length(x)-by-6-by-6, where `H[k, , ]` is the 6-by-6 Hessian matrix at
#'   `x[k]`.
logistic6_gradient_hessian <- function(x, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  b <- exp(-eta * (x - phi))

  f <- xi + nu * b
  g <- f^(-1 / nu)
  h <- (beta - alpha) * g

  q <- (x - phi) * b
  r <- -eta * b

  s <- h / f
  t <- q / f
  u <- r / f
  v <- u / eta + log(f) / nu

  gradient <- matrix(1, nrow = length(x), ncol = 6)
  hessian <- array(0, dim = c(length(x), 6, 6))

  gradient[, 2] <- h
  gradient[, 3] <- q * s
  gradient[, 4] <- r * s
  gradient[, 5] <- v * h
  gradient[, 6] <- -xi * s / nu

  hessian[, 2, 2] <- gradient[, 2]
  hessian[, 3, 2] <- gradient[, 3]
  hessian[, 4, 2] <- gradient[, 4]
  hessian[, 5, 2] <- gradient[, 5]
  hessian[, 6, 2] <- gradient[, 6]

  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- ((nu + 1) * t - (x - phi)) * gradient[, 3]
  hessian[, 4, 3] <- ((nu + 1) * t - (x - phi) + 1 / eta) * gradient[, 4]
  hessian[, 5, 3] <- (nu * u / eta + v) * gradient[, 3]
  hessian[, 6, 3] <- (nu + 1) * t * gradient[, 6]

  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- ((nu + 1) * u + eta) * gradient[, 4]
  hessian[, 5, 4] <- (nu * u / eta + v) * gradient[, 4]
  hessian[, 6, 4] <- (nu + 1) * u * gradient[, 6]

  hessian[, 2, 5] <- hessian[, 5, 2]
  hessian[, 3, 5] <- hessian[, 5, 3]
  hessian[, 4, 5] <- hessian[, 5, 4]
  hessian[, 5, 5] <- (nu * (u / eta)^2 / v + v - 1) * gradient[, 5]
  hessian[, 6, 5] <- (nu * u / eta + v - 1) * gradient[, 6]

  hessian[, 2, 6] <- hessian[, 6, 2]
  hessian[, 3, 6] <- hessian[, 6, 3]
  hessian[, 4, 6] <- hessian[, 6, 4]
  hessian[, 5, 6] <- hessian[, 6, 5]
  hessian[, 6, 6] <- -((nu + 1) * xi / (nu * f) - 1) * gradient[, 6]

  list(G = gradient, H = hessian)
}

#' Residual sum of squares
#'
#' Evaluate the residual sum of squares (RSS) against the mean of a
#' 6-parameter logistic model.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (exp(xi') + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)`, `nu = exp(nu')`,
#' `xi = exp(xi')`.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
logistic6_rss <- function(stats) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- logistic6_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' @rdname logistic6_rss
logistic6_rss_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- logistic6_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 3] - mu)^2)
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient and Hessian of the residual sum of squares (RSS)
#' against the mean of a 6-parameter logistic model.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (exp(xi') + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)`, `nu = exp(nu')`,
#' `xi = exp(xi')`.
#'
#' @param stats matrix of sufficient statistics.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
logistic6_rss_gradient_hessian <- function(stats) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- logistic6_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic6_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 6, 6))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- stats[, 2] * (r * H[, , 6] + G[, 6] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

#' @rdname logistic6_rss_gradient_hessian
logistic6_rss_gradient_hessian_fixed <- function(stats, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- logistic6_function(stats[, 1], theta)
    mu_gradient_hessian <- logistic6_gradient_hessian(stats[, 1], theta)

    r <- mu - stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(stats), 6, 6))
    hessian[, , 1] <- stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- stats[, 2] * (r * H[, , 6] + G[, 6] * G)

    list(
      G = apply(gradient[, idx, drop = FALSE], 2, sum),
      H = apply(hessian[, idx, idx, drop = FALSE], 2:3, sum)
    )
  }
}

#' Initialize vector of parameters
#'
#' Given the sufficient statistics, try to guess a good approximation to the
#' Maximum Likelihood estimator of the six parameters of the logistic function.
#'
#' @param stats numeric matrix of sufficient statistics.
#'
#' @return Numeric vector of length 6 with a (hopefully) good starting point.
#'
#' @importFrom stats median
logistic6_init <- function(stats) {
  k <- nrow(stats)
  delta <- mean(diff(stats[, 1]))

  rss <- logistic6_rss(stats)

  w <- stats[, 2]

  eta_set <- seq(-2, -0.01, length.out = 5)
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[k, 1] + 0.5 * delta, length.out = 5
  )
  log_nu_set <- seq(-1, 0.5, length.out = 5)
  log_xi_set <- seq(-1, 0.5, length.out = 5)

  theta <- c(
    min(stats[, 3]),
    log(max(stats[, 3]) - min(stats[, 3])),
    -1,
    stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1],
    0,
    0
  )

  if (stats[k, 3] > stats[1, 3]) {
    theta[3] <- 1
  }

  best_rss <- rss(theta)

  for (eta in eta_set) {
    for (phi in phi_set) {
      for (log_nu in log_nu_set) {
        for (log_xi in log_xi_set) {
          xi <- exp(log_xi)
          nu <- exp(log_nu)

          g <- 1 / (xi + nu * exp(-eta * (stats[, 1] - phi)))^(1 / nu)

          t1 <- 0
          t2 <- 0
          t3 <- 0
          t4 <- 0
          t5 <- 0

          for (i in seq_len(k)) {
            t1 <- t1 + w[i] * g[i]
            t2 <- t2 + w[i] * g[i] * stats[i, 3]
            t3 <- t3 + w[i] * stats[i, 3]
            t4 <- t4 + w[i] * g[i]^2
            t5 <- t5 + w[i]
          }

          alpha <- (t1 * t2 - t3 * t4) / (t1^2 - t4 * t5)
          log_beta <- log((t1 * t3 - t2 * t5) / (t1^2 - t4 * t5))

          current_par <- c(alpha, log_beta, eta, phi, log_nu, log_xi)
          current_rss <- rss(current_par)

          if (!is.nan(current_rss) && (current_rss < best_rss)) {
            theta <- current_par
            best_rss <- current_rss
          }
        }
      }
    }
  }

  names(theta) <- NULL

  theta
}

#' 6-parameter logistic fit
#'
#' Evaluate the Fisher information matrix at the maximum likelihood estimate.
#'
#' @details
#' Let `mu(x; theta)` be the 6-parameter logistic function. We assume that our
#' observations `y` are independent and such that
#' `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
#' distribution `N(0, 1)`.
#'
#' The 6-by-6 (symmetric) Fisher information matrix is the expected value of
#' the negative Hessian matrix of the log-likelihood function.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param n effective sample size.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
logistic6_fisher_info_normal <- function(stats, n, theta, sigma) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  b <- exp(-eta * (stats[, 1] - phi))

  f <- xi + nu * b
  g <- f^(-1 / nu)
  h <- (beta - alpha) * g

  q <- (stats[, 1] - phi) * b
  r <- -eta * b

  t <- q / f

  u <- r / f
  v <- u / eta + log(f) / nu

  gradient <- matrix(1, nrow = nrow(stats), ncol = 6)

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- t * h
  gradient[, 4] <- u * h
  gradient[, 5] <- v * h / nu
  gradient[, 6] <- -h / (nu * f)

  tmp <- array(0, dim = c(nrow(stats), 6, 6))
  tmp[, , 1] <- stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- stats[, 2] * gradient[, 2] * gradient
  tmp[, , 3] <- stats[, 2] * gradient[, 3] * gradient
  tmp[, , 4] <- stats[, 2] * gradient[, 4] * gradient
  tmp[, , 5] <- stats[, 2] * gradient[, 5] * gradient
  tmp[, , 6] <- stats[, 2] * gradient[, 6] * gradient

  fim <- matrix(0, nrow = 7, ncol = 7)
  fim[1:6, 1:6] <- apply(tmp, 2:3, sum)
  fim[7, 7] <- n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' Fit a 6-parameter logistic function
#'
#' Use a Newton trust-region method to fit a 6-parameter logistic function to
#' observed data.
#'
#' @param stats numeric matrix of sufficient statistics.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param constraint boolean matrix with a representation of the constraints.
#' @param lower_bound numeric vector of length 6 with the minimum admissible
#'   values.
#' @param upper_bound numeric vector of length 6 with the maximum admissible
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
logistic6_ntrm <- function(stats, start, max_iter) {
  init <- if (!is.null(start)) {
    start
  } else {
    logistic6_init(stats)
  }

  fn <- logistic6_rss(stats)
  gh <- logistic6_rss_gradient_hessian(stats)

  ntrm(fn, gh, init, max_iter)
}

#' @rdname logistic6_ntrm
logistic6_ntrm_constrained <- function(
  stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
) {
  init <- if (!is.null(start)) {
    # equality constraints have the priority over the provided starting values
    ifelse(is.na(known_param), start, known_param)
  } else {
    ifelse(is.na(known_param), logistic6_init(stats), known_param)
  }

  if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and Hessian
    idx <- which(!constraint[, 2])

    fn <- logistic6_rss_fixed(stats, known_param)
    gh <- logistic6_rss_gradient_hessian_fixed(stats, known_param)

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
    fn <- logistic6_rss(stats)
    gh <- logistic6_rss_gradient_hessian(stats)

    ntrm_constrained(
      fn, gh, init, max_iter, lower_bound, upper_bound
    )
  }
}

#' 6-parameter logistic fit
#'
#' Fit a 6-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (exp(xi') + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)`, `nu = exp(nu')`,
#' `xi = exp(xi')`.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of optional weights.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param lower_bound numeric vector of length 6 with the minimum admissible
#'   value of `alpha`, `omega`, `eta`, `phi`, `nu'`, and `xi'` respectively.
#'   Values can be equal to `-Inf`.
#' @param upper_bound numeric vector of length 6 with the maximum admissible
#'   value of `alpha`, `omega`, `eta`, `phi`, `nu'`, and `xi'` respectively.
#'   Values can be equal to `Inf`.
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
logistic6_fit <- function(x, y, start, max_iter) {
  stats <- suff_stats(x, y)

  solution <- logistic6_ntrm(stats, start, max_iter)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])
  theta[5] <- exp(theta[5])
  theta[6] <- exp(theta[6])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 6),
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 6,
    fitted.values = logistic6_function(x, theta)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic6_fit
logistic6_fit_constrained <- function(
  x, y, start, max_iter, lower_bound, upper_bound
) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 6, 3)

  for (i in seq_len(6)) {
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

  solution <- logistic6_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  # bring the parameters back to their natural scale
  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])
  theta[5] <- exp(theta[5])
  theta[6] <- exp(theta[6])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic6_function(x, theta)
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic6_fit
logistic6_fit_weighted <- function(x, y, w, start, max_iter) {
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
          estimated = rep(NA_real_, 6),
          coefficients = rep(NA_real_, 6),
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

  solution <- logistic6_ntrm(stats, start, max_iter)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])
  theta[5] <- exp(theta[5])
  theta[6] <- exp(theta[6])

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 6),
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - 6,
    fitted.values = logistic6_function(x, theta),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}

#' @rdname logistic6_fit
logistic6_fit_weighted_constrained <- function(
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
          estimated = rep(NA_real_, 6),
          coefficients = rep(NA_real_, 6),
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
  constraint <- matrix(FALSE, 6, 3)

  for (i in seq_len(6)) {
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

  solution <- logistic6_ntrm_constrained(
    stats, start, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  # bring the parameters back to their natural scale
  theta <- lower_bound
  theta[!constraint[, 2]] <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])
  theta[5] <- exp(theta[5])
  theta[6] <- exp(theta[6])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(stats[, 2] * stats[, 4]) + solution$minimum,
    df.residual = length(y) - sum(estimated),
    fitted.values = logistic6_function(x, theta),
    weights = w
  )

  result$residuals <- y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result
}
