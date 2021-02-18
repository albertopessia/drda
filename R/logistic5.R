#' @rdname logistic6_new
logistic5_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 5) {
      stop("'start' must be of length 5", call. = FALSE)
    }

    if (start[2] <= start[1]) {
      stop("parameter 'beta' is smaller than 'alpha'", call. = FALSE)
    }

    if (start[3] == 0) {
      stop("parameter 'eta' cannot be initialized to zero", call. = FALSE)
    }

    if (start[5] <= 0) {
      stop("parameter 'nu' cannot be negative nor zero", call. = FALSE)
    }

    start[2] <- log(start[2] - start[1])
    start[5] <- log(start[5])
  }

  if (!is.null(lower_bound)) {
    if (length(lower_bound) != 5) {
      stop("'lower_bound' must be of length 5", call. = FALSE)
    }

    if (lower_bound[2] < lower_bound[1]) {
      stop(
        "'lower_bound[2]' cannot be smaller than 'lower_bound[1]'",
        call. = FALSE
      )
    }
  }

  if (!is.null(upper_bound)) {
    if (length(upper_bound) != 5) {
      stop("'upper_bound' must be of length 5", call. = FALSE)
    }

    if (upper_bound[2] < upper_bound[1]) {
      stop(
        "'upper_bound[2]' cannot be smaller than 'upper_bound[1]'",
        call. = FALSE
      )
    }
  }

  if (!is.null(lower_bound) && !is.null(upper_bound)) {
    if (lower_bound[2] < upper_bound[1]) {
      stop(
        "'lower_bound[2]' cannot be smaller than 'upper_bound[1]'",
        call. = FALSE
      )
    }

    if (!is.infinite(lower_bound[2])) {
      lower_bound[2] <- log(lower_bound[2] - upper_bound[1])
    }

    lower_bound[5] <- if (lower_bound[5] > 0) {
      log(lower_bound[5])
    } else {
      -Inf
    }

    if (!is.infinite(upper_bound[2])) {
      upper_bound[2] <- log(upper_bound[2] - lower_bound[1])
    }

    upper_bound[5] <- if (upper_bound[5] > 0) {
      log(upper_bound[5])
    } else {
      -Inf
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
    class = "logistic5"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    object$lower_bound <- if (!is.null(lower_bound)) {
      lower_bound
    } else {
      rep(-Inf, 5)
    }

    object$upper_bound <- if (!is.null(upper_bound)) {
      upper_bound
    } else {
      rep(Inf, 5)
    }
  }

  object
}

#' 5-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
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

#' 5-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' @param object object of class `logistic5`.
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, beta, eta, phi, nu)`.
#'
#' @return Numeric vector with the values of the logistic function.
fn.logistic5 <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

#' @rdname fn.logistic5
fn.logistic5_fit <- function(object, x, theta) {
  logistic5_fn(x, theta)
}

#' 5-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 5-parameter logistic function.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)` and `nu = exp(nu')`.
#'
#' NOTE: This function expects the input vector
#' `theta = c(alpha, beta, eta, phi, nu)` but computes the gradient and
#' Hessian of the alternative parametrization
#' `theta' = (alpha, omega, eta, phi, nu')`.
#' This duality is unfortunately a trade-off between what users expect (first
#' version is easier to interpret) and the optimization stability of the second
#' version.
#'
#' @param object object of class `logistic5`.
#' @param theta numeric vector with the five parameters in the form
#'   `c(alpha, beta, eta, phi, nu)`.
#'
#' @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic5 <- function(object, theta) {
  x <- object$stats[, 1]

  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  b <- exp(-eta * (x - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)
  h <- (beta - alpha) * g

  q <- (x - phi) * b
  r <- -eta * b

  s <- h / f
  t <- q / f
  u <- r / f
  v <- u / eta + log(f) / nu

  gradient <- matrix(1, nrow = length(x), ncol = 5)
  hessian <- array(0, dim = c(length(x), 5, 5))

  gradient[, 2] <- h
  gradient[, 3] <- q * s
  gradient[, 4] <- r * s
  gradient[, 5] <- v * h

  hessian[, 2, 2] <- gradient[, 2]
  hessian[, 3, 2] <- gradient[, 3]
  hessian[, 4, 2] <- gradient[, 4]
  hessian[, 5, 2] <- gradient[, 5]

  hessian[, 2, 3] <- hessian[, 3, 2]
  hessian[, 3, 3] <- ((nu + 1) * t - (x - phi)) * gradient[, 3]
  hessian[, 4, 3] <- ((nu + 1) * t - (x - phi) + 1 / eta) * gradient[, 4]
  hessian[, 5, 3] <- (nu * u / eta + v) * gradient[, 3]

  hessian[, 2, 4] <- hessian[, 4, 2]
  hessian[, 3, 4] <- hessian[, 4, 3]
  hessian[, 4, 4] <- ((nu + 1) * u + eta) * gradient[, 4]
  hessian[, 5, 4] <- (nu * u / eta + v) * gradient[, 4]

  hessian[, 2, 5] <- hessian[, 5, 2]
  hessian[, 3, 5] <- hessian[, 5, 3]
  hessian[, 4, 5] <- hessian[, 5, 4]
  hessian[, 5, 5] <- (nu * (u / eta)^2 / v + v - 1) * gradient[, 5]

  # When `b` is infinite, gradient and Hessian show NaNs
  # In the limit b -> Inf, both gradient and Hessian converge to zero
  if (any(is.nan(gradient))) {
    gradient[is.nan(gradient)] <- 0
  }

  if (any(is.nan(hessian))) {
    hessian[is.nan(hessian)] <- 0
  }

  list(G = gradient, H = hessian)
}

#' Residual sum of squares
#'
#' Evaluate the residual sum of squares (RSS) against the mean of a
#' 5-parameter logistic model.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)` and `nu = exp(nu')`.
#'
#' @param object object of class `logistic5`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
rss.logistic5 <- function(object) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

#' @rdname rss.logistic5
rss_fixed.logistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient and Hessian of the residual sum of squares (RSS)
#' against the mean of a 5-parameter logistic model.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)` and `nu = exp(nu')`.
#'
#' @param object object of class `logistic5`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic5 <- function(object) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
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

#' @rdname rss_gradient_hessian.logistic5
rss_gradient_hessian_fixed.logistic5 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 5)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
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

#' Maximum likelihood estimators
#'
#' Given a set of parameters, compute the maximum likelihood estimates of the
#' lower and upper horizontal asymptotes.
#'
#' @param object object of class `logistic5`.
#' @param theta vector of parameters.
#'
#' @return Numeric vector of length 2 with the MLE of the two asymptotes.
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
    t1 <- t1 + w[i] * g[i]
    t2 <- t2 + w[i] * g[i] * y[i]
    t3 <- t3 + w[i] * y[i]
    t4 <- t4 + w[i] * g[i]^2
    t5 <- t5 + w[i]
  }

  alpha <- (t1 * t2 - t3 * t4) / (t1^2 - t4 * t5)
  omega <- (t1 * t3 - t2 * t5) / (t1^2 - t4 * t5)

  if (omega > 0) {
    theta[1] <- alpha
    theta[2] <- log(omega)
  }

  theta
}

#' Initialize vector of parameters
#'
#' Given the sufficient statistics, try to guess a good approximation to the
#' Maximum Likelihood estimator of the five parameters of the logistic function.
#'
#' @param object object of class `logistic5`.
#'
#' @return Numeric vector of length 5 with a (hopefully) good starting point.
#'
#' @importFrom stats coefficients lm median nls nls.control
init.logistic5 <- function(object) {
  m <- object$m
  stats <- object$stats

  linear_fit <- summary(lm(stats[, 3] ~ stats[, 1], weights = stats[, 3]))
  linear_coef <- linear_fit$coefficients

  if (linear_coef[2, 4] > 0.2) {
    # we are in big problems as a flat horizontal line is likely the best model
    tiny <- 1.0e-30
    return(
      c(
        # alpha
        sum(object$w * object$y) / sum(object$w),
        # log(omega)
        log(tiny),
        # eta
        if (linear_coef[2, 1] <= 0) -tiny else tiny,
        # phi
        object$stats[m, 1] + 1000,
        # log(nu)
        0
      )
    )
  }

  delta <- mean(diff(stats[, 1]))

  rss_fn <- rss(object)

  eta_set <- seq(-2, -0.01, length.out = 5)
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[m, 1] + 0.5 * delta, length.out = 5
  )
  log_nu_set <- seq(-1, 0.5, length.out = 5)

  theta <- c(
    min(stats[, 3]),
    log(max(stats[, 3]) - min(stats[, 3])),
    -1,
    stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1],
    0
  )

  if (linear_coef[2, 1] > 0) {
    eta_set <- seq(0.01, 2, length.out = 5)
    theta[3] <- 1
  }

  best_rss <- rss_fn(theta)

  for (eta in eta_set) {
    for (phi in phi_set) {
      for (log_nu in log_nu_set) {
        current_par <- mle_asy(object, c(theta[1], theta[2], eta, phi, log_nu))
        current_rss <- rss_fn(current_par)

        if (!is.nan(current_rss) && (current_rss < best_rss)) {
          theta <- current_par
          best_rss <- current_rss
        }
      }
    }
  }

  D <- data.frame(y = object$y, x = object$x)
  frm <- y ~ alpha + exp(omega) / (1 + exp(psi - eta * (x - phi)))^exp(-psi)
  start <- c(
    "alpha" = theta[1], "omega" = theta[2], "eta" = theta[3], "phi" = theta[4],
    "psi" = theta[5]
  )
  ctrl <- nls.control(
    tol = sqrt(.Machine$double.eps), minFactor = 1.0e-5, warnOnly = TRUE
  )

  fit_nls <- tryCatch(
    {
      suppressWarnings(
        if (!object$constrained) {
          nls(
            formula = frm, data = D, start = start, control = ctrl,
            weights = object$w
          )
        } else {
          nls(
            formula = frm, data = D, start = start, control = ctrl,
            algorithm = "port", weights = object$w, lower = object$lower_bound,
            upper = object$upper_bound
          )
        }
      )
    },
    error = function(e) NULL
  )

  if (!is.null(fit_nls)) {
    current_par <- coefficients(fit_nls)
    current_rss <- rss_fn(current_par)

    if (!is.nan(current_rss) && (current_rss < best_rss)) {
      theta <- current_par
    }
  }

  names(theta) <- NULL

  theta
}

#' 5-parameter logistic fit
#'
#' Fit a 5-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' In our optimization algorithm, however, we consider instead the equivalent
#' function `f(x; theta')`
#'
#' `alpha + exp(omega) / (1 + exp(nu' - eta * (x - phi)))^(-exp(-nu'))`
#'
#' Note that it is simply `beta = alpha + exp(omega)` and `nu = exp(nu')`.
#'
#' @param object object of class `logistic5`.
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
fit.logistic5 <- function(object) {
  solution <- find_optimum(object)

  # bring the parameters back to their natural scale
  theta <- solution$optimum
  theta[2] <- theta[1] + exp(theta[2])
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

#' @rdname fit.logistic5
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
  theta[2] <- theta[1] + exp(theta[2])
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

#' 5-parameter logistic fit
#'
#' Evaluate the Fisher information matrix at the maximum likelihood estimate.
#'
#' @details
#' Let `mu(x; theta)` be the 5-parameter logistic function. We assume that our
#' observations `y` are independent and such that
#' `y = mu(x; theta) + sigma * epsilon`, where `epsilon` has a standard Normal
#' distribution `N(0, 1)`.
#'
#' The 5-by-5 (symmetric) Fisher information matrix is the expected value of
#' the negative Hessian matrix of the log-likelihood function.
#'
#' @param object object of class `logistic5`.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic5 <- function(object, theta, sigma) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]

  b <- exp(-eta * (object$stats[, 1] - phi))

  f <- 1 + nu * b
  g <- f^(-1 / nu)
  h <- (beta - alpha) * g

  q <- (object$stats[, 1] - phi) * b
  r <- -eta * b

  t <- q / f

  u <- r / f
  v <- u / eta + log(f) / nu

  gradient <- matrix(1, nrow = nrow(object$stats), ncol = 5)

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- t * h
  gradient[, 4] <- u * h
  gradient[, 5] <- v * h / nu

  tmp <- array(0, dim = c(nrow(object$stats), 5, 5))
  tmp[, , 1] <- object$stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- object$stats[, 2] * gradient[, 2] * gradient
  tmp[, , 3] <- object$stats[, 2] * gradient[, 3] * gradient
  tmp[, , 4] <- object$stats[, 2] * gradient[, 4] * gradient
  tmp[, , 5] <- object$stats[, 2] * gradient[, 5] * gradient

  fim <- matrix(0, nrow = 6, ncol = 6)
  fim[1:5, 1:5] <- apply(tmp, 2:3, sum)
  fim[6, 6] <- object$n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' 5-parameter logistic fit
#'
#' Evaluate the normalized area under the curve (AUC) and area above the curve
#' (AAC).
#'
#' @details
#' The 5-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu)`.
#'
#' The area under the curve (AUC) is simply the integral of `f(x; theta)`
#' between `lower_bound` and `upper_bound` with respect to `x`.
#'
#' When the interval of integration is fixed, since the the curve ranges between
#' `alpha` and `beta`, the curve `f(x; theta)` is contained into the rectangle
#' of height `beta - alpha` and width `upper_bound - lower_bound`. The maximum
#' area the curve can have is obviously
#' `(upper_bound - lower_bound) * (beta - alpha)`.
#'
#' We first shift the curve by `alpha` to set the minimum to 0. We then
#' integrate the curve and define the normalized AUC (NAUC) by dividing its
#' value by the maximum area. As a consequence, the normalized area above the
#' curve is simply `NAAC = 1 - NAUC`.
#'
#' Default values of `lower_bound` and `upper_bound` were chosen based on common
#' dose ranges used in the literature. They are also symmetric around zero
#' so that `NAUC` and `NAAC` are equal to `0.5` in the standard logistic model.
#'
#' *Note*: Integral of a 5-parameter logistic function involves Gaussian
#' hypergeometric functions, which are prone to numerical errors and not part
#' of base R. For this reason, we opt here for a brute force approach and simply
#' use the `integrate` function.
#'
#' @param object object of class `logistic5_fit`.
#' @param lower_bound numeric value with the lower bound of the integration
#'   interval.
#' @param upper_bound numeric value with the upper bound of the integration
#'   interval.
#'
#' @return Numeric value with the requested area.
#'
#' @importFrom stats integrate
nauc.logistic5_fit <- function(object, lower_bound = -10, upper_bound = 10) {
  alpha <- object$coefficients[1]
  beta <- object$coefficients[2]
  eta <- object$coefficients[3]
  phi <- object$coefficients[4]
  nu <- object$coefficients[5]

  f <- function(z) {
    1 / (1 + nu * exp(-eta * (z - phi)))^(1 / nu)
  }

  I <- integrate(f, lower = lower_bound, upper = upper_bound)

  nauc <- I$value / ((beta - alpha) * (upper_bound - lower_bound))
  names(nauc) <- NULL

  nauc
}

#' @rdname nauc.logistic5_fit
naac_logistic5_fit <- function(object, lower_bound = -10, upper_bound = 10) {
  1 - nauc(object, lower_bound, upper_bound)
}
