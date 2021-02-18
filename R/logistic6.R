#' Fit dose-response data
#'
#' Use a Newton trust-region method to fit a logistic function to dose-response
#' data.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w an optional vector of weights to be used in the fitting
#'   process.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#' @param lower_bound numeric vector with the minimum admissible values of the
#'   parameters.
#' @param upper_bound numeric vector with the maximum admissible values of the
#'   parameters.
#'
#' @return An object of class `logistic*`.
logistic6_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 6) {
      stop("'start' must be of length 6", call. = FALSE)
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

    if (start[6] <= 0) {
      stop("parameter 'xi' cannot be negative nor zero", call. = FALSE)
    }

    start[2] <- log(start[2] - start[1])
    start[5] <- log(start[5])
    start[6] <- log(start[6])
  }

  if (!is.null(lower_bound)) {
    if (length(lower_bound) != 6) {
      stop("'lower_bound' must be of length 6", call. = FALSE)
    }

    if (lower_bound[2] < lower_bound[1]) {
      stop(
        "'lower_bound[2]' cannot be smaller than 'lower_bound[1]'",
        call. = FALSE
      )
    }
  }

  if (!is.null(upper_bound)) {
    if (length(upper_bound) != 6) {
      stop("'upper_bound' must be of length 6", call. = FALSE)
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

    lower_bound[6] <- if (lower_bound[6] > 0) {
      log(lower_bound[6])
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

    upper_bound[6] <- if (upper_bound[6] > 0) {
      log(upper_bound[6])
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
    class = "logistic6"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    object$lower_bound <- if (!is.null(lower_bound)) {
      lower_bound
    } else {
      rep(-Inf, 6)
    }

    object$upper_bound <- if (!is.null(upper_bound)) {
      upper_bound
    } else {
      rep(Inf, 6)
    }
  }

  object
}

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
logistic6_fn <- function(x, theta) {
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
#' Evaluate at a particular set of parameters the 6-parameter logistic function.
#'
#' @details
#' The 6-parameter logistic function `f(x; theta)` is defined in this package as
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `theta = c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @param object object of class `logistic6`.
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @return Numeric vector with the values of the logistic function.
fn.logistic6 <- function(object, x, theta) {
  logistic6_fn(x, theta)
}

#' @rdname fn.logistic6
fn.logistic6_fit <- function(object, x, theta) {
  logistic6_fn(x, theta)
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
#' @param object object of class `logistic6`.
#' @param theta numeric vector with the six parameters in the form
#'   `c(alpha, beta, eta, phi, nu, xi)`.
#'
#' @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic6 <- function(object, theta) {
  x <- object$stats[, 1]

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
#' @param object object of class `logistic6`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
rss.logistic6 <- function(object) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

#' @rdname rss.logistic6
rss_fixed.logistic6 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
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
#' @param object object of class `logistic6`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic6 <- function(object) {
  function(theta) {
    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 6, 6))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- object$stats[, 2] * (r * H[, , 6] + G[, 6] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

#' @rdname rss_gradient_hessian.logistic6
rss_gradient_hessian_fixed.logistic6 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 6)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    theta[2] <- theta[1] + exp(theta[2])
    theta[5] <- exp(theta[5])
    theta[6] <- exp(theta[6])

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(nrow(object$stats), 6, 6))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)
    hessian[, , 3] <- object$stats[, 2] * (r * H[, , 3] + G[, 3] * G)
    hessian[, , 4] <- object$stats[, 2] * (r * H[, , 4] + G[, 4] * G)
    hessian[, , 5] <- object$stats[, 2] * (r * H[, , 5] + G[, 5] * G)
    hessian[, , 6] <- object$stats[, 2] * (r * H[, , 6] + G[, 6] * G)

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
#' @param object object of class `logistic6`.
#' @param theta vector of parameters.
#'
#' @return Numeric vector of length 2 with the MLE of the two asymptotes.
mle_asy.logistic6 <- function(object, theta) {
  m <- object$m

  x <- object$stats[, 1]
  y <- object$stats[, 3]
  w <- object$stats[, 2]

  eta <- theta[3]
  phi <- theta[4]
  log_nu <- theta[5]
  log_xi <- theta[6]

  g <- 1 / (exp(log_xi) + exp(log_nu - eta * (x - phi)))^exp(-log_nu)

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
#' Maximum Likelihood estimator of the six parameters of the logistic function.
#'
#' @param object object of class `logistic6`.
#'
#' @return Numeric vector of length 6 with a (hopefully) good starting point.
#'
#' @importFrom stats coefficients lm median nls nls.control
init.logistic6 <- function(object) {
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
        median(object$stats[, 1]),
        # log(nu)
        0,
        # log(xi)
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
  log_xi_set <- seq(-1, 0.5, length.out = 5)

  theta <- c(
    min(stats[, 3]),
    log(max(stats[, 3]) - min(stats[, 3])),
    -1,
    stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1],
    0,
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
        for (log_xi in log_xi_set) {
          current_par <- mle_asy(
            object, c(theta[1], theta[2], eta, phi, log_nu, log_xi)
          )
          current_rss <- rss_fn(current_par)

          if (!is.nan(current_rss) && (current_rss < best_rss)) {
            theta <- current_par
            best_rss <- current_rss
          }
        }
      }
    }
  }

  D <- data.frame(y = object$y, x = object$x)
  frm <- y ~ alpha + exp(omega) / (
    exp(lambda) + exp(psi - eta * (x - phi))
  )^exp(-psi)
  start <- c(
    "alpha" = theta[1], "omega" = theta[2], "eta" = theta[3], "phi" = theta[4],
    "psi" = theta[5], "lambda" = theta[6]
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
#' @param object object of class `logistic6`.
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
fit.logistic6 <- function(object) {
  solution <- find_optimum(object)

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
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - 6,
    fitted.values = logistic6_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic6_fit"

  result
}

#' @rdname fit.logistic6
fit_constrained.logistic6 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 6, 3)

  for (i in seq_len(6)) {
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
  theta[6] <- exp(theta[6])

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = object$n - sum(estimated),
    fitted.values = logistic6_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("alpha", "beta", "eta", "phi", "nu", "xi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic6_fit"

  result
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
#' @param object object of class `logistic6`.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic6 <- function(object, theta, sigma) {
  alpha <- theta[1]
  beta <- theta[2]
  eta <- theta[3]
  phi <- theta[4]
  nu <- theta[5]
  xi <- theta[6]

  b <- exp(-eta * (object$stats[, 1] - phi))

  f <- xi + nu * b
  g <- f^(-1 / nu)
  h <- (beta - alpha) * g

  q <- (object$stats[, 1] - phi) * b
  r <- -eta * b

  t <- q / f

  u <- r / f
  v <- u / eta + log(f) / nu

  gradient <- matrix(1, nrow = nrow(object$stats), ncol = 6)

  gradient[, 1] <- 1 - g
  gradient[, 2] <- g
  gradient[, 3] <- t * h
  gradient[, 4] <- u * h
  gradient[, 5] <- v * h / nu
  gradient[, 6] <- -h / (nu * f)

  tmp <- array(0, dim = c(nrow(object$stats), 6, 6))
  tmp[, , 1] <- object$stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- object$stats[, 2] * gradient[, 2] * gradient
  tmp[, , 3] <- object$stats[, 2] * gradient[, 3] * gradient
  tmp[, , 4] <- object$stats[, 2] * gradient[, 4] * gradient
  tmp[, , 5] <- object$stats[, 2] * gradient[, 5] * gradient
  tmp[, , 6] <- object$stats[, 2] * gradient[, 6] * gradient

  fim <- matrix(0, nrow = 7, ncol = 7)
  fim[1:6, 1:6] <- apply(tmp, 2:3, sum)
  fim[7, 7] <- object$n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}
