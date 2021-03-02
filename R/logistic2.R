#' @rdname logistic6_new
logistic2_new <-  function(
  x, y, w, start, max_iter, lower_bound, upper_bound
) {
  if (!is.null(start)) {
    if (length(start) != 2) {
      stop("'start' must be of length 2", call. = FALSE)
    }

    if (start[1] == 0) {
      stop("parameter 'eta' cannot be initialized to zero", call. = FALSE)
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
    class = "logistic2"
  )

  object$m <- nrow(object$stats)

  if (!is.null(lower_bound) || !is.null(upper_bound)) {
    object$constrained <- TRUE

    if (is.null(lower_bound)) {
      rep(-Inf, 2)
    } else {
      if (length(lower_bound) != 2) {
        stop("'lower_bound' must be of length 2", call. = FALSE)
      }
    }

    if (is.null(upper_bound)) {
      rep(Inf, 2)
    } else {
      if (length(upper_bound) != 2) {
        stop("'upper_bound' must be of length 2", call. = FALSE)
      }
    }

    object$lower_bound <- lower_bound
    object$upper_bound <- upper_bound
  }

  object
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the parameters in the form `c(eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
#'
#' @export
logistic2_fn <- function(x, theta) {
  eta <- theta[1]
  phi <- theta[2]

  1 / (1 + exp(-eta * (x - phi)))
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param object object of class `logistic2`.
#' @param x numeric vector at which the logistic function is to be evaluated.
#' @param theta numeric vector with the parameters in the form `c(eta, phi)`.
#'
#' @return Numeric vector of the same length of `x` with the values of the
#'   logistic function.
fn.logistic2 <- function(object, x, theta) {
  logistic2_fn(x, theta)
}

#' @rdname fn.logistic2
fn.logistic2_fit <- function(object, x, theta) {
  logistic2_fn(x, theta)
}

#' 2-parameter logistic function
#'
#' Evaluate at a particular set of parameters the gradient and Hessian of the
#' 2-parameter logistic function.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param object object of class `logistic2`.
#' @param theta numeric vector with the parameters in the form `c(eta, phi)`.
#'
#' @return List of two elements: `G` the gradient and `H` the Hessian.
gradient_hessian.logistic2 <- function(object, theta) {
  x <- object$stats[, 1]

  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  s <- 1 / f^2
  t <- q * s
  u <- r * s

  gradient <- matrix(0, nrow = length(x), ncol = 2)
  hessian <- array(0, dim = c(length(x), 2, 2))

  gradient[, 1] <- t
  gradient[, 2] <- u

  hessian[, 1, 1] <- q * t * (2 / f - 1 / b)
  hessian[, 2, 1] <- (1 / eta + (2 - f / b) * t * f) * u

  hessian[, 1, 2] <- hessian[, 2, 1]
  hessian[, 2, 2] <- (2 / f - 1 / b) * r * u

  # When `b` is infinite, gradient and Hessian show NaNs
  # these are the limits for b -> Inf
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
#' 2-parameter logistic model.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param object object of class `logistic2`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(p)` to evaluate the RSS associated to a particular
#'   parameter choice `p`.
rss.logistic2 <- function(object) {
  function(theta) {
    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

#' @rdname rss.logistic2
rss_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    mu <- fn(object, object$stats[, 1], theta)
    sum(object$stats[, 2] * (object$stats[, 3] - mu)^2)
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient and Hessian of the residual sum of squares (RSS)
#' against the mean of a 2-parameter logistic model.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param object object of class `logistic2`.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any.
#'
#' @return Function handle `f(theta)` to evaluate the gradient and Hessian of
#'   the RSS associated to a particular parameter choice `theta`.
rss_gradient_hessian.logistic2 <- function(object) {
  function(theta) {
    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

    list(G = apply(gradient, 2, sum), H = apply(hessian, 2:3, sum))
  }
}

#' @rdname rss_gradient_hessian.logistic2
rss_gradient_hessian_fixed.logistic2 <- function(object, known_param) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 2)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    mu <- fn(object, object$stats[, 1], theta)
    mu_gradient_hessian <- gradient_hessian(object, theta)

    r <- mu - object$stats[, 3]

    G <- mu_gradient_hessian$G
    H <- mu_gradient_hessian$H

    gradient <- object$stats[, 2] * r * G

    hessian <- array(0, dim = c(object$m, 2, 2))
    hessian[, , 1] <- object$stats[, 2] * (r * H[, , 1] + G[, 1] * G)
    hessian[, , 2] <- object$stats[, 2] * (r * H[, , 2] + G[, 2] * G)

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
#' @param object object of class `logistic2`.
#' @param theta vector of parameters.
#'
#' @return Numeric vector `theta`.
mle_asy.logistic2 <- function(object, theta) {
  theta
}

#' Initialize vector of parameters
#'
#' Given the sufficient statistics, try to guess a good approximation to the
#' Maximum Likelihood estimator of the four parameters of the logistic function.
#'
#' @param object object of class `logistic2`.
#'
#' @return Numeric vector of length 2 with a (hopefully) good starting point.
#'
#' @importFrom stats coefficients lm median nls nls.control optim
init.logistic2 <- function(object) {
  m <- object$m
  stats <- object$stats
  rss_fn <- rss(object)

  linear_fit <- summary(lm(stats[, 3] ~ stats[, 1], weights = stats[, 2]))
  linear_coef <- linear_fit$coefficients

  theta <- c(
    if (linear_coef[2, 1] < 0) -1 else 1,
    stats[which.min(abs(stats[, 3] - mean(range(stats[, 3])))), 1]
  )

  best_rss <- rss_fn(theta)

  if (linear_coef[2, 4] > 0.2) {
    # we are in big problems as a flat horizontal line is likely the best model
    current_par <- c(-10, stats[m, 1] + 100)
    current_rss <- rss_fn(current_par)

    if (!is.nan(current_rss) && (current_rss < best_rss)) {
      theta <- current_par
      best_rss <- current_rss
    }
  }

  delta <- mean(diff(stats[, 1]))

  eta_set <- if (linear_coef[2, 1] < 0) {
    seq(-2, -0.01, length.out = 15)
  } else {
    seq(0.01, 2, length.out = 15)
  }
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[m, 1] + 0.5 * delta, length.out = 15
  )

  for (phi in phi_set) {
    for (eta in eta_set) {
      current_par <- c(eta, phi)
      current_rss <- rss_fn(current_par)

      if (!is.nan(current_rss) && (current_rss < best_rss)) {
        theta <- current_par
        best_rss <- current_rss
      }
    }
  }

  fit_optim <- tryCatch(
    {
      rss_gh <- rss_gradient_hessian(object)
      gr <- function(par) {
        rss_gh(par)$G
      }

      suppressWarnings(
        if (!object$constrained) {
          optim(
            theta, rss_fn, gr, method = "BFGS",
            control = list(trace = 0, maxit = 10000, reltol = 1.0e-10)
          )
        } else {
          optim(
            theta, rss_fn, gr, method = "L-BFGS-B",
            lower = object$lower_bound, upper = object$upper_bound,
            control = list(trace = 0, maxit = 10000, reltol = 1.0e-10)
          )
        }
      )
    },
    error = function(e) NULL
  )

  if (!is.null(fit_optim)) {
    current_par <- mle_asy(object, fit_optim$par)
    current_rss <- rss_fn(current_par)

    if (!is.nan(current_rss) && (current_rss < best_rss)) {
      theta <- current_par
      best_rss <- current_rss
    }
  }

  D <- data.frame(y = object$y, x = object$x)
  frm <- y ~ 1 / (1 + exp(-eta * (x - phi)))
  start <- c(eta = theta[1], phi = theta[2])
  ctrl <- nls.control(tol = 1.0e-10, minFactor = 1.0e-5, warnOnly = TRUE)

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
    current_par <- mle_asy(object, coefficients(fit_nls))
    current_rss <- rss_fn(current_par)

    if (!is.nan(current_rss) && (current_rss < best_rss)) {
      theta <- current_par
      best_rss <- current_rss
    }
  }

  names(theta) <- NULL

  theta
}

#' 2-parameter logistic fit
#'
#' Fit a 2-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' @param object object of class `logistic2`.
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
fit.logistic2 <- function(object) {
  solution <- find_optimum(object)

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = FALSE,
    estimated = rep(TRUE, 2),
    coefficients = solution$optimum,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - 2,
    fitted.values = logistic2_fn(object$x, solution$optimum),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic2_fit"

  result
}

#' @rdname fit.logistic2
fit_constrained.logistic2 <- function(object) {
  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 2, 3)

  for (i in seq_len(2)) {
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

  theta <- object$lower_bound
  theta[!constraint[, 2]] <- solution$optimum

  estimated <- !constraint[, 2]

  result <- list(
    converged = solution$converged,
    iterations = solution$iterations,
    constrained = !all(constraint[estimated, 1]),
    estimated = estimated,
    coefficients = theta,
    rss = sum(object$stats[, 2] * object$stats[, 4]) + solution$minimum,
    df.residual = length(object$y) - sum(estimated),
    fitted.values = logistic2_fn(object$x, theta),
    weights = object$w
  )

  result$residuals <- object$y - result$fitted.values

  param_names <- c("eta", "phi")

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  class(result) <- "logistic2_fit"

  result
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
#' @param object object of class `logistic2`.
#' @param theta numeric vector with the model parameters.
#' @param sigma estimate of the standard deviation.
#'
#' @return Fisher information matrix evaluated at `theta`.
fisher_info.logistic2 <- function(object, theta, sigma) {
  w <- object$stats[, 2]
  d <- fn(object, object$stats[, 1], theta) - object$stats[, 3]

  gh <- gradient_hessian(object, theta)

  # in case of theta being the maximum likelihood estimator, this gradient G
  # should be zero. We compute it anyway because we likely have rounding errors
  # in our estimate.
  G <- matrix(0, nrow = object$m, ncol = 2)
  G[, 1] <- w * d * gh$G[, 1]
  G[, 2] <- w * d * gh$G[, 2]

  G <- apply(G, 2, sum)

  H <- array(0, dim = c(object$m, 2, 2))

  H[, , 1] <- w * (d * gh$H[, , 1] + gh$G[, 1] * gh$G)
  H[, , 2] <- w * (d * gh$H[, , 2] + gh$G[, 2] * gh$G)

  H <- apply(H, 2:3, sum)

  mu <- fn(object, object$x, theta)
  z <- 3 * sum(object$w * (object$y - mu)^2) / sigma^2 - object$n

  fim <- rbind(cbind(H, -2 * G / sigma), c(-2 * G / sigma, z)) / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' 2-parameter logistic fit
#'
#' Evaluate the variance of the maximum likelihood curve at different predictor
#' values.
#'
#' @param object object of class `logistic2_fit`.
#' @param x numeric vector at which to evaluate the variance.
#'
#' @return Numeric vector with the variances of the maximum likelihood curve.
curve_variance.logistic2_fit <- function(object, x) {
  m <- length(x)

  V <- object$vcov[1:2, 1:2]

  if (any(is.na(V))) {
    return(rep(NA_real_, m))
  }

  eta <- object$coefficients[1]
  phi <- object$coefficients[2]

  b <- exp(-eta * (x - phi))

  f <- 1 + b

  q <- (x - phi) * b
  r <- -eta * b

  s <- 1 / f^2
  t <- q * s
  u <- r * s

  G <- matrix(0, nrow = m, ncol = 2)

  G[, 1] <- t
  G[, 2] <- u

  # When `b` is infinite, gradient shows NaNs
  if (any(is.nan(G))) {
    # these are the limits for b -> Inf
    G[is.nan(G)] <- 0
  }

  variance <- rep(NA_real_, m)

  for (i in seq_len(m)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, ], V), G[i, ]))
  }

  variance
}

#' 2-parameter logistic fit
#'
#' Evaluate the normalized area under the curve (AUC) and area above the curve
#' (AAC).
#'
#' @details
#' The 2-parameter logistic function `f(x; theta)` is defined here as
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(eta, phi)`, `eta` is the steepness of the curve or growth
#' rate (also known as the Hill coefficient), and `phi` is the value of `x` at
#' which the curve is equal to its mid-point, i.e. 1 / 2.
#'
#' The area under the curve (AUC) is simply the integral of `f(x; theta)`
#' between `lower_bound` and `upper_bound` with respect to `x`.
#'
#' When the interval of integration is fixed, since the maximum value the curve
#' can attain is 1, the curve `f(x; theta)` is contained into the rectangle
#' of height 1 and width `upper_bound - lower_bound`. The maximum area the curve
#' can have is obviously `upper_bound - lower_bound`.
#'
#' The normalized AUC is simply `NAUC = AUC / (upper_bound - lower_bound)`. As a
#' consequence, the normalized area above the curve is `NAAC = 1 - NAUC`.
#'
#' Default values of `lower_bound` and `upper_bound` were chosen based on common
#' dose ranges used in the literature. They are also symmetric around zero
#' so that `NAUC` and `NAAC` are equal to `0.5` in the standard logistic model.
#'
#' @param object object of class `logistic2_fit`.
#' @param lower_bound numeric value with the lower bound of the integration
#'   interval.
#' @param upper_bound numeric value with the upper bound of the integration
#'   interval.
#'
#' @return Numeric value with the requested area.
#'
#' @export
nauc.logistic2_fit <- function(object, lower_bound = -10, upper_bound = 10) {
  eta <- object$coefficients[1]
  phi <- object$coefficients[2]

  t1 <- 1 + exp(-eta * (upper_bound - phi))
  t2 <- 1 + exp(-eta * (lower_bound - phi))

  nauc <- 1 + log(t1 / t2) / (eta * (upper_bound - lower_bound))
  names(nauc) <- NULL

  nauc
}

#' @rdname nauc.logistic2_fit
#'
#' @export
naac.logistic2_fit <- function(object, lower_bound = -10, upper_bound = 10) {
  1 - nauc(object, lower_bound, upper_bound)
}
