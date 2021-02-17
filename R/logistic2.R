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

  if (!is.null(lower_bound)) {
    if (length(lower_bound) != 2) {
      stop("'lower_bound' must be of length 2", call. = FALSE)
    }
  }

  if (!is.null(upper_bound)) {
    if (length(upper_bound) != 2) {
      stop("''upper_bound' must be of length 2", call. = FALSE)
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

    object$lower_bound <- if (!is.null(lower_bound)) {
      lower_bound
    } else {
      rep(-Inf, 2)
    }

    object$upper_bound <- if (!is.null(upper_bound)) {
      upper_bound
    } else {
      rep(Inf, 2)
    }
  }

  object
}

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
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
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
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
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
#' The 2-parameter logistic function `f(x; theta)` is defined in this package as
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
#' @importFrom stats lm median
init.logistic2 <- function(object) {
  m <- object$m
  stats <- object$stats

  linear_fit <- summary(lm(stats[, 3] ~ stats[, 1], weights = stats[, 3]))
  linear_coef <- linear_fit$coefficients

  if (linear_coef[2, 4] > 0.2) {
    # we are in big problems as a flat horizontal line is likely the best model
    return(c(-10, stats[m, 1] + 1000))
  }

  delta <- mean(diff(stats[, 1]))

  rss_fn <- rss(object)

  eta_set <- seq(-2, -0.01, length.out = 15)
  phi_set <- seq(
    stats[1, 1] - 0.5 * delta, stats[m, 1] + 0.5 * delta, length.out = 15
  )

  theta <- c(-1, stats[which.min(abs(stats[, 3] - median(stats[, 3]))), 1])

  if (linear_coef[2, 1] > 0) {
    eta_set <- seq(0.01, 2, length.out = 15)
    theta[1] <- 1
  }

  best_rss <- rss_fn(theta)

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

  names(theta) <- NULL

  theta
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
  eta <- theta[1]
  phi <- theta[2]

  b <- exp(-eta * (object$stats[, 1] - phi))

  f <- 1 + b

  q <- (object$stats[, 1] - phi) * b
  r <- -eta * b

  gradient <- matrix(1, nrow = object$m, ncol = 2)

  gradient[, 1] <- q / f^2
  gradient[, 2] <- r / f^2

  tmp <- array(0, dim = c(object$m, 2, 2))
  tmp[, , 1] <- object$stats[, 2] * gradient[, 1] * gradient
  tmp[, , 2] <- object$stats[, 2] * gradient[, 2] * gradient

  fim <- matrix(0, nrow = 3, ncol = 3)
  fim[1:2, 1:2] <- apply(tmp, 2:3, sum)
  fim[3, 3] <- object$n - 3
  fim <- fim / sigma^2

  lab <- c(names(theta), "sigma")
  rownames(fim) <- lab
  colnames(fim) <- lab

  fim
}

#' 2-parameter logistic fit
#'
#' Evaluate the normalized area under the curve (AUC) and area above the curve
#' (AAC).
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
