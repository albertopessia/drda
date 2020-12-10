#' Fit a four parameter logistic function
#'
#' Use a Newton trust-region method to fit a four parameter logistic function to
#' observed data, when each observation is given its own weight.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of observation weights.
#' @param init numeric vector of length 4 with an initial guess of the solution.
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
#'     \item{converged}{boolean value assessing if the optimization algorithm
#'       converged or not.}
#'     \item{iterations}{total number of iterations performed by the
#'       optimization algorithm}
#'     \item{constrained}{boolean value set to `TRUE` if optimization was
#'       constrained.}
#'     \item{coefficients}{maximum likelihood estimates of the model
#'       parameters.}
#'     \item{sigma}{corrected maximum likelihood estimate of the standard
#'       deviation.}
#'     \item{rss}{value of the residual sum of squares.}
#'     \item{loglik}{maximum value of the log-likelihood function.}
#'     \item{df.residual}{residual degrees of freedom.}
#'   }
logistic4_weighted_ntrm_constrained <- function(
  x,
  y,
  w,
  init,
  max_iter,
  constraint,
  lower_bound,
  upper_bound,
  known_param
) {
  constrained <- TRUE

  n <- length(y)
  log_w <- sum(log(w))

  if (any(constraint[, 2])) {
    # there are equality constraints, so we must subset the gradient and hessian
    idx <- which(!constraint[, 2])
    k <- length(idx)

    fn <- logistic4_weighted_rss_fixed(x, y, w, known_param)
    gradient <- logistic4_weighted_rss_gradient_fixed(x, y, w, known_param)
    hessian <- logistic4_weighted_rss_hessian_fixed(x, y, w, known_param)

    result <- if (all(constraint[idx, 1])) {
      # we only have equality constraints, so after fixing the parameters what
      # remains is an unconstrained optimization
      constrained <- FALSE
      ntrm_unconstrained(fn, gradient, hessian, init[idx], max_iter)
    } else {
      ntrm_constrained(
        fn, gradient, hessian, init[idx], max_iter, lower_bound[idx],
        upper_bound[idx]
      )
    }

    coef <- lower_bound
    coef[!constraint[, 2]] <- result$optimum

    v <- logistic4_variance(result$minimum, n, k)
    loglik <- logistic4_loglik(v, n, k, log_w = log_w)

    list(
      converged = result$converged,
      iterations = result$iterations,
      constrained = constrained,
      coefficients = coef,
      sigma = sqrt(v),
      rss = result$minimum,
      loglik = loglik,
      df.residual = n - k
    )
  } else {
    fn <- logistic4_weighted_rss(x, y, w)
    gradient <- logistic4_weighted_rss_gradient(x, y, w)
    hessian <- logistic4_weighted_rss_hessian(x, y, w)
    result <- ntrm_constrained(
      fn, gradient, hessian, init, max_iter, lower_bound, upper_bound
    )

    # we estimated all 4 parameters (k = 4)
    v <- logistic4_variance(result$minimum, n)
    loglik <- logistic4_loglik(v, n, log_w = log_w)

    list(
      converged = result$converged,
      iterations = result$iterations,
      constrained = constrained,
      coefficients = result$optimum,
      sigma = sqrt(v),
      rss = result$minimum,
      loglik = loglik,
      df.residual = n - 4
    )
  }
}

#' @rdname logistic4_weighted_ntrm_constrained
logistic4_weighted_ntrm_unconstrained <- function(
  x,
  y,
  w,
  init,
  max_iter
) {
  fn <- logistic4_weighted_rss(x, y, w)
  gradient <- logistic4_weighted_rss_gradient(x, y, w)
  hessian <- logistic4_weighted_rss_hessian(x, y, w)
  result <- ntrm_unconstrained(fn, gradient, hessian, init, max_iter)

  n <- length(y)

  # we estimated 4 parameters => k - 4 degrees of freedom
  v <- logistic4_variance(result$minimum, n)
  loglik <- logistic4_loglik(v, n, log_w = sum(log(w)))

  list(
    converged = result$converged,
    iterations = result$iterations,
    constrained = FALSE,
    coefficients = result$optimum,
    sigma = sqrt(v),
    rss = result$minimum,
    loglik = loglik,
    df.residual = n - 4
  )
}

#' 4-parameter logistic fit
#'
#' Fit a 4-parameter logistic function to observed data with a Maximum
#' Likelihood approach.
#'
#' @details
#' The 4-parameter logistic function is defined in this package as
#'
#' `f(x; theta) = alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' where `theta = c(alpha, beta, eta, phi)`, `alpha` is the lower horizontal
#' asymptote, `beta` is the upper horizontal asymptote, `eta` is the steepness
#' of the curve or growth rate (also known as the Hill coefficient), and `phi`
#' is the value of `x` at which the curve is equal to its mid-point.
#'
#' The model assumes the observed data to be normally distributed around the
#' logistic function, with a possibly unequal standard deviation at each and
#' every predictor `x[k]`, i.e.
#'
#' `y[i] ~ N(f(x[k]; theta), sigma^2 / w[i])`
#'
#' for any predictor `x[k]` and for any statistical unit `i` observed at
#' `x[k]`.
#'
#' @param x numeric vector representing the fixed predictor variable.
#' @param y numeric vector of observed values.
#' @param w numeric vector of observation weights.
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
#'     \item{converged}{boolean value assessing if the optimization algorithm
#'       converged or not.}
#'     \item{iterations}{total number of iterations performed by the
#'       optimization algorithm}
#'     \item{constrained}{boolean value set to `TRUE` if optimization was
#'       constrained.}
#'     \item{coefficients}{maximum likelihood estimates of the model
#'       parameters.}
#'     \item{sigma}{corrected maximum likelihood estimate of the standard
#'       deviation.}
#'     \item{rss}{value of the residual sum of squares.}
#'     \item{loglik}{maximum value of the log-likelihood function.}
#'     \item{df.residual}{residual degrees of freedom.}
#'     \item{estimated}{boolean vector indicating which parameters were
#'       estimated from the data.}
#'     \item{fitted.values}{the fitted mean values.}
#'     \item{residuals}{the residuals, that is response minus fitted values.}
#'     \item{weights}{the specified weights.}
#'     \item{fisher_info}{Fisher information matrix evaluated at the maximum
#'       likelihood estimator.}
#'   }
logistic4_weighted_fit_constrained <- function(
  x,
  y,
  w,
  start,
  max_iter,
  lower_bound,
  upper_bound
) {
  if (!is.numeric(lower_bound) || !is.null(dim(lower_bound))) {
    stop("'lower_bound' must be a numeric vector", call. = FALSE)
  }

  if (length(lower_bound) != 4) {
    stop("'lower_bound' must be of length 4", call. = FALSE)
  }

  if (!is.numeric(upper_bound) || !is.null(dim(upper_bound))) {
    stop("'upper_bound' must be a numeric vector", call. = FALSE)
  }

  if (length(upper_bound) != 4) {
    stop("'upper_bound' must be of length 4", call. = FALSE)
  }

  if (any(lower_bound > upper_bound)) {
    stop("'lower_bound' cannot be larger than 'upper_bound'", call. = FALSE)
  }

  # when the lower bound is fixed, we don't care about the constraints on the
  # upper bound
  if ((lower_bound[1] != upper_bound[1]) && (lower_bound[2] < lower_bound[1])) {
    stop("'lower_bound[2]' is smaller than 'lower_bound[1]'", call. = FALSE)
  }

  # when the upper bound is fixed, we don't care about the constraints on the
  # lower bound
  if ((lower_bound[2] != upper_bound[2]) && (upper_bound[1] > upper_bound[2])) {
    stop("'upper_bound[1]' is greater than 'upper_bound[2]'", call. = FALSE)
  }

  # process constraints
  # first column is for unconstrained parameters
  # second column is for equality parameters
  # third column is for inequality parameters
  constraint <- matrix(FALSE, 4, 3)

  for (i in 1:4) {
    lb_is_inf <- is.infinite(lower_bound[i])
    ub_is_inf <- is.infinite(upper_bound[i])

    if (lb_is_inf && (lower_bound[i] > 0)) {
      stop("'lower_bound' cannot be equal to infinity", call. = FALSE)
    }

    if (ub_is_inf && (upper_bound[i] < 0)) {
      stop("'upper_bound' cannot be equal to -infinity", call. = FALSE)
    }

    if (lower_bound[i] == upper_bound[i]) {
      constraint[i, 2] <- TRUE
    } else if (!lb_is_inf || !ub_is_inf) {
      constraint[i, 3] <- TRUE
    } else {
      constraint[i, 1] <- TRUE
    }
  }

  known_param <- ifelse(constraint[, 2], lower_bound, NA_real_)

  n <- length(y)

  w_complete <- w
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
          coefficients = rep(NA_real_, 4),
          sigma = NA_real_,
          loglik = NA_real_,
          df.residual = 0L,
          fitted.values = numeric(0),
          residuals = numeric(0)
        )
      )
    }
  }

  stats <- logistic4_suff_stat(x, y)

  init <- if (!is.null(start)) {
    if (!is.numeric(start) || !is.null(dim(start))) {
      stop("'start' must be a numeric vector", call. = FALSE)
    }

    if (length(start) != 4) {
      stop("'start' must be of length 4", call. = FALSE)
    }

    if (any(is.infinite(start) | is.na(start))) {
      stop("'start' must be finite", call. = FALSE)
    }

    # equality constraints have the priority over the provided starting values
    ifelse(is.na(known_param), start, known_param)
  } else {
    ifelse(is.na(known_param), logistic4_init(stats), known_param)
  }

  result <- logistic4_weighted_ntrm_constrained(
    x, y, w, init, max_iter, constraint, lower_bound, upper_bound, known_param
  )

  if (result$coefficients[1] > result$coefficients[2]) {
    # logistic function has the symmetric property
    # f(x; alpha, beta, eta, phi) = f(x; beta, alpha, -eta, phi)
    #
    # we want to force our interpretation that the first parameter must always
    # be the lower horizontal asymptote
    tmp <- result$coefficients[1]
    result$coefficients[1] <- result$coefficients[2]
    result$coefficients[2] <- tmp
    result$coefficients[3] <- -result$coefficients[3]
  }

  result$estimated <- !constraint[, 2]

  param_names <- c(
    "minimum", "maximum", "growth_rate", "x_midpoint"
  )

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result$fitted.values <- logistic4_function(x, result$coefficients)
  result$residuals <- y - result$fitted.values

  result$weights <- w_complete

  # evaluate the gradient at the maximum likelihood estimator
  g <- logistic4_gradient(x, result$coefficients)

  tmp <- array(0, dim = c(length(x), 4, 4))
  tmp[, , 1] <- w * g[, 1] * g
  tmp[, , 2] <- w * g[, 2] * g
  tmp[, , 3] <- w * g[, 3] * g
  tmp[, , 4] <- w * g[, 4] * g

  fim <- matrix(0, nrow = 5, ncol = 5)
  fim[1:4, 1:4] <- apply(tmp, 2:3, sum)
  fim[5, 5] <- 3 * sum(w) - n
  fim <- fim / result$sigma^2

  # if there are equality constraints, we set the corresponding elements to zero
  if (any(constraint[, 2])) {
    fim[constraint[, 2], ] <- 0
    fim[, constraint[, 2]] <- 0
  }

  result$fisher.info <- fim

  result
}

#' @rdname logistic4_weighted_fit_constrained
logistic4_weighted_fit_unconstrained <- function(
  x,
  y,
  w,
  start,
  max_iter
) {
  n <- length(y)

  w_complete <- w
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
          coefficients = rep(NA_real_, 4),
          sigma = NA_real_,
          loglik = NA_real_,
          df.residual = 0L,
          fitted.values = numeric(0),
          residuals = numeric(0)
        )
      )
    }
  }

  stats <- logistic4_suff_stat(x, y)

  init <- if (!is.null(start)) {
    if (!is.numeric(start) || !is.null(dim(start))) {
      stop("'start' must be a numeric vector", call. = FALSE)
    }

    if (length(start) != 4) {
      stop("'start' must be of length 4", call. = FALSE)
    }

    if (any(is.infinite(start) | is.na(start))) {
      stop("'start' must be finite", call. = FALSE)
    }

    start
  } else {
    logistic4_init(stats)
  }

  result <- logistic4_weighted_ntrm_unconstrained(x, y, w, init, max_iter)

  if (result$coefficients[1] > result$coefficients[2]) {
    tmp <- result$coefficients[1]
    result$coefficients[1] <- result$coefficients[2]
    result$coefficients[2] <- tmp
    result$coefficients[3] <- -result$coefficients[3]
  }

  result$estimated <- rep(TRUE, 4)

  param_names <- c(
    "minimum", "maximum", "growth_rate", "x_midpoint"
  )

  names(result$coefficients) <- param_names
  names(result$estimated) <- param_names

  result$fitted.values <- logistic4_function(x, result$coefficients)
  result$residuals <- y - result$fitted.values

  result$weights <- w_complete

  g <- logistic4_gradient(x, result$coefficients)

  tmp <- array(0, dim = c(length(x), 4, 4))
  tmp[, , 1] <- w * g[, 1] * g
  tmp[, , 2] <- w * g[, 2] * g
  tmp[, , 3] <- w * g[, 3] * g
  tmp[, , 4] <- w * g[, 4] * g

  fim <- matrix(0, nrow = 5, ncol = 5)
  fim[1:4, 1:4] <- apply(tmp, 2:3, sum)
  fim[5, 5] <- 3 * sum(w) - n
  fim <- fim / result$sigma^2

  result$fisher.info <- fim

  result
}
