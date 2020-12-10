#' Residual sum of squares
#'
#' Evaluate the residual sum of squares (RSS) of a 4-parameter logistic model.
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
#' @param stats matrix of sufficient statistics. First column are the unique
#'   elements of `x`, second column are the corresponding sample sizes,
#'   third column are the sample means, while the fourth column are the
#'   (uncorrected) sample variances.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any. Vector is of the form `c(alpha, beta, eta, phi)`
#'   where `NA` represents an unknown parameter.
#'
#' @return Function handle `f(theta)` to evaluate the RSS associated to a
#'   particular parameter choice `theta`.
logistic4_rss <- function(
  stats
) {
  function(z) {
    mu <- logistic4_function(stats[, 1], z)
    sum(stats[, 2] * (stats[, 4] + (mu - stats[, 3])^2))
  }
}

#' @rdname logistic4_rss
logistic4_rss_fixed <- function(
  stats,
  known_param
) {
  function(z) {
    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    mu <- logistic4_function(stats[, 1], theta)
    sum(stats[, 2] * (stats[, 4] + (mu - stats[, 3])^2))
  }
}

#' Residual sum of squares
#'
#' Evaluate the gradient of the residual sum of squares (RSS) of a 4-parameter
#' logistic model.
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
#' @param stats matrix of sufficient statistics. First column are the unique
#'   elements of `x`, second column are the corresponding sample sizes,
#'   third column are the sample means, while the fourth column are the
#'   (uncorrected) sample variances.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any. Vector is of the form `c(alpha, beta, eta, phi)`
#'   where `NA` represents an unknown parameter.
#'
#' @return Function handle `f(theta)` to evaluate the gradient of the RSS
#'   associated to a particular parameter choice `theta`.
logistic4_rss_gradient <- function(
  stats
) {
  function(z) {
    n <- stats[, 2]

    d <- z[2] - z[1]
    u <- stats[, 1] - z[4]
    v <- -z[3] * u

    log_denom <- log1p(exp(v))

    t1 <- exp(v - log_denom)
    t2 <- exp(-log_denom)
    t3 <- t1 * t2
    t4 <- d * t3

    predicted <- z[1] + d * t2
    r <- predicted - stats[, 3]

    G <- matrix(0, nrow = nrow(stats), ncol = 4)
    G[, 1] <- n * t1 * r
    G[, 2] <- n * t2 * r
    G[, 3] <- n * u * t4 * r
    G[, 4] <- -n * z[3] * t4 * r

    apply(G, 2, sum)
  }
}

#' @rdname logistic4_rss_gradient
logistic4_rss_gradient_fixed <- function(
  stats,
  known_param
) {
  function(z) {
    n <- stats[, 2]

    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    d <- theta[2] - theta[1]
    u <- stats[, 1] - theta[4]
    v <- -theta[3] * u

    log_denom <- log1p(exp(v))

    t1 <- exp(v - log_denom)
    t2 <- exp(-log_denom)
    t3 <- t1 * t2
    t4 <- d * t3

    predicted <- theta[1] + d * t2
    r <- predicted - stats[, 3]

    G <- matrix(0, nrow = nrow(stats), ncol = 4)
    G[, 1] <- n * t1 * r
    G[, 2] <- n * t2 * r
    G[, 3] <- n * u * t4 * r
    G[, 4] <- -n * theta[3] * t4 * r

    apply(G[, idx, drop = FALSE], 2, sum)
  }
}

#' Residual sum of squares
#'
#' Evaluate the Hessian matrix of the residual sum of squares (RSS) of a
#' 4-parameter logistic model.
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
#' @param stats matrix of sufficient statistics. First column are the unique
#'   elements of `x`, second column are the corresponding sample sizes,
#'   third column are the sample means, while the fourth column are the
#'   (uncorrected) sample variances.
#' @param known_param numeric vector with the known fixed values of the model
#'   parameters, if any. Vector is of the form `c(alpha, beta, eta, phi)`
#'   where `NA` represents an unknown parameter.
#'
#' @return Function handle `f(theta)` to evaluate the Hessian matrix of the RSS
#'   associated to a particular parameter choice `theta`.
logistic4_rss_hessian <- function(
  stats
) {
  function(z) {
    n <- stats[, 2]

    d <- z[2] - z[1]
    u <- stats[, 1] - z[4]
    v <- -z[3] * u

    log_denom <- log1p(exp(v))

    t1 <- exp(v - log_denom)
    t2 <- exp(-log_denom)
    t3 <- t1 * t2
    t4 <- d * t3

    predicted <- z[1] + d * t2
    r <- predicted - stats[, 3]

    t5 <- d * t1 - r
    t6 <- d * t2 + r

    t7 <- 2 * t1 - 1
    t8 <- t4 + t7 * r

    H <- array(0, dim = c(nrow(stats), 4, 4))

    H[, 1, 1] <- n * t1^2
    H[, 2, 1] <- n * t3
    H[, 3, 1] <- n * u * t3 * t5
    H[, 4, 1] <- -n * z[3] * t3 * t5

    H[, 1, 2] <- H[, 2, 1]
    H[, 2, 2] <- n * t2^2
    H[, 3, 2] <- n * u * t3 * t6
    H[, 4, 2] <- -n * z[3] * t3 * t6

    H[, 1, 3] <- H[, 3, 1]
    H[, 2, 3] <- H[, 3, 2]
    H[, 3, 3] <- n * u^2 * t4 * t8
    H[, 4, 3] <- -n * t4 * (r - v * t8)

    H[, 1, 4] <- H[, 4, 1]
    H[, 2, 4] <- H[, 4, 2]
    H[, 3, 4] <- H[, 4, 3]
    H[, 4, 4] <- n * z[3]^2 * t4 * t8

    apply(H, 2:3, sum)
  }
}

#' @rdname logistic4_rss_hessian
logistic4_rss_hessian_fixed <- function(
  stats,
  known_param
) {
  function(z) {
    n <- stats[, 2]

    idx <- is.na(known_param)

    theta <- rep(0, 4)
    theta[ idx] <- z
    theta[!idx] <- known_param[!idx]

    d <- theta[2] - theta[1]
    u <- stats[, 1] - theta[4]
    v <- -theta[3] * u

    log_denom <- log1p(exp(v))

    t1 <- exp(v - log_denom)
    t2 <- exp(-log_denom)
    t3 <- t1 * t2
    t4 <- d * t3

    predicted <- theta[1] + d * t2
    r <- predicted - stats[, 3]

    t5 <- d * t1 - r
    t6 <- d * t2 + r

    t7 <- 2 * t1 - 1
    t8 <- t4 + t7 * r

    H <- array(0, dim = c(nrow(stats), 4, 4))

    H[, 1, 1] <- n * t1^2
    H[, 2, 1] <- n * t3
    H[, 3, 1] <- n * u * t3 * t5
    H[, 4, 1] <- -n * theta[3] * t3 * t5

    H[, 1, 2] <- H[, 2, 1]
    H[, 2, 2] <- n * t2^2
    H[, 3, 2] <- n * u * t3 * t6
    H[, 4, 2] <- -n * theta[3] * t3 * t6

    H[, 1, 3] <- H[, 3, 1]
    H[, 2, 3] <- H[, 3, 2]
    H[, 3, 3] <- n * u^2 * t4 * t8
    H[, 4, 3] <- -n * t4 * (r * (1 - v * t7) - v * t4)

    H[, 1, 4] <- H[, 4, 1]
    H[, 2, 4] <- H[, 4, 2]
    H[, 3, 4] <- H[, 4, 3]
    H[, 4, 4] <- n * theta[3]^2 * t4 * t8

    apply(H[, idx, idx, drop = FALSE], 2:3, sum)
  }
}
