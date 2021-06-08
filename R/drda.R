#' Fit non-linear growth curves
#'
#' Use the Newton's with a trust-region method to fit non-linear growth curves
#' to observed data.
#'
#' @param formula an object of class \code{link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the model to be
#'   fitted. Currently supports only formulas of the type `y ~ x`.
#' @param data an optional data frame, list or environment (or object coercible
#'   by \code{\link[base]{as.data.frame}} to a data frame) containing the
#'   variables in the model. If not found in `data`, the variables are taken
#'   from `environment(formula)`, typically the environment from which `drda`
#'   is called.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'   process. If provided, weighted least squares is used with weights `weights`
#'   (that is, minimizing `sum(weights * residuals^2)`), otherwise ordinary
#'   least squares is used.
#' @param na.action a function which indicates what should happen when the data
#'   contain `NA`s. The default is set by the `na.action` setting of
#'   \code{\link[base]{options}}, and is \code{\link[stats]{na.fail}} if that is
#'   unset. The ‘factory-fresh’ default is \code{na.omit}. Another
#'   possible value is `NULL`, no action. Value \code{na.exclude} can be useful.
#' @param mean_function the model to be fitted. See `details` for available
#'   models.
#' @param is_log a logical value indicating whether the predictor variable `x`
#'   is already log-transformed. Default to `TRUE`. Set to `FALSE` if `x` is
#'   on its natural scale, i.e. strictly positive.
#' @param lower_bound numeric vector with the minimum admissible values of the
#'   parameters. Use `-Inf` to specify an unbounded parameter.
#' @param upper_bound numeric vector with the maximum admissible values of the
#'   parameters. Use `Inf` to specify an unbounded parameter.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#'
#' @details
#'
#' ## Available models
#'
#' ### Generalized logistic function
#'
#' The most general model in this package is the generalized logistic function
#' selected by setting `mean_function = "logistic6"`. It is defined in this
#' package as the 6-parameter function
#'
#' `alpha + (beta - alpha) / (xi + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `eta != 0`, `nu > 0`, and `xi > 0`. Although `beta` can be any real
#' value, we use the convention `beta > alpha` to avoid identifiability
#' problems: when `beta < alpha` it is always possible to adjust the other
#' parameters to obtain the same exact curve. When `beta > alpha` and `eta > 0`
#' the curve is monotonically increasing. If `beta > alpha` and `eta < 0` the
#' curve is monotonically decreasing.
#'
#' Parameter `alpha` represents the lower horizontal asymptote of the curve.
#' Parameter `beta` is related to the upper horizontal asymptote of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the value of the function at `x = 0`.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#' Parameter `xi` affects the value of the upper asymptote.
#'
#' **Note**: the 6-parameter logistic function is non-identifiable from data and
#' should not be used in real applications. It is available only for theoretical
#' research convenience.
#'
#' ### 5-parameter logistic function
#'
#' The 5-parameter logistic function can be selected by choosing
#' `mean_function = "logistic5"`. The function is obtained by setting `xi = 1`
#' in the generalized logistic function, that is
#'
#' `alpha + (beta - alpha) / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' Parameter `alpha` represents the lower horizontal asymptote of the curve.
#' Parameter `beta` represents the upper horizontal asymptote of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the value of the function at `x = 0`.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#'
#' ### 4-parameter logistic function
#'
#' The 4-parameter logistic function is the default model of `drda`. It can be
#' explicitly selected by choosing `mean_function = "logistic4"`. The function
#' is obtained by setting `xi = 1` and `nu = 1` in the generalized logistic
#' function, that is
#'
#' `alpha + (beta - alpha) / (1 + exp(-eta * (x - phi)))`
#'
#' Parameter `alpha` represents the lower horizontal asymptote of the curve.
#' Parameter `beta` represents the upper horizontal asymptote of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; alpha, beta, eta, phi) = (alpha + beta) / 2`.
#'
#' ### 2-parameter logistic function
#'
#' The 2-parameter logistic function can be selected by choosing
#' `mean_function = "logistic2"`. The function is obtained by setting `xi = 1`,
#' `nu = 1`, `beta = 1`, and `alpha = 0` in the generalized logistic function,
#' that is
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; eta, phi) = 1 / 2`.
#'
#' ### Gompertz function
#'
#' The Gompertz function is the limit for `nu -> 0` of the 5-parameter logistic
#' function. It can be selected by choosing `mean_function = "gompertz"`. The
#' function is defined in this package as
#'
#' `alpha + (beta - alpha) * exp(-exp(-eta * (x - phi)))`
#'
#' where `eta != 0`.
#'
#' Parameter `alpha` represents the lower horizontal asymptote of the curve.
#' Parameter `beta` represents the upper horizontal asymptote of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the value of the function at `x = 0`.
#'
#' ### Constrained optimization
#'
#' It is possible to search for the maximum likelihood estimates within
#' pre-specified interval regions. Since the upper horizontal asymptote `beta`
#' must be greater than the lower horizontal asymptote `alpha`, intervals are
#' adjusted to satisfy this constraint.
#'
#' *Note*: Hypothesis testing is not available for constrained estimates
#' because asymptotic approximations might not be valid
#'
#' @return An object of class `drda` and `model_fit`, where `model` is the
#' chosen mean function. It is a list containing the following components:
#'   \describe{
#'     \item{converged}{boolean value assessing if the optimization algorithm
#'       converged or not.}
#'     \item{iterations}{total number of iterations performed by the
#'       optimization algorithm}
#'     \item{constrained}{boolean value set to `TRUE` if optimization was
#'       constrained.}
#'     \item{estimated}{boolean vector indicating which parameters were
#'       estimated from the data.}
#'     \item{coefficients}{maximum likelihood estimates of the model
#'       parameters.}
#'     \item{rss}{minimum value (found) of the residual sum of squares.}
#'     \item{df.residuals}{residual degrees of freedom.}
#'     \item{fitted.values}{fitted mean values.}
#'     \item{residuals}{residuals, that is response minus fitted values.}
#'     \item{weights}{(only for weighted fits) the specified weights.}
#'     \item{mean_function}{model that was used for fitting.}
#'     \item{n}{effective sample size.}
#'     \item{sigma}{corrected maximum likelihood estimate of the standard
#'       deviation.}
#'     \item{loglik}{maximum value (found) of the log-likelihood function.}
#'     \item{fisher.info}{observed Fisher information matrix evaluated at the
#'       maximum likelihood estimator.}
#'     \item{vcov}{approximate variance-covariance matrix of the model
#'       parameters.}
#'     \item{call}{the matched call.}
#'     \item{terms}{the \code{\link[stats]{terms}} object used.}
#'     \item{model}{the model frame used.}
#'     \item{na.action}{(where relevant) information returned by
#'       \code{\link[stats]{model.frame}} on the special handling of `NA`s.}
#'     \item{is_log}{boolean value. It is `TRUE` if the predictor variable was
#'       given on the log scale.}
#'   }
#'
#' @importFrom stats model.frame model.matrix model.response model.weights terms
#'
#' @export
drda <- function(
  formula, data, subset, weights, na.action, mean_function = "logistic4",
  is_log = TRUE, lower_bound = NULL, upper_bound = NULL, start = NULL,
  max_iter = 10000
) {
  # first, we expand the call to this function
  model_frame <- match.call(expand.dots = FALSE)

  # select arguments to pass to the stats::model.frame function
  arg_name <- names(model_frame)
  arg_idx <- c(
    1,
    match(c("formula", "data", "subset", "weights", "na.action"), arg_name, 0)
  )

  # the first element of `model_frame` is now the name of this function.
  # Instead, overwrite it with `stats::model.frame` and evaluate the call
  model_frame <- model_frame[arg_idx]
  model_frame[[1]] <- quote(model.frame)
  model_frame <- eval(model_frame, parent.frame())

  model_terms <- attr(model_frame, "terms")
  model_matrix <- model.matrix(model_terms, model_frame)

  # we expect a `y ~ x` formula (we don't care about the actual names used)
  # by default we set our `x` to be the first non-intercept variable in the
  # formula, i.e. if the user gives `response ~ a + b + c + d` we simply select
  # `a` as our `x` variable
  x <- model_matrix[, attr(model_matrix, "assign") == 1, drop = TRUE]
  n <- length(x)

  if (n == 0) {
    stop("0 (non-NA) cases", call. = FALSE)
  }

  y <- model.response(model_frame, "numeric")

  if (is.matrix(y)) {
    stop("response variable must be a vector", call. = FALSE)
  }

  if (length(y) != n) {
    stop("incompatible dimensions", call. = FALSE)
  }

  w <- as.vector(model.weights(model_frame))

  if (is.null(w)) {
    w <- rep(1, length(y))
  } else {
    if (!is.numeric(w)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    }

    if (length(w) != n) {
      stop("incompatible dimensions", call. = FALSE)
    }

    if (any(w < 0 | is.na(w))) {
      stop("missing or negative weights not allowed", call. = FALSE)
    }
  }

  w_zero <- w == 0

  if (any(w_zero)) {
    w_positive <- !w_zero

    x <- x[w_positive]
    y <- y[w_positive]
    w <- w[w_positive]

    if (length(y) == 0) {
      stop("weights cannot be all zero", call. = FALSE)
    }
  }

  if (!is_log) {
    if (any(x <= 0)) {
      stop(
        paste(
          "predictor variable `x` is not strictly positive",
          "and cannot be log-transformed"
        ),
        call. = FALSE
      )
    }

    x <- log(x)
  }

  max_iter <- ceiling(max_iter[1])

  if (max_iter <= 0) {
    stop("maximum number of iterations must be positive", call. = FALSE)
  }

  if (!is.null(lower_bound)) {
    if (!is.numeric(lower_bound) || !is.null(dim(lower_bound))) {
      stop("'lower_bound' must be a numeric vector", call. = FALSE)
    }
  }

  if (!is.null(upper_bound)) {
    if (!is.numeric(upper_bound) || !is.null(dim(upper_bound))) {
      stop("'upper_bound' must be a numeric vector", call. = FALSE)
    }
  }

  if (!is.null(lower_bound) && !is.null(upper_bound)) {
    if (length(lower_bound) != length(upper_bound)) {
      stop(
        "'lower_bound' and 'upper_bound' must have the same length",
        call. = FALSE
      )
    }

    if (any(lower_bound > upper_bound)) {
      stop("'lower_bound' cannot be larger than 'upper_bound'", call. = FALSE)
    }

    if (any(is.infinite(lower_bound) & (lower_bound > 0))) {
      stop("'lower_bound' cannot be equal to infinity", call. = FALSE)
    }

    if (any(is.infinite(upper_bound) & (upper_bound < 0))) {
      stop("'upper_bound' cannot be equal to -infinity", call. = FALSE)
    }
  }

  if (!is.null(start)) {
    if (!is.numeric(start) || !is.null(dim(start))) {
      stop("'start' must be a numeric vector", call. = FALSE)
    }

    if (any(is.infinite(start) | is.na(start))) {
      stop("'start' must be finite", call. = FALSE)
    }
  }

  object <- switch(mean_function,
    logistic2 = logistic2_new(
      x, y, w, start, max_iter, lower_bound, upper_bound
    ),
    logistic4 = logistic4_new(
      x, y, w, start, max_iter, lower_bound, upper_bound
    ),
    logistic5 = logistic5_new(
      x, y, w, start, max_iter, lower_bound, upper_bound
    ),
    logistic6 = logistic6_new(
      x, y, w, start, max_iter, lower_bound, upper_bound
    ),
    gompertz = gompertz_new(
      x, y, w, start, max_iter, lower_bound, upper_bound
    )
  )

  if (is.null(object)) {
    stop(
      "chosen 'mean_function' is wrongly typed or not yet available",
      call. = FALSE
    )
  }

  result <- if (!object$constrained) {
    fit(object)
  } else {
    fit_constrained(object)
  }

  result$mean_function <- mean_function

  log_w <- sum(log(result$weights))
  v <- variance_normal(result$rss, result$df.residual)

  result$n <- object$n
  result$sigma <- sqrt(v)
  result$loglik <- loglik_normal(result$rss, object$n, log_w)

  result$fisher.info <- fisher_info(object, result$coefficients, result$sigma)
  result$vcov <- if (!object$constrained) {
    approx_vcov(result$fisher.info)
  } else {
    k <- nrow(result$fisher.info)
    idx <- which(!result$estimated)
    vcov <- matrix(NA, nrow = k, ncol = k)
    vcov[-idx, -idx] <- approx_vcov(result$fisher.info[-idx, -idx])
    vcov
  }

  result$call <- match.call()
  result$terms <- model_terms
  result$model <- model_frame
  result$na.action <- attr(model_frame, "na.action")
  result$is_log <- is_log

  class(result) <- c("drda", class(result))

  result
}

#' @importFrom stats anova pchisq var
#'
#' @export
anova.drda <- function(object, ...) {
  ## check for multiple objects
  dotargs <- list(...)

  named <- if (is.null(names(dotargs))) {
    rep_len(FALSE, length(dotargs))
  } else {
    names(dotargs) != ""
  }

  if (any(named)) {
    warning(
      "the following arguments to 'anova.drda' are invalid and dropped: ",
      paste(deparse(dotargs[named]), collapse = ", ")
    )
  }

  dotargs <- dotargs[!named]

  is_drda <- vapply(dotargs, function(x) inherits(x, "drda"), FALSE)
  dotargs <- dotargs[is_drda]

  if (length(dotargs)) {
    return(anova.drdalist(c(list(object), dotargs)))
  }

  if (object$constrained) {
    # solution to the constrained problem is unlikely the maximum likelihood
    # estimator, therefore the asymptotic approximation might not hold
    stop("hypothesis testing is not available for constrained optimization")
  }

  y <- object$model[, 1]
  x <- object$model[, 2]
  w <- object$weights

  idx <- !is.na(y) & !is.na(x) & !is.na(w) & !(w == 0)

  if (sum(idx) != length(y)) {
    y <- y[idx]
    x <- x[idx]
    w <- w[idx]
  }

  n <- length(y)
  log_n <- log(n)
  log_w <- sum(log(w))
  k <- sum(object$estimated)

  l <- if (k >= 5) {
    # we compare the full model against a flat horizontal line
    2
  } else {
    # we compare the estimated model against the baseline and the full model
    3
  }

  deviance_df <- rep(-1, l)
  deviance_value <- rep(-1, l)
  loglik <- rep(-1, l)

  # constant model: horizontal line
  weighted_mean <- sum(w * y) / sum(w)
  deviance_df[1] <- n - 1
  deviance_value[1] <- sum(w * (y - weighted_mean)^2)

  # fitted model
  deviance_df[2] <- object$df.residual
  deviance_value[2] <- object$rss

  if (k < 5) {
    # at least a parameter was considered fixed, so we now fit the full model
    fit <- drda(y ~ x, weights = w, mean_function = "logistic5")
    deviance_df[3] <- fit$df.residual
    deviance_value[3] <- fit$rss
  }

  loglik <- loglik_normal(deviance_value, n, log_w)

  df <- n - deviance_df

  aic <- 2 * (df - loglik)
  bic <- log_n * df - 2 * loglik
  lrt <- -2 * (loglik - loglik[l])
  lrt[l] <- NA_real_

  table <- data.frame(deviance_df, deviance_value, df, aic, bic, lrt)
  table$pvalue <- pchisq(lrt, -diff(deviance_df), lower.tail = FALSE)

  rownames(table)[1:2] <- c("Constant model", "Estimated model")
  if (k < 5) {
    rownames(table)[3] <- "Full model (logistic5)"
  }
  colnames(table) <- c(
    "Resid. Df", "Resid. Dev", "Df", "AIC", "BIC", "LRT", "p-value"
  )

  model <- switch(object$mean_function,
    logistic2 = "2-parameter logistic",
    logistic4 = "4-parameter logistic",
    logistic5 = "5-parameter logistic",
    logistic6 = "6-parameter logistic"
  )

  title <- paste("Analysis of Deviance Table\n\nModel: ", model, "\n", sep = "")

  structure(table, heading = title, class = c("anova", "data.frame"))
}

#' @importFrom stats anova pchisq var
#'
#' @export
anova.drdalist <- function(object, ...) {
  n_models <- length(object)

  if (n_models == 1) {
    return(anova.drda(object[[1L]]))
  }

  is_constrained <- any(vapply(object, function(x) x$constrained, FALSE))

  if (any(is_constrained)) {
    # solution to the constrained problem is unlikely the maximum likelihood
    # estimator, therefore the asymptotic approximation might not hold
    stop("hypothesis testing is not available for constrained optimization")
  }

  n_residuals <- vapply(object, function(x) length(x$residuals), 0)

  if (any(n_residuals != n_residuals[1L])) {
    stop(
      "models were not all fitted to the same size of dataset", call. = FALSE
    )
  }

  Y <- vapply(object, function(z) z$model[, 1], numeric(n_residuals[1L]))
  y <- Y[, 1]

  if (!all(Y[, -1] == y)) {
    stop("models were not all fitted on the same data", call. = FALSE)
  }

  X <- vapply(object, function(z) z$model[, 2], numeric(n_residuals[1L]))
  x <- X[, 1]

  if (!all(X[, -1] == x)) {
    stop("models were not all fitted on the same data", call. = FALSE)
  }

  W <- vapply(object, function(x) x$weights, numeric(n_residuals[1L]))
  w <- W[, 1]

  if (!all(W[, -1] == w)) {
    stop("models were not all fitted with the same weights", call. = FALSE)
  }

  idx <- !is.na(y) & !is.na(x) & !is.na(w) & !(w == 0)

  if (sum(idx) != length(y)) {
    y <- y[idx]
    x <- x[idx]
    w <- w[idx]
  }

  n_obs <- length(y)
  log_n <- log(n_obs)
  log_w <- sum(log(w))

  n_params <- vapply(object, function(x) sum(x$estimated), 0)

  tmp_dev <- vapply(object, function(x) x$rss, 0)

  ord <- order(n_params, -tmp_dev)

  object <- object[ord]
  n_params <- n_params[ord]
  tmp_dev <- tmp_dev[ord]

  tmp_df <- vapply(object, function(x) x$df.residual, 0)

  k <- max(n_params)

  df <- if (k >= 5) {
    c(1, n_params)
  } else {
    c(1, n_params, 5)
  }

  deviance_df <- if (k >= 5) {
    c(-1, tmp_df)
  } else {
    c(-1, tmp_df, -1)
  }

  deviance_value <- if (k >= 5) {
    c(-1, tmp_dev)
  } else {
    c(-1, tmp_dev, -1)
  }

  l <- length(deviance_df)

  loglik <- rep(-1, l)

  # constant model: horizontal line
  weighted_mean <- sum(w * y) / sum(w)
  deviance_df[1] <- n_obs - 1
  deviance_value[1] <- sum(w * (y - weighted_mean)^2)

  if (k < 5) {
    fit <- drda(y ~ x, weights = w, mean_function = "logistic5")
    deviance_df[l] <- fit$df.residual
    deviance_value[l] <- fit$rss
  }

  loglik <- loglik_normal(deviance_value, n_obs, log_w)

  aic <- 2 * (df - loglik)
  bic <- log_n * df - 2 * loglik
  lrt <- -2 * (loglik - loglik[l])
  lrt[l] <- NA_real_

  table <- data.frame(deviance_df, deviance_value, df, aic, bic, lrt)
  table$pvalue <- pchisq(lrt, -diff(deviance_df), lower.tail = FALSE)

  colnames(table) <- c(
    "Resid. Df", "Resid. Dev", "Df", "AIC", "BIC", "LRT", "p-value"
  )
  rownames(table) <- paste("Model", seq_len(l))

  title <- "Analysis of Deviance Table\n"

  str <- vapply(object, function(z) z$mean_function, "string")
  str <- paste("Model ", 2:(n_models + 1), ": ", str, sep = "")

  topnote <- if (k >= 5) {
    str[n_models] <- paste(str[n_models], "(Full)")
    paste(c("Model 1: Constant", str, "\n"), collapse = "\n")
  } else {
    paste(
      c(
        "Model 1: Constant", str,
        paste("Model ", l, ": logistic5 (Full)\n", sep = "")
      ),
      collapse = "\n"
    )
  }

  comment <- paste(
    "Model", which.min(aic),
    "is the best model according to the Akaike Information Criterion.\n"
  )

  structure(
    table, heading = c(title, topnote, comment),
    class = c("anova", "data.frame")
  )
}

#' @importFrom stats deviance
#'
#' @export
deviance.drda <- function(object, ...) {
  object$rss
}

#' @importFrom stats logLik
#'
#' @export
logLik.drda <- function(object, ...) {
  object$loglik
}

#' @importFrom stats predict
#'
#' @export
predict.drda <- function(object, x, ...) {
  if (missing(x) || is.null(x)) {
    return(object$fitted.values)
  }

  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("variable `x` is not a numeric vector")
  }

  fn(object, x, object$coefficients)
}

#' @export
print.drda <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  if (length(x$coefficients)) {
    cat("Coefficients:\n")

    print.default(
      format(x$coefficients, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else {
    cat("No coefficients\n")
  }

  cat("\n")

  invisible(x)
}

#' @importFrom stats naprint printCoefmat quantile setNames
#'
#' @export
print.summary.drda <- function(
  x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
  signif.stars = getOption("show.signif.stars"), ...
) {
  cat(
    "\nCall: ",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n",
    sep = ""
  )

  cat("\nPearson Residuals:\n")
  pearson_resid <- setNames(
    quantile(x$pearson_resid, na.rm = TRUE),
    c("Min", "1Q", "Median", "3Q", "Max")
  )
  print.default(
    zapsmall(pearson_resid, digits + 1L),
    digits = digits,
    na.print = "",
    print.gap = 2L
  )

  cat("\nParameters:\n")
  printCoefmat(
    x$param, digits = digits, cs.ind = numeric(0), P.values = FALSE,
    has.Pvalue = FALSE
  )

  cat(
    "\nResidual standard error on", x$df.residual, "degrees of freedom\n"
  )

  msg <- naprint(x$na.action)
  if (nzchar(msg)) {
    cat("  (", msg, ")\n", sep = "")
  }

  cat("\nLog-likelihood:", format(x$loglik, digits = max(4L, digits + 1L)))
  cat("\nAIC:", format(x$aic, digits = max(4L, digits + 1L)))
  cat("\nBIC:", format(x$bic, digits = max(4L, digits + 1L)))
  cat("\n")

  if (x$converged) {
    cat(
      "\nOptimization algorithm converged in",
      x$iterations,
      "iterations\n"
    )
  } else {
    cat(
      "\nOptimization algorithm DID NOT converge in",
      x$iterations,
      "iterations\n"
    )
  }

  invisible(x)
}

#' @importFrom stats coef naresid residuals
#'
#' @export
residuals.drda <- function(
  object, type = c("response", "weighted", "pearson"), ...
) {
  r <- object$residuals

  type <- match.arg(type)
  if (type != "response") {
    if (!is.null(object$weights)) {
      r <- sqrt(object$weights) * r
    }

    if (type == "pearson") {
      r <- r / object$sigma
    }
  }

  naresid(object$na.action, r)
}

#' @importFrom stats sigma
#'
#' @export
sigma.drda <- function(object, ...) {
  object$sigma
}

#' @importFrom stats qnorm
#'
#' @export
summary.drda <- function(object, ...) {
  std_err <- sqrt(diag(object$vcov))

  object$pearson_resid <- residuals(object, type = "pearson")

  object$param <- c(object$coefficients, sigma = object$sigma)

  if (inherits(object, "logistic2_fit")) {
    names(object$param) <- c(
      "Growth rate", "Midpoint at", "Residual std err."
    )
  } else if (inherits(object, "logistic4_fit")) {
    names(object$param) <- c(
      "Minimum", "Maximum", "Growth rate", "Midpoint at", "Residual std err."
    )
  }

  object$param <- matrix(
    c(
      object$param,
      object$param + qnorm(0.025) * std_err,
      object$param + qnorm(0.975) * std_err
    ),
    ncol = 3,
    dimnames = list(
      names(object$param),
      c("Estimate", "Lower .95", "Upper .95")
    )
  )

  if (!object$is_log) {
    # give the user summaries on the same scale they provided
    if (inherits(object, "logistic2_fit")) {
      object$param[2, ] <- exp(object$param[2, ])
    } else if (inherits(object, "logistic4_fit")) {
      object$param[4, ] <- exp(object$param[4, ])
    }
  }

  k <- sum(object$estimated)

  object$aic <- 2 * (k - object$loglik)
  object$bic <- log(object$n) * k - 2 * object$loglik

  class(object) <- "summary.drda"

  object
}

#' @importFrom stats naresid weights
#'
#' @export
weights.drda <- function(object, ...) {
  if (is.null(object$na.action)) {
    object$weights
  } else {
    naresid(object$na.action, object$weights)
  }
}
