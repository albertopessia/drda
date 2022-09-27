#' Fit non-linear growth curves
#'
#' Use the Newton's with a trust-region method to fit non-linear growth curves
#' to observed data.
#'
#' @param formula an object of class \code{\link[stats]{formula}} (or one that
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
#'   unset. The 'factory-fresh' default is \code{na.omit}. Another
#'   possible value is `NULL`, no action. Value \code{na.exclude} can be useful.
#' @param mean_function the model to be fitted. See `details` for available
#'   models.
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
#' ### Generalized (5-parameter) logistic function
#'
#' The 5-parameter logistic function can be selected by choosing
#' `mean_function = "logistic5"` or `mean_function = "l5"`. The function is
#' defined here as
#'
#' `alpha + delta / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)`
#'
#' where `eta > 0` and `nu > 0`. When `delta` is positive (negative) the curve
#' is monotonically increasing (decreasing).
#'
#' Parameter `alpha` is the value of the function when `x -> -Inf`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the mid-value of the function.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' ### 4-parameter logistic function
#'
#' The 4-parameter logistic function is the default model of `drda`. It can be
#' explicitly selected by choosing `mean_function = "logistic4"` or
#' `mean_function = "l4"`. The function is obtained by setting `nu = 1` in the
#' generalized logistic function, that is
#'
#' `alpha + delta / (1 + exp(-eta * (x - phi)))`
#'
#' where `eta > 0`. When `delta` is positive (negative) the curve is
#' monotonically increasing (decreasing).
#'
#' Parameter `alpha` is the value of the function when `x -> -Inf`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; alpha, delta, eta, phi) = alpha + delta / 2`.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' ### 2-parameter logistic function
#'
#' The 2-parameter logistic function can be selected by choosing
#' `mean_function = "logistic2"` or `mean_function = "l2"`. For a monotonically
#' increasing curve set `nu = 1`, `alpha = 0`, and `delta = 1`:
#'
#' `1 / (1 + exp(-eta * (x - phi)))`
#'
#' For a monotonically decreasing curve set `nu = 1`, `alpha = 1`, and
#' `delta = -1`:
#'
#' `1 - 1 / (1 + exp(-eta * (x - phi)))`
#'
#' where `eta > 0`. The lower bound of the curve is zero while the upper bound
#' of the curve is one.
#'
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; eta, phi) = 1 / 2`.
#'
#' ### Gompertz function
#'
#' The Gompertz function is the limit for `nu -> 0` of the 5-parameter logistic
#' function. It can be selected by choosing `mean_function = "gompertz"` or
#' `mean_function = "gz"`. The function is defined in this package as
#'
#' `alpha + delta * exp(-exp(-eta * (x - phi)))`
#'
#' where `eta > 0`.
#'
#' Parameter `alpha` is the value of the function when `x -> -Inf`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` sets the displacement along the `x`-axis.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' The mid-point of the function, that is `alpha + delta / 2`, is achieved at
#' `x = phi - log(log(2)) / eta`.
#'
#' ### Generalized (5-parameter) log-logistic function
#'
#' The 5-parameter log-logistic function is selected by setting
#' `mean_function = "loglogistic5"` or `mean_function = "ll5"`. The function is
#' defined here as
#'
#' `alpha + delta * (x^eta / (x^eta + nu * phi^eta))^(1 / nu)`
#'
#' where `x >= 0`, `eta > 0`, `phi > 0`, and `nu > 0`. When `delta` is
#' positive (negative) the curve is monotonically increasing (decreasing). The
#' function is defined only for positive values of the predictor variable `x`.
#'
#' Parameter `alpha` is the value of the function at `x = 0`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` is related to the mid-value of the function.
#' Parameter `nu` affects near which asymptote maximum growth occurs.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' ### 4-parameter log-logistic function
#'
#' The 4-parameter log-logistic function is selected by setting
#' `mean_function = "loglogistic4"` or `mean_function = "ll4"`. The function is
#' obtained by setting `nu = 1` in the generalized log-logistic function, that
#' is
#'
#' `alpha + delta * x^eta / (x^eta + phi^eta)`
#'
#' where `x >= 0` and `eta > 0`. When `delta` is positive (negative) the curve
#' is monotonically increasing (decreasing). The function is defined only for
#' positive values of the predictor variable `x`.
#'
#' Parameter `alpha` is the value of the function at `x = 0`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; alpha, delta, eta, phi) = alpha + delta / 2`.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' ### 2-parameter log-logistic function
#'
#' The 2-parameter log-logistic function is selected by setting
#' `mean_function = "loglogistic2"` or `mean_function = "ll2"`. For a
#' monotonically increasing curve set `nu = 1`, `alpha = 0`, and `delta = 1`:
#'
#' `x^eta / (x^eta + phi^eta)`
#'
#' For a monotonically decreasing curve set `nu = 1`, `alpha = 1`, and
#' `delta = -1`:
#'
#' `1 - x^eta / (x^eta + phi^eta)`
#'
#' where `x >= 0`, `eta > 0`, and `phi > 0`. The lower bound of the curve is
#' zero while the upper bound of the curve is one.
#'
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` represents the `x` value at which the curve is equal to its
#' mid-point, i.e. `f(phi; eta, phi) = 1 / 2`.
#'
#' ### log-Gompertz function
#'
#' The log-Gompertz function is the limit for `nu -> 0` of the 5-parameter
#' log-logistic function. It can be selected by choosing
#' `mean_function = "loggompertz"` or `mean_function = "lgz"`. The function is
#' defined in this package as
#'
#' `alpha + delta * exp(-(phi / x)^eta)`
#'
#' where `x > 0`, `eta > 0`, and `phi > 0`. Note that the limit for `x -> 0` is
#' `alpha`. When `delta` is positive (negative) the curve is monotonically
#' increasing (decreasing). The function is defined only for positive values of
#' the predictor variable `x`.
#'
#' Parameter `alpha` is the value of the function at `x = 0`.
#' Parameter `delta` is the (signed) height of the curve.
#' Parameter `eta` represents the steepness (growth rate) of the curve.
#' Parameter `phi` sets the displacement along the `x`-axis.
#'
#' The value of the function when `x -> Inf` is `alpha + delta`. In
#' dose-response studies `delta` can be interpreted as the maximum theoretical
#' achievable effect.
#'
#' ### Constrained optimization
#'
#' It is possible to search for the maximum likelihood estimates within
#' pre-specified interval regions.
#'
#' *Note*: Hypothesis testing is not available for constrained estimates
#' because asymptotic approximations might not be valid.
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
#'   }
#'
#' @importFrom stats model.frame model.matrix model.response model.weights terms
#'
#' @export
#'
#' @examples
#' # by default `drda` uses a 4-parameter logistic function for model fitting
#' fit_l4 <- drda(response ~ log_dose, data = voropm2)
#'
#' # get a general overview of the results
#' summary(fit_l4)
#'
#' # compare the model against a flat horizontal line and the full model
#' anova(fit_l4)
#'
#' # 5-parameter logistic curve appears to be a better model
#' fit_l5 <- drda(response ~ log_dose, data = voropm2, mean_function = "l5")
#' plot(fit_l4, fit_l5)
#'
#' # fit a 2-parameter logistic function
#' fit_l2 <- drda(response ~ log_dose, data = voropm2, mean_function = "l2")
#'
#' # compare our models
#' anova(fit_l2, fit_l4)
#'
#' # use log-logistic functions when utilizing doses (instead of log-doses)
#' # here we show the use of other arguments as well
#' fit_ll5 <- drda(
#'   response ~ dose, weights = weight, data = voropm2,
#'   mean_function = "loglogistic5", lower_bound = c(0.5, -1.5, 0, -Inf, 0.25),
#'   upper_bound = c(1.5, 0.5, 5, Inf, 3), start = c(1, -1, 3, 100, 1),
#'   max_iter = 10000
#' )
#'
#' # note that the maximum likelihood estimate is outside the region of
#' # optimization: not only the variance-covariance matrix is now singular but
#' # asymptotic assumptions do not hold anymore.
drda <- function(
  formula, data, subset, weights, na.action, mean_function = "logistic4",
  lower_bound = NULL, upper_bound = NULL, start = NULL, max_iter = 1000
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

  object <- if (mean_function == "logistic4" || mean_function == "l4") {
    # we want to make sure it is the full name instead of the abbreviation
    mean_function <- "logistic4"
    logistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)
  } else if (mean_function == "logistic2" || mean_function == "l2") {
    mean_function <- "logistic2"
    logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)
  } else if (mean_function == "logistic5" || mean_function == "l5") {
    mean_function <- "logistic5"
    logistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)
  } else if (mean_function == "gompertz" || mean_function == "gz") {
    mean_function <- "gompertz"
    gompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)
  } else if (mean_function == "logistic6" || mean_function == "l6") {
    mean_function <- "logistic6"
    logistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)
  } else {
    if (any(x < 0)) {
      stop("predictor variable 'x' is negative", call. = FALSE)
    }

    if (mean_function == "loglogistic4" || mean_function == "ll4") {
      mean_function <- "loglogistic4"
      loglogistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)
    } else if (mean_function == "loglogistic2" || mean_function == "ll2") {
      mean_function <- "loglogistic2"
      loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)
    } else if (mean_function == "loglogistic5" || mean_function == "ll5") {
      mean_function <- "loglogistic5"
      loglogistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)
    } else if (mean_function == "loggompertz" || mean_function == "lgz") {
      mean_function <- "loggompertz"
      loggompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)
    } else if (mean_function == "loglogistic6" || mean_function == "ll6") {
      mean_function <- "loglogistic6"
      loglogistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)
    } else {
      stop(
        "chosen 'mean_function' is wrongly typed or not yet available",
        call. = FALSE
      )
    }
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
    vcov <- matrix(NA_real_, nrow = k, ncol = k)
    vcov[-idx, -idx] <- approx_vcov(result$fisher.info[-idx, -idx])
    vcov
  }

  result$call <- match.call()
  result$terms <- model_terms
  result$model <- model_frame
  result$na.action <- attr(model_frame, "na.action")

  if (any(w_zero)) {
    # fitting was done with only positive weights but we want to report all
    # of them
    result$weights <- model_frame[, 3]
  }

  class(result) <- c(class(result), "drda")

  result
}

#' @importFrom stats anova pchisq var
#'
#' @export
anova.drda <- function(object, ...) {
  # check for multiple objects
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
    stop(
      "hypothesis testing is not available for constrained optimization",
      call. = FALSE
    )
  }

  s <- substr(object$mean_function, 1, 8)

  model_type <- if (s == "logistic" || s == "gompertz") {
    1
  } else if (s == "loglogis" || s == "loggompe") {
    2
  } else {
    stop("model not supported", call. = FALSE)
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
    s <- substr(object$mean_function, 1, 8)
    fit <- if (s == "logistic" || s == "gompertz") {
      drda(y ~ x, weights = w, mean_function = "logistic5")
    } else {
      drda(y ~ x, weights = w, mean_function = "loglogistic5")
    }

    deviance_df[3] <- fit$df.residual
    deviance_value[3] <- fit$rss
  }

  loglik <- loglik_normal(deviance_value, n, log_w)

  df <- n - deviance_df + 1

  aic <- 2 * (df - loglik)
  bic <- log_n * df - 2 * loglik
  lrt <- 2 * diff(loglik)
  dvn <- c(NA_real_, diff(deviance_value))

  table <- data.frame(
    deviance_df, deviance_value, df - 1, aic, bic, dvn, c(NA_real_, lrt)
  )

  pvalue <- pchisq(lrt, diff(df), lower.tail = FALSE)
  pvalue[pvalue == 0] <- NA_real_

  table$pvalue <- c(NA_real_, pvalue)

  colnames(table) <- c(
    "Resid. Df", "Resid. Dev", "Df", "AIC", "BIC", "Deviance", "LRT", "Pr(>Chi)"
  )
  rownames(table) <- paste("Model", seq_len(l))

  title <- "Analysis of Deviance Table\n"

  str <- switch(object$mean_function,
    logistic2 = if (object$coefficients[2] >= 0) {
      "1 / (1 + exp(-e * (x - p)))"
    } else {
      "1 - 1 / (1 + exp(-e * (x - p)))"
    },
    logistic4 = "a + d / (1 + exp(-e * (x - p)))",
    logistic5 = "a + d / (1 + n * exp(-e * (x - p)))^(1 / n)",
    logistic6 = "a + d / (w + n * exp(-e * (x - p)))^(1 / n)",
    gompertz = "a + d * exp(-exp(-e * (x - p)))",
    loglogistic2 = if (object$coefficients[2] >= 0) {
      "x^e / (x^e + p^e)"
    } else {
      "1 - x^e / (x^e + p^e)"
    },
    loglogistic4 = "a + d * x^e / (x^e + p^e)",
    loglogistic5 = "a + d * (x^e / (x^e + n * p^e))^(1 / n)",
    loglogistic6 = "a + d * (x^e / (w * x^e + n * p^e))^(1 / n)",
    loggompertz = "a + d * exp(-(p / x)^e)"
  )

  topnote <- if (k >= 5) {
    paste(
      paste0(
        c(
          "Model 1: a", "\n",
          "Model 2: ", str, " (Full)", "\n"
        )
      ),
      collapse = ""
    )
  } else {
    tmp <- if (model_type == 1) {
      "a + d / (1 + n * exp(-e * (x - p)))^(1 / n)"
    } else {
      "a + d * (x^e / (x^e + n * p^e))^(1 / n)"
    }

    paste(
      paste0(
        c(
          "Model 1: a", "\n",
          "Model 2: ", str, " (Fit)", "\n",
          "Model 3: ", tmp, " (Full)", "\n"
        )
      ),
      collapse = ""
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
    stop(
      "hypothesis testing is not available for constrained optimization",
      call. = FALSE
    )
  }

  model_type <- NULL
  for (i in seq_len(n_models)) {
    s <- substr(object[[i]]$mean_function, 1, 8)

    if (s == "logistic" || s == "gompertz") {
      if (!is.null(model_type) && model_type != 1) {
        stop("curves defined on different domains", call. = FALSE)
      }

      model_type <- 1
    } else if (s == "loglogis" || s == "loggompe") {
      if (!is.null(model_type) && model_type != 2) {
        stop("curves defined on different domains", call. = FALSE)
      }

      model_type <- 2
    } else {
      stop("model not supported", call. = FALSE)
    }
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
    c(1, n_params) + 1
  } else {
    c(1, n_params, 5) + 1
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
    fit <- if (model_type == 1) {
      drda(y ~ x, weights = w, mean_function = "logistic5")
    } else {
      drda(y ~ x, weights = w, mean_function = "loglogistic5")
    }

    deviance_df[l] <- fit$df.residual
    deviance_value[l] <- fit$rss
  }

  loglik <- loglik_normal(deviance_value, n_obs, log_w)

  aic <- 2 * (df - loglik)
  bic <- log_n * df - 2 * loglik
  df <- diff(df)
  lrt <- 2 * diff(loglik)
  dvn <- c(NA_real_, diff(deviance_value))

  table <- data.frame(
    deviance_df, deviance_value, c(NA_real_, df), aic, bic, dvn,
    c(NA_real_, lrt)
  )

  pvalue <- pchisq(lrt, df, lower.tail = FALSE)
  pvalue[pvalue == 0] <- NA_real_

  table$pvalue <- c(NA_real_, pvalue)

  colnames(table) <- c(
    "Resid. Df", "Resid. Dev", "Df", "AIC", "BIC", "Deviance", "LRT", "Pr(>Chi)"
  )
  rownames(table) <- paste("Model", seq_len(l))

  title <- "Analysis of Deviance Table\n"

  f <- function(x) {
    switch(x$mean_function,
      logistic2 = if (x$coefficients[2] >= 0) {
        "1 / (1 + exp(-e * (x - p)))"
      } else {
        "1 - 1 / (1 + exp(-e * (x - p)))"
      },
      logistic4 = "a + d / (1 + exp(-e * (x - p)))",
      logistic5 = "a + d / (1 + n * exp(-e * (x - p)))^(1 / n)",
      logistic6 = "a + d / (w + n * exp(-e * (x - p)))^(1 / n)",
      gompertz = "a + d * exp(-exp(-e * (x - p)))",
      loglogistic2 = if (x$coefficients[2] >= 0) {
        "x^e / (x^e + p^e)"
      } else {
        "1 - x^e / (x^e + p^e)"
      },
      loglogistic4 = "a + d * x^e / (x^e + p^e)",
      loglogistic5 = "a + d * (x^e / (x^e + n * p^e))^(1 / n)",
      loglogistic6 = "a + d * (x^e / (w * x^e + n * p^e))^(1 / n)",
      loggompertz = "a + d * exp(-(p / x)^e)"
    )
  }

  str <- vapply(object, f, "a")
  str <- paste0("Model ", 2:(n_models + 1), ": ", str)

  topnote <- if (k >= 5) {
    str[n_models] <- paste(str[n_models], "(Full)\n")
    paste(c("Model 1: a", str), collapse = "\n")
  } else {
    tmp <- if (model_type == 1) {
      "a + d / (1 + n * exp(-e * (x - p)))^(1 / n)"
    } else {
      "a + d * (x^e / (x^e + n * p^e))^(1 / n)"
    }

    paste(
      c("Model 1: a", str, paste0("Model ", l, ": ", tmp, " (Full)\n")),
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
  structure(
    object$loglik,
    nobs = object$n,
    df = sum(object$estimated) + 1,
    class = "logLik"
  )
}

#' @importFrom stats predict
#'
#' @export
predict.drda <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    # did the user provide arguments?
    dotargs <- list(...)

    if (length(dotargs) == 0) {
      return(object$fitted.values)
    } else {
      # we only consider the first argument
      newdata <- dotargs[[1]]
    }
  }

  if (is.data.frame(newdata) || is.matrix(newdata)) {
    if (ncol(newdata) == 1) {
      newdata <- newdata[, 1]
    } else {
      # in case of multiple columns, pick the one corresponding to our predictor
      idx <- match(colnames(object$model)[2], colnames(newdata))

      if (!is.na(idx)) {
        newdata <- newdata[, idx]
      } else {
        stop(
          "cannot find the predictor variable in 'newdata'",
          call. = FALSE
        )
      }
    }
  } else if (!is.numeric(newdata) || !is.null(dim(newdata))) {
    stop(
      "variable `newdata` is not a data.frame nor a numeric vector",
      call. = FALSE
    )
  }

  res <- fn(object, newdata, object$coefficients)
  names(res) <- names(newdata)

  res
}

#' @export
print.drda <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  if (length(x$coefficients) > 0) {
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

#' @importFrom stats AIC BIC qnorm
#'
#' @export
summary.drda <- function(object, level = 0.95, ...) {
  if (level <= 0 || level >= 1) {
    stop("Confidence level must be in the interval (0, 1)", call. = FALSE)
  }

  is_2 <- inherits(object, "logistic2_fit") ||
    inherits(object, "loglogistic2_fit")

  is_4 <- inherits(object, "logistic4_fit") ||
    inherits(object, "loglogistic4_fit")

  std_err <- if (is_2) {
    c(alpha = NA_real_, delta = NA_real_, sqrt(diag(object$vcov)))
  } else {
    sqrt(diag(object$vcov))
  }

  object$pearson_resid <- residuals(object, type = "pearson")

  object$param <- c(object$coefficients, sigma = object$sigma)

  if (is_2 || is_4) {
    names(object$param) <-  {
      c("Maximum", "Height", "Growth rate", "Midpoint at", "Residual std err.")
    }

    if (object$coefficients[2] > 0) {
      names(object$param)[1] <- "Minimum"
    }
  }

  q <- qnorm((1 - level) / 2)
  l <- round(level * 100)

  object$param <- matrix(
    c(
      object$param,
      std_err,
      object$param + q * std_err,
      object$param - q * std_err
    ),
    ncol = 4,
    dimnames = list(
      names(object$param),
      c(
        "Estimate",
        "Std. Error",
        paste0(c("Lower .", "Upper ."), c(l, l))
      )
    )
  )

  object$aic <- AIC(object)
  object$bic <- BIC(object)

  class(object) <- "summary.drda"

  object
}

#' @importFrom stats vcov
#'
#' @export
vcov.drda <- function(object, ...) {
  p <- nrow(object$vcov)
  object$vcov[-p, -p]
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
