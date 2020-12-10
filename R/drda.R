#' Fit a parametric model to dose-response data
#'
#' Use a Newton trust-region method to fit a curve to dose-response data.
#'
#' @param formula an object of class \code{link[stats]{formula}} (or one that
#'   can be coerced to that class): a symbolic description of the model to be
#'   fitted. Currently supports only formulas of the type `response ~ log_dose`.
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
#'   unset. The ‘factory-fresh’ default is \code{\link[stats]{na.omit}}. Another
#'   possible value is `NULL`, no action. Value \code{\link[stats]{na.exclude}}
#'   can be useful.
#' @param mean_function the model to be fitted. Currently only "logistic4" is
#'   supported.
#' @param lower_bound numeric vector of length 4 with the minimum admissible
#'   values of the parameters.
#' @param upper_bound numeric vector of length 4 with the maximum admissible
#'   values of the parameters.
#' @param start starting values for the parameters.
#' @param max_iter maximum number of iterations in the optimization algorithm.
#'
#' @return An object of class "drda", that is a list containing the following
#'   components:
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
#'     \item{loglik}{maximum value of the log-likelihood function.}
#'     \item{df.residuals}{residual degrees of freedom.}
#'     \item{estimated}{boolean vector indicating which parameters were
#'       estimated from the data.}
#'     \item{fitted.values}{fitted mean values.}
#'     \item{residuals}{residuals, that is response minus fitted values.}
#'     \item{weights}{(only for weighted fits) the specified weights.}
#'     \item{fisher.info}{Fisher information matrix evaluated at the maximum
#'       likelihood estimator.}
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
drda <- function(
  formula,
  data,
  subset,
  weights,
  na.action,
  mean_function = "logistic4",
  lower_bound = NULL,
  upper_bound = NULL,
  start = NULL,
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
    stop("0 (non-NA) cases")
  }

  y <- model.response(model_frame, "numeric")

  if (is.matrix(y)) {
    stop("response variable must be a vector")
  }

  if (length(y) != n) {
    stop("incompatible dimensions")
  }

  w <- as.vector(model.weights(model_frame))

  max_iter <- ceiling(max_iter[1])

  if (max_iter <= 0) {
    stop("maximum number of iterations must be positive")
  }

  result <- if (is.null(w)) {
    if (is.null(lower_bound) && is.null(upper_bound)) {
      logistic4_fit_unconstrained(x, y, start, max_iter)
    } else {
      # the user might provide only one of either lower_bound and upper_bound
      # in this case, we must initialize the NULL variable to the proper
      # unconstrained bounds
      if (is.null(lower_bound)) {
        lower_bound <- rep(-Inf, 4)
      }

      if (is.null(upper_bound)) {
        upper_bound <- rep(Inf, 4)
      }

      logistic4_fit_constrained(x, y, start, max_iter, lower_bound, upper_bound)
    }
  } else {
    if(!is.numeric(w)) {
      stop("'weights' must be a numeric vector")
    }

    if (length(w) != n) {
      stop("incompatible dimensions")
    }

    if (any(w < 0 | is.na(w))) {
      stop("missing or negative weights not allowed")
    }

    if (is.null(lower_bound) && is.null(upper_bound)) {
      logistic4_weighted_fit_unconstrained(x, y, w, start, max_iter)
    } else {
      if (is.null(lower_bound)) {
        lower_bound <- rep(-Inf, 4)
      }

      if (is.null(upper_bound)) {
        upper_bound <- rep(Inf, 4)
      }

      logistic4_weighted_fit_constrained(
        x, y, w, start, max_iter, lower_bound, upper_bound
      )
    }
  }

  result$call <- match.call()
  result$terms <- model_terms
  result$model <- model_frame
  result$na.action <- attr(model_frame, "na.action")

  class(result) <- "drda"

  result
}

#' @importFrom stats anova pchisq
#'
#' @export
anova.drda <- function(
  object,
  ...
) {
  if (object$constrained) {
    # solution to the constrained problem is unlikely the maximum likelihood
    # estimator, therefore the asymptotic approximation might not hold
    stop("hypothesis testing is not available for constrained optimization")
  }

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

  y <- object$model[, 1]
  x <- object$model[, 2]
  w <- object$weights

  idx <- if (is.null(w)) {
    !is.na(y) & !is.na(x)
  } else {
    !is.na(y) & !is.na(x) & !is.na(w) & !(w == 0)
  }

  if (sum(idx) != length(y)) {
    y <- y[idx]
    x <- x[idx]
    w <- w[idx]
  }

  n <- length(y)
  log_n <- log(n)
  k <- sum(object$estimated)

  l <- if (k == 4) {
    # we compare the full model against a flat horizontal line
    2
  } else {
    # we compare the estimated model against the baseline and the full model
    3
  }

  deviance_df <- rep(-1, l)
  deviance_value <- rep(-1, l)
  loglik <- rep(-1, l)

  if (is.null(w)) {
    # constant model: horizontal line
    deviance_df[1] <- n - 1
    deviance_value[1] <- (n - 1) * var(y)

    # fitted model
    deviance_df[2] <- object$df.residual
    deviance_value[2] <- object$rss

    if (k < 4) {
      # at least a parameter was fixed, so we now fit the full model
      fit <- drda(y ~ x)
      deviance_df[3] <- fit$df.residual
      deviance_value[3] <- fit$rss
    }

    loglik <- - n * (1 + log(2 * pi * deviance_value / n)) / 2
  } else {
    log_w <- sum(log(w))

    # constant model: horizontal line
    weighted_mean <- sum(w * y) / sum(w)
    deviance_df[1] <- n - 1
    deviance_value[1] <- sum(w * (y - weighted_mean)^2)

    # fitted model
    deviance_df[2] <- object$df.residual
    deviance_value[2] <- object$rss

    if (k < 4) {
      # at least a parameter was considered fixed, so we now fit the full model
      fit <- drda(y ~ x, weights = w)
      deviance_df[3] <- fit$df.residual
      deviance_value[3] <- fit$rss
    }

    loglik <- -(n * (1 + log(2 * pi * deviance_value / n)) - log_w) / 2
  }

  df <- n - deviance_df

  aic <- 2 * (df - loglik)
  bic <- log_n * df - 2 * loglik
  lrt <- -2 * (loglik - loglik[l])
  lrt[l] <- NA_real_

  table <- data.frame(deviance_df, deviance_value, df, aic, bic, lrt)
  table$pvalue <- pchisq(lrt, -diff(deviance_df), lower.tail = FALSE)

  rownames(table)[1:2] <- c("Constant model", "Estimated model")
  if (k < 4) {
    rownames(table)[3] <- "Full model"
  }
  colnames(table) <- c(
    "Resid. Df", "Resid. Dev", "Df", "AIC", "BIC", "LRT", "p-value"
  )

  title <- "Analysis of Deviance Table\n\nModel: 4 parameter logistic\n"

  structure(table, heading = title, class = c("anova", "data.frame"))
}

#' @importFrom stats anova pchisq
#'
#' @export
anova.drdalist <- function(
  object,
  ...
) {
  n_models <- length(object)

  if (n_models == 1) {
    return(anova.drda(object[[1L]]))
  }

  n_residuals <- vapply(object, function(x) length(x$residuals), 0)

  if (any(n_residuals != n_residuals[1L])) {
    stop("models were not all fitted to the same size of dataset")
  }

  y <- object[[1L]]$model[, 1]
  n_obs <- length(y)
  log_n <- log(n_obs)

  weights <- unlist(sapply(object, function(x) x$weights))
  w <- NULL

  if (!is.null(weights)) {
    if (length(weights) != n_models * n_obs) {
      stop("models were not all fitted with weights")
    }

    w <- weights[, 1]

    # there are models with weights, check that they are all the same
    if (!all(weights[, 2:n_models] == w)) {
      stop("models were not all fitted with the same weights")
    }
  }

  n_params <- vapply(object, function(x) sum(x$estimated), 0)

  tmp_dev <- vapply(object, function(x) x$rss, 0)

  ord <- order(n_params, -tmp_dev)

  object <- object[ord]
  n_params <- n_params[ord]
  tmp_dev <- tmp_dev[ord]

  tmp_df <- vapply(object, function(x) x$df.residual, 0)

  k <- max(n_params)

  df <- if (k == 4) {
    c(1, n_params)
  } else {
    c(1, n_params, 4)
  }

  deviance_df <- if (k == 4) {
    c(-1, tmp_df)
  } else {
    c(-1, tmp_df, -1)
  }

  deviance_value <- if (k == 4) {
    c(-1, tmp_dev)
  } else {
    c(-1, tmp_dev, -1)
  }

  l <- length(deviance_df)

  loglik <- rep(-1, l)

  if (is.null(w)) {
    # constant model: horizontal line
    deviance_df[1] <- n_obs - 1
    deviance_value[1] <- (n_obs - 1) * var(y)

    if (k < 4) {
      fit <- drda(y ~ x)
      deviance_df[l] <- fit$df.residual
      deviance_value[l] <- fit$rss
    }

    loglik <- - n_obs * (1 + log(2 * pi * deviance_value / n_obs)) / 2
  } else {
    log_w <- sum(log(w))

    # constant model: horizontal line
    weighted_mean <- sum(w * y) / sum(w)
    deviance_df[1] <- n_obs - 1
    deviance_value[1] <- sum(w * (y - weighted_mean)^2)

    if (k < 4) {
      fit <- drda(y ~ x, weights = w)
      deviance_df[l] <- fit$df.residual
      deviance_value[l] <- fit$rss
    }

    loglik <- -(n_obs * (1 + log(2 * pi * deviance_value / n_obs)) - log_w) / 2
  }

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

  str <- vapply(
    object,
    function(x) {
      paste(names(x$estimated)[!x$estimated], collapse = ", ")
    },
    "string"
  )

  str <- paste("Model ", 2:(n_models + 1), ": Fixed ", str, sep = "")

  topnote <- if (k == 4) {
    str[n_models] <- paste("Model ", l, ": Complete", sep = "")
    paste(
      c("Model 1: Constant", str, "\n"), collapse = "\n"
    )
  } else {
    paste(
      c("Model 1: Constant", str, paste("Model ", l, ": Complete\n", sep = "")),
      collapse = "\n"
    )
  }

  structure(
    table, heading = c(title, topnote), class = c("anova", "data.frame")
  )
}

#' @importFrom stats deviance
#'
#' @export
deviance.drda <- function(
  object,
  ...
) {
  object$rss
}

#' @importFrom stats logLik
#'
#' @export
logLik.drda <- function(
  object,
  ...
) {
  object$loglik
}

#' @importFrom stats predict
#'
#' @export
predict.drda <- function(
  object,
  x,
  ...
) {
  if (missing(x) || is.null(x)) {
    return(object$fitted.values)
  }

  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("variable `x` is not a numeric vector")
  }

  logistic4_function(x, object$coefficients)
}

#' @export
print.drda <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call),sep = "\n", collapse = "\n"),
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
  x,
  digits = max(3L, getOption("digits") - 3L),
  symbolic.cor = x$symbolic.cor,
  signif.stars = getOption("show.signif.stars"),
  ...
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
  printCoefmat(x$param, digits = digits, signif.stars = signif.stars, ...)

  cat(
    "\nResidual standard error:",
    format(signif(x$sigma, digits)),
    "on",
    x$df.residual,
    "degrees of freedom\n"
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
  object,
  type = c("response", "weighted", "pearson"),
  ...
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
sigma.drda <- function(
  object,
  ...
) {
  object$sigma
}

#' @export
summary.drda <- function(
  object,
  ...
) {
  object$pearson_resid <- residuals(object, type = "pearson")

  object$param <- object$coefficients
  object$param <- matrix(
    object$param,
    ncol = 1,
    dimnames = list(
      names(object$param),
      "Estimate"
    )
  )

  k <- sum(object$estimated)
  n <- object$df.residual + k

  object$aic <- 2 * (k - object$loglik)
  object$bic <- log(n) * k - 2 * object$loglik

  class(object) <- "summary.drda"

  object
}
