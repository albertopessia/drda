test_that("drda: default arguments", {
  result <- drda(y ~ x, data = ltd$D)

  expect_s3_class(result, "drda")
  expect_s3_class(result, "logistic")
  expect_s3_class(result, "logistic4_fit")

  expect_named(
    result,
    c(
      "converged", "iterations", "constrained", "estimated", "coefficients",
      "rss", "df.residual", "fitted.values", "weights", "residuals",
      "mean_function", "n", "sigma", "loglik", "fisher.info", "vcov", "call",
      "terms", "model"
    )
  )
})

test_that("drda: weights", {
  result <- drda(
    y ~ x, data = lltd$D, weights = w, mean_function = "loggompertz"
  )

  expect_s3_class(result, "drda")
  expect_s3_class(result, "loglogistic")
  expect_s3_class(result, "loggompertz_fit")

  expect_named(
    result,
    c(
      "converged", "iterations", "constrained", "estimated", "coefficients",
      "rss", "df.residual", "fitted.values", "weights", "residuals",
      "mean_function", "n", "sigma", "loglik", "fisher.info", "vcov", "call",
      "terms", "model"
    )
  )

  expect_length(result$fitted.values, nrow(lltd$D))
  expect_length(result$weights, nrow(lltd$D))
  expect_length(result$residuals, nrow(lltd$D))

  expect_identical(result$weights, lltd$D$w)

  expect_true(all(residuals(result) != 0))
  expect_identical(
    residuals(result, type = "weighted")[lltd$D$w == 0],
    rep(0, sum(lltd$D$w == 0))
  )
})

test_that("drda: predict", {
  fit <- drda(y ~ x, data = ltd$D, mean_function = "logistic5")

  expect_error(predict(fit, "a"))
  expect_error(predict(fit, matrix(c("a", "b"), ncol = 1)))
  expect_error(predict(fit, 0, level = 0))
  expect_error(predict(fit, 0, level = 1))
  expect_error(predict(fit, 0, level = -1))
  expect_error(predict(fit, 0, level = 2))
  expect_error(predict(fit, data.frame(a = c(-1, 0, 1), b = c(-2, 0, 2))))

  expect_warning(predict(fit, a = c(-1, 0, 1), b = c(-2, 0, 2)))
  expect_warning(predict(fit, interval = TRUE))

  expect_identical(fit$fitted.values, predict(fit))

  newdata <- 0

  result <- predict(fit, newdata)
  expect_type(result, "double")
  expect_length(result, length(newdata))
  expect_identical(as.numeric(logistic5_fn(newdata, fit$coefficients)), result)

  result <- predict(fit, newdata, se.fit = TRUE)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("fit", "se.fit", "df"))
  expect_length(result$fit, length(newdata))
  expect_length(result$se.fit, length(newdata))
  expect_length(result$df, 1)
  expect_identical(
    as.numeric(logistic5_fn(newdata, fit$coefficients)), result$fit
  )

  result <- predict(fit, newdata, interval = TRUE, level = 0.99)
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_identical(dim(result), as.integer(c(length(newdata), 3)))
  expect_identical(colnames(result), c("fit", "lwr", "upr"))
  expect_identical(
    as.numeric(logistic5_fn(newdata, fit$coefficients)), as.numeric(result[, 1])
  )

  result <- predict(fit, newdata, se.fit = TRUE, interval = TRUE, level = 0.99)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_named(result, c("fit", "se.fit", "lwr", "upr", "df"))
  expect_length(result$fit, length(newdata))
  expect_length(result$se.fit, length(newdata))
  expect_length(result$lwr, length(newdata))
  expect_length(result$upr, length(newdata))
  expect_length(result$df, 1)
  expect_identical(
    as.numeric(logistic5_fn(newdata, fit$coefficients)), result$fit
  )

  newdata <- c(-2, -1, 0, 1, 2)

  result <- predict(fit, newdata)
  expect_type(result, "double")
  expect_length(result, length(newdata))
  expect_identical(logistic5_fn(newdata, fit$coefficients), result)

  result <- predict(fit, newdata, se.fit = TRUE)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("fit", "se.fit", "df"))
  expect_length(result$fit, length(newdata))
  expect_length(result$se.fit, length(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata, fit$coefficients), result$fit)

  result <- predict(fit, newdata, interval = TRUE, level = 0.99)
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_identical(dim(result), as.integer(c(length(newdata), 3)))
  expect_identical(colnames(result), c("fit", "lwr", "upr"))
  expect_identical(logistic5_fn(newdata, fit$coefficients), result[, 1])

  result <- predict(fit, newdata, se.fit = TRUE, interval = TRUE, level = 0.99)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_named(result, c("fit", "se.fit", "lwr", "upr", "df"))
  expect_length(result$fit, length(newdata))
  expect_length(result$se.fit, length(newdata))
  expect_length(result$lwr, length(newdata))
  expect_length(result$upr, length(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata, fit$coefficients), result$fit)

  newdata <- data.frame(a = c(-2, -1, 0, 1, 2))

  result <- predict(fit, newdata)
  expect_type(result, "double")
  expect_length(result, nrow(newdata))
  expect_identical(logistic5_fn(newdata$a, fit$coefficients), result)

  result <- predict(fit, newdata, se.fit = TRUE)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("fit", "se.fit", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata$a, fit$coefficients), result$fit)

  result <- predict(fit, newdata, interval = TRUE, level = 0.99)
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_identical(dim(result), as.integer(c(nrow(newdata), 3)))
  expect_identical(colnames(result), c("fit", "lwr", "upr"))
  expect_identical(logistic5_fn(newdata$a, fit$coefficients), result[, 1])

  result <- predict(fit, newdata, se.fit = TRUE, interval = TRUE, level = 0.99)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_named(result, c("fit", "se.fit", "lwr", "upr", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$lwr, nrow(newdata))
  expect_length(result$upr, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata$a, fit$coefficients), result$fit)

  newdata <- data.frame(a = c(-3, -2, 0, 2, 3), x = c(-2, -1, 0, 1, 2))

  result <- predict(fit, newdata)
  expect_type(result, "double")
  expect_length(result, nrow(newdata))
  expect_identical(logistic5_fn(newdata$x, fit$coefficients), result)

  result <- predict(fit, newdata, se.fit = TRUE)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("fit", "se.fit", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata$x, fit$coefficients), result$fit)

  result <- predict(fit, newdata, interval = TRUE, level = 0.99)
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_identical(dim(result), as.integer(c(nrow(newdata), 3)))
  expect_identical(colnames(result), c("fit", "lwr", "upr"))
  expect_identical(logistic5_fn(newdata$x, fit$coefficients), result[, 1])

  result <- predict(fit, newdata, se.fit = TRUE, interval = TRUE, level = 0.99)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_named(result, c("fit", "se.fit", "lwr", "upr", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$lwr, nrow(newdata))
  expect_length(result$upr, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata$x, fit$coefficients), result$fit)

  newdata <- matrix(
    c(-2, 1, 0, 1, 2, -3, -2, 0, 2, 3), ncol = 2,
    dimnames = list(NULL, c("a", "x"))
  )

  result <- predict(fit, newdata)
  expect_type(result, "double")
  expect_length(result, nrow(newdata))
  expect_identical(logistic5_fn(newdata[, "x"], fit$coefficients), result)

  result <- predict(fit, newdata, se.fit = TRUE)
  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("fit", "se.fit", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata[, "x"], fit$coefficients), result$fit)

  result <- predict(fit, newdata, interval = TRUE, level = 0.99)
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_identical(dim(result), as.integer(c(nrow(newdata), 3)))
  expect_identical(colnames(result), c("fit", "lwr", "upr"))
  expect_identical(logistic5_fn(newdata[, "x"], fit$coefficients), result[, 1])

  result <- predict(fit, newdata, se.fit = TRUE, interval = TRUE, level = 0.99)
  expect_type(result, "list")
  expect_length(result, 5)
  expect_named(result, c("fit", "se.fit", "lwr", "upr", "df"))
  expect_length(result$fit, nrow(newdata))
  expect_length(result$se.fit, nrow(newdata))
  expect_length(result$lwr, nrow(newdata))
  expect_length(result$upr, nrow(newdata))
  expect_length(result$df, 1)
  expect_identical(logistic5_fn(newdata[, "x"], fit$coefficients), result$fit)
})
