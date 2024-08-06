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
