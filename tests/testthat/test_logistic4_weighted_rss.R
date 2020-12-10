context("Weighted RSS of the four-parameter logistic model")

test_that("Value of the weighted RSS", {
  x <- -rep(log(c(1000, 100, 10, 1, 0.1)), each = 3)

  y <- c(
    10204 / 10000, 9501 / 10000, 10375 / 10000,
    8026 / 10000, 8272 / 10000, 8430 / 10000,
    7046 / 10000, 8214 / 10000, 9764 / 10000,
    611 / 10000, 2143 / 10000, 1181 / 10000,
    1828 / 10000, 108 / 10000, 814 / 10000
  )

  w <- c(
    1 / 2, 1 / 3, 3 / 2,
    2 / 3, 1 / 2, 1 / 2,
    2 / 3, 2 / 3, 1,
    1, 1, 1,
    1 / 2, 1 / 2, 1 / 2
  )

  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  true_value <- 0.49979565465993548

  f <- logistic4_weighted_rss(x, y, w)

  expect_type(f, "closure")

  value <- f(param)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_weighted_rss_fixed(x, y, w, known_param)

  expect_type(f, "closure")

  value <- f(c(9 / 10, -2))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient of the weighted RSS", {
  x <- -rep(log(c(1000, 100, 10, 1, 0.1)), each = 3)

  y <- c(
    10204 / 10000, 9501 / 10000, 10375 / 10000,
    8026 / 10000, 8272 / 10000, 8430 / 10000,
    7046 / 10000, 8214 / 10000, 9764 / 10000,
    611 / 10000, 2143 / 10000, 1181 / 10000,
    1828 / 10000, 108 / 10000, 814 / 10000
  )

  w <- c(
    1 / 2, 1 / 3, 3 / 2,
    2 / 3, 1 / 2, 1 / 2,
    2 / 3, 2 / 3, 1,
    1, 1, 1,
    1 / 2, 1 / 2, 1 / 2
  )

  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  true_value <- c(
    1.2840456550254584, -0.18741445858679887, 0.034758912990703695,
    0.027751482494782981
  )

  f <- logistic4_weighted_rss_gradient(x, y, w)

  expect_type(f, "closure")

  gradient <- f(param)

  expect_type(gradient, "double")
  expect_length(gradient, 4)
  expect_equal(gradient, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_weighted_rss_gradient_fixed(x, y, w, known_param)

  expect_type(f, "closure")

  gradient <- f(c(9 / 10, -2))

  expect_type(gradient, "double")
  expect_length(gradient, 2)
  expect_equal(gradient, true_value[2:3])
})

test_that("Hessian of the RSS", {
  x <- -rep(log(c(1000, 100, 10, 1, 0.1)), each = 3)

  y <- c(
    10204 / 10000, 9501 / 10000, 10375 / 10000,
    8026 / 10000, 8272 / 10000, 8430 / 10000,
    7046 / 10000, 8214 / 10000, 9764 / 10000,
    611 / 10000, 2143 / 10000, 1181 / 10000,
    1828 / 10000, 108 / 10000, 814 / 10000
  )

  w <- c(
    1 / 2, 1 / 3, 3 / 2,
    2 / 3, 1 / 2, 1 / 2,
    2 / 3, 2 / 3, 1,
    1, 1, 1,
    1 / 2, 1 / 2, 1 / 2
  )

  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  true_value <- matrix(
    c(
      # alpha
      4.2859839145917578, 0.46465369253092314, 0.0069024697646605563,
      0.12871077587714882,
      # beta
      0.46465369253092314, 5.6180420336797293, -0.039558533363869521,
      0.33594291665377432,
      # eta
      0.0069024697646605563, -0.039558533363869521, 0.054293021835232003,
      0.021398969711958365,
      # phi
      0.12871077587714882, 0.33594291665377432, 0.021398969711958365,
      0.13960836417931997
    ),
    nrow = 4,
    ncol = 4
  )

  f <- logistic4_weighted_rss_hessian(x, y, w)

  expect_type(f, "closure")

  hessian <- f(param)

  expect_type(hessian, "double")
  expect_length(hessian, 4 * 4)
  expect_equal(hessian, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_weighted_rss_hessian_fixed(x, y, w, known_param)

  expect_type(f, "closure")

  hessian <- f(c(9 / 10, -2))

  expect_type(hessian, "double")
  expect_length(hessian, 2 * 2)
  expect_equal(hessian, true_value[2:3, 2:3])
})
