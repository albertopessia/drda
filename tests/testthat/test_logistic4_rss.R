context("RSS of the four-parameter logistic model")

test_that("Value of the RSS", {
  # actual sample
  # y <- c(
  #   10204 / 10000, 9501 / 10000, 10375 / 10000,
  #   8026 / 10000, 8272 / 10000, 8430 / 10000,
  #   7046 / 10000, 8214 / 10000, 9764 / 10000,
  #   611 / 10000, 2143 / 10000, 1181 / 10000,
  #   1828 / 10000, 108 / 10000, 814 / 10000
  # )

  x <- -log(c(1000, 100, 10, 1, 0.1))
  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  n <- rep(3, 5)

  # mean by dose
  m <- c(
    376 / 375, 3091 / 3750, 1564 / 1875, 787 / 6000, 11 / 120
  )

  # sample variance by dose
  v <- c(
    643663 / 450000000, 31087 / 112500000, 1394281 / 112500000,
    449671 / 112500000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_value <- 0.66098935323236375

  f <- logistic4_rss(stats)

  expect_type(f, "closure")

  value <- f(param)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_rss_fixed(stats, known_param)

  expect_type(f, "closure")

  value <- f(c(9 / 10, -2))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient of the RSS", {
  # actual sample
  # y <- c(
  #   10204 / 10000, 9501 / 10000, 10375 / 10000,
  #   8026 / 10000, 8272 / 10000, 8430 / 10000,
  #   7046 / 10000, 8214 / 10000, 9764 / 10000,
  #   611 / 10000, 2143 / 10000, 1181 / 10000,
  #   1828 / 10000, 108 / 10000, 814 / 10000
  # )

  x <- -log(c(1000, 100, 10, 1, 0.1))
  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  n <- rep(3, 5)

  # mean by dose
  m <- c(
    376 / 375, 3091 / 3750, 1564 / 1875, 787 / 6000, 11 / 120
  )

  # sample variance by dose
  v <- c(
    643663 / 450000000, 31087 / 112500000, 1394281 / 112500000,
    449671 / 112500000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_value <- c(
    1.7528316378106113, -0.086573687534061886, 0.032911020845508688,
    0.033129036826382914
  )

  f <- logistic4_rss_gradient(stats)

  expect_type(f, "closure")

  gradient <- f(param)

  expect_type(gradient, "double")
  expect_length(gradient, 4)
  expect_equal(gradient, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_rss_gradient_fixed(stats, known_param)

  expect_type(f, "closure")

  gradient <- f(c(9 / 10, -2))

  expect_type(gradient, "double")
  expect_length(gradient, 2)
  expect_equal(gradient, true_value[2:3])
})

test_that("Hessian of the RSS", {
  # actual sample
  # y <- c(
  #   10204 / 10000, 9501 / 10000, 10375 / 10000,
  #   8026 / 10000, 8272 / 10000, 8430 / 10000,
  #   7046 / 10000, 8214 / 10000, 9764 / 10000,
  #   611 / 10000, 2143 / 10000, 1181 / 10000,
  #   1828 / 10000, 108 / 10000, 814 / 10000
  # )

  x <- -log(c(1000, 100, 10, 1, 0.1))
  param <- c(4 / 10, 9 / 10, -2, -3 / 2)

  n <- rep(3, 5)

  # mean by dose
  m <- c(
    376 / 375, 3091 / 3750, 1564 / 1875, 787 / 6000, 11 / 120
  )

  # sample variance by dose
  v <- c(
    643663 / 450000000, 31087 / 112500000, 1394281 / 112500000,
    449671 / 112500000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_value <- matrix(
    c(
      # alpha
      5.8031474358280998, 0.56093666361880140, 0.0057751967022874164,
      0.13423785257079926,
      # beta
      0.56093666361880140, 8.0749792369342974, -0.078452752070457148,
      0.42669881104800214,
      # eta
      0.0057751967022874164, -0.078452752070457148, 0.056075449207447959,
      0.017625413119076090,
      # phi
      0.13423785257079926, 0.42669881104800214, 0.017625413119076090,
      0.14602991636978138
    ),
    nrow = 4,
    ncol = 4
  )

  f <- logistic4_rss_hessian(stats)

  expect_type(f, "closure")

  hessian <- f(param)

  expect_type(hessian, "double")
  expect_length(hessian, 4 * 4)
  expect_equal(hessian, true_value)

  known_param <- c(4 / 10, NA, NA, -3 / 2)
  f <- logistic4_rss_hessian_fixed(stats, known_param)

  expect_type(f, "closure")

  hessian <- f(c(9 / 10, -2))

  expect_type(hessian, "double")
  expect_length(hessian, 2 * 2)
  expect_equal(hessian, true_value[2:3, 2:3])
})
