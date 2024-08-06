test_that("curve_variance: call", {
  result <- drda(y ~ x, data = ltd$D)

  expect_length(curve_variance(result, c(-1, 0, 1)), 3)
  expect_length(curve_variance(result, ltd$D$x), length(ltd$D$x))
})

test_that("curve_variance: fixed params", {
  result <- drda(
    y ~ x, data = ltd$D, lower_bound = c(0.85, -0.75, -Inf, -Inf),
    upper_bound = c(0.85, -0.75, Inf, Inf)
  )

  expect_length(curve_variance(result, c(-1, 0, 1)), 3)
  expect_length(curve_variance(result, ltd$D$x), length(ltd$D$x))
})
