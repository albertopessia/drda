test_that("plot show=FALSE with level=0 does not crash", {
  result <- drda(y ~ x, data = ltd$D, mean_function = "logistic4")

  expect_silent({
    params <- plot(result, show = FALSE, level = 0)
  })
  expect_true(is.list(params))
  expect_true(all(is.na(params$confidence_interval$y_lower)))
  expect_true(all(is.na(params$confidence_interval$y_upper)))
})
