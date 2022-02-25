test_that("Constructor", {
  x <- rep(c(0, 10^seq(-3, 2)), times = c(3, 2, 2, 5, 3, 4, 1))

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- matrix(
    c(
      0, 0.001, 0.01, 0.1, 1, 10, 100, 3, 2, 2, 5, 3, 4, 1, 0.932, 0.902, 0.89,
      0.5542, 0.2556666667, 0.16425, 0.092, 0.0014186667, 0.002116, 0.000049,
      0.00160656, 0.0000862222, 0.0014676875, 0
    ),
    nrow = 7,
    ncol = 4
  )
  colnames(stats) <- c("x", "n", "m", "v")

  start <- c(1, 1)

  lower_bound <- c(0.5, 1)
  upper_bound <- c(2, 5)

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, 0, 0))
  expect_equal(object$lower_bound, c(log(0.5), 0))
  expect_equal(object$upper_bound, c(log(2), log(5)))

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  stats <- matrix(
    c(
      0, 0.001, 0.01, 0.1, 1, 10, 100, 4.549, 0.96, 1.126, 3.756, 1.797, 2.61,
      0.157, 0.9353000659, 0.948, 0.8836838366, 0.55221459, 0.2606149137,
      0.1807233716, 0.092, 0.0014467345, 0, 0.0000091061, 0.0007707846,
      0.0000597738, 0.0014230308, 0
    ),
    nrow = 7,
    ncol = 4
  )
  colnames(stats) <- c("x", "n", "m", "v")

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, 0, 0))
  expect_equal(object$lower_bound, c(log(0.5), 0))
  expect_equal(object$upper_bound, c(log(2), log(5)))
})

test_that("Constructor: errors", {
  x <- rep(c(0, 10^seq(-3, 2)), times = c(3, 2, 2, 5, 3, 4, 1))

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  max_iter <- 10000

  expect_error(
    loglogistic2_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(-1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, -Inf, Inf),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, -Inf, rep(Inf, 2)),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), Inf),
    "'upper_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(
      x, y, w, NULL, max_iter, rep(-Inf, 2), c(0, Inf)
    ),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(
      x, y, w, NULL, max_iter, rep(-Inf, 2), c(-1, Inf)
    ),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(
      x, y, w, NULL, max_iter, rep(-Inf, 2), c(Inf, 0)
    ),
    "'upper_bound[2]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(
      x, y, w, NULL, max_iter, rep(-Inf, 2), c(Inf, -1)
    ),
    "'upper_bound[2]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value (increasing)", {
  x <- c(0, 2, 4, 6, 8, 10)

  true_value <- c(
    0, 0.5, 0.88888888888888889, 0.96428571428571429, 0.98461538461538462,
    0.99206349206349206
  )

  value <- loglogistic2_fn(x, c(0, 1, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  value <- loglogistic2_fn(x, c(-0.5, 3, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(
      stats = matrix(x, nrow = 6, ncol = 1),
      start = c(0, 1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  value <- fn(object, object$stats[, 1], c(3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic2_fit"
  )

  value <- fn(object, object$stats[, 1], c(0, 1, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Function value (decreasing)", {
  x <- c(0, 2, 4, 6, 8, 10)

  true_value <- c(
    1, 0.5, 0.11111111111111111, 0.035714285714285714, 0.015384615384615385,
    0.0079365079365079365
  )

  value <- loglogistic2_fn(x, c(1, -1, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  value <- loglogistic2_fn(x, c(0.5, -3, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(
      stats = matrix(x, nrow = 6, ncol = 1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  value <- fn(object, object$stats[, 1], c(3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic2_fit"
  )

  value <- fn(object, object$stats[, 1], c(1, -1, 3, 2))

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(3, 2)

  true_gradient <- matrix(
    c(
      # log_eta
      0, 0, 0.20537694238813194, 0.11350458594657766, 0.062998465641424615,
      0.038015823706398818,
      # log_phi
      0, -0.75, -0.29629629629629630, -0.10331632653061224,
      -0.045443786982248521, -0.023620559334845049
    ),
    nrow = 6,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (log_eta, log_eta)
      0, 0, -0.12678810427136534, -0.23386711296076076, -0.19094314883742754,
      -0.14262297118028387,
      # (log_eta, log_phi)
      0, -0.75000000000000000, 0.18291656927601157, 0.21287502003485409,
      0.13773636695973998, 0.088616634464999081,
      # (log_phi, log_eta)
      0, -0.75000000000000000, 0.18291656927601157, 0.21287502003485409,
      0.13773636695973998, 0.088616634464999081,
      # (log_phi, log_phi)
      0, 0, -0.69135802469135802, -0.28780976676384840, -0.13213654984069185,
      -0.069736889464780621
    ),
    dim = c(6, 2, 2)
  )

  object <- structure(
    list(
      stats = matrix(x, nrow = 6, ncol = 1),
      start = c(0, 1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 2)
  expect_length(gradient_hessian$H, 6 * 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  object <- structure(
    list(
      stats = matrix(x, nrow = 6, ncol = 1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 2)
  expect_length(gradient_hessian$H, 6 * 2 * 2)

  expect_equal(gradient_hessian$G, -true_gradient)
  expect_equal(gradient_hessian$H, -true_hessian)
})

test_that("Value of the RSS", {
  x <- c(0, 2, 4, 6, 8)
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(log(3), log(2))

  true_value <- 1.6216590112971088

  object <- structure(
    list(
      stats = cbind(x, n, m, v), m = 5, start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(NA_real_, log(2))
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(log(3))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  x <- c(0, 2, 4, 6, 8)
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(log(3), log(2))

  true_gradient <- c(0.38748581655130477, -1.2518775105266555)

  true_hessian <- matrix(
    c(
      # eta
      -0.19761486113458578, -0.49427996475403962,
      # phi
      -0.49427996475403962, 0.66697956359255129
    ),
    nrow = 2,
    ncol = 2
  )

  object <- structure(
    list(
      stats = cbind(x, n, m, v), m = 5, start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(NA, log(2))
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(log(3))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 1)
  expect_length(gradient_hessian$H, 1)

  expect_equal(gradient_hessian$G, true_gradient[1])
  expect_equal(gradient_hessian$H, true_hessian[1, 1, drop = FALSE])
})

test_that("mle_asy", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  theta <- c(log(3), log(2))

  object <- loglogistic2_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_equal(result, theta)
})

test_that("fit", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3774037382542154),
    phi = exp(1.8461172076558907)
  )

  rss_value <- 0.066796495711301507

  fitted_values <- c(
    rep(1, 3), rep(0.98975897104553856, 2),
    rep(0.8609270809729283, 2), rep(0.553669054087477, 5),
    rep(0.283932781514331, 3), rep(0.140673135711566, 4),
    0.000017760285330
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04175897104553856, -0.13375897104553856,
    0.0360729190270717, 0.0220729190270717, -0.065669054087477,
    -0.021669054087477, 0.032330945912523, 0.012330945912523, 0.045330945912523,
    -0.024932781514331, -0.018932781514331, -0.040932781514331,
    -0.023673135711566, 0.002326864288434, 0.037326864288434, 0.078326864288434,
    0.091982239714670
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic2_new(x, y, w, c(1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3774037382542154),
    phi = exp(1.8461172076558907)
  )

  rss_value <- 0.066796495711301507

  fitted_values <- c(
    rep(1, 3), rep(0.98975897104553856, 2),
    rep(0.8609270809729283, 2), rep(0.553669054087477, 5),
    rep(0.283932781514331, 3), rep(0.140673135711566, 4),
    0.000017760285330
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04175897104553856, -0.13375897104553856,
    0.0360729190270717, 0.0220729190270717, -0.065669054087477,
    -0.021669054087477, 0.032330945912523, 0.012330945912523, 0.045330945912523,
    -0.024932781514331, -0.018932781514331, -0.040932781514331,
    -0.023673135711566, 0.002326864288434, 0.037326864288434, 0.078326864288434,
    0.091982239714670
  )

  object <- loglogistic2_new(
    x, y, w, NULL, 10000, c(1, 1), c(10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(x, y, w, c(2, 2), 10000, c(1, 1), c(10, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(x, y, w, c(0.5, 20), 10000, c(1, 1), c(10, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = 4,
    phi = exp(1.8465108937874798)
  )

  rss_value <- 0.066824532927887730

  fitted_values <- c(
    rep(1, 3), rep(0.990179895281417676, 2), rep(0.86305111097279304, 2),
    rep(0.5545336299750062, 5), rep(0.2825753687406050, 3),
    rep(0.13891908975957001, 4), 0.00001613284500299034
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.042179895281417676, -0.134179895281417676,
    0.03394888902720696, 0.01994888902720696, -0.0665336299750062,
    -0.0225336299750062, 0.0314663700249938, 0.0114663700249938,
    0.0444663700249938, -0.0235753687406050, -0.0175753687406050,
    -0.0395753687406050, -0.02191908975957001, 0.00408091024042999,
    0.03908091024042999, 0.08008091024042999, 0.09198386715499700966
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic2_new(x, y, w, c(4, 1), 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic2_new(x, y, w, c(2, 1), 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities and inequalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = 4,
    phi = exp(1.8465108937874798)
  )

  rss_value <- 0.066824532927887730

  fitted_values <- c(
    rep(1, 3), rep(0.990179895281417676, 2), rep(0.86305111097279304, 2),
    rep(0.5545336299750062, 5), rep(0.2825753687406050, 3),
    rep(0.13891908975957001, 4), 0.00001613284500299034
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.042179895281417676, -0.134179895281417676,
    0.03394888902720696, 0.01994888902720696, -0.0665336299750062,
    -0.0225336299750062, 0.0314663700249938, 0.0114663700249938,
    0.0444663700249938, -0.0235753687406050, -0.0175753687406050,
    -0.0395753687406050, -0.02191908975957001, 0.00408091024042999,
    0.03908091024042999, 0.08008091024042999, 0.09198386715499700966
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(x, y, w, c(4, 2), 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(x, y, w, c(0.5, 0.5), 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit (weighted)", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3311583628023397),
    phi = exp(1.8484421002737211)
  )

  rss_value <- 0.040523243758802340

  fitted_values <- c(
    rep(1, 3), rep(0.98754725333008165, 2),
    rep(0.8518803828461078, 2), rep(0.5534371122345799, 5),
    rep(0.294333091961414, 3), rep(0.151983958075273, 4),
    0.000029373464997
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03954725333008165, -0.13154725333008165,
    0.0451196171538922, 0.0311196171538922, -0.0654371122345799,
    -0.0214371122345799, 0.0325628877654201, 0.0125628877654201,
    0.0455628877654201, -0.035333091961414, -0.029333091961414,
    -0.051333091961414, -0.034983958075273, -0.008983958075273,
    0.026016041924727, 0.067016041924727, 0.091970626535003
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic2_new(x, y, w, c(1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): inequalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3311583628023397),
    phi = exp(1.8484421002737211)
  )

  rss_value <- 0.040523243758802340

  fitted_values <- c(
    rep(1, 3), rep(0.98754725333008165, 2),
    rep(0.8518803828461078, 2), rep(0.5534371122345799, 5),
    rep(0.294333091961414, 3), rep(0.151983958075273, 4),
    0.000029373464997
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03954725333008165, -0.13154725333008165,
    0.0451196171538922, 0.0311196171538922, -0.0654371122345799,
    -0.0214371122345799, 0.0325628877654201, 0.0125628877654201,
    0.0455628877654201, -0.035333091961414, -0.029333091961414,
    -0.051333091961414, -0.034983958075273, -0.008983958075273,
    0.026016041924727, 0.067016041924727, 0.091970626535003
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, c(1, 1), c(10, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(x, y, w, c(2, 3), 10000, c(1, 1), c(10, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(x, y, w, c(0.5, 20), 10000, c(1, 1), c(10, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = 4,
    phi = exp(1.8502402740371774)
  )

  rss_value <- 0.041235903961182446

  fitted_values <- c(
    rep(1, 3), rep(0.990323892718425514, 2), rep(0.86480474095238314, 2),
    rep(0.5582155831860413, 5), rep(0.2856093302034027, 3),
    rep(0.14071316259945311, 4), 0.00001637530708649895
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.042323892718425514, -0.134323892718425514,
    0.03219525904761686, 0.01819525904761686, -0.0702155831860413,
    -0.0262155831860413, 0.0277844168139587, 0.0077844168139587,
    0.0407844168139587, -0.0266093302034027, -0.0206093302034027,
    -0.0426093302034027, -0.02371316259945311, 0.00228683740054689,
    0.03728683740054689, 0.07828683740054689, 0.09198362469291350105
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic2_new(x, y, w, c(4, 1), 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic2_new(x, y, w, c(1, 1), 10000, c(4, -Inf), c(4, Inf))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities and inequalities", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = 4,
    phi = exp(1.8502402740371774)
  )

  rss_value <- 0.041235903961182446

  fitted_values <- c(
    rep(1, 3), rep(0.990323892718425514, 2), rep(0.86480474095238314, 2),
    rep(0.5582155831860413, 5), rep(0.2856093302034027, 3),
    rep(0.14071316259945311, 4), 0.00001637530708649895
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.042323892718425514, -0.134323892718425514,
    0.03219525904761686, 0.01819525904761686, -0.0702155831860413,
    -0.0262155831860413, 0.0277844168139587, 0.0077844168139587,
    0.0407844168139587, -0.0266093302034027, -0.0206093302034027,
    -0.0426093302034027, -0.02371316259945311, 0.00228683740054689,
    0.03728683740054689, 0.07828683740054689, 0.09198362469291350105
  )

  object <- loglogistic2_new(x, y, w, NULL, 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(x, y, w, c(4, 2), 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(x, y, w, c(4, 0.5), 10000, c(4, 1), c(4, 10))

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  theta <- c(alpha = 1, delta = -1, eta = 3, phi = 2)

  sigma <- 0.05

  true_value <- matrix(c(
      # eta
      -46.642067026144702, 15.343013377087543, -2366.7441624635119,
      # phi
      15.343013377087543, 24.608783839908594, 6491.3837695465929,
      # sigma
      -2366.7441624635119, 6491.3837695465929, 993480.12792657081
    ),
    nrow = 3,
    ncol = 3
  )

  rownames(true_value) <- colnames(true_value) <- c("eta", "phi", "sigma")

  object <- loglogistic2_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 3 * 3)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", lower_bound = c("a", "b")),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = matrix(-Inf, nrow = 2, ncol = 2), upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3), upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c(2, -Inf), upper_bound = c(1, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c(Inf, -Inf), upper_bound = c(Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3), upper_bound = rep(Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", upper_bound = c("a", "b")),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 2), upper_bound = matrix(Inf, nrow = 2, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c(-Inf, -Inf), upper_bound = c(-Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3), upper_bound = rep(Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'start' argument errors", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c("a", "b")),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(Inf, 1)),
    "'start' must be finite"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(-Inf, 1)),
    "'start' must be finite"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(1, 1, 1)),
    "'start' must be of length 2"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(-1, 1)),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(0, 1)),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(1, -1)),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(y ~ x, mean_function = "loglogistic2", start = c(1, 0)),
    "parameter 'phi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "loglogistic2")

  expect_equal(nauc(result), 0.65330597222122364)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.99792774825879220)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.64500935253752706)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.47352740160136782
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})

test_that("naac: decreasing", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "loglogistic2")

  expect_equal(naac(result), 1.0 - 0.65330597222122364)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.99792774825879220)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.64500935253752706)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.47352740160136782
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 1.0)
})

test_that("nauc: increasing", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "loglogistic2")

  expect_equal(nauc(result), 0.39437658083386260)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.0057169814918619080)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.41147487335299511)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.65461731829532533
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 1.0)
})

test_that("naac: increasing", {
  x <- rep(
    c(0, 2, 4, 6, 8, 10, 100), times = c(3, 2, 2, 5, 3, 4, 1)
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "loglogistic2")

  expect_equal(naac(result), 1.0 - 0.39437658083386260)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.0057169814918619080)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.41147487335299511)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.65461731829532533
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})
