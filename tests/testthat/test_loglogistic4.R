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

  start <- c(0, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1)
  upper_bound <- c(3, 2, 2, 5)

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_null(object$start)
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- loglogistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5)))

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

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_null(object$start)
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- loglogistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5)))
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
    loglogistic4_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, 0)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, -1)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2)

  true_value <- c(
    0.04, -0.41, -0.76, -0.82785714285714286, -0.84615384615384615,
    -0.85285714285714286
  )

  value <- loglogistic4_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic4"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic4_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # delta
      0, 0.5, 0.88888888888888889, 0.96428571428571429, 0.98461538461538462,
      0.99206349206349206,
      # log_eta
      0, 0, -0.18483924814931875, -0.10215412735191989, -0.056698619077282154,
      -0.034214241335758937,
      # log_phi
      0, 0.675, 0.26666666666666667, 0.092984693877551020, 0.040899408284023669,
      0.021258503401360544
    ),
    nrow = 6,
    ncol = 4
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, 6),
      # (alpha, delta)
      rep(0, 6),
      # (alpha, log_eta)
      rep(0, 6),
      # (alpha, log_phi)
      rep(0, 6),
      # (delta, alpha)
      rep(0, 6),
      # (delta, delta)
      rep(0, 6),
      # (delta, log_eta)
      0, 0, 0.20537694238813194, 0.11350458594657766, 0.062998465641424615,
      0.038015823706398818,
      # (delta, log_phi)
      0, -0.75000000000000000, -0.29629629629629630, -0.10331632653061224,
      -0.045443786982248521, -0.023620559334845049,
      # (log_eta, alpha)
      rep(0, 6),
      # (log_eta, delta)
      0, 0, 0.20537694238813194, 0.11350458594657766, 0.062998465641424615,
      0.038015823706398818,
      # (log_eta, log_eta)
      0, 0, 0.11410929384422880, 0.21048040166468468, 0.17184883395368479,
      0.12836067406225549,
      # (log_eta, log_phi)
      0, 0.67500000000000000, -0.16462491234841041, -0.19158751803136868,
      -0.12396273026376598, -0.079754971018499173,
      # (log_phi, alpha)
      rep(0, 6),
      # (log_phi, delta)
      0, -0.75000000000000000, -0.29629629629629630, -0.10331632653061224,
      -0.045443786982248521, -0.023620559334845049,
      # (log_phi, log_eta)
      0, 0.67500000000000000, -0.16462491234841041, -0.19158751803136868,
      -0.12396273026376598, -0.079754971018499173,
      # (log_phi, log_phi)
      0, 0, 0.62222222222222222, 0.25902879008746356, 0.11892289485662267,
      0.062763200518302559
    ),
    dim = c(6, 4, 4)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic4"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 4)
  expect_length(gradient_hessian$H, 6 * 4 * 4)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

test_that("Value of the RSS", {
  x <- c(0, 2, 4, 6, 8)
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, -9 / 10, log(3), log(2))

  true_value <- 19.276313893957252

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic4"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, log(2))
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(-9 / 10, log(3)))

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

  theta <- c(4 / 100, -9 / 10, log(3), log(2))

  true_gradient <- c(
    -16.612290109890110, -11.322024907083149, 1.1701819464814164,
    -3.8609372916475908
  )

  true_hessian <- matrix(
    c(
      # alpha
      15, 10.088766788766789, -0.94839086293816353, 3.0529703336956084,
      # beta
      10.088766788766789, 8.9580370355461931, -2.1903073311110781,
      6.2559702057947107,
      # eta
      -0.94839086293816353, -2.1903073311110781, -1.5611816136875038,
      -1.0026466887351753,
      # phi
      3.0529703336956084, 6.2559702057947107, -1.0026466887351753,
      -1.8579749594331474
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic4"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 4)
  expect_length(gradient_hessian$H, 4 * 4)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, log(2))
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(c(-9 / 10, log(3)))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[c(2, 3)])
  expect_equal(gradient_hessian$H, true_hessian[c(2, 3), c(2, 3)])
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

  theta <- c(0, 1, 1.8115476511704276, 1.8210577595704281)

  true_value <- c(
    0.92551260060548019, -0.80876869456840958, 1.8115476511704276,
    1.8210577595704281
  )

  object <- loglogistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(result, true_value)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.92551260060548019,
    delta = -0.80876869456840958,
    eta = exp(1.8115476511704276),
    phi = exp(1.8210577595704281)
  )

  rss_value <- 0.025373136628166469

  fitted_values <- c(
    rep(0.92551260060548019, 3), rep(0.92470055974557833, 2),
    rep(0.8726743286848762, 2), rep(0.557285088913756, 5),
    rep(0.254730718198142, 3), rep(0.157087720976285, 4),
    0.116743938254306
  )

  residuals <- c(
    0.00248739939451981, -0.03751260060548019, 0.05448739939451981,
    0.02329944025442167, -0.06870055974557833, 0.0243256713151238,
    0.0103256713151238, -0.069285088913756, -0.025285088913756,
    0.028714911086244, 0.008714911086244, 0.041714911086244, 0.004269281801858,
    0.010269281801858, -0.011730718198142, -0.040087720976285,
    -0.014087720976285, 0.020912279023715, 0.061912279023715, -0.024743938254306
  )

  object <- loglogistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic4_new(x, y, w, c(0, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.92551260060548019,
    delta = -0.80876869456840958,
    eta = exp(1.8115476511704276),
    phi = exp(1.8210577595704281)
  )

  rss_value <- 0.025373136628166469

  fitted_values <- c(
    rep(0.92551260060548019, 3), rep(0.92470055974557833, 2),
    rep(0.8726743286848762, 2), rep(0.557285088913756, 5),
    rep(0.254730718198142, 3), rep(0.157087720976285, 4),
    0.116743938254306
  )

  residuals <- c(
    0.00248739939451981, -0.03751260060548019, 0.05448739939451981,
    0.02329944025442167, -0.06870055974557833, 0.0243256713151238,
    0.0103256713151238, -0.069285088913756, -0.025285088913756,
    0.028714911086244, 0.008714911086244, 0.041714911086244, 0.004269281801858,
    0.010269281801858, -0.011730718198142, -0.040087720976285,
    -0.014087720976285, 0.020912279023715, 0.061912279023715, -0.024743938254306
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic4_new(
    x, y, w, c(0, 0, 2, 2), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic4_new(
    x, y, w, c(-2, -5, 0.5, 20), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3774037382542154),
    phi = exp(1.8461172076558907)
  )

  rss_value <- 0.066796495711301507

  fitted_values <- c(
    rep(1, 3), rep(0.98975897104553855, 2), rep(0.8609270809729283, 2),
    rep(0.553669054087477, 5), rep(0.283932781514331, 3),
    rep(0.140673135711566, 4), 0.000017760285330
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04175897104553855, -0.13375897104553855,
    0.0360729190270717, 0.0220729190270717, -0.065669054087477,
    -0.021669054087477, 0.032330945912523, 0.012330945912523, 0.045330945912523,
    -0.024932781514331, -0.018932781514331, -0.040932781514331,
    -0.023673135711566, 0.002326864288434, 0.037326864288434, 0.078326864288434,
    0.091982239714670
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  # initial values with same equalities
  object <- loglogistic4_new(
    x, y, w, c(1, -1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  # initial values with different equalities
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3774037382542154),
    phi = exp(1.8461172076558907)
  )

  rss_value <- 0.066796495711301507

  fitted_values <- c(
    rep(1, 3), rep(0.98975897104553855, 2), rep(0.8609270809729283, 2),
    rep(0.553669054087477, 5), rep(0.283932781514331, 3),
    rep(0.140673135711566, 4), 0.000017760285330
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04175897104553855, -0.13375897104553855,
    0.0360729190270717, 0.0220729190270717, -0.065669054087477,
    -0.021669054087477, 0.032330945912523, 0.012330945912523, 0.045330945912523,
    -0.024932781514331, -0.018932781514331, -0.040932781514331,
    -0.023673135711566, 0.002326864288434, 0.037326864288434, 0.078326864288434,
    0.091982239714670
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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
  object <- loglogistic4_new(
    x, y, w, c(1, -1, 2, 2), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 0.5, 0.5), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.93806056213457730,
    delta = -0.79670860565100278,
    eta = exp(1.8199023333421624),
    phi = exp(1.8020330741835448)
  )

  rss_value <- 0.013912488651588458

  fitted_values <- c(
    rep(0.93806056213457730, 3), rep(0.93721164651371502, 2),
    rep(0.8811887559772096, 2), rep(0.552330070097719, 5),
    rep(0.263176550555172, 3), rep(0.176058052614237, 4),
    0.141351980945642
  )

  residuals <- c(
    -0.01006056213457730, -0.05006056213457730, 0.04193943786542270,
    0.01078835348628498, -0.08121164651371502, 0.0158112440227904,
    0.0018112440227904, -0.064330070097719, -0.020330070097719,
    0.033669929902281, 0.013669929902281, 0.046669929902281,
    -0.004176550555172, 0.001823449444828, -0.020176550555172,
    -0.059058052614237, -0.033058052614237, 0.001941947385763,
    0.042941947385763, -0.049351980945642
  )

  object <- loglogistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic4_new(x, y, w, c(1, -1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.93806056213457730,
    delta = -0.79670860565100278,
    eta = exp(1.8199023333421624),
    phi = exp(1.8020330741835448)
  )

  rss_value <- 0.013912488651588458

  fitted_values <- c(
    rep(0.93806056213457730, 3), rep(0.93721164651371502, 2),
    rep(0.8811887559772096, 2), rep(0.552330070097719, 5),
    rep(0.263176550555172, 3), rep(0.176058052614237, 4),
    0.141351980945642
  )

  residuals <- c(
    -0.01006056213457730, -0.05006056213457730, 0.04193943786542270,
    0.01078835348628498, -0.08121164651371502, 0.0158112440227904,
    0.0018112440227904, -0.064330070097719, -0.020330070097719,
    0.033669929902281, 0.013669929902281, 0.046669929902281,
    -0.004176550555172, 0.001823449444828, -0.020176550555172,
    -0.059058052614237, -0.033058052614237, 0.001941947385763,
    0.042941947385763, -0.049351980945642
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic4_new(
    x, y, w, c(0.3, 0.6, 2, 3), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic4_new(
    x, y, w, c(-2, -5, 0.5, 20), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3311583629301975),
    phi = exp(1.8484421002700265)
  )

  rss_value <- 0.040523243758802340

  fitted_values <- c(
    rep(1, 3), rep(0.98754725333678601, 2), rep(0.8518803828725668, 2),
    rep(0.5534371122379037, 5), rep(0.294333091935288, 3),
    rep(0.151983958045141, 4), 0.000029373464958
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03954725333678601, -0.13154725333678601,
    0.0451196171274332, 0.0311196171274332, -0.0654371122379037,
    -0.0214371122379037, 0.0325628877620963, 0.0125628877620963,
    0.0455628877620963, -0.035333091935288, -0.029333091935288,
    -0.051333091935288, -0.034983958045141, -0.008983958045141,
    0.026016041954859, 0.067016041954859, 0.091970626535042
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  # initial values with same equalities
  object <- loglogistic4_new(
    x, y, w, c(1, -1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  # initial values with different equalities
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 1,
    delta = -1,
    eta = exp(1.3311583629301975),
    phi = exp(1.8484421002700265)
  )

  rss_value <- 0.040523243758802340

  fitted_values <- c(
    rep(1, 3), rep(0.98754725333678601, 2), rep(0.8518803828725668, 2),
    rep(0.5534371122379037, 5), rep(0.294333091935288, 3),
    rep(0.151983958045141, 4), 0.000029373464958
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03954725333678601, -0.13154725333678601,
    0.0451196171274332, 0.0311196171274332, -0.0654371122379037,
    -0.0214371122379037, 0.0325628877620963, 0.0125628877620963,
    0.0455628877620963, -0.035333091935288, -0.029333091935288,
    -0.051333091935288, -0.034983958045141, -0.008983958045141,
    0.026016041954859, 0.067016041954859, 0.091970626535042
  )

  object <- loglogistic4_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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
  object <- loglogistic4_new(
    x, y, w, c(1, -1, 2, 2), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 0.5, 0.5), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
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

  theta <- c(alpha = 4 / 100, delta = -9 / 10, eta = 3, phi = 2)

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5982, 3847.3537344782560, -104.40263338110749, 285.30029971608626,
      275947.34244121722,
      # delta
      3847.3537344782560, 3636.0201791330070, -258.69909747088174,
      654.17722274224129, 193269.92078353599,
      # eta
      -104.40263338110749, -258.69909747088174, -132.47182391898161,
      52.152763193549342, -5742.3762879391817,
      # phi
      285.30029971608626, 654.17722274224129, 52.152763193549342,
      -35.178966118180261, 15953.701512060766,
      # sigma
      275947.34244121722, 193269.92078353599, -5742.3762879391817,
      15953.701512060766, 9955612.1461136607
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- loglogistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
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
    drda(
      y ~ x, mean_function = "loglogistic4", lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2), upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5), upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c( 0, -Inf, -Inf, -Inf), upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5), upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
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
    drda(
      y ~ x, mean_function = "loglogistic4", upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 4), upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5), upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
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
    drda(
      y ~ x, mean_function = "loglogistic4", start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(1, 1, 1, 1, 1)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(0, 1, 1, -1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4", start = c(0, 1, 1, 0)
    ),
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

  result <- drda(y ~ x, mean_function = "loglogistic4")

  expect_equal(nauc(result), 0.63097122559352374)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.92539849554119428)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.64001041004420353)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.45969717209770064
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

  result <- drda(y ~ x, mean_function = "loglogistic4")

  expect_equal(naac(result), 1.0 - 0.63097122559352374)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.92539849554119428)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.64001041004420353)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.45969717209770064
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

  result <- drda(y ~ x, mean_function = "loglogistic4")

  expect_equal(nauc(result), 0.43871949389022659)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.15536599581899617)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.40467437230093021)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.63244069235812449
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

  result <- drda(y ~ x, mean_function = "loglogistic4")

  expect_equal(naac(result), 1.0 - 0.43871949389022659)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.15536599581899617)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.40467437230093021)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.63244069235812449
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})
