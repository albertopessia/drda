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

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loggompertz"))
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
    loggompertz_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(
      x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, 0)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(
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
    0.04, -0.29109149705429809, -0.75424721232613586, -0.82727639987115761,
    -0.84604679330486757, -0.85282872335335457
  )

  value <- loggompertz_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loggompertz"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loggompertz_fit"
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
      0, 0.36787944117144232, 0.88249690258459540, 0.96364044430128623,
      0.98449643700540841, 0.99203191483706063,
      # log_eta
      0, 0, -0.20644883095929649, -0.10586672339669901, -0.057577578433448398,
      -0.034486857520200622,
      # log_phi
      0, 0.99327449116289427, 0.29784270462230095, 0.096364044430128623,
      0.041533443436165667, 0.021427889360480510
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
      0, 0, 0.22938758995477388, 0.11762969266299891, 0.063975087148275998,
      0.038318730578000691,
      # (delta, log_phi)
      0, -1.1036383235143270, -0.33093633846922328, -0.10707116047792069,
      -0.046148270484628519, -0.023808765956089455,
      # (log_eta, alpha)
      rep(0, 6),
      # (log_eta, delta)
      0, 0, 0.22938758995477388, 0.11762969266299891, 0.063975087148275998,
      0.038318730578000691,
      # (log_eta, log_eta)
      0, 0, 0.16918715995270293, 0.23012978387004229, 0.17813930072973977,
      0.13069440345760184,
      # (log_eta, log_phi)
      0, 0.99327449116289427, -0.24408547664585235, -0.20947315649366853,
      -0.12850034287511163, -0.081204998619636540,
      # (log_phi, alpha)
      rep(0, 6),
      # (log_phi, delta)
      0, -1.1036383235143270, -0.33093633846922328, -0.10707116047792069,
      -0.046148270484628519, -0.023808765956089455,
      # (log_phi, log_eta)
      0, 0.99327449116289427, -0.24408547664585235, -0.20947315649366853,
      -0.12850034287511163, -0.081204998619636540,
      # (log_phi, log_phi)
      0, 0, 0.78183709963353999, 0.27838501724259380, 0.12265345014742674,
      0.063769398736789997
    ),
    dim = c(6, 4, 4)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loggompertz"
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

  true_value <- 18.394916331375033

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loggompertz"
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
    -16.241414895214399, -10.664816362826166, 1.2561542563386518,
    -4.7998198260950907
  )

  true_hessian <- matrix(
    c(
      # alpha
      15, 9.6766832169048879, -1.0090972908057342, 4.0855653907622962,
      # delta
      9.6766832169048879, 8.5857187428544014, -2.3382324471552245,
      7.4491492064911019,
      # log_eta
      -1.0090972908057342, -2.3382324471552245, -1.8152174744087680,
      -1.5116103497794037,
      # log_phi
      4.0855653907622962, 7.4491492064911019, -1.5116103497794037,
      -0.83283728941545387
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loggompertz"
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

  theta <- c(0, 1, 1.4803633738747118, 1.7465438898240248)

  true_value <- c(
    0.91370245143275957, -0.82104372701784045, 1.4803633738747118,
    1.7465438898240248
  )

  object <- loggompertz_new(x, y, w, NULL, 10000, NULL, NULL)

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
    alpha = 0.91370245143275957,
    delta = -0.82104372701784045,
    eta = exp(1.4803633738747118),
    phi = exp(1.7465438898240248)
  )

  rss_value <- 0.024834712402586844

  fitted_values <- c(
    rep(0.91370245143275957, 3), rep(0.91370245143275957, 2),
    rep(0.9074034438609566, 2), rep(0.552014427473422, 5),
    rep(0.2623699947389637, 3), rep(0.1609592835840543, 4),
    0.0926615991940691
  )

  residuals <- c(
    0.01429754856724043, -0.02570245143275957, 0.06629754856724043,
    0.03429754856724043, -0.05770245143275957, -0.0104034438609566,
    -0.0244034438609566, -0.064014427473422, -0.020014427473422,
    0.033985572526578, 0.013985572526578, 0.046985572526578,
    -0.0033699947389637, 0.0026300052610363, -0.0193699947389637,
    -0.0439592835840543, -0.0179592835840543, 0.0170407164159457,
    0.0580407164159457, -0.0006615991940691
  )

  object <- loggompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  object <- loggompertz_new(x, y, w, c(0, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    alpha = 0.91370245143275957,
    delta = -0.82104372701784045,
    eta = exp(1.4803633738747118),
    phi = exp(1.7465438898240248)
  )

  rss_value <- 0.024834712402586844

  fitted_values <- c(
    rep(0.91370245143275957, 3), rep(0.91370245143275957, 2),
    rep(0.9074034438609566, 2), rep(0.552014427473422, 5),
    rep(0.2623699947389637, 3), rep(0.1609592835840543, 4),
    0.0926615991940691
  )

  residuals <- c(
    0.01429754856724043, -0.02570245143275957, 0.06629754856724043,
    0.03429754856724043, -0.05770245143275957, -0.0104034438609566,
    -0.0244034438609566, -0.064014427473422, -0.020014427473422,
    0.033985572526578, 0.013985572526578, 0.046985572526578,
    -0.0033699947389637, 0.0026300052610363, -0.0193699947389637,
    -0.0439592835840543, -0.0179592835840543, 0.0170407164159457,
    0.0580407164159457, -0.0006615991940691
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 0, 2, 2), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(-2, -5, 0.5, 20), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    eta = exp(1.0663677629450972),
    phi = exp(1.7022792018028630)
  )

  rss_value <- 0.069741177998016680

  fitted_values <- c(
    rep(1, 3), rep(0.9999999928250493173528, 2), rep(0.9182398745646568, 2),
    rep(0.5375009924562357, 5), rep(0.2841925412906962, 3),
    rep(0.1604247796868143, 4), 0.0002176866854686176
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.0519999928250493173528,
    -0.1439999928250493173528, -0.0212398745646568, -0.0352398745646568,
    -0.0495009924562357, -0.0055009924562357, 0.0484990075437643,
    0.0284990075437643, 0.0614990075437643, -0.0251925412906962,
    -0.0191925412906962, -0.0411925412906962, -0.0434247796868143,
    -0.0174247796868143, 0.0175752203131857, 0.0585752203131857,
    0.0917823133145313824
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(1, -1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    eta = exp(1.0663677629450972),
    phi = exp(1.7022792018028630)
  )

  rss_value <- 0.069741177998016680

  fitted_values <- c(
    rep(1, 3), rep(0.9999999928250493173528, 2), rep(0.9182398745646568, 2),
    rep(0.5375009924562357, 5), rep(0.2841925412906962, 3),
    rep(0.1604247796868143, 4), 0.0002176866854686176
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.0519999928250493173528,
    -0.1439999928250493173528, -0.0212398745646568, -0.0352398745646568,
    -0.0495009924562357, -0.0055009924562357, 0.0484990075437643,
    0.0284990075437643, 0.0614990075437643, -0.0251925412906962,
    -0.0191925412906962, -0.0411925412906962, -0.0434247796868143,
    -0.0174247796868143, 0.0175752203131857, 0.0585752203131857,
    0.0917823133145313824
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(1, -1, 2, 2), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 0.5, 0.5), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    alpha = 0.93332069218009519,
    delta = -0.85244786873956615,
    eta = exp(1.3247858062455946),
    phi = exp(1.7286245728947466)
  )

  rss_value <- 0.015372205748979479

  fitted_values <- c(
    rep(0.93332069218009519, 3), rep(0.93332069218009519, 2),
    rep(0.9105855700196264, 2), rep(0.5459062976415841, 5),
    rep(0.2807900729515784, 3), rep(0.1738223869604868, 4),
    0.0808898721981491
  )

  residuals <- c(
    -0.00532069218009519, -0.04532069218009519, 0.04667930781990481,
    0.01467930781990481, -0.07732069218009519, -0.0135855700196264,
    -0.0275855700196264, -0.0579062976415841, -0.0139062976415841,
    0.0400937023584159, 0.0200937023584159, 0.0530937023584159,
    -0.0217900729515784, -0.0157900729515784, -0.0377900729515784,
    -0.0568223869604868, -0.0308223869604868, 0.0041776130395132,
    0.0451776130395132, 0.0111101278018509
  )

  object <- loggompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  object <- loggompertz_new(x, y, w, c(1, -1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    alpha = 0.93332069218009519,
    delta = -0.85244786873956615,
    eta = exp(1.3247858062455946),
    phi = exp(1.7286245728947466)
  )

  rss_value <- 0.015372205748979479

  fitted_values <- c(
    rep(0.93332069218009519, 3), rep(0.93332069218009519, 2),
    rep(0.9105855700196264, 2), rep(0.5459062976415841, 5),
    rep(0.2807900729515784, 3), rep(0.1738223869604868, 4),
    0.0808898721981491
  )

  residuals <- c(
    -0.00532069218009519, -0.04532069218009519, 0.04667930781990481,
    0.01467930781990481, -0.07732069218009519, -0.0135855700196264,
    -0.0275855700196264, -0.0579062976415841, -0.0139062976415841,
    0.0400937023584159, 0.0200937023584159, 0.0530937023584159,
    -0.0217900729515784, -0.0157900729515784, -0.0377900729515784,
    -0.0568223869604868, -0.0308223869604868, 0.0041776130395132,
    0.0451776130395132, 0.0111101278018509
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0.3, 0.6, 2, 3), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(-2, -5, 0.5, 20), 10000, c(-1, -3, 1, 1), c(1, 3, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    eta = exp(1.0192530558433254),
    phi = exp(1.6997766671166766)
  )

  rss_value <- 0.039862989913205782

  fitted_values <- c(
    rep(1, 3), rep(0.999999914333498004037, 2), rep(0.9078010022598841, 2),
    rep(0.5392949652230244, 5), rep(0.2947508809799500, 3),
    rep(0.1715187409889579, 4), 0.0003186665179221630
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.051999914333498004037,
    -0.143999914333498004037, -0.0108010022598841, -0.0248010022598841,
    -0.0512949652230244, -0.0072949652230244, 0.0467050347769756,
    0.0267050347769756, 0.0597050347769756, -0.0357508809799500,
    -0.0297508809799500, -0.0517508809799500, -0.0545187409889579,
    -0.0285187409889579, 0.0064812590110421, 0.0474812590110421,
    0.0916813334820778370
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(1, -1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 2)), c(1, -1, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
    eta = exp(1.0192530558433254),
    phi = exp(1.6997766671166766)
  )

  rss_value <- 0.039862989913205782

  fitted_values <- c(
    rep(1, 3), rep(0.999999914333498004037, 2), rep(0.9078010022598841, 2),
    rep(0.5392949652230244, 5), rep(0.2947508809799500, 3),
    rep(0.1715187409889579, 4), 0.0003186665179221630
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.051999914333498004037,
    -0.143999914333498004037, -0.0108010022598841, -0.0248010022598841,
    -0.0512949652230244, -0.0072949652230244, 0.0467050347769756,
    0.0267050347769756, 0.0597050347769756, -0.0357508809799500,
    -0.0297508809799500, -0.0517508809799500, -0.0545187409889579,
    -0.0285187409889579, 0.0064812590110421, 0.0474812590110421,
    0.0916813334820778370
  )

  object <- loggompertz_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(1, -1, 2, 2), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 0.5, 0.5), 10000, c(1, -1, 1, 1), c(1, -1, 10, 10)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
      5982, 3792.6525694635765, -109.81168928761127, 356.28470561183336,
      273978.10050068875,
      # delta
      3792.6525694635765, 3584.7863935594863, -272.35268776124911,
      751.68703757415076, 189464.87516122746,
      # eta
      -109.81168928761127, -272.35268776124911, -144.40272840427726,
      47.656954802768353, -6063.1276169219224,
      # phi
      356.28470561183336, 751.68703757415076, 47.656954802768353,
      15.642222155343835, 18964.339787363599,
      # sigma
      273978.10050068875, 189464.87516122746, -6063.1276169219224,
      18964.339787363599, 9799938.9220864622
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- loggompertz_new(x, y, w, NULL, 10000, NULL, NULL)

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
      y ~ x, mean_function = "loggompertz", lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2), upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = rep(-Inf, 5), upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c( 0, -Inf, -Inf, -Inf), upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
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
      y ~ x, mean_function = "loggompertz", upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = rep(-Inf, 4), upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
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
      y ~ x, mean_function = "loggompertz", start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(1, 1, 1, 1, 1)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(0, 1, 1, -1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz", start = c(0, 1, 1, 0)
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

  result <- drda(y ~ x, mean_function = "loggompertz")

  expect_equal(nauc(result), 0.63548959655000862)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.91370245143275946)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.64363690292771473)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.46392060394389913
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

  result <- drda(y ~ x, mean_function = "loggompertz")

  expect_equal(naac(result), 1.0 - 0.63548959655000862)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.91370245143275946)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.64363690292771473)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.46392060394389913
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

  result <- drda(y ~ x, mean_function = "loggompertz")

  expect_equal(nauc(result), 0.43541376016353408)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.17004252069434817)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.39856406583794775)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.63687322049437206
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

  result <- drda(y ~ x, mean_function = "loggompertz")

  expect_equal(naac(result), 1.0 - 0.43541376016353408)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.17004252069434817)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.39856406583794775)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.63687322049437206
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})
