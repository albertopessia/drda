context("2-parameter logistic - core functions")

test_that("Constructor", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  stats <- matrix(
    c(
      -6.908, -4.605, -2.303, 0, 2.303, 4.605, 6.908, 3, 2, 2, 5, 3, 4, 1,
      0.932, 0.902, 0.89, 0.5542, 0.2556666667, 0.16425, 0.092, 0.0014186667,
      0.002116, 0.000049, 0.00160656, 0.0000862222, 0.0014676875, 0
    ),
    nrow = 7,
    ncol = 4
  )
  colnames(stats) <- c("x", "n", "m", "v")

  start <- c(-1, 0)

  max_iter <- 10000

  lower_bound <- c(-Inf, -10)
  upper_bound <- c(0, Inf)

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
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

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, start)
  expect_equal(object$lower_bound, lower_bound)
  expect_equal(object$upper_bound, upper_bound)

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0, 0.055, 1.071, 0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269, 0, 1.294, 0.038, 1.278, 0.157
  )

  stats <- matrix(
    c(
      -6.908, -4.605, -2.303, 0.0, 2.303, 4.605, 6.908, 4.549, 0.96, 1.126,
      3.756, 1.797, 2.61, 0.157, 0.9353000659, 0.948, 0.8836838366, 0.55221459,
      0.2606149137, 0.1807233716, 0.092, 0.0014467345, 0, 0.0000091061,
      0.0007707846, 0.0000597738, 0.0014230308, 0
    ),
    nrow = 7,
    ncol = 4
  )
  colnames(stats) <- c("x", "n", "m", "v")

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
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

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, start)
  expect_equal(object$lower_bound, lower_bound)
  expect_equal(object$upper_bound, upper_bound)
})

test_that("Constructor: errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

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
    logistic2_new(x, y, w, c(-1, 0, 1), max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, c(0, 0), max_iter, NULL, NULL),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 2)),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), rep(Inf, 3)),
    "'upper_bound' must be of length 2"
  )
})

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(-2, -3 / 2)

  true_value <- c(
    0.99997991486649750, 0.99799547250877509, 0.83273975003305088,
    0.047425873177566781, 0.00049762293180936533, 4.9786820493880368e-06
  )

  value <- logistic2_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic2"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(-2, -3 / 2)

  true_gradient <- matrix(
    c(
      # eta
      -0.00010861330517063830, -0.0062119220238277483, -0.11178746975982441,
      0.067764989596368199, 0.0018913119136747743, 0.000030395549882905923,
      # phi
      0.000040169460179818566, 0.0040010187215236721, 0.27856851749588564,
      0.090353319461824265, 0.00099475060645420557, 9.9573145242261759e-06
    ),
    nrow = 6,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (eta, eta)
      -0.00058733058023020478, -0.019211744103572242, -0.059706126526206196,
      0.092006042987122827, 0.0071847168058047548, 0.00018556815714372614,
      # (eta, phi)
      0.00019713315420042931, 0.010373526813013228, 0.0095002802309020117,
      0.077498064251918303, 0.0032814838834046505, 0.000055811837184584424,
      # (phi, eta)
      0.00019713315420042931, 0.010373526813013228, 0.0095002802309020117,
      0.077498064251918303, 0.0032814838834046505, 0.000055811837184584424,
      # (phi, phi)
      -0.000080335693123755393, -0.0079699568349665452, -0.37076327551463419,
      0.16356629864377391, 0.0019875211700555996, 0.000019914430751240024
    ),
    dim = c(6, 2, 2)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic2"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 2)
  expect_length(gradient_hessian$H, 6 * 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

context("2-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(-2, -3 / 2)

  true_value <- 0.16210551379791279

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic2"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(NA, -3 / 2)
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(-2)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(-2, -3 / 2)

  true_gradient <- c(-0.015329328113800766, -0.070203605072346703)

  true_hessian <- matrix(
    c(
      # eta
      0.0036156597529564309, -0.064762568286274501,
      # phi
      -0.064762568286274501, 0.16862644463302931
    ),
    nrow = 2,
    ncol = 2
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic2"
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

  known_param <- c(NA, -3 / 2)
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(-2)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 1)
  expect_length(gradient_hessian$H, 1)

  expect_equal(gradient_hessian$G, true_gradient[1])
  expect_equal(gradient_hessian$H, true_hessian[1, 1, drop = FALSE])
})

context("2-parameter logistic - fit")

test_that("fit", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(eta = TRUE, phi = TRUE)

  theta <- c(eta = -0.48361565993858155, phi = 0.55079105348629500)

  rss_value <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552140, 3), rep(0.92367933936267508, 2),
    rep(0.79901316975086140, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936759, 3), rep(0.12339358998664537, 4),
    0.04417373430922940
  )

  residuals <- c(
    -0.045588472981552140, -0.085588472981552140, 0.006411527018447860,
    0.02432066063732492, -0.06767933936267508, 0.09798683024913860,
    0.08398683024913860, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936759, -0.03497945701936759, -0.05697945701936759,
    -0.00639358998664537, 0.01960641001335463, 0.05460641001335463,
    0.09560641001335463, 0.04782626569077060
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic2_new(x, y, w, c(-1, 0), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(eta = TRUE, phi = TRUE)

  theta <- c(eta = -0.48361565993858155, phi = 0.55079105348629500)

  rss_value <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552140, 3), rep(0.92367933936267508, 2),
    rep(0.79901316975086140, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936759, 3), rep(0.12339358998664537, 4),
    0.04417373430922940
  )

  residuals <- c(
    -0.045588472981552140, -0.085588472981552140, 0.006411527018447860,
    0.02432066063732492, -0.06767933936267508, 0.09798683024913860,
    0.08398683024913860, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936759, -0.03497945701936759, -0.05697945701936759,
    -0.00639358998664537, 0.01960641001335463, 0.05460641001335463,
    0.09560641001335463, 0.04782626569077060
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-1, 0), c(0, 1))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(x, y, w, c(-0.5, 0.5), max_iter, c(-1, 0), c(0, 1))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(x, y, w, c(-3, -1), max_iter, c(-1, 0), c(0, 1))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(eta = TRUE, phi = FALSE)

  theta <- c(eta = -0.46193439818028055, phi = 0)

  rss_value <- 0.099093708890039360

  fitted_values <- c(
    rep(0.96049580791205290, 3), rep(0.89351965363258953, 2),
    rep(0.74342272097556632, 2), rep(0.5, 5),
    rep(0.25657727902443368, 3), rep(0.10648034636741047, 4),
    0.03950419208794710
  )

  residuals <- c(
    -0.03249580791205290, -0.07249580791205290, 0.01950419208794710,
    0.05448034636741047, -0.03751965363258953, 0.15357727902443368,
    0.13957727902443368, -0.012, 0.032, 0.086, 0.066, 0.099,
    0.00242272097556632, 0.00842272097556632, -0.01357727902443368,
    0.01051965363258953, 0.03651965363258953, 0.07151965363258953,
    0.11251965363258953, 0.05249580791205290
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-Inf, 0), c(Inf, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  object <- logistic2_new(x, y, w, c(-1, 0), max_iter, c(-Inf, 0), c(Inf, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with different equalities
  object <- logistic2_new(x, y, w, c(-1, 1), max_iter, c(-Inf, 0), c(Inf, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("fit_constrained: equalities and inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  n <- length(y)

  w <- rep(1, n)

  estimated <- c(eta = TRUE, phi = FALSE)

  theta <- c(eta = -0.46193439818028055, phi = 0)

  rss_value <- 0.099093708890039360

  fitted_values <- c(
    rep(0.96049580791205290, 3), rep(0.89351965363258953, 2),
    rep(0.74342272097556632, 2), rep(0.5, 5),
    rep(0.25657727902443368, 3), rep(0.10648034636741047, 4),
    0.03950419208794710
  )

  residuals <- c(
    -0.03249580791205290, -0.07249580791205290, 0.01950419208794710,
    0.05448034636741047, -0.03751965363258953, 0.15357727902443368,
    0.13957727902443368, -0.012, 0.032, 0.086, 0.066, 0.099,
    0.00242272097556632, 0.00842272097556632, -0.01357727902443368,
    0.01051965363258953, 0.03651965363258953, 0.07151965363258953,
    0.11251965363258953, 0.05249580791205290
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-1, 0), c(0, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  object <- logistic2_new(x, y, w, c(-0.5, 0), max_iter, c(-1, 0), c(0, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  object <- logistic2_new(x, y, w, c(-5, -1), max_iter, c(-1, 0), c(0, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("2-parameter logistic - weighted fit")

test_that("fit (weighted)", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 1, 2, 4, 3, 3, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.897, 0.883, 0.488, 0.532, 0.566, 0.599, 0.259,
    0.265, 0.243, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0.055, 1.071, 0.134, 1.825, 1.169, 0.628, 0.327,
    1.201, 0.269, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(eta = TRUE, phi = TRUE)

  theta <- c(eta = -0.45801428680160555, phi = 0.56172839498959744)

  rss_value <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275445, 3), rep(0.91423176703152782, 1),
    rep(0.78786209884225020, 2), rep(0.56396744745889191, 4),
    rep(0.31055470762448179, 3), rep(0.13565276472766797, 3),
    0.05182418247944723
  )

  residuals <- c(
    -0.04036185505275445, -0.08036185505275445, 0.01163814494724555,
    0.03376823296847218, 0.10913790115774980, 0.09513790115774980,
    -0.07596744745889191, -0.03196744745889191, 0.00203255254110809,
    0.03503255254110809, -0.05155470762448179, -0.04555470762448179,
    -0.06755470762448179, 0.00734723527233203, 0.04234723527233203,
    0.08334723527233203, 0.04017581752055277
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic2_new(x, y, w, c(-1, 0), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 1, 2, 4, 3, 3, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.897, 0.883, 0.488, 0.532, 0.566, 0.599, 0.259,
    0.265, 0.243, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0.055, 1.071, 0.134, 1.825, 1.169, 0.628, 0.327,
    1.201, 0.269, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(eta = TRUE, phi = TRUE)

  theta <- c(eta = -0.45801428680160555, phi = 0.56172839498959744)

  rss_value <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275445, 3), rep(0.91423176703152782, 1),
    rep(0.78786209884225020, 2), rep(0.56396744745889191, 4),
    rep(0.31055470762448179, 3), rep(0.13565276472766797, 3),
    0.05182418247944723
  )

  residuals <- c(
    -0.04036185505275445, -0.08036185505275445, 0.01163814494724555,
    0.03376823296847218, 0.10913790115774980, 0.09513790115774980,
    -0.07596744745889191, -0.03196744745889191, 0.00203255254110809,
    0.03503255254110809, -0.05155470762448179, -0.04555470762448179,
    -0.06755470762448179, 0.00734723527233203, 0.04234723527233203,
    0.08334723527233203, 0.04017581752055277
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-1, 0), c(0, 1))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(x, y, w, c(-0.5, 0.5), max_iter, c(-1, 0), c(0, 1))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(x, y, w, c(-3, -1), max_iter, c(-1, 0), c(0, 1))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 1, 2, 4, 3, 3, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.897, 0.883, 0.488, 0.532, 0.566, 0.599, 0.259,
    0.265, 0.243, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0.055, 1.071, 0.134, 1.825, 1.169, 0.628, 0.327,
    1.201, 0.269, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(eta = TRUE, phi = FALSE)

  theta <- c(eta = -0.44146360310995570, phi = 0)

  rss_value <- 0.065691040002891072

  fitted_values <- c(
    rep(0.95476657450173949, 3), rep(0.88421240668361712, 1),
    rep(0.73432748259669405, 2), rep(0.5, 4),
    rep(0.26567251740330595, 3), rep(0.11578759331638288, 3),
    0.04523342549826051
  )

  residuals <- c(
    -0.02676657450173949, -0.06676657450173949, 0.02523342549826051,
    0.06378759331638288, 0.16267251740330595, 0.14867251740330595,
    -0.012, 0.032, 0.066, 0.099, -0.00667251740330595, -0.00067251740330595,
    -0.02267251740330595, 0.02721240668361712, 0.06221240668361712,
    0.10321240668361712, 0.04676657450173949
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-Inf, 0), c(Inf, 0))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic2_new(x, y, w, c(-1, 0), max_iter, c(-Inf, 0), c(Inf, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic2_new(x, y, w, c(-1, 1), max_iter, c(-Inf, 0), c(Inf, 0))

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities and inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 1, 2, 4, 3, 3, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.897, 0.883, 0.488, 0.532, 0.566, 0.599, 0.259,
    0.265, 0.243, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0.055, 1.071, 0.134, 1.825, 1.169, 0.628, 0.327,
    1.201, 0.269, 1.294, 0.038, 1.278, 0.157
  )

  estimated <- c(eta = TRUE, phi = FALSE)

  theta <- c(eta = -0.44146360310995570, phi = 0)

  rss_value <- 0.065691040002891072

  fitted_values <- c(
    rep(0.95476657450173949, 3), rep(0.88421240668361712, 1),
    rep(0.73432748259669405, 2), rep(0.5, 4),
    rep(0.26567251740330595, 3), rep(0.11578759331638288, 3),
    0.04523342549826051
  )

  residuals <- c(
    -0.02676657450173949, -0.06676657450173949, 0.02523342549826051,
    0.06378759331638288, 0.16267251740330595, 0.14867251740330595,
    -0.012, 0.032, 0.066, 0.099, -0.00667251740330595, -0.00067251740330595,
    -0.02267251740330595, 0.02721240668361712, 0.06221240668361712,
    0.10321240668361712, 0.04676657450173949
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, c(-1, 0), c(0, 0))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(x, y, w, c(-0.8, 0), max_iter, c(-1, 0), c(0, 0))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(x, y, w, c(-5, -1), max_iter, c(-1, 0), c(0, 0))
  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

context("2-parameter logistic - general functions")

test_that("fisher_info", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 1, 2, 4, 3, 3, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.897, 0.883, 0.488, 0.532, 0.566, 0.599, 0.259,
    0.265, 0.243, 0.143, 0.178, 0.219, 0.092
  )

  w <- c(
    1.46, 1.385, 1.704, 0.96, 0.055, 1.071, 0.134, 1.825, 1.169, 0.628, 0.327,
    1.201, 0.269, 1.294, 0.038, 1.278, 0.157
  )

  theta <- c(eta = -2, phi = -3 / 2)
  sigma <- 0.05

  true_value <- matrix(c(
      # eta
      -57.679888255261216, -64.222434378395683, 1972.9878566880451,
      # phi
      -64.222434378395683, -68.914476021821838, 3000.1366707555258,
      # sigma
      1972.9878566880451, 3000.1366707555258, 570580.59195114280
    ),
    nrow = 3,
    ncol = 3
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "eta", "phi", "sigma"
  )

  object <- logistic2_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 3 * 3)
  expect_equal(fim, true_value)
})

context("2-parameter logistic - drda fit")

test_that("drda: 'lower_bound' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c("a", "b")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = matrix(-Inf, nrow = 2, ncol = 2),
      upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c( 0, -Inf),
      upper_bound = c(-1,  Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c(Inf, -Inf),
      upper_bound = c(Inf,  Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep( Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      upper_bound = c("a", "b")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = rep(-Inf, 2),
      upper_bound = matrix(Inf, nrow = 2, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c(-Inf, -Inf),
      upper_bound = c(-Inf,  Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep(Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'start' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c("a", "b")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(-1, Inf)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(-Inf, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(-1, 0, 1)
    ),
    "'start' must be of length 2"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(0, 0)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )
})

context("2-parameter logistic - Area under and above the curve")

test_that("nauc: decreasing", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "logistic2")

  expect_equal(nauc(result), 0.52710078010808987)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.50588473625160243)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.52753955267431475)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.50980789375267071
  )
})

test_that("naac: decreasing", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "logistic2")

  expect_equal(naac(result), 1 - 0.52710078010808987)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.50588473625160243)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.52753955267431475)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.50980789375267071
  )
})

test_that("nauc: increasing", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "logistic2")

  expect_equal(nauc(result), 0.51082835161586119)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.57318877340106192)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.51115463147639835)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.62198128900176987
  )
})

test_that("naac: increasing", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "logistic2")

  expect_equal(naac(result), 1 - 0.51082835161586119)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.57318877340106192)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.51115463147639835)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.62198128900176987
  )
})
