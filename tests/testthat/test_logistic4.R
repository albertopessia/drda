context("4-parameter logistic - core functions")

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

  max_iter <- 10000

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

  start <- c(0, 1, -1, 0)

  lower_bound <- c(0, -1, -Inf, -10)
  upper_bound <- c(3, 2, 0, 5)

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic4"))
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

  object <- logistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, -1, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10))
  expect_equal(object$upper_bound, c(2, 2, 0, 5))

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

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic4"))
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

  object <- logistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, -1, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10))
  expect_equal(object$upper_bound, c(2, 2, 0, 5))
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
    logistic4_new(x, y, w, c(0, 1, -1, 0, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, -1, -1, 0), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, 0, -1, 0), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, 1, 0, 0), max_iter, NULL, NULL),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 5)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 5)),
    "'upper_bound' must be of length 4"
  )
})

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_value <- c(
    0.89998272678518785, 0.89827610635754658, 0.75615618502842375,
    0.080786250932707432, 0.040427955721356054, 0.040004281666562474
  )

  value <- logistic4_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic4"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic4_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      0.000020085133502497096, 0.0020045274912249125, 0.16726024996694912,
      0.95257412682243322, 0.99950237706819063, 0.99999502131795061,
      # beta
      0.99997991486649750, 0.99799547250877509, 0.83273975003305088,
      0.047425873177566781, 0.00049762293180936533, 4.9786820493880368e-06,
      # eta
      -0.000093407442446748936, -0.0053422529404918635, -0.096137223993448991,
      0.058277891052876651, 0.0016265282457603059, 0.000026140172899299094,
      # phi
      0.000034545735754643967, 0.0034408761005103580, 0.23956892504646165,
      0.077703854737168868, 0.00085548552155061679, 8.5632904908345112e-06
    ),
    nrow = 6,
    ncol = 4
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, 6),
      # (alpha, beta)
      rep(0, 6),
      # (alpha, eta)
      0.00010861330517063830, 0.0062119220238277483, 0.11178746975982441,
      -0.067764989596368199, -0.0018913119136747743, -0.000030395549882905923,
      # (alpha, phi)
      -0.000040169460179818566, -0.0040010187215236721, -0.27856851749588564,
      -0.090353319461824265, -0.00099475060645420557, -9.9573145242261759e-06,
      # (beta, alpha)
      rep(0, 6),
      # (beta, beta)
      rep(0, 6),
      # (beta, eta)
      -0.00010861330517063830, -0.0062119220238277483, -0.11178746975982441,
      0.067764989596368199, 0.0018913119136747743, 0.000030395549882905923,
      # (beta, phi)
      0.000040169460179818566, 0.0040010187215236721, 0.27856851749588564,
      0.090353319461824265, 0.00099475060645420557, 9.9573145242261759e-06,
      # (eta, alpha)
      0.00010861330517063830, 0.0062119220238277483, 0.11178746975982441,
      -0.067764989596368199, -0.0018913119136747743, -0.000030395549882905923,
      # (eta, beta)
      -0.00010861330517063830, -0.0062119220238277483, -0.11178746975982441,
      0.067764989596368199, 0.0018913119136747743, 0.000030395549882905923,
      # (eta, eta)
      -0.00050510429899797611, -0.016522099929072128, -0.051347268812537329,
      0.079125196968925631, 0.0061788564529920891, 0.00015958861514360448,
      # (eta, phi)
      0.00016953451261236921, 0.0089212330591913757, 0.0081702409985757300,
      0.066648335256649741, 0.0028220761397279995, 0.000047998179978742605,
      # (phi, alpha)
      -0.000040169460179818566, -0.0040010187215236721, -0.27856851749588564,
      -0.090353319461824265, -0.00099475060645420557, -9.9573145242261759e-06,
      # (phi, beta)
      0.000040169460179818566, 0.0040010187215236721, 0.27856851749588564,
      0.090353319461824265, 0.00099475060645420557, 9.9573145242261759e-06,
      # (phi, eta)
      0.00016953451261236921, 0.0089212330591913757, 0.0081702409985757300,
      0.066648335256649741, 0.0028220761397279995, 0.000047998179978742605,
      # (phi, phi)
      -0.000069088696086429638, -0.0068541628780712289, -0.31885641694258540,
      0.14066701683364557, 0.0017092682062478156, 0.000017126410446066421
    ),
    dim = c(6, 4, 4)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic4"
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

context("4-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_value <- 0.11303184522146127

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic4"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, -3 / 2)
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(9 / 10, -2))

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

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_gradient <- c(
    -0.44448183274141616, -0.33640042687863516, 0.011139573467794527,
    -0.087637515178570666
  )

  true_hessian <- matrix(
    c(
      # alpha
      6.6825689117006145, 0.46682906460177071, 0.18178820493109173,
      0.48070540381088593,
      # beta
      0.46682906460177071, 7.3837729590958440, -0.15237848511801803,
      0.32224058730415969,
      # eta
      0.18178820493109173, -0.15237848511801803, 0.022131256419080663,
      -0.045877037358683589,
      # phi
      0.48070540381088593, 0.32224058730415969, -0.045877037358683589,
      0.19227987353030561
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic4"
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

  known_param <- c(4 / 100, NA, NA, -3 / 2)
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(c(9 / 10, -2))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[2:3])
  expect_equal(gradient_hessian$H, true_hessian[2:3, 2:3])
})

context("4-parameter logistic - support functions")

test_that("mle_asy", {
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

  theta <- c(0, 1, -0.89523277708641287, 0.13714073896752251)

  true_value <- c(
    0.14236056369991001, 0.93473303092200540, -0.89523277708641287,
    0.13714073896752251
  )

  object <- logistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(result, true_value)
})

context("4-parameter logistic - fit")

test_that("fit", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.14236056369991001,
    beta = 0.93473303092200540,
    eta = -0.89523277708641287,
    phi = 0.13714073896752251
  )

  rss_value <- 0.030080292750857811

  fitted_values <- c(
    rep(0.9332908321140832, 3), rep(0.9235378622597109, 2),
    rep(0.8545832854276004, 2), rep(0.56283675776095411, 5),
    rep(0.24201206136017809, 3), rep(0.15661550382373066, 4),
    0.14420322010290027
  )

  residuals <- c(
    -0.0052908321140832, -0.0452908321140832, 0.0467091678859168,
    0.0244621377402891, -0.0675378622597109, 0.0424167145723996,
    0.0284167145723996, -0.07483675776095411, -0.03083675776095411,
    0.02316324223904589, 0.00316324223904589, 0.03616324223904589,
    0.01698793863982191, 0.02298793863982191, 0.00098793863982191,
    -0.03961550382373066, -0.01361550382373066, 0.02138449617626934,
    0.06238449617626934, -0.05220322010290027
  )

  object <- logistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic4_new(x, y, w, c(0, 1, -1, 0), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.14236056369991001,
    beta = 0.93473303092200540,
    eta = -0.89523277708641287,
    phi = 0.13714073896752251
  )

  rss_value <- 0.030080292750857811

  fitted_values <- c(
    rep(0.9332908321140832, 3), rep(0.9235378622597109, 2),
    rep(0.8545832854276004, 2), rep(0.56283675776095411, 5),
    rep(0.24201206136017809, 3), rep(0.15661550382373066, 4),
    0.14420322010290027
  )

  residuals <- c(
    -0.0052908321140832, -0.0452908321140832, 0.0467091678859168,
    0.0244621377402891, -0.0675378622597109, 0.0424167145723996,
    0.0284167145723996, -0.07483675776095411, -0.03083675776095411,
    0.02316324223904589, 0.00316324223904589, 0.03616324223904589,
    0.01698793863982191, 0.02298793863982191, 0.00098793863982191,
    -0.03961550382373066, -0.01361550382373066, 0.02138449617626934,
    0.06238449617626934, -0.05220322010290027
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic4_new(
    x, y, w, c(0, 1, -0.5, 0.5), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic4_new(
    x, y, w, c(-1, 2, 3, -2), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.48361565993858148,
    phi = 0.55079105348629506
  )

  rss_value <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552128, 3), rep(0.92367933936267506, 2),
    rep(0.79901316975086137, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936762, 3), rep(0.12339358998664541, 4),
    0.04417373430922942
  )

  residuals <- c(
    -0.045588472981552128, -0.085588472981552128, 0.006411527018447872,
    0.02432066063732494, -0.06767933936267506, 0.09798683024913863,
    0.08398683024913863, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936762, -0.03497945701936762, -0.05697945701936762,
    -0.00639358998664541, 0.01960641001335459, 0.05460641001335459,
    0.09560641001335459, 0.04782626569077058
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic4_new(
    x, y, w, c(0, 1, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic4_new(
    x, y, w, c(1, 3, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fit_constrained: equalities and inequalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.48361565993858148,
    phi = 0.55079105348629506
  )

  rss_value <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552128, 3), rep(0.92367933936267506, 2),
    rep(0.79901316975086137, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936762, 3), rep(0.12339358998664541, 4),
    0.04417373430922942
  )

  residuals <- c(
    -0.045588472981552128, -0.085588472981552128, 0.006411527018447872,
    0.02432066063732494, -0.06767933936267506, 0.09798683024913863,
    0.08398683024913863, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936762, -0.03497945701936762, -0.05697945701936762,
    -0.00639358998664541, 0.01960641001335459, 0.05460641001335459,
    0.09560641001335459, 0.04782626569077058
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(0, 1, -1, 0), c(0, 1, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0, 1, -0.8, 0.5), 10000, c(0, 1, -1, 0), c(0, 1, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(1, 2, -5, -1), 10000, c(0, 1, -1, 0), c(0, 1, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

context("4-parameter logistic - weighted fit")

test_that("fit (weighted)", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.17286261579329597,
    beta = 0.94203946208670380,
    eta = -0.96536781199483224,
    phi = -0.0061161452508499899
  )

  rss_value <- 0.015426894657174924

  fitted_values <- c(
    rep(0.9410580985201252, 3), rep(0.9330686602484605, 1),
    rep(0.8665029710263503, 2), rep(0.55631567380252894, 4),
    rep(0.24759852669347845, 3), rep(0.18172932146179907, 3),
    0.17383247345211994
  )

  residuals <- c(
    -0.0130580985201252, -0.0530580985201252, 0.0389419014798748,
    0.0149313397515395, 0.0304970289736497, 0.0164970289736497,
    -0.06831567380252894, -0.02431567380252894, 0.00968432619747106,
    0.04268432619747106, 0.01140147330652155, 0.01740147330652155,
    -0.00459852669347845, -0.03872932146179907, -0.00372932146179907,
    0.03727067853820093, -0.08183247345211994
  )

  object <- logistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic4_new(x, y, w, c(0, 1, -1, 0), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): inequalities", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.17286261579329597,
    beta = 0.94203946208670380,
    eta = -0.96536781199483224,
    phi = -0.0061161452508499899
  )

  rss_value <- 0.015426894657174924

  fitted_values <- c(
    rep(0.9410580985201252, 3), rep(0.9330686602484605, 1),
    rep(0.8665029710263503, 2), rep(0.55631567380252894, 4),
    rep(0.24759852669347845, 3), rep(0.18172932146179907, 3),
    0.17383247345211994
  )

  residuals <- c(
    -0.0130580985201252, -0.0530580985201252, 0.0389419014798748,
    0.0149313397515395, 0.0304970289736497, 0.0164970289736497,
    -0.06831567380252894, -0.02431567380252894, 0.00968432619747106,
    0.04268432619747106, 0.01140147330652155, 0.01740147330652155,
    -0.00459852669347845, -0.03872932146179907, -0.00372932146179907,
    0.03727067853820093, -0.08183247345211994
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic4_new(
    x, y, w, c(0, 1, -0.5, 0.5), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic4_new(
    x, y, w, c(-1, 2, 3, -2), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.45801428680160602,
    phi = 0.56172839498959692
  )

  rss_value <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275455, 3), rep(0.91423176703152800, 1),
    rep(0.78786209884225038, 2), rep(0.56396744745889192, 4),
    rep(0.31055470762448156, 3), rep(0.13565276472766772, 3),
    0.05182418247944707
  )

  residuals <- c(
    -0.04036185505275455, -0.08036185505275455, 0.01163814494724545,
    0.03376823296847200, 0.10913790115774962, 0.09513790115774962,
    -0.07596744745889192, -0.03196744745889192, 0.00203255254110808,
    0.03503255254110808, -0.05155470762448156, -0.04555470762448156,
    -0.06755470762448156, 0.00734723527233228, 0.04234723527233228,
    0.08334723527233228, 0.04017581752055293
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic4_new(
    x, y, w, c(0, 1, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic4_new(
    x, y, w, c(1, 3, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf,Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fit_constrained (weighted): equalities and inequalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.45801428680160602,
    phi = 0.56172839498959692
  )

  rss_value <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275455, 3), rep(0.91423176703152800, 1),
    rep(0.78786209884225038, 2), rep(0.56396744745889192, 4),
    rep(0.31055470762448156, 3), rep(0.13565276472766772, 3),
    0.05182418247944707
  )

  residuals <- c(
    -0.04036185505275455, -0.08036185505275455, 0.01163814494724545,
    0.03376823296847200, 0.10913790115774962, 0.09513790115774962,
    -0.07596744745889192, -0.03196744745889192, 0.00203255254110808,
    0.03503255254110808, -0.05155470762448156, -0.04555470762448156,
    -0.06755470762448156, 0.00734723527233228, 0.04234723527233228,
    0.08334723527233228, 0.04017581752055293
  )

  object <- logistic4_new(
    x, y, w, NULL, 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0, 1, -0.8, 0.5), 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(1, 2, -5, -1), 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

context("4-parameter logistic - general functions")

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

  theta <- c(
    alpha = 4 / 100,
    beta = 9 / 10,
    eta = -2,
    phi = -3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      3200.7334651687024, 131.73980832572013, 119.11651396790572,
      210.06856397497571, 39707.343906259385,
      # beta
      131.73980832572013, 2517.7869181798573, -75.884747727720387,
      16.523906345262924, 6591.5447019487505,
      # eta
      119.11651396790572, -75.884747727720387, -44.467697395836343,
      -51.884823618362492, 1436.4916158378913,
      # phi
      210.06856397497571, 16.523906345262924, -51.884823618362492,
      -46.573859172100720, 2759.2902437346177,
      # sigma
      39707.343906259385, 6591.5447019487505, 1436.4916158378913,
      2759.2902437346177, 479738.08567645750
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "sigma"
  )

  object <- logistic4_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

context("4-parameter logistic - drda fit")

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
      y ~ x, mean_function = "logistic4",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
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
      y ~ x, mean_function = "logistic4",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
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
      y ~ x, mean_function = "logistic4",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, Inf, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(-Inf, 0, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 0, -1, 0, 1)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, -1, -1, 0, 1, 1)
    ),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 1, 0, 0)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )
})

context("4-parameter logistic - Area under and above the curve")

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

  result <- drda(y ~ x, mean_function = "logistic4")

  expect_equal(nauc(result), 0.54397871471228701)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.48271527244189229)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.52125978065641581)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.47119212073648715
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

  result <- drda(y ~ x, mean_function = "logistic4")

  expect_equal(naac(result), 1 - 0.54397871471228701)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.48271527244189229)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.52125978065641581)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.47119212073648715
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

  result <- drda(y ~ x, mean_function = "logistic4")

  expect_equal(nauc(result), 0.53061051681313758)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.56420053421792491)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50464326593339348)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.60700089036320819
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

  result <- drda(y ~ x, mean_function = "logistic4")

  expect_equal(naac(result), 1 - 0.53061051681313758)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.56420053421792491)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50464326593339348)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.60700089036320819
  )
})
