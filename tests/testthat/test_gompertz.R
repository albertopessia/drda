context("Gompertz - core functions")

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

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "gompertz"))
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
    gompertz_new(x, y, w, c(0, 1, -1, 0, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, -1, -1, 0), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, 0, -1, 0), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, 1, 0, 0), max_iter, NULL, NULL),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 5)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 5)),
    "'upper_bound' must be of length 4"
  )
})

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_value <- c(
    0.89998272661171928, 0.89827437740755882, 0.74350643245438481,
    0.040000001627273678, 0.04, 0.04
  )

  value <- gompertz_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "gompertz"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "gompertz_fit"
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
      0.000020085335210141429, 0.0020065378981874184, 0.18196926458792463,
      0.99999999810782131, 1.0, 1.0,
      # beta
      0.99997991466478986, 0.99799346210181258, 0.81803073541207537,
      1.8921786948382926e-09, 0.0, 0.0,
      # eta
      -0.000093409318566541349, -0.0053529723590708084, -0.11340771690782248,
      4.9026998302172214e-08, 0.0, 0.0,
      # phi
      0.000034546429617326587, 0.0034477803395290860, 0.28260608849525154,
      6.5369331069562953e-08, 0.0,0.0
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
      0.00010861548670528064, 0.0062243864640358238, 0.13186943826490986,
      -5.7008137560665366e-08, 0.0, 0.0,
      # (alpha, phi)
      -0.000040170266996891380, -0.0040090469064291698, -0.32861173080843203,
      -7.6010850080887154e-08, 0.0, 0.0,
      # (beta, alpha)
      rep(0, 6),
      # (beta, beta)
      rep(0, 6),
      # (beta, eta)
      -0.00010861548670528064, -0.0062243864640358238, -0.13186943826490986,
      5.7008137560665366e-08, 0.0, 0.0,
      # (beta, phi)
      0.000040170266996891380, 0.0040090469064291698, 0.32861173080843203,
      7.6010850080887154e-08, 0.0, 0.0,
      # (eta, alpha)
      0.00010861548670528064, 0.0062243864640358238, 0.13186943826490986,
      -5.7008137560665366e-08, 0.0, 0.0,
      # (eta, beta)
      -0.00010861548670528064, -0.0062243864640358238, -0.13186943826490986,
      5.7008137560665366e-08, 0.0, 0.0,
      # (eta, eta)
      -0.00050512458968195536, -0.016588504216919084, -0.072737619271040149,
      1.4035598794937503e-06, 0.0, 0.0,
      # (eta, phi)
      0.00016954166997178533, 0.0089605510835836899, 0.039955291861489536,
      1.8387285071235523e-06, 0.0, 0.0,
      # (phi, alpha)
      -0.000040170266996891380, -0.0040090469064291698, -0.32861173080843203,
      -7.6010850080887154e-08, 0.0, 0.0,
      # (phi, beta)
      0.000040170266996891380, 0.0040090469064291698, 0.32861173080843203,
      7.6010850080887154e-08, 0.0, 0.0,
      # (phi, eta)
      0.00016954166997178533, 0.0089605510835836899, 0.039955291861489536,
      1.8387285071235523e-06, 0.0, 0.0,
      # (phi, phi)
      -0.000069091471467477887, -0.0068817105751956417, -0.45168627648672276,
      2.4952175635444451e-06, 0.0, 0.0
    ),
    dim = c(6, 4, 4)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1), m = 6),
    class = "gompertz"
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

context("Gompertz - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_value <- 0.14821441202238303

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "gompertz"
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
    -0.62991438860939107, -0.34070142791491030, 0.034085930123255337,
    -0.087075519079779177
  )

  true_hessian <- matrix(
    c(
      # alpha
      7.0662376911651834, 0.30378070014214381, -0.080940272503724425,
      0.20412286562894163,
      # beta
      0.30378070014214381, 7.3262009085505290, -0.16221411023683938,
      0.37153655314632497,
      # eta
      -0.080940272503724425, -0.16221411023683938, 0.044886511003499472,
      -0.074635938611291648,
      # phi
      0.20412286562894163, 0.37153655314632497, -0.074635938611291648,
      0.29863869130618845
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "gompertz"
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

context("Gompertz - support functions")

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

  theta <- c(0, 1, -0.55813290400196051, 0.77264827378630365)

  true_value <- c(
    0.15662111821394763, 0.95551649247615962, -0.55813290400196051,
    0.77264827378630365
  )

  object <- gompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(result, true_value)
})

context("Gompertz - fit")

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
    alpha = 0.15662111821394763,
    beta = 0.95551649247615962,
    eta = -0.55813290400196051,
    phi = 0.77264827378630365
  )

  rss_value <- 0.040219212122823052

  fitted_values <- c(
    rep(0.9446080059051645, 3), rep(0.9167702852554140, 2),
    rep(0.8241333174784071, 2), rep(0.57380509005971934, 5),
    rep(0.23286046858289257, 3), rep(0.15678520063890135, 4),
    0.15662111821398466
  )

  residuals <- c(
    -0.0166080059051645, -0.0566080059051645, 0.0353919940948355,
    0.0312297147445860, -0.0607702852554140, 0.0728666825215929,
    0.0588666825215929, -0.08580509005971934, -0.04180509005971934,
    0.01219490994028066, -0.00780509005971934, 0.02519490994028066,
    0.02613953141710743, 0.03213953141710743, 0.01013953141710743,
    -0.03978520063890135, -0.01378520063890135, 0.02121479936109865,
    0.06221479936109865, -0.06462111821398466
  )

  object <- gompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- gompertz_new(x, y, w, c(0, 1, -1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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
    alpha = 0.15662111821394763,
    beta = 0.95551649247615962,
    eta = -0.55813290400196051,
    phi = 0.77264827378630365
  )

  rss_value <- 0.040219212122823052

  fitted_values <- c(
    rep(0.9446080059051645, 3), rep(0.9167702852554140, 2),
    rep(0.8241333174784071, 2), rep(0.57380509005971934, 5),
    rep(0.23286046858289257, 3), rep(0.15678520063890135, 4),
    0.15662111821398466
  )

  residuals <- c(
    -0.0166080059051645, -0.0566080059051645, 0.0353919940948355,
    0.0312297147445860, -0.0607702852554140, 0.0728666825215929,
    0.0588666825215929, -0.08580509005971934, -0.04180509005971934,
    0.01219490994028066, -0.00780509005971934, 0.02519490994028066,
    0.02613953141710743, 0.03213953141710743, 0.01013953141710743,
    -0.03978520063890135, -0.01378520063890135, 0.02121479936109865,
    0.06221479936109865, -0.06462111821398466
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- gompertz_new(
    x, y, w, c(0, 1, -0.5, 0.5), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- gompertz_new(
    x, y, w, c(-1, 2, 3, -2), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
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
    eta = -0.31497925998576551,
    phi = 1.9292985149742921
  )

  rss_value <- 0.095762067353560977

  fitted_values <- c(
    rep(0.94005425461042752, 3), rep(0.88012827563932673, 2),
    rep(0.76823294694391275, 2), rep(0.58006912883470518, 5),
    rep(0.32467979372885109, 3), rep(0.09799490761511883, 4),
    0.008246676243829067
  )

  residuals <- c(
    -0.01205425461042752, -0.05205425461042752, 0.03994574538957248,
    0.06787172436067327, -0.02412827563932673, 0.12876705305608725,
    0.11476705305608725, -0.09206912883470518, -0.04806912883470518,
    0.00593087116529482, -0.01406912883470518, 0.01893087116529482,
    -0.06567979372885109, -0.05967979372885109, -0.08167979372885109,
    0.01900509238488117, 0.04500509238488117, 0.08000509238488117,
    0.12100509238488117, 0.083753323756170933
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(1, 3, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
    eta = -0.31497925998576551,
    phi = 1.9292985149742921
  )

  rss_value <- 0.095762067353560977

  fitted_values <- c(
    rep(0.94005425461042752, 3), rep(0.88012827563932673, 2),
    rep(0.76823294694391275, 2), rep(0.58006912883470518, 5),
    rep(0.32467979372885109, 3), rep(0.09799490761511883, 4),
    0.008246676243829067
  )

  residuals <- c(
    -0.01205425461042752, -0.05205425461042752, 0.03994574538957248,
    0.06787172436067327, -0.02412827563932673, 0.12876705305608725,
    0.11476705305608725, -0.09206912883470518, -0.04806912883470518,
    0.00593087116529482, -0.01406912883470518, 0.01893087116529482,
    -0.06567979372885109, -0.05967979372885109, -0.08167979372885109,
    0.01900509238488117, 0.04500509238488117, 0.08000509238488117,
    0.12100509238488117, 0.083753323756170933
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(0, 1, -1, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, -0.8, 0.5), 10000, c(0, 1, -1, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(1, 2, -5, -1), 10000, c(0, 1, -1, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

context("Gompertz - weighted fit")

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
    alpha = 0.18884106587251166,
    beta = 0.95399962011994593,
    eta = -0.61538799796323045,
    phi = 0.55914376657765553
  )

  rss_value <- 0.021335475841876475

  fitted_values <- c(
    rep(0.9463101911623085, 3), rep(0.9227699728144719, 1),
    rep(0.8332069394171869, 2), rep(0.56545409142841065, 4),
    rep(0.22992063611644200, 3), rep(0.18884549994968278, 3),
    0.18884106587251166
  )

  residuals <- c(
    -0.0183101911623085, -0.0583101911623085, 0.0336898088376915,
    0.0252300271855281, 0.0637930605828131, 0.0497930605828131,
    -0.07745409142841065, -0.03345409142841065, 0.00054590857158935,
    0.03354590857158935, 0.02907936388355800, 0.03507936388355800,
    0.01307936388355800, -0.04584549994968278, -0.01084549994968278,
    0.03015450005031722, -0.09684106587251166
  )

  object <- gompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- gompertz_new(x, y, w, c(0, 1, -1, 0), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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
    alpha = 0.18884106587251166,
    beta = 0.95399962011994593,
    eta = -0.61538799796323045,
    phi = 0.55914376657765553
  )

  rss_value <- 0.021335475841876475

  fitted_values <- c(
    rep(0.9463101911623085, 3), rep(0.9227699728144719, 1),
    rep(0.8332069394171869, 2), rep(0.56545409142841065, 4),
    rep(0.22992063611644200, 3), rep(0.18884549994968278, 3),
    0.18884106587251166
  )

  residuals <- c(
    -0.0183101911623085, -0.0583101911623085, 0.0336898088376915,
    0.0252300271855281, 0.0637930605828131, 0.0497930605828131,
    -0.07745409142841065, -0.03345409142841065, 0.00054590857158935,
    0.03354590857158935, 0.02907936388355800, 0.03507936388355800,
    0.01307936388355800, -0.04584549994968278, -0.01084549994968278,
    0.03015450005031722, -0.09684106587251166
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- gompertz_new(
    x, y, w, c(0, 1, -0.5, 0.5), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- gompertz_new(
    x, y, w, c(-1, 2, 3, -2), 10000, c(-0.5, 0.9, -2, -1), c(0.5, 1.5, 0, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 4)
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
    eta = -0.30463172419202541,
    phi = 1.9984790796325234
  )

  rss_value <- 0.061582287492551908

  fitted_values <- c(
    rep(0.93582714909417189, 3), rep(0.87479075734980850, 1),
    rep(0.76359246595046023, 2), rep(0.58042005084694914, 4),
    rep(0.33380250491484725, 3), rep(0.10944953162467914, 3),
    0.01153932721903692
  )

  residuals <- c(
    -0.00782714909417189, -0.04782714909417189, 0.04417285090582811,
    0.07320924265019150, 0.13340753404953977, 0.11940753404953977,
    -0.09242005084694914, -0.04842005084694914, -0.01442005084694914,
    0.01857994915305086, -0.07480250491484725, -0.06880250491484725,
    -0.09080250491484725, 0.03355046837532086, 0.06855046837532086,
    0.10955046837532086, 0.08046067278096308
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(1, 3, -1, 0), 10000, c(0, 1, -Inf, -Inf), c(0, 1, Inf,Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
    eta = -0.30463172419202541,
    phi = 1.9984790796325234
  )

  rss_value <- 0.061582287492551908

  fitted_values <- c(
    rep(0.93582714909417189, 3), rep(0.87479075734980850, 1),
    rep(0.76359246595046023, 2), rep(0.58042005084694914, 4),
    rep(0.33380250491484725, 3), rep(0.10944953162467914, 3),
    0.01153932721903692
  )

  residuals <- c(
    -0.00782714909417189, -0.04782714909417189, 0.04417285090582811,
    0.07320924265019150, 0.13340753404953977, 0.11940753404953977,
    -0.09242005084694914, -0.04842005084694914, -0.01442005084694914,
    0.01857994915305086, -0.07480250491484725, -0.06880250491484725,
    -0.09080250491484725, 0.03355046837532086, 0.06855046837532086,
    0.10955046837532086, 0.08046067278096308
  )

  object <- gompertz_new(
    x, y, w, NULL, 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- gompertz_new(
    x, y, w, c(0, 1, -0.8, 0.5), 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- gompertz_new(
    x, y, w, c(1, 2, -5, -1), 10000, c(0, 1, -2, 0), c(0, 1, 0, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

context("Gompertz - general functions")

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
      3342.8932199135800, 67.811647332874071, -17.733850799910351,
      43.927487173638108, 43592.808103392150,
      # beta
      67.811647332874071, 2503.4834854206719, -35.563453987617197,
      84.661231390286641, 5397.2473298858682,
      # eta
      -17.733850799910351, -35.563453987617197, 10.741431979171633,
      -17.141181890023916, -290.45079925888223,
      # phi
      43.927487173638108, 84.661231390286641, -17.141181890023916,
      64.543412202089920, 715.35910155304622,
      # sigma
      43592.808103392150, 5397.2473298858682, -290.45079925888223,
      715.35910155304622, 554058.41271733448
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "sigma"
  )

  object <- gompertz_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

context("Gompertz - drda fit")

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
      y ~ x, mean_function = "gompertz",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
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
      y ~ x, mean_function = "gompertz",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
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
      y ~ x, mean_function = "gompertz",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(0, Inf, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(-Inf, 0, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
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
      y ~ x, mean_function = "gompertz",
      start = c(0, 1, 0, 0)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )
})

context("Gompertz - Area under and above the curve")

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

  result <- drda(y ~ x, mean_function = "gompertz")

  expect_equal(nauc(result), 0.54579664519347487)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.49159337345681468)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.52113889669997744)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.48598895576135779
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

  result <- drda(y ~ x, mean_function = "gompertz")

  expect_equal(naac(result), 1 - 0.54579664519347487)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.49159337345681468)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.52113889669997744)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.48598895576135779
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

  result <- drda(y ~ x, mean_function = "gompertz")

  expect_equal(nauc(result), 0.53560914488319814)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.56876475421652998)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50298998768397183)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.61460792369421663
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

  result <- drda(y ~ x, mean_function = "gompertz")

  expect_equal(naac(result), 1 - 0.53560914488319814)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.56876475421652998)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50298998768397183)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.61460792369421663
  )
})
