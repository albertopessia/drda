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

  start <- c(0, 1, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1, 0)
  upper_bound <- c(3, 2, 2, 5, 2)

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic5"))
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

  object <- loglogistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0, -Inf))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5), log(2)))

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

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic5"))
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

  object <- loglogistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0, -Inf))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5), log(2)))
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
    loglogistic5_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 5"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, -1, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    loglogistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 4)),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    loglogistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 5)),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    loglogistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 4)),
    "'upper_bound' must be of length 5"
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, 0, rep(Inf, 2))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, -1, rep(Inf, 2))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, Inf, 0, Inf)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, Inf, -1, Inf)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, Inf, Inf, Inf, 0)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, Inf, Inf, Inf, -1)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2, 1 / 2)

  true_value <- c(
    0.040000000000000000, -0.36000000000000000, -0.75723183391003460,
    -0.82757024793388430, -0.84610059491617090, -0.85284297074649609
  )

  value <- loglogistic5_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic5"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic5_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2, 1 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # delta
      0, 0.44444444444444444, 0.88581314878892734, 0.96396694214876033,
      0.98455621657352323, 0.99204774527388454,
      # log_eta
      0, 0, -0.19503494044495591, -0.10397709116323680,
      -0.057134709679872937, -0.034350007816008794,
      # log_phi
      0, 0.80000000000000000, 0.28137594138001221, 0.094644027047332832,
      0.041213981158891670, 0.021342859858481978,
      # log_nu
      0, -0.057705419819864839, -0.0028717764016332145, -0.00029032455170337206,
      -0.000053524831221610322, -0.000014209639660043913
    ),
    nrow = 6,
    ncol = 5
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
      # (alpha, log_nu)
      rep(0, 6),
      # (delta, alpha)
      rep(0, 6),
      # (delta, delta)
      rep(0, 6),
      # (delta, log_eta)
      0, 0, 0.21670548938328435, 0.11553010129248533, 0.063483010755414374,
      0.038166675351120882,
      # (delta, log_phi)
      0, -0.88888888888888889, -0.31263993486668024, -0.10516003005259204,
      -0.045793312398768522, -0.023714288731646642,
      # (delta, log_nu)
      0, 0.064117133133183154, 0.0031908626684813494, 0.00032258283522596896,
      0.000059472034690678136, 0.000015788488511159903,
      # (log_eta, alpha)
      rep(0, 6),
      # (log_eta, delta)
      0, 0, 0.21670548938328435, 0.11553010129248533, 0.063483010755414374,
      0.038166675351120882,
      # (log_eta, log_eta)
      0, 0, 0.13895874198822747, 0.22002217382468715, 0.17495590096040998,
      0.12952030463176099,
      # (log_eta, log_phi)
      0, 0.80000000000000000, -0.20047508795458475, -0.20027281334293882,
      -0.12620400534492205, -0.080475490002516200,
      # (log_eta, log_nu)
      0, 0, 0.010770091655824237, 0.0018556975716730271, 0.00043945350891208111,
      0.00013630593872216487,
      # (log_phi, alpha)
      rep(0, 6),
      # (log_phi, delta)
      0, -0.88888888888888889, -0.31263993486668024, -0.10516003005259204,
      -0.045793312398768522, -0.023714288731646642,
      # (log_phi, log_eta)
      0, 0.80000000000000000, -0.20047508795458475, -0.20027281334293882,
      -0.12620400534492205, -0.080475490002516200,
      # (log_phi, log_phi)
      0, 0, 0.69516409046826547, 0.26844487671607131, 0.12076654944233373,
      0.063263297747850963,
      # (log_phi, log_nu)
      0, -0.15125582702693699, -0.015537957821777231, -0.0016891287224929564,
      -0.00031699869900435665, -0.000084691641516022768,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, delta)
      0, 0.064117133133183154, 0.0031908626684813494, 0.00032258283522596896,
      0.000059472034690678136, 0.000015788488511159903,
      # (log_nu, log_eta)
      0, 0, 0.010770091655824237, 0.0018556975716730271, 0.00043945350891208111,
      0.00013630593872216487,
      # (log_nu, log_phi)
      0, -0.15125582702693699, -0.015537957821777231, -0.0016891287224929564,
      -0.00031699869900435665, -0.000084691641516022768,
      # (log_nu, log_nu)
      0, -0.039508257760491174, -0.0026557435888732806, -0.00028337276671025439,
      -0.000052974477287249613, -0.000014134357879474393
    ),
    dim = c(6, 5, 5)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic5"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 5)
  expect_length(gradient_hessian$H, 6 * 5 * 5)

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

  theta <- c(4 / 100, -9 / 10, log(3), log(2), -log(2))

  true_value <- 18.892649300006485

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic5"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, log(2), -log(2))
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

  theta <- c(4 / 100, -9 / 10, log(3), log(2), -log(2))

  true_gradient <- c(
    -16.455446444304119, -11.031868204076790, 1.2111196567236504,
    -4.2582569609536717, 0.21580757337191640
  )

  true_hessian <- matrix(
    c(
      # alpha
      15, 9.9144960492267990, -0.97738237458247783,
      3.4649699344260308, -0.18018168496333926,
      # beta
      9.9144960492267990, 8.7869043547305381, -2.2608964521904574,
      6.7832236046641192, -0.32309201468648505,
      # eta
      -0.97738237458247783, -2.2608964521904574, -1.6790404717333001,
      -1.2004702738499023, -0.042876540518495818,
      # phi
      3.4649699344260308, 6.7832236046641192, -1.2004702738499023,
      -1.5670527915537013, 0.45608285384179623,
      # log_nu
      -0.18018168496333926, -0.32309201468648505, -0.042876540518495818,
      0.45608285384179623, 0.16041914038129725
    ),
    nrow = 5,
    ncol = 5
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic5"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 5)
  expect_length(gradient_hessian$H, 5 * 5)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, log(2), -log(2))
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

  theta <- c(0, 1, 1.5969577079412751, 1.7805526373684888, -0.93859189577840407)

  true_value <- c(
    0.91962793079794114, -0.81961827766883428, 1.5969577079412751,
    1.7805526373684888, -0.93859189577840407
  )

  object <- loglogistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.024098464516913181

  fitted_values <- c(
    rep(0.91962793079794114, 3), rep(0.91961835704678826, 2),
    rep(0.8915156044334301, 2), rep(0.553185547862475, 5),
    rep(0.261160239132016, 3), rep(0.15910736421012, 4),
    0.10001037109684
  )

  residuals <- c(
    0.00837206920205886, -0.03162793079794114, 0.06037206920205886,
    0.02838164295321174, -0.06361835704678826, 0.0054843955665699,
    -0.0085156044334301, -0.065185547862475, -0.021185547862475,
    0.032814452137525, 0.012814452137525, 0.045814452137525, -0.002160239132016,
    0.003839760867984, -0.018160239132016, -0.04210736421012, -0.01610736421012,
    0.01889263578988, 0.05989263578988, -0.00801037109684
  )

  object <- loglogistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic5_new(x, y, w, c(0, 1, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.024098464516913181

  fitted_values <- c(
    rep(0.91962793079794114, 3), rep(0.91961835704678826, 2),
    rep(0.8915156044334301, 2), rep(0.553185547862475, 5),
    rep(0.261160239132016, 3), rep(0.15910736421012, 4),
    0.10001037109684
  )

  residuals <- c(
    0.00837206920205886, -0.03162793079794114, 0.06037206920205886,
    0.02838164295321174, -0.06361835704678826, 0.0054843955665699,
    -0.0085156044334301, -0.065185547862475, -0.021185547862475,
    0.032814452137525, 0.012814452137525, 0.045814452137525, -0.002160239132016,
    0.003839760867984, -0.018160239132016, -0.04210736421012, -0.01610736421012,
    0.01889263578988, 0.05989263578988, -0.00801037109684
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1, 0.2), c(1, 3, 5, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 0, 2, 2, 2), 10000,
    c(-1, -3, 1, 1, 0.2), c(1, 3, 5, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, -5, 0.5, 20, 0.1), 10000,
    c(-1, -3, 1, 1, 0.2), c(1, 3, 5, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.066645539419639160

  fitted_values <- c(
    rep(1, 3), rep(0.99402574079604513, 2), rep(0.8720114280687310, 2),
    rep(0.551659663201459, 5), rep(0.282854004623884, 3),
    rep(0.143438974880828, 4), 0.000029593250303
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04602574079604513, -0.13802574079604513,
    0.0249885719312690, 0.0109885719312690, -0.063659663201459,
    -0.019659663201459, 0.034340336798541, 0.014340336798541, 0.047340336798541,
    -0.023854004623884, -0.017854004623884, -0.039854004623884,
    -0.026438974880828, -0.000438974880828, 0.034561025119172,
    0.075561025119172, 0.091970406749697
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(1, -1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.066645539419639160

  fitted_values <- c(
    rep(1, 3), rep(0.99402574079604513, 2), rep(0.8720114280687310, 2),
    rep(0.551659663201459, 5), rep(0.282854004623884, 3),
    rep(0.143438974880828, 4), 0.000029593250303
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04602574079604513, -0.13802574079604513,
    0.0249885719312690, 0.0109885719312690, -0.063659663201459,
    -0.019659663201459, 0.034340336798541, 0.014340336798541, 0.047340336798541,
    -0.023854004623884, -0.017854004623884, -0.039854004623884,
    -0.026438974880828, -0.000438974880828, 0.034561025119172,
    0.075561025119172, 0.091970406749697
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(1, -1, 2, 2, 2), 10000,
    c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10), 10000,
    c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.013857640726518951

  fitted_values <- c(
    rep(0.93735105801003028, 3), rep(0.93691762149598806, 2),
    rep(0.8856728142843277, 2), rep(0.551303815077571, 5),
    rep(0.266468096226581, 3), rep(0.175184022887416, 4),
    0.13295499048021
  )

  residuals <- c(
    -0.00935105801003028, -0.04935105801003028, 0.04264894198996972,
    0.01108237850401194, -0.08091762149598806, 0.0113271857156723,
    -0.0026728142843277, -0.063303815077571, -0.019303815077571,
    0.034696184922429, 0.014696184922429, 0.047696184922429, -0.007468096226581,
    -0.001468096226581, -0.023468096226581, -0.058184022887416,
    -0.032184022887416, 0.002815977112584, 0.043815977112584, -0.04095499048021
  )

  object <- loglogistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic5_new(x, y, w, c(1, -1, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.013857640726518951

  fitted_values <- c(
    rep(0.93735105801003028, 3), rep(0.93691762149598806, 2),
    rep(0.8856728142843277, 2), rep(0.551303815077571, 5),
    rep(0.266468096226581, 3), rep(0.175184022887416, 4),
    0.13295499048021
  )

  residuals <- c(
    -0.00935105801003028, -0.04935105801003028, 0.04264894198996972,
    0.01108237850401194, -0.08091762149598806, 0.0113271857156723,
    -0.0026728142843277, -0.063303815077571, -0.019303815077571,
    0.034696184922429, 0.014696184922429, 0.047696184922429, -0.007468096226581,
    -0.001468096226581, -0.023468096226581, -0.058184022887416,
    -0.032184022887416, 0.002815977112584, 0.043815977112584, -0.04095499048021
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1, 0.5), c(1, 3, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.3, 0.6, 2, 3, 0.7), 10000,
    c(-1, -3, 1, 1, 0.5), c(1, 3, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, -5, 0.5, 20, 0.1), 10000,
    c(-1, -3, 1, 1, 0.5), c(1, 3, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.039067731194279039

  fitted_values <- c(
    rep(1, 3), rep(0.998950082417103934, 2), rep(0.8860150020889389, 2),
    rep(0.546177768205627, 5), rep(0.292942675583851, 3),
    rep(0.162371071685693, 4), 0.00013155549984
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.050950082417103934, -0.142950082417103934,
    0.0109849979110611, -0.0030150020889389, -0.058177768205627,
    -0.014177768205627, 0.039822231794373, 0.019822231794373, 0.052822231794373,
    -0.033942675583851, -0.027942675583851, -0.049942675583851,
    -0.045371071685693, -0.019371071685693, 0.015628928314307,
    0.056628928314307, 0.09186844450016
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(1, -1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 3)), c(1, -1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.039067731194279039

  fitted_values <- c(
    rep(1, 3), rep(0.998950082417103934, 2), rep(0.8860150020889389, 2),
    rep(0.546177768205627, 5), rep(0.292942675583851, 3),
    rep(0.162371071685693, 4), 0.00013155549984
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.050950082417103934, -0.142950082417103934,
    0.0109849979110611, -0.0030150020889389, -0.058177768205627,
    -0.014177768205627, 0.039822231794373, 0.019822231794373, 0.052822231794373,
    -0.033942675583851, -0.027942675583851, -0.049942675583851,
    -0.045371071685693, -0.019371071685693, 0.015628928314307,
    0.056628928314307, 0.09186844450016
  )

  object <- loglogistic5_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(1, -1, 2, 2, 2), 10000,
    c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10), 10000,
    c(1, -1, 1, 1, 0.2), c(1, -1, 10, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  theta <- c(
    alpha = 4 / 100, delta = -9 / 10, eta = 3, phi = 2, nu = 1 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5982, 3824.0971949074954, -106.99802097838297, 314.01641102322962,
      -47.883642843988663, 275110.10701666983,
      # delta
      3824.0971949074954, 3612.3736605916052, -265.25273235105052,
      697.33881993895583, -93.525543818853014, 191586.78848347920,
      # eta
      -106.99802097838297, -265.25273235105052, -138.11661261765622,
      51.132858116483664, -8.0020098371613302, -5896.1169068345876,
      # phi
      314.01641102322962, 697.33881993895583, 51.132858116483664,
      -19.058759215163788, 73.166051563367763, 17235.736993621231,
      # nu
      -47.883642843988663, -93.525543818853014, -8.0020098371613302,
      73.166051563367763, -32.127948336840477, -2541.2806074379009,
      # sigma
      275110.10701666983, 191586.78848347920, -5896.1169068345876,
      17235.736993621231, -2541.2806074379009, 9887707.4179328366
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "sigma"
  )

  object <- loglogistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
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
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c("a", "b", "c", "d", "e")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = matrix(-Inf, nrow = 5, ncol = 2),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
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
      y ~ x, mean_function = "loglogistic5",
      upper_bound = c("a", "b", "c", "d", "e")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 5),
      upper_bound = matrix(Inf, nrow = 5, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
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
      y ~ x, mean_function = "loglogistic5",
      start = c("a", "b", "c", "d", "e")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(-Inf, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(1, 1, 1, 1, 1, 1)
    ),
    "'start' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, -1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 0, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, -1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 0, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, -1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, 0)
    ),
    "parameter 'nu' cannot be negative nor zero"
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

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.63385859186512101)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.91962722237824168)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.64213201569311209)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.46158084161667038
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

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(naac(result), 1.0 - 0.63385859186512101)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.91962722237824168)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.64213201569311209)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.46158084161667038
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

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.44217263241239896)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.14940394882918324)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.40787801655210269)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.62476097122666880
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

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(naac(result), 1.0 - 0.44217263241239896)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.14940394882918324)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.40787801655210269)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.62476097122666880
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})
