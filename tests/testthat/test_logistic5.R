context("5-parameter logistic - core functions")

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

  start <- c(0, 1, -1, 0, 1)

  max_iter <- 10000

  lower_bound <- c(-Inf, 1, -Inf, -10, 1)
  upper_bound <- c(1, Inf, 0, Inf, Inf)

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic5"))
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

  object <- logistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 0, -1, 0, 0))
  expect_equal(object$lower_bound, c(-Inf, -Inf, -Inf, -10, 0))
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

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic5"))
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

  object <- logistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 0, -1, 0, 0))
  expect_equal(object$lower_bound, c(-Inf, -Inf, -Inf, -10, 0))
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
    logistic5_new(x, y, w, c(0, 1, -1, 0), max_iter, NULL, NULL),
    "'start' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, c(1, 0, -1, 0, 1), max_iter, NULL, NULL),
    "parameter 'beta' is smaller than 'alpha'"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 0, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, -1, 0, 0), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, -1, 0, -1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 4)),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 5)),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, c(0, -1, rep(-Inf, 3)), rep(Inf, 5)),
    "'lower_bound[2]' cannot be smaller than 'lower_bound[1]'",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 4)),
    "'upper_bound' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 0, rep(Inf, 3))),
    "'upper_bound[2]' cannot be smaller than 'upper_bound[1]'",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(
      x, y, w, NULL, max_iter, c(-0.5, 0, -1, 0, 0), c(0.5, 1, 0, 1, 1)
    ),
    "'lower_bound[2]' cannot be smaller than 'upper_bound[1]'",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2)

  true_value <- c(
    0.89998272669845415, 0.89827524246036657, 0.75019144346736942,
    0.047052490645792413, 0.040000850995162837, 0.040000000085267377
  )

  value <- logistic5_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic5"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # log_omega
      0.85998272669845415, 0.85827524246036657, 0.71019144346736942,
      0.0070524906457924134, 8.5099516283666758e-07, 8.5267376783781469e-11,
      # eta
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # phi
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # log_nu
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09
    ),
    nrow = 6,
    ncol = 5
  )

  true_hessian <- array(
    c(
      # (alpha, ...)
      rep(0, 30),
      # (log_omega, alpha)
      rep(0, 6),
      # (log_omega, log_omega)
      0.85998272669845415, 0.85827524246036657, 0.71019144346736942,
      0.0070524906457924134, 8.5099516283666758e-07, 8.5267376783781469e-11,
      # (log_omega, eta)
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # (log_omega, phi)
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # (log_omega, log_nu)
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09,
      # (eta, alpha)
      rep(0, 6),
      # (eta, log_omega)
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # (eta, eta)
      -0.00050511444418713689, -0.016555252126485643, -0.060637799043344759,
      0.049883501712141417, 0.000049098048382061985, 1.2712402410105175e-08,
      # (eta, phi)
      0.00016953809123729721, 0.0089408616320290739, 0.021478649997785072,
      0.053683659136403524, 0.000024123213411596926, 3.9939380477898473e-09,
      # (eta, log_nu)
      9.3805989611171948e-10, 5.3597039009804577e-06, 0.0085715606355158565,
      0.039930425361417561, 0.000070015304961175913, 2.0858519166902048e-08,
      # (phi, alpha)
      rep(0, 6),
      # (phi, log_omega)
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # (phi, eta)
      0.00016953809123729721, 0.0089408616320290739, 0.021478649997785072,
      0.053683659136403524, 0.000024123213411596926, 3.9939380477898473e-09,
      # (phi, phi)
      -0.000069090083756049672, -0.0068679160064153061, -0.37654877818008748,
      0.088681780821584741, 0.000013582081688859556, 1.3642440673798294e-09,
      # (phi, log_nu)
      -3.4693134127485286e-10, -3.4521160387059138e-06, -0.021359879993633144,
      0.053240567148556747, 0.000036825108839864453, 6.8330672303858801e-09,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, log_omega)
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09,
      # (log_nu, eta)
      9.3805989611171948e-10, 5.3597039009804577e-06, 0.0085715606355158565,
      0.039930425361417561, 0.000070015304961175913, 2.0858519166902048e-08,
      # (log_nu, phi)
      -3.4693134127485286e-10, -3.4521160387059138e-06, -0.021359879993633144,
      0.053240567148556747, 0.000036825108839864453, 6.8330672303858801e-09,
      # (log_nu, log_nu)
      8.6733125674846536e-11, 8.6331893228632990e-07, 0.0055845141256710526,
      0.053441912635589904, 0.00011068910773858665, 3.6103283407515775e-08
    ),
    dim = c(6, 5, 5)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic5"
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

context("5-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2, -log(2), log(3 / 2))

  true_value <- 0.13844046588658472

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic5"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2))
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(log(43 / 50), -2))

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

  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2, -log(2), log(3 / 2))

  true_gradient <- c(
    -0.92903069002014085, -0.28833791237257394, 0.022267352061200450,
    -0.086374079772717690, -0.010097206230551137
  )

  true_hessian <- matrix(
    c(
      # alpha
      15.000000000000000, 6.6033693099798592, -0.14741189907774185,
      0.63157849846052230, 0.096833141608282843,
      # log_omega
      6.6033693099798592, 5.1492248569956934, -0.13897258392352741,
      0.29154887863609586, -0.00055060290012352080,
      # eta
      -0.14741189907774185, -0.13897258392352741, 0.018237239090142242,
      -0.077452292586773564, -0.017846532854813774,
      # phi
      0.63157849846052230, 0.29154887863609586, -0.077452292586773564,
      0.21294298805991181, -0.0090213900412146364,
      # log_nu
      0.096833141608282843, -0.00055060290012352080, -0.017846532854813774,
      -0.0090213900412146364, -0.020700058407993170
    ),
    nrow = 5,
    ncol = 5
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic5"
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

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2))
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(c(log(43 / 50), -2))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[2:3])
  expect_equal(gradient_hessian$H, true_hessian[2:3, 2:3])
})

context("5-parameter logistic - general functions")

test_that("fisher_info", {
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

  theta <- c(
    alpha = 4 / 100,
    beta = 9 / 10,
    eta = -2,
    phi = -3 / 2,
    nu = 1 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5191.5736239671087, 132.90987565500003, 23.679974158167395,
      86.980812626015735, 85.283559876424630, 0,
      # beta
      132.90987565500003, 2542.6066247228912, -72.803115749724175,
      174.41094047950246, 9.0075218161025453, 0,
      # eta
      23.679974158167395, -72.803115749724175, 9.4208692315701512,
      -20.590332645351424, 0.57269249241964757, 0,
      # phi
      86.980812626015735, 174.41094047950246, -20.590332645351424,
      55.031696576965175, 4.7687936450855367, 0,
      # nu
      85.283559876424630, 9.0075218161025453, 0.57269249241964757,
      4.7687936450855367, 3.6713059622893367, 0,
      # sigma
      rep(0, 5), 6800
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "nu", "sigma"
  )

  object <- structure(
    list(stats = suff_stats(x, y, w), n = n, m = 7),
    class = "logistic5"
  )

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
  expect_equal(fim, true_value)
})

context("5-parameter logistic - fit")

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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.093212121358460132,
    beta = 0.093212121358460132 + exp(-0.18985173265400376),
    eta = -1.7617462932355769,
    phi = -0.47836972294568200,
    nu = exp(1.3765334489390746)
  )

  rss_value <- 0.024883087882351184

  fitted_values <- c(
    rep(0.92028391866853350, 3), rep(0.91971916938223621, 2),
    rep(0.89002742082078101, 2), rep(0.55337929341255560, 5),
    rep(0.26271803120343567, 3), rep(0.15412982230606348, 4),
    0.11508521369102614
  )

  residuals <- c(
    0.00771608133146650, -0.03228391866853350, 0.05971608133146650,
    0.02828083061776379, -0.06371916938223621, 0.00697257917921899,
    -0.00702742082078101, -0.06537929341255560, -0.02137929341255560,
    0.03262070658744440, 0.01262070658744440, 0.04562070658744440,
    -0.00371803120343567, 0.00228196879656433, -0.01971803120343567,
    -0.03712982230606348, -0.01112982230606348, 0.02387017769393652,
    0.06487017769393652, -0.02308521369102614
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic5_new(x, y, w, c(0, 1, -1, 0, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.093212121358461020,
    beta = 0.093212121358461020 + exp(-0.18985173265400473),
    eta = -1.7617462932355648,
    phi = -0.47836972294567659,
    nu = exp(1.3765334489390619)
  )

  rss_value <- 0.024883087882351184

  fitted_values <- c(
    rep(0.92028391866853359, 3), rep(0.91971916938223627, 2),
    rep(0.89002742082078076, 2), rep(0.55337929341255556, 5),
    rep(0.26271803120343537, 3), rep(0.15412982230606358, 4),
    0.11508521369102661
  )

  residuals <- c(
    0.00771608133146641, -0.03228391866853359, 0.05971608133146641,
    0.02828083061776373, -0.06371916938223627, 0.00697257917921924,
    -0.00702742082078076, -0.06537929341255556, -0.02137929341255556,
    0.03262070658744444, 0.01262070658744444, 0.04562070658744444,
    -0.00371803120343537, 0.00228196879656463, -0.01971803120343537,
    -0.03712982230606358, -0.01112982230606358, 0.02387017769393642,
    0.06487017769393642, -0.02308521369102661
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(-1, 0.5, -2, -1, 0), c(0.2, 1, -1, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 0.75, -1.5, -0.5, 1), max_iter, c(-1, 0.5, -2, -1, 0),
    c(0.2, 1, -1, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(-2, 2, 10, -10, 7), max_iter, c(-1, 0.5, -2, -1, 0),
    c(0.2, 1, -1, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605569,
    phi = -0.30897155961772666,
    nu = exp(0.85096464960136723)
  )

  rss_value <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706217, 3), rep(0.96390948874813836, 2),
    rep(0.83729075842125315, 2), rep(0.55558345438725122, 5),
    rep(0.29108572543250873, 3), rep(0.14085756611083844, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706217, -0.105387436096706217, -0.013387436096706217,
    -0.01590948874813836, -0.10790948874813836, 0.05970924157874685,
    0.04570924157874685, -0.06758345438725122, -0.02358345438725122,
    0.03041654561274878, 0.01041654561274878, 0.04341654561274878,
    -0.03208572543250873, -0.02608572543250873, -0.04808572543250873,
    -0.02385756611083844, 0.00214243388916156, 0.03714243388916156,
    0.07814243388916156, 0.02497284742870274
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(0, 1, -Inf, -Inf, -Inf), c(0, 1, Inf, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic5_new(
    x, y, w, c(0, 1, -1, 0, 1), max_iter, c(0, 1, -Inf, -Inf, -Inf),
    c(0, 1, Inf, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic5_new(
    x, y, w, c(1, 3, -1, 0, 1), max_iter, c(0, 1, -Inf, -Inf, -Inf),
    c(0, 1, Inf, Inf, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605569,
    phi = -0.30897155961772666,
    nu = exp(0.85096464960136723)
  )

  rss_value <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706217, 3), rep(0.96390948874813836, 2),
    rep(0.83729075842125315, 2), rep(0.55558345438725122, 5),
    rep(0.29108572543250873, 3), rep(0.14085756611083844, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706217, -0.105387436096706217, -0.013387436096706217,
    -0.01590948874813836, -0.10790948874813836, 0.05970924157874685,
    0.04570924157874685, -0.06758345438725122, -0.02358345438725122,
    0.03041654561274878, 0.01041654561274878, 0.04341654561274878,
    -0.03208572543250873, -0.02608572543250873, -0.04808572543250873,
    -0.02385756611083844, 0.00214243388916156, 0.03714243388916156,
    0.07814243388916156, 0.02497284742870274
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(0, 1, -1, -1, 0), c(0, 1, 0, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, -0.8, -0.5, 1), max_iter, c(0, 1, -1, -1, 0),
    c(0, 1, 0, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(1, 2, -5, -2, 8), max_iter, c(0, 1, -1, -1, 0),
    c(0, 1, 0, 0, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

context("5-parameter logistic - weighted fit")

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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.14021510699415423,
    beta = 0.14021510699415423 + exp(-0.22629299711936212),
    eta = -1.3532016035342649,
    phi = -0.36746911119363781,
    nu = exp(0.88825453177528800)
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.93758527225933412, 3), rep(0.93513517257012652, 1),
    rep(0.88595597201927951, 2), rep(0.55164531553540003, 4),
    rep(0.26479534106274697, 3), rep(0.17495286982311324, 3),
    0.14985594300774285
  )

  residuals <- c(
    -0.00958527225933412, -0.04958527225933412, 0.04241472774066588,
    0.01286482742987348, 0.01104402798072049, -0.00295597201927951,
    -0.06364531553540003, -0.01964531553540003, 0.01435468446459997,
    0.04735468446459997, -0.00579534106274697, 0.00020465893725303,
    -0.02179534106274697, -0.03195286982311324, 0.00304713017688676,
    0.04404713017688676, -0.05785594300774285
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic5_new(x, y, w, c(0, 1, -1, 0, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.14021510699416073,
    beta = 0.14021510699416073 + exp(-0.22629299711936962),
    eta = -1.3532016035341854,
    phi = -0.36746911119357900,
    nu = exp(0.88825453177516012)
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.93758527225933452, 3), rep(0.93513517257012633, 1),
    rep(0.88595597201927642, 2), rep(0.55164531553540060, 4),
    rep(0.26479534106274417, 3), rep(0.17495286982311410, 3),
    0.14985594300774688
  )

  residuals <- c(
    -0.00958527225933452, -0.04958527225933452, 0.04241472774066548,
    0.01286482742987367, 0.01104402798072358, -0.00295597201927642,
    -0.06364531553540060, -0.01964531553540060, 0.01435468446459940,
    0.04735468446459940, -0.00579534106274417, 0.00020465893725583,
    -0.02179534106274417, -0.03195286982311410, 0.00304713017688590,
    0.04404713017688590, -0.05785594300774688
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(-0.5, 0.5, -2, -1, 0.1), c(0.5, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0.1, 0.6, -0.5, 0.5, 2), max_iter, c(-0.5, 0.5, -2, -1, 0.1),
    c(0.5, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(2, 3, -3, 2, 4), max_iter, c(-0.5, 0.5, -2, -1, 0.1),
    c(0.5, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011306,
    phi = -0.33982624766042563,
    nu = exp(0.86063297591301602)
  )

  rss_value <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829322, 3), rep(0.95779593343449219, 1),
    rep(0.82640795793649349, 2), rep(0.55492427795777896, 4),
    rep(0.30125409492277097, 3), rep(0.15182795902080914, 3),
    0.07523599271751830
  )

  residuals <- c(
    -0.063572996129829322, -0.103572996129829322, -0.011572996129829322,
    -0.00979593343449219, 0.07059204206350651, 0.05659204206350651,
    -0.06692427795777896, -0.02292427795777896, 0.01107572204222104,
    0.04407572204222104, -0.04225409492277097, -0.03625409492277097,
    -0.05825409492277097, -0.00882795902080914, 0.02617204097919086,
    0.06717204097919086, 0.01676400728248170
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(0, 1, rep(-Inf, 3)), c(0, 1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic5_new(
    x, y, w, c(0, 1, -1, 0, 1), max_iter, c(0, 1, rep(-Inf, 3)),
    c(0, 1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic5_new(
    x, y, w, c(1, 2, -1, 0, 1), max_iter, c(0, 1, rep(-Inf, 3)),
    c(0, 1, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011306,
    phi = -0.33982624766042563,
    nu = exp(0.86063297591301602)
  )

  rss_value <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829322, 3), rep(0.95779593343449219, 1),
    rep(0.82640795793649349, 2), rep(0.55492427795777896, 4),
    rep(0.30125409492277097, 3), rep(0.15182795902080914, 3),
    0.07523599271751830
  )

  residuals <- c(
    -0.063572996129829322, -0.103572996129829322, -0.011572996129829322,
    -0.00979593343449219, 0.07059204206350651, 0.05659204206350651,
    -0.06692427795777896, -0.02292427795777896, 0.01107572204222104,
    0.04407572204222104, -0.04225409492277097, -0.03625409492277097,
    -0.05825409492277097, -0.00882795902080914, 0.02617204097919086,
    0.06717204097919086, 0.01676400728248170
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter, c(0, 1, -2, -1, 0.1), c(0, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, -1, 0.1, 2), max_iter, c(0, 1, -2, -1, 0.1),
    c(0, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(1, 2, -5, 2, 5), max_iter, c(0, 1, -2, -1, 0.1),
    c(0, 1, 0, 1, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

context("5-parameter logistic - drda fit")

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
      y ~ x, mean_function = "logistic5",
      lower_bound = c("a", "b", "c", "d", "e")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = matrix(-Inf, nrow = 5, ncol = 2),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = c(0, -1, -Inf, -Inf, 0),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound[2]' cannot be smaller than 'lower_bound[1]'",
    fixed = TRUE
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
      y ~ x, mean_function = "logistic5",
      upper_bound = c("a", "b", "c", "d", "e")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 5),
      upper_bound = matrix(Inf, nrow = 5, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 5),
      upper_bound = c(1, 0, Inf, Inf, 0)
    ),
    "'upper_bound[2]' cannot be smaller than 'upper_bound[1]'",
    fixed = TRUE
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
      y ~ x, mean_function = "logistic5",
      start = c("a", "b", "c", "d", "e")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, Inf, -1, 0, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(-Inf, 0, -1, 0, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 0, -1, 0, 1, 1)
    ),
    "'start' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, -1, -1, 0, 1)
    ),
    "parameter 'beta' is smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, 0, 0, 1)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, -1, 0, 0)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, -1, 0, -1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )
})
