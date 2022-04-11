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

  start <- c(0, -1, 1, 0, 1, 1)

  lower_bound <- c(0, -1, 0, -10, 0, 0.5)
  upper_bound <- c(3, 2, Inf, 5, 2, 1)

  object <- logistic6_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic6"))
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

  object <- logistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic6"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, -1, 0, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, -Inf, -10, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(3, 2, Inf, 5, log(2), 0))

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

  object <- logistic6_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic6"))
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

  object <- logistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic6"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, -1, 0, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, -Inf, -10, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(3, 2, Inf, 5, log(2), 0))
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
    logistic6_new(x, y, w, c(0, -1, 1, 0, 1), max_iter, NULL, NULL),
    "'start' must be of length 6"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, 0, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, -1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, 1, 0, 0, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, 1, 0, -1, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, 1, 0, 1, 0), max_iter, NULL, NULL),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, 1, 0, 1, -1), max_iter, NULL, NULL),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 5)),
    "'lower_bound' must be of length 6"
  )

  expect_error(
    logistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 6)),
    "'lower_bound' must be of length 6"
  )

  expect_error(
    logistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 6), rep(Inf, 5)),
    "'upper_bound' must be of length 6"
  )

  expect_error(
    logistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, 0, Inf, Inf, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, -1, Inf, Inf, Inf)
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, Inf, Inf, 0, Inf)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, Inf, Inf, -1, Inf)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, rep(Inf, 4), 0)),
    "'upper_bound[6]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, rep(Inf, 4), -1)),
    "'upper_bound[6]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(19 / 45, -43 / 50, 2, -3 / 2, 1 / 2, 3 / 2)

  true_value <- c(
    0.42222222083459441, 0.42220851001976607, 0.36818494892151421,
    0.052377517227330391, 0.040126833251911149, 0.040001268645102526
  )

  value <- logistic6_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic6"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic6_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- c(-500, -log(c(100, 10, 1, 0.1)), 500)
  theta <- c(19 / 45, -43 / 50, 2, -3 / 2, 1 / 2, 3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # beta
      0.0, 0.000015944421460645258, 0.062834038721753500, 0.43005198255219980,
      0.44429696391896636, 0.44444444444444444,
      # log_eta
      0.0, 0.00016929477606169280, 0.10825016156008817, -0.036225775505554982,
      -0.00096434839234998967, 0.0,
      # phi
      0.0, 0.000054520289041040683, 0.13487686540035359, 0.024150517003703322,
      0.00025360336948848910, 0.0,
      # log_nu
      0.0, -0.00012421034982980512, -0.082096662214012938, -0.30001789547377784,
      -0.30985270691418151, -0.30995554930935233,
      # log_xi
      0.0, 1.6426039178950299e-07, 0.040636113901239228, 0.72761415148793200,
      0.76406397625587790, 0.76444444444444444
    ),
    nrow = 6,
    ncol = 6
  )

  true_hessian <- array(
    c(
      # (alpha, *)
      rep(0, 36),
      # (beta, alpha)
      rep(0, 6),
      # (beta, beta)
      rep(0, 6),
      # (beta, log_eta)
      0.0, -0.00019685439076941023, -0.12587228088382345, 0.042122994773901142,
      0.0011213353399418484, 0.0,
      # (beta, phi)
      0.0, -0.000063395684931442655, -0.15683356441901580,
      -0.028081996515934095, -0.00029488763894010361, 0.0,
      # (beta, log_nu)
      0.0, 0.00014443063933698270, 0.095461235132573184, 0.34885801799276493,
      0.36029384524904826, 0.36041342942947945,
      # (beta, log_xi)
      0.0, -1.9100045556918952e-07, -0.047251295233999102, -0.84606296684643256,
      -0.88844648401846268, -0.88888888888888889,
      # (log_eta, alpha)
      rep(0, 6),
      # (log_eta, beta)
      0.0, -0.00019685439076941023, -0.12587228088382345, 0.042122994773901142,
      0.0011213353399418484, 0.0,
      # (log_eta, log_eta)
      0.0, -0.0019145696794883463, -0.043268092124344362, 0.067129154367644463,
      0.0063660344555865443, 0.0,
      # (log_eta, phi)
      0.0, -0.00061657479777686132, -0.053910909263132009,
      -0.044752769578429642, -0.0016741333329569523, 0.0,
      # (log_eta, log_nu)
      0.0, 0.0013652557558011332, 0.096912112894594106, -0.028794959028545277,
      -0.00078185926316106587, 0.0,
      # (log_eta, log_xi)
      0.0, -3.0420087144345961e-06, -0.12210643578843353, 0.10690319430217645,
      0.0028925651359227001, 0.0,
      # (phi, alpha)
      rep(0, 6),
      # (phi, beta)
      0.0, -0.000063395684931442655, -0.15683356441901580,
      -0.028081996515934095, -0.00029488763894010361, 0.0,
      # (phi, log_eta)
      0.0, -0.00061657479777686132, -0.053910909263132009,
      -0.044752769578429642, -0.0016741333329569523, 0.0,
      # (phi, phi)
      0.0, -0.00021612183765198489, -0.23522462142825546, 0.045935524388088642,
      0.00050695425751211717, 0.0,
      # (phi, log_nu)
      0.0, 0.00043967179704409577, 0.12074995379376314, 0.019196639352363518,
      0.00020561256199146683, 0.0,
      # (phi, log_xi)
      0.0, -9.7965925608892296e-07, -0.15214142008657945, -0.071268796201450965,
      -0.00076068386773303679, 0.0,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, beta)
      0.0, 0.00014443063933698270, 0.095461235132573184, 0.34885801799276493,
      0.36029384524904826, 0.36041342942947945,
      # (log_nu, log_eta)
      0.0, 0.0013652557558011332, 0.096912112894594106, -0.028794959028545277,
      -0.00078185926316106587, 0.0,
      # (log_nu, phi)
      0.0, 0.00043967179704409577, 0.12074995379376314, 0.019196639352363518,
      0.00020561256199146683, 0.0,
      # (log_nu, log_nu)
      0.0, -0.0010280311514900970, -0.084711033634174417, 0.056446357685499487,
      0.058583754728948002, 0.058603228690468293,
      # (log_nu, log_xi)
      0.0, 1.1603962252525160e-06, -0.0042561990702717135,
      -0.14925193871840975, -0.14458817039325487, -0.14453334582573979,
      # (log_xi, alpha)
      rep(0, 6),
      # (log_xi, beta)
      0.0, -1.9100045556918952e-07, -0.047251295233999102,
      -0.84606296684643256, -0.88844648401846268, -0.88888888888888889,
      # (log_xi, log_eta)
      0.0, -3.0420087144345961e-06, -0.12210643578843353, 0.10690319430217645,
      0.0028925651359227001, 0.0,
      # (log_xi, phi)
      0.0, -9.7965925608892296e-07, -0.15214142008657945, -0.071268796201450965,
      -0.00076068386773303679, 0.0,
      # (log_xi, log_nu)
      0.0, 1.1603962252525160e-06, -0.0042561990702717135, -0.14925193871840975,
      -0.14458817039325487, -0.14453334582573979,
      # (log_xi, log_xi)
      0.0, 1.6130884446545550e-07, -0.0052015177591887325, -1.4195939048751385,
      -1.5277476105778893, -1.5288888888888889
    ),
    dim = c(6, 6, 6)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic6"
  )

  gh <- gradient_hessian(object, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 6 * 6)
  expect_length(gh$H, 6 * 6 * 6)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(19 / 45, -43 / 50, log(2), -3 / 2, -log(2), log(3 / 2))

  true_value <- 2.1010793055397967

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic6"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(19 / 45, NA, NA, -3 / 2, -log(2), NA)
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(-43 / 50, log(2), log(3 / 2)))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  x <- c(-500, -log(c(10, 1, 0.1)), 500)
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(19 / 45, -43 / 50, log(2), -3 / 2, -log(2), log(3 / 2))

  true_gradient <- c(
    -5.3759161191064853, -1.0688047637852019, -0.086377512187944453,
    -0.22553860902584356, 0.79792437296940233, -1.7255748891976254
  )

  true_hessian <- matrix(
    c(
      # alpha
      15.0, 4.1591272702788589, 0.24844154009975457, 0.45394604368642136,
      -3.0156032531743775, 6.9267258830364266,
      # delta
      4.1591272702788589, 1.7639255258434534, 0.087972590786873860,
      0.30890142735400208, -2.1652815479152137, 5.0171107183441425,
      # log_eta
      0.24844154009975457, 0.087972590786873860, -0.019331559108568847,
      0.19228267597978567, -0.087250267492429096, -0.057597786349391573,
      # phi
      0.45394604368642136, 0.30890142735400208, 0.19228267597978567,
      0.29960404329994302, -0.24582663358150991, 0.38150843339891208,
      # log_nu
      -3.0156032531743775, -2.1652815479152137, -0.087250267492429096,
      -0.24582663358150991, 0.85924377265794259, -1.7630267519128655,
      # log_xi
      6.9267258830364266, 5.0171107183441425, -0.057597786349391573,
      0.38150843339891208, -1.7630267519128655, 8.4386762188524690
    ),
    nrow = 6,
    ncol = 6
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "logistic6"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta)

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 6)
  expect_length(gh$H, 6 * 6)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  known_param <- c(19 / 45, NA, NA, -3 / 2, -log(2), NA)
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(c(-43 / 50, log(2), log(3 / 2)))

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 3)
  expect_length(gh$H, 3 * 3)

  expect_equal(gh$G, true_gradient[c(2, 3, 6)])
  expect_equal(gh$H, true_hessian[c(2, 3, 6), c(2, 3, 6)])
})

test_that("mle_asy", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  true_value <- c(
    0.92000003919821951, -0.98167504806989653, 2.6045134462190666,
    0.16970988944039289, 1.3142240447570685, 0.90296118014410684
  )
  theta <- c(0, 1, true_value[-c(1, 2)])

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 6)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.028004265218879379

  fitted_values <- c(
    rep(0.92000003919821842, 3), rep(0.92000001908978876, 2),
    rep(0.8899999155532745, 2), rep(0.5542000029977990, 5),
    rep(0.2556666657219331, 3), rep(0.1498000001928848, 4),
    0.1498000001928848
  )

  residuals <- c(
    0.00799996080178158, -0.03200003919821842, 0.05999996080178158,
    0.02799998091021124, -0.06400001908978876, 0.0070000844467255,
    -0.0069999155532745, -0.0662000029977990, -0.0222000029977990,
    0.0317999970022010, 0.0117999970022010, 0.0447999970022010,
    0.0033333342780669, 0.0093333342780669, -0.0126666657219331,
    -0.0328000001928848, -0.0068000001928848, 0.0281999998071152,
    0.0691999998071152, -0.0578000001928848
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- logistic6_new(x, y, w, c(1, -1, 1, 0, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.028004265218879379

  fitted_values <- c(
    rep(0.92000003919821842, 3), rep(0.92000001908978876, 2),
    rep(0.8899999155532745, 2), rep(0.5542000029977990, 5),
    rep(0.2556666657219331, 3), rep(0.1498000001928848, 4),
    0.1498000001928848
  )

  residuals <- c(
    0.00799996080178158, -0.03200003919821842, 0.05999996080178158,
    0.02799998091021124, -0.06400001908978876, 0.0070000844467255,
    -0.0069999155532745, -0.0662000029977990, -0.0222000029977990,
    0.0317999970022010, 0.0117999970022010, 0.0447999970022010,
    0.0033333342780669, 0.0093333342780669, -0.0126666657219331,
    -0.0328000001928848, -0.0068000001928848, 0.0281999998071152,
    0.0691999998071152, -0.0578000001928848
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000,
    c(0.5, -0.5, 5, -1, 0.5, 0), c(1.0, 0.0, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0.6, -0.1, 10.0, 0.0, 2.0, 0.3), 10000,
    c(0.5, -0.5, 5, -1, 0.5, 0), c(1.0, 0.0, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, -2, 1, -2, 0.1, 2), 10000,
    c(0.5, -0.5, 5, -1, 0.5, 0), c(1.0, 0.0, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.030004266666749054

  fitted_values <- c(
    rep(0.90000000000000002, 3), rep(0.89999999998970159, 2),
    rep(0.89000000001556041, 2), rep(0.55419999999971517, 5),
    rep(0.25566666666689792, 3), rep(0.14979999999994013, 4),
    0.14979999999994013
  )

  residuals <- c(
    0.02800000000000000000027, -0.01199999999999999999973,
    0.08000000000000000000027, 0.048000000010298387944347,
    -0.04399999998970161205565, 0.00699999998443959, -0.00700000001556041,
    -0.06619999999971517, -0.02219999999971517,
    0.03180000000028483, 0.01180000000028483, 0.04480000000028483,
    0.00333333333310208, 0.00933333333310208, -0.01266666666689792,
    -0.03279999999994013, -0.00679999999994013, 0.02820000000005987,
    0.06920000000005987, -0.05779999999994013
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0.9, -0.5, rep(-Inf, 4)), c(0.9, -0.5, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic6_new(
    x, y, w, c(0.9, -0.5, 1, 0, 1, 1), 10000,
    c(0.9, -0.5, rep(-Inf, 4)), c(0.9, -0.5, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic6_new(
    x, y, w, c(1, 3, 1, 0, 1, 1), 10000,
    c(0.9, -0.5, rep(-Inf, 4)), c(0.9, -0.5, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities and inequalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.030004266666749054

  fitted_values <- c(
    rep(0.90000000000000002, 3), rep(0.89999999998970159, 2),
    rep(0.89000000001556041, 2), rep(0.55419999999971517, 5),
    rep(0.25566666666689792, 3), rep(0.14979999999994013, 4),
    0.14979999999994013
  )

  residuals <- c(
    0.02800000000000000000027, -0.01199999999999999999973,
    0.08000000000000000000027, 0.048000000010298387944347,
    -0.04399999998970161205565, 0.00699999998443959, -0.00700000001556041,
    -0.06619999999971517, -0.02219999999971517,
    0.03180000000028483, 0.01180000000028483, 0.04480000000028483,
    0.00333333333310208, 0.00933333333310208, -0.01266666666689792,
    -0.03279999999994013, -0.00679999999994013, 0.02820000000005987,
    0.06920000000005987, -0.05779999999994013
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0.9, -0.5, 5, -1, 0.5, 0), c(0.9, -0.5, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0.9, -0.5, 10, 0.0, 2, 0.5), 10000,
    c(0.9, -0.5, 5, -1, 0.5, 0), c(0.9, -0.5, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, -2, 10, 3), 10000,
    c(0.9, -0.5, 5, -1, 0.5, 0), c(0.9, -0.5, 20, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit (weighted)", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.014601665143055133

  fitted_values <- c(
    rep(0.93751313951088907, 3), rep(0.93751233830256449, 2),
    rep(0.8836847290952300, 2), rep(0.5522145569269031, 5),
    rep(0.2606149158850656, 3), rep(0.1756891937008905, 4),
    0.1756891937008905
  )

  residuals <- c(
    -0.00951313951088907, -0.04951313951088907, 0.04248686048911093,
    0.01048766169743551, -0.08151233830256449, 0.0133152709047700,
    -0.0006847290952300, -0.0642145569269031, -0.0202145569269031,
    0.0337854430730969, 0.0137854430730969, 0.0467854430730969,
    -0.0016149158850656, 0.0043850841149344, -0.0176149158850656,
    -0.0586891937008905, -0.0326891937008905, 0.0023108062991095,
    0.0433108062991095, -0.0836891937008905
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- logistic6_new(x, y, w, c(1, -1, 1, 0, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): inequalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.014601665143055133

  fitted_values <- c(
    rep(0.93751313951088907, 3), rep(0.93751233830256449, 2),
    rep(0.8836847290952300, 2), rep(0.5522145569269031, 5),
    rep(0.2606149158850656, 3), rep(0.1756891937008905, 4),
    0.1756891937008905
  )

  residuals <- c(
    -0.00951313951088907, -0.04951313951088907, 0.04248686048911093,
    0.01048766169743551, -0.08151233830256449, 0.0133152709047700,
    -0.0006847290952300, -0.0642145569269031, -0.0202145569269031,
    0.0337854430730969, 0.0137854430730969, 0.0467854430730969,
    -0.0016149158850656, 0.0043850841149344, -0.0176149158850656,
    -0.0586891937008905, -0.0326891937008905, 0.0023108062991095,
    0.0433108062991095, -0.0836891937008905
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0.5, -1.0, 5, -1, 1, 0), c(1.0, -0.5, 25, 1, 9, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0.7, -0.9, 10, 0, 3, 1.5), 10000,
    c(0.5, -1.0, 5, -1, 1, 0), c(1.0, -0.5, 25, 1, 9, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, -2, 0.5, 5), 10000,
    c(0.5, -1.0, 5, -1, 1, 0), c(1.0, -0.5, 25, 1, 9, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022354118749285909

  fitted_values <- c(
    rep(0.9, 3), rep(0.899999999611796, 2),
    rep(0.88368384359209389, 2), rep(0.55221458986513316, 5),
    rep(0.26061491383540554, 3), rep(0.17568919404694744, 4),
    0.17568919404694744
  )

  residuals <- c(
    0.028, -0.012, 0.08, 0.04800000038820447, -0.04399999961179553,
    0.01331615640790611, -0.00068384359209389, -0.06421458986513316,
    -0.02021458986513316, 0.03378541013486684, 0.01378541013486684,
    0.04678541013486684, -0.00161491383540554, 0.00438508616459446,
    -0.01761491383540554, -0.05868919404694744, -0.03268919404694744,
    0.00231080595305256, 0.04331080595305256, -0.08368919404694744
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0.9, -0.8, rep(-Inf, 4)), c(0.9, -0.8, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic6_new(
    x, y, w, c(0.9, -0.8, 1, 0, 1, 1), 10000,
    c(0.9, -0.8, rep(-Inf, 4)), c(0.9, -0.8, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic6_new(
    x, y, w, c(1, 3, 1, 0, 1, 1), 10000,
    c(0.9, -0.8, rep(-Inf, 4)), c(0.9, -0.8, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities and inequalities", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022354118749285909

  fitted_values <- c(
    rep(0.9, 3), rep(0.899999999611796, 2),
    rep(0.88368384359209389, 2), rep(0.55221458986513316, 5),
    rep(0.26061491383540554, 3), rep(0.17568919404694744, 4),
    0.17568919404694744
  )

  residuals <- c(
    0.028, -0.012, 0.08, 0.04800000038820447, -0.04399999961179553,
    0.01331615640790611, -0.00068384359209389, -0.06421458986513316,
    -0.02021458986513316, 0.03378541013486684, 0.01378541013486684,
    0.04678541013486684, -0.00161491383540554, 0.00438508616459446,
    -0.01761491383540554, -0.05868919404694744, -0.03268919404694744,
    0.00231080595305256, 0.04331080595305256, -0.08368919404694744
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0.9, -0.8, 3, -1, 2, 0), c(0.9, -0.8, 20, 1, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0.9, -0.8, 5, -0.1, 6, 0.5), 10000,
    c(0.9, -0.8, 3, -1, 2, 0), c(0.9, -0.8, 20, 1, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 1, 1, -2, 20, 3), 10000,
    c(0.9, -0.8, 3, -1, 2, 0), c(0.9, -0.8, 20, 1, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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

  theta <- c(
    alpha = 9 / 10,
    delta = -4 / 5,
    eta = 2,
    phi = 1 / 10,
    nu = 1 / 2,
    xi = 3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5982.0000000000000, 1059.1500861068465, 30.211576920320743,
      502.69044249645634, -1478.7444026902667, 962.19661101515078,
      -34445.389244580913,
      # delta
      1059.1500861068465, 352.40713784936091, 4.6512167334071279,
      -21.251700020328802, 62.374462295725370, -44.560403416923947,
      -12826.357513890556,
      # eta
      30.211576920320743, 4.6512167334071279, 2.1514973627787064,
      57.014121258535741, -8.0872587583127677, 3.5872015251675816,
      63.018884697913377,
      # phi
      502.69044249645634, -21.251700020328802, 57.014121258535741,
      156.14784701809923, -28.281873927206431, -74.717762361050467,
      -4231.3465829708541,
      # nu
      -1478.7444026902667, 62.374462295725370, -8.0872587583127677,
      -28.281873927206431, 1596.7868783260568, -608.11234964965376,
      17241.948008527057,
      # xi
      962.19661101515078, -44.560403416923947, 3.5872015251675816,
      -74.717762361050467, -608.11234964965376, -227.14182875480719,
      -12270.999153826308,
      # sigma
      -34445.389244580913, -12826.357513890556, 63.018884697913377,
      -4231.3465829708541, 17241.948008527057, -12270.999153826308,
      385240.26258494415
    ),
    nrow = 7,
    ncol = 7
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "xi", "sigma"
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 7 * 7)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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
      y ~ x, mean_function = "logistic6",
      lower_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = matrix(-Inf, nrow = 6, ncol = 2),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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
      y ~ x, mean_function = "logistic6",
      upper_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 6),
      upper_bound = matrix(Inf, nrow = 6, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
  )
})

test_that("drda: 'start' argument errors", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
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
      y ~ x, mean_function = "logistic6",
      start = c("a", "b", "c", "d", "e", "f")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, Inf, -1, 0, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(-Inf, 0, -1, 0, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 0, -1, 0, 1, 1, 1)
    ),
    "'start' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 0, 0, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 1, 0, 0, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 1, 0, -1, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 1, 0, 1, 0)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 1, 0, 1, -1)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(nauc(result), 0.53306791881552503)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.39522120881029673)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50125330577908701)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.34168870519391313
  )
})

test_that("naac: decreasing", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  )

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(naac(result), 1 - 0.53306791881552503)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.39522120881029673)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50125330577908701)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.34168870519391313
  )
})

test_that("nauc: increasing", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(nauc(result), 0.54157960218789312)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.68754939436736376)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50492003174203455)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.69946687828023036
  )
})

test_that("naac: increasing", {
  x <- round(
    rep(
      -log(c(10000, 100, 2, 1, 0.8, 0.001, 0.0001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- rev(c(
    0.928, 0.888, 0.98, 0.948, 0.856, 0.897, 0.883, 0.488, 0.532, 0.586, 0.566,
    0.599, 0.259, 0.265, 0.243, 0.117, 0.143, 0.178, 0.219, 0.092
  ))

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(naac(result), 1 - 0.54157960218789312)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.68754939436736376)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50492003174203455)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.69946687828023036
  )
})
