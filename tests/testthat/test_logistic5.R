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

  start <- c(0, 1, -1, 0, 1)

  lower_bound <- c(0, -1, -Inf, -10, 0)
  upper_bound <- c(3, 2, 0, 5, 2)

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
  expect_equal(object$start, c(0, 1, -1, 0, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10, -Inf))
  expect_equal(object$upper_bound, c(2, 2, 0, 5, log(2)))

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
  expect_equal(object$start, c(0, 1, -1, 0, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10, -Inf))
  expect_equal(object$upper_bound, c(2, 2, 0, 5, log(2)))
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
    logistic5_new(x, y, w, c(0, -1, -1, 0, 1), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 0, -1, 0, 1), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
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
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 4)),
    "'upper_bound' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, rep(Inf, 3), 0)),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, rep(Inf, 3), -1)),
    "'upper_bound[5]' cannot be negative nor zero",
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

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic5_fit"
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
      0.000020085234355644039, 0.0020055320228295701, 0.17419599596817509,
      0.99179942948163673, 0.99999901047074089, 0.99999999990085189,
      # beta
      0.99997991476564436, 0.99799446797717043, 0.82580400403182491,
      0.0082005705183632714, 9.8952925911240417e-07, 9.9148112539280778e-11,
      # eta
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # phi
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # log_nu
      8.6734287058557399e-11, 8.6447455933317324e-7, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09
    ),
    nrow = 6,
    ncol = 5
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, 6),
      # (alpha, beta)
      rep(0, 6),
      # (alpha, eta)
      0.00010861439592700533, 0.0062181479955229955, 0.12097343462887282,
      -0.022373854325206450, -7.5180523838040284e-06, -1.2106201466980158e-09,
      # (alpha, phi)
      -0.000040169863584303702, -0.0040050287894570444, -0.30145946064754610,
      -0.029831805766941934, -3.9541796961521937e-06, -3.9658850116135887e-10,
      # (alpha, log_nu)
      -1.0085382216111325e-10, -1.0052029759688061e-06, -0.0073273532118546466,
      -0.024475960355018867, -0.000011704177808949019, -2.0855236572853840e-09,
      # (beta, alpha)
      rep(0, 6),
      # (beta, beta)
      rep(0, 6),
      # (beta, eta)
      -0.00010861439592700533, -0.0062181479955229955, -0.12097343462887282,
      0.022373854325206450, 7.5180523838040284e-06, 1.2106201466980158e-09,
      # (beta, phi)
      0.000040169863584303702, 0.0040050287894570444, 0.30145946064754610,
      0.029831805766941934, 3.9541796961521937e-06, 3.9658850116135887e-10,
      # (beta, log_nu)
      1.0085382216111325e-10, 1.0052029759688061e-06, 0.0073273532118546466,
      0.024475960355018867, 0.000011704177808949019, 2.0855236572853840e-09,
      # (eta, alpha)
      0.00010861439592700533, 0.0062181479955229955, 0.12097343462887282,
      -0.022373854325206450, -7.5180523838040284e-06, -1.2106201466980158e-09,
      # (eta, beta)
      -0.00010861439592700533, -0.0062181479955229955, -0.12097343462887282,
      0.022373854325206450, 7.5180523838040284e-06, 1.2106201466980158e-09,
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
      -0.000040169863584303702, -0.0040050287894570444, -0.30145946064754610,
      -0.029831805766941934, -3.9541796961521937e-06, -3.9658850116135887e-10,
      # (phi, beta)
      0.000040169863584303702, 0.0040050287894570444, 0.30145946064754610,
      0.029831805766941934, 3.9541796961521937e-06, 3.9658850116135887e-10,
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
      -1.0085382216111325e-10, -1.0052029759688061e-06, -0.0073273532118546466,
      -0.024475960355018867, -0.000011704177808949019, -2.0855236572853840e-09,
      # (log_nu, beta)
      1.0085382216111325e-10, 1.0052029759688061e-06, 0.0073273532118546466,
      0.024475960355018867, 0.000011704177808949019, 2.0855236572853840e-09,
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

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, -log(2))

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

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, -log(2))

  true_gradient <- c(
    -0.59375404772645021, -0.33527664229369063, 0.022267352061200450,
    -0.086374079772717690, -0.010097206230551137
  )

  true_hessian <- matrix(
    c(
      # alpha
      6.9953590538168059, 0.32630453922986950, 0.014184128740313277,
      0.29256817446506200, 0.097473377538659030,
      # beta
      0.32630453922986950, 7.3520318677234551, -0.16159602781805513,
      0.33901032399546030, -0.00064023593037618698,
      # eta
      0.014184128740313277, -0.16159602781805513, 0.018237239090142242,
      -0.077452292586773564, -0.017846532854813774,
      # phi
      0.29256817446506200, 0.33901032399546030, -0.077452292586773564,
      0.21294298805991181, -0.0090213900412146364,
      # log_nu
      0.097473377538659030, -0.00064023593037618698, -0.017846532854813774,
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

  gradient_hessian <- rss_gh(c(9 / 10, -2))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[2:3])
  expect_equal(gradient_hessian$H, true_hessian[2:3, 2:3])
})

context("5-parameter logistic - support functions")

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

  theta <- c(
    0, 1, -1.7617462932355768, -0.47836972294568214, 1.3765334489390748
  )

  true_value <- c(
    0.093212121358460102, 0.92029387542791528, -1.7617462932355768,
    -0.47836972294568214, 1.3765334489390748
  )

  object <- logistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_equal(result, true_value)
})

context("5-parameter logistic - fit")

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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.093212121358460102,
    beta = 0.92029387542791528,
    eta = -1.7617462932355768,
    phi = -0.47836972294568214,
    nu = 3.96114628361664067
  )

  rss_value <- 0.024883087882351184

  fitted_values <- c(
    rep(0.9202839186685335, 3), rep(0.9197191693822362, 2),
    rep(0.8900274208207810, 2), rep(0.5533792934125556, 5),
    rep(0.26271803120343568, 3), rep(0.15412982230606348, 4),
    0.11508521369102612
  )

  residuals <- c(
    0.0077160813314665, -0.0322839186685335, 0.0597160813314665,
    0.0282808306177638, -0.0637191693822362, 0.0069725791792190,
    -0.0070274208207810, -0.0653792934125556, -0.0213792934125556,
    0.0326207065874444, 0.0126207065874444, 0.0456207065874444,
    -0.00371803120343568, 0.00228196879656432, -0.01971803120343568,
    -0.03712982230606348, -0.01112982230606348, 0.02387017769393652,
    0.06487017769393652, -0.02308521369102612
  )

  object <- logistic5_new(x, y, w, NULL, 10000, NULL, NULL)

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

  object <- logistic5_new(x, y, w, c(0, 1, -1, 0, 1), 10000, NULL, NULL)

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
    alpha = 0.093212121358460102,
    beta = 0.92029387542791528,
    eta = -1.7617462932355768,
    phi = -0.47836972294568214,
    nu = 3.96114628361664067
  )

  rss_value <- 0.024883087882351184

  fitted_values <- c(
    rep(0.9202839186685335, 3), rep(0.9197191693822362, 2),
    rep(0.8900274208207810, 2), rep(0.5533792934125556, 5),
    rep(0.26271803120343568, 3), rep(0.15412982230606348, 4),
    0.11508521369102612
  )

  residuals <- c(
    0.0077160813314665, -0.0322839186685335, 0.0597160813314665,
    0.0282808306177638, -0.0637191693822362, 0.0069725791792190,
    -0.0070274208207810, -0.0653792934125556, -0.0213792934125556,
    0.0326207065874444, 0.0126207065874444, 0.0456207065874444,
    -0.00371803120343568, 0.00228196879656432, -0.01971803120343568,
    -0.03712982230606348, -0.01112982230606348, 0.02387017769393652,
    0.06487017769393652, -0.02308521369102612
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000,
    c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(-0.1, 1.2, -1.3, 0.3, 2), 10000,
    c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(-1, 2, 1, -2, 0.1), 10000,
    c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605610,
    phi = -0.30897155961772727,
    nu = 2.34190488025700727
  )

  rss_value <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706232, 3), rep(0.96390948874813840, 2),
    rep(0.83729075842125321, 2), rep(0.55558345438725121, 5),
    rep(0.29108572543250870, 3), rep(0.14085756611083843, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706232, -0.105387436096706232, -0.013387436096706232,
    -0.01590948874813840, -0.10790948874813840, 0.05970924157874679,
    0.04570924157874679, -0.06758345438725121, -0.02358345438725121,
    0.03041654561274879, 0.01041654561274879, 0.04341654561274879,
    -0.03208572543250870, -0.02608572543250870, -0.04808572543250870,
    -0.02385756611083843, 0.00214243388916157, 0.03714243388916157,
    0.07814243388916157, 0.02497284742870274
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000, c(0, 1, -Inf, -Inf, -Inf), c(0, 1, Inf, Inf, Inf)
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
    x, y, w, c(0, 1, -1, 0, 1), 10000,
    c(0, 1, -Inf, -Inf, -Inf), c(0, 1, Inf, Inf, Inf)
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
    x, y, w, c(1, 3, -1, 0, 1), 10000,
    c(0, 1, -Inf, -Inf, -Inf), c(0, 1, Inf, Inf, Inf)
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605610,
    phi = -0.30897155961772727,
    nu = 2.34190488025700727
  )

  rss_value <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706232, 3), rep(0.96390948874813840, 2),
    rep(0.83729075842125321, 2), rep(0.55558345438725121, 5),
    rep(0.29108572543250870, 3), rep(0.14085756611083843, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706232, -0.105387436096706232, -0.013387436096706232,
    -0.01590948874813840, -0.10790948874813840, 0.05970924157874679,
    0.04570924157874679, -0.06758345438725121, -0.02358345438725121,
    0.03041654561274879, 0.01041654561274879, 0.04341654561274879,
    -0.03208572543250870, -0.02608572543250870, -0.04808572543250870,
    -0.02385756611083843, 0.00214243388916157, 0.03714243388916157,
    0.07814243388916157, 0.02497284742870274
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000, c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, -1.2, -0.3, 2), 10000,
    c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(-1, 2, 1, 3, 5), 10000,
    c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

context("5-parameter logistic - weighted fit")

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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.14021510699415424,
    beta = 0.93769951379161088,
    eta = -1.3532016035342649,
    phi = -0.36746911119363776,
    nu = 2.43088291720838878
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.9375852722593340, 3), rep(0.9351351725701264, 1),
    rep(0.8859559720192794, 2), rep(0.5516453155354000, 4),
    rep(0.26479534106274689, 3), rep(0.17495286982311317, 3),
    0.14985594300774278
  )

  residuals <- c(
    -0.0095852722593340, -0.0495852722593340, 0.0424147277406660,
    0.0128648274298736, 0.0110440279807206, -0.0029559720192794,
    -0.0636453155354000, -0.0196453155354000, 0.0143546844646000,
    0.0473546844646000, -0.00579534106274689, 0.00020465893725311,
    -0.02179534106274689, -0.03195286982311317, 0.00304713017688683,
    0.04404713017688683, -0.05785594300774278
  )

  object <- logistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic5_new(x, y, w, c(0, 1, -1, 0, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.14021510699415424,
    beta = 0.93769951379161088,
    eta = -1.3532016035342649,
    phi = -0.36746911119363776,
    nu = 2.43088291720838878
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.9375852722593340, 3), rep(0.9351351725701264, 1),
    rep(0.8859559720192794, 2), rep(0.5516453155354000, 4),
    rep(0.26479534106274689, 3), rep(0.17495286982311317, 3),
    0.14985594300774278
  )

  residuals <- c(
    -0.0095852722593340, -0.0495852722593340, 0.0424147277406660,
    0.0128648274298736, 0.0110440279807206, -0.0029559720192794,
    -0.0636453155354000, -0.0196453155354000, 0.0143546844646000,
    0.0473546844646000, -0.00579534106274689, 0.00020465893725311,
    -0.02179534106274689, -0.03195286982311317, 0.00304713017688683,
    0.04404713017688683, -0.05785594300774278
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000, c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0.1, 1.2, -0.5, 0.5, 2), 10000,
    c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(2, 3, -3, 2, 6), 10000,
    c(-0.5, 0.9, -2, -1, 0.5), c(0.5, 1.5, 0, 1, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 5)
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011087,
    phi = -0.33982624766042181,
    nu = 2.36465699101360594
  )

  rss_value <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829225, 3), rep(0.95779593343449193, 1),
    rep(0.82640795793649317, 2), rep(0.55492427795777899, 4),
    rep(0.30125409492277104, 3), rep(0.15182795902080912, 3),
    0.07523599271751825
  )

  residuals <- c(
    -0.063572996129829225, -0.103572996129829225, -0.011572996129829225,
    -0.00979593343449193, 0.07059204206350683, 0.05659204206350683,
    -0.06692427795777899, -0.02292427795777899, 0.01107572204222101,
    0.04407572204222101, -0.04225409492277104, -0.03625409492277104,
    -0.05825409492277104, -0.00882795902080912, 0.02617204097919088,
    0.06717204097919088, 0.01676400728248175
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000, c(0, 1, rep(-Inf, 3)), c(0, 1, rep(Inf, 3))
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
    x, y, w, c(0, 1, -1, 0, 1), 10000,
    c(0, 1, rep(-Inf, 3)), c(0, 1, rep(Inf, 3))
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
    x, y, w, c(1, 2, -1, 0, 1), 10000,
    c(0, 1, rep(-Inf, 3)), c(0, 1, rep(Inf, 3))
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011087,
    phi = -0.33982624766042181,
    nu = 2.36465699101360594
  )

  rss_value <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829225, 3), rep(0.95779593343449193, 1),
    rep(0.82640795793649317, 2), rep(0.55492427795777899, 4),
    rep(0.30125409492277104, 3), rep(0.15182795902080912, 3),
    0.07523599271751825
  )

  residuals <- c(
    -0.063572996129829225, -0.103572996129829225, -0.011572996129829225,
    -0.00979593343449193, 0.07059204206350683, 0.05659204206350683,
    -0.06692427795777899, -0.02292427795777899, 0.01107572204222101,
    0.04407572204222101, -0.04225409492277104, -0.03625409492277104,
    -0.05825409492277104, -0.00882795902080912, 0.02617204097919088,
    0.06717204097919088, 0.01676400728248175
  )

  object <- logistic5_new(
    x, y, w, NULL, 10000, c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, -1, 0.1, 2), 10000,
    c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(1, 2, -5, 2, 5), 10000,
    c(0, 1, -2, -2, 1), c(0, 1, 0, 2, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

context("5-parameter logistic - general functions")

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
    phi = -3 / 2,
    nu = 1 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      3317.1075212912263, 77.779705234498998, 30.105094962347675,
      99.367820533674209, 101.76674312410457, 42879.475233843036,
      # beta
      77.779705234498998, 2509.3330682397757, -50.270268072025062,
      57.263278623851573, -32.834177409666645, 5566.4533586419151,
      # eta
      30.105094962347675, -50.270268072025062, -28.431570455178267,
      -53.631172444195592, -61.035018201080491, 329.86733099061305,
      # phi
      99.367820533674209, 57.263278623851573, -53.631172444195592,
      -13.327710758643254, -75.171417594972841, 1404.2234763310936,
      # nu
      101.76674312410457, -32.834177409666645, -61.035018201080491,
      -75.171417594972841, -95.495915267527460, 1308.4123739466697,
      # sigma
      42879.475233843036, 5566.4533586419151, 329.86733099061305,
      1404.2234763310936, 1308.4123739466697, 540135.75988991146
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "nu", "sigma"
  )

  object <- logistic5_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
  expect_equal(fim, true_value)
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
      y ~ x, mean_function = "logistic6",
      start = c(0, -1, -1, 0, 1, 1)
    ),
    "parameter 'beta' cannot be smaller than 'alpha'"
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

context("5-parameter logistic - Area under and above the curve")

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

  result <- drda(y ~ x, mean_function = "logistic5")

  expect_equal(nauc(result), 0.53873944885547561)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.48525519463583185)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.52493778668672264)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.47542532439305309
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

  result <- drda(y ~ x, mean_function = "logistic5")

  expect_equal(naac(result), 1 - 0.53873944885547561)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.48525519463583185)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.52493778668672264)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.47542532439305309
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

  result <- drda(y ~ x, mean_function = "logistic5")

  expect_equal(nauc(result), 0.53465218185475628)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.56806493427090424)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50308844614941488)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.61344155711817374
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

  result <- drda(y ~ x, mean_function = "logistic5")

  expect_equal(naac(result), 1 - 0.53465218185475628)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.56806493427090424)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50308844614941488)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.61344155711817374
  )
})
