context("6-parameter logistic - core functions")

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

  start <- c(0, 1, -1, 0, 1, 1)

  lower_bound <- c(0, -1, -Inf, -10, 0, 0.5)
  upper_bound <- c(3, 2, 0, 5, 2, 1)

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
  expect_equal(object$start, c(0, 1, -1, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(2, 2, 0, 5, log(2), 0))

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
  expect_equal(object$start, c(0, 1, -1, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, 0, -Inf, -10, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(2, 2, 0, 5, log(2), 0))
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
    logistic6_new(x, y, w, c(0, 1, -1, 0, 1), max_iter, NULL, NULL),
    "'start' must be of length 6"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, -1, -1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 0, -1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 1, 0, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 1, -1, 0, 0, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 1, -1, 0, -1, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 1, -1, 0, 1, 0), max_iter, NULL, NULL),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    logistic6_new(x, y, w, c(0, 1, -1, 0, 1, -1), max_iter, NULL, NULL),
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
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2, 3 / 2)

  true_value <- c(
    0.42221710418125004, 0.42171092652477733, 0.37575797785755273,
    0.046454735980600713, 0.040000850149265818, 0.040000000085266528
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
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2, 3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      0.55556150676598832, 0.55615008543630543, 0.60958374667726427,
      0.99249449304581312, 0.99999901145434207, 0.99999999990085287,
      # beta
      0.44443849323401168, 0.44384991456369457, 0.39041625332273573,
      0.0075055069541868754, 9.8854565792745283e-07, 9.9147125297712429e-11,
      # eta
      -0.000027676835132687997, -0.0015860669501348013, -0.033819316025770067,
      0.016847800200509503, 6.4558872590255383e-06, 1.0411177759771462e-09,
      # phi
      0.000010235979146562789, 0.0010215652316203735, 0.084275963560716097,
      0.022463733600679338, 3.3955254655155444e-06, 3.4106101689568101e-10,
      # log_nu
      0.30995139895241304, 0.30954109512977705, 0.27365642866850772,
      0.020345465001078442, 0.000010057277525685746, 1.7935341844903725e-09,
      # log_xi
      -0.76442909037292680, -0.76291107043374447, -0.62937797393474741,
      -0.0016776051608617570, -2.5357988774466817e-09, -2.5470642248697705e-15
    ),
    nrow = 6,
    ncol = 6
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, 6),
      # (alpha, beta)
      rep(0, 6),
      # (alpha, eta)
      0.000032182366433358136, 0.0018442638955055830, 0.039324786076476822,
      -0.019590465349429655, -7.5068456500296957e-06, -1.2106020650897049e-09,
      # (alpha, phi)
      -0.000011902301333212546, -0.0011878665483957831, -0.097995306465948950,
      -0.026120620465906206, -3.9482854250180749e-06, -3.9658257778567560e-10,
      # (alpha, log_nu)
      -0.36040860343303842, -0.35993150596485703, -0.31820514961454386,
      -0.023657517443114468, -0.000011694508750797379, -2.0855048656864797e-09,
      # (alpha, log_xi)
      0.88887103531735675, 0.88710589585319124, 0.73183485341249698,
      0.0019507036754206476, 2.9486033458682346e-09, 2.9617025870578727e-15,
      # (beta, alpha)
      rep(0, 6),
      # (beta, beta)
      rep(0, 6),
      # (beta, eta)
      -0.000032182366433358136, -0.0018442638955055830, -0.039324786076476822,
      0.019590465349429655, 7.5068456500296957e-06, 1.2106020650897049e-09,
      # (beta, phi)
      0.000011902301333212546, 0.0011878665483957831, 0.097995306465948950,
      0.026120620465906206, 3.9482854250180749e-06, 3.9658257778567560e-10,
      # (beta, log_nu)
      0.36040860343303842, 0.35993150596485703, 0.31820514961454386,
      0.023657517443114468, 0.000011694508750797379, 2.0855048656864797e-09,
      # (beta, log_xi)
      -0.88887103531735675, -0.88710589585319124, -0.73183485341249698,
      -0.0019507036754206476, -2.9486033458682346e-09, -2.9617025870578727e-15,
      # (eta, alpha)
      0.000032182366433358136, 0.0018442638955055830, 0.039324786076476822,
      -0.019590465349429655, -7.5068456500296957e-06, -1.2106020650897049e-09,
      # (eta, beta)
      -0.000032182366433358136, -0.0018442638955055830, -0.039324786076476822,
      0.019590465349429655, 7.5068456500296957e-06, 1.2106020650897049e-09,
      # (eta, eta)
      -0.00014966654512113986, -0.0049151222824604186, -0.022033188829639366,
      0.040691115014078006, 0.000048988285040711751, 1.2712117605288108e-08,
      # (eta, phi)
      0.000050234568891350354, 0.0026549841461145773, 0.012767570346070067,
      0.043022953218431006, 0.000024068017042147825, 3.9938472952908821e-09,
      # (eta, log_nu)
      -0.000022443796599450935, -0.0012851291354092528, -0.025441940164379255,
      0.038446220915158790, 0.000069926976157246812, 2.0858233544335177e-08,
      # (eta, log_xi)
      0.000083029949497691868, 0.0047550172812319847, 0.095091408114754192,
      -0.0065681903916336688, -2.8884630598278052e-08, -4.6650084931125022e-14,
      # (phi, alpha)
      -0.000011902301333212546, -0.0011878665483957831, -0.097995306465948950,
      -0.026120620465906206, -3.9482854250180749e-06, -3.9658257778567560e-10,
      # (phi, beta)
      0.000011902301333212546, 0.0011878665483957831, 0.097995306465948950,
      0.026120620465906206, 3.9482854250180749e-06, 3.9658257778567560e-10,
      # (phi, eta)
      0.000050234568891350354, 0.0026549841461145773, 0.012767570346070067,
      0.043022953218431006, 0.000024068017042147825, 3.9938472952908821e-09,
      # (phi, phi)
      -0.000020471547105604366, -0.0020390294716921548, -0.13682175910245932,
      0.072339760025027567, 0.000013551717657746544, 1.3642135032685379e-09,
      # (phi, log_nu)
      8.3005962517132876e-06, 0.00082773507307800837, 0.063399981849819896,
      0.051261627886878386, 0.000036778651599976862, 6.8329736629477352e-09,
      # (phi, log_xi)
      -0.000030707731845927761, -0.0030626451990868243, -0.23696280667266186,
      -0.0087575871888448917, -1.5192102157816607e-08, -1.5282157093078623e-14,
      # (log_nu, alpha)
      -0.36040860343303842, -0.35993150596485703, -0.31820514961454386,
      -0.023657517443114468, -0.000011694508750797379, -2.0855048656864797e-09,
      # (log_nu, beta)
      0.36040860343303842, 0.35993150596485703, 0.31820514961454386,
      0.023657517443114468, 0.000011694508750797379, 2.0855048656864797e-09,
      # (log_nu, eta)
      -0.000022443796599450935, -0.0012851291354092528, -0.025441940164379255,
      0.038446220915158790, 0.000069926976157246812, 2.0858233544335177e-08,
      # (log_nu, phi)
      8.3005962517132876e-06, 0.00082773507307800837, 0.063399981849819896,
      0.051261627886878386, 0.000036778651599976862, 6.8329736629477352e-09,
      # (log_nu, log_nu)
      -0.058602443935219952, -0.058524387549721785, -0.047971111555920137,
      0.053556142044358048, 0.00011061568486094975, 3.6102991956890671e-08,
      # (log_nu, log_xi)
      0.14453556075579985, 0.14475353032653957, 0.15590308930299996,
      -0.0021506440975162657, -2.4930719747841768e-08, -4.8481994812079120e-14,
      # (log_xi, alpha)
      0.88887103531735675, 0.88710589585319124, 0.73183485341249698,
      0.0019507036754206476, 2.9486033458682346e-09, 2.9617025870578727e-15,
      # (log_xi, beta)
      -0.88887103531735675, -0.88710589585319124, -0.73183485341249698,
      -0.0019507036754206476, -2.9486033458682346e-09, -2.9617025870578727e-15,
      # (log_xi, eta)
      0.000083029949497691868, 0.0047550172812319847, 0.095091408114754192,
      -0.0065681903916336688, -2.8884630598278052e-08, -4.6650084931125022e-14,
      # (log_xi, phi)
      -0.000030707731845927761, -0.0030626451990868243, -0.23696280667266186,
      -0.0087575871888448917, -1.5192102157816607e-08, -1.5282157093078623e-14,
      # (log_xi, log_nu)
      0.14453556075579985, 0.14475353032653957, 0.15590308930299996,
      -0.0021506440975162657, -2.4930719747841768e-08, -4.8481994812079120e-14,
      # (log_xi, log_xi)
      1.5288428268799306, 1.5242908182679455, 1.1402745445331639,
      -0.0010235832726989319, -2.5244533240149399e-09, -2.5469500967997703e-15
    ),
    dim = c(6, 6, 6)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "logistic6"
  )

  gradient_hessian <- gradient_hessian(object, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 6)
  expect_length(gradient_hessian$H, 6 * 6 * 6)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

context("6-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, -log(2), log(3 / 2))

  true_value <- 2.0908902035928604

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

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2), NA)
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(9 / 10, -2, log(3 / 2)))

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

  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, -log(2), log(3 / 2))

  true_gradient <- c(
    -2.8218971065888476, -1.7213813512077646, 0.030726379870560891,
    -0.098256479294546431, -1.2078753918106346, 2.9116476844300818
  )

  true_hessian <- matrix(
    c(
      # alpha
      9.5372145480529441, 1.9870627284547445, -0.012747398476564869,
      0.30791044040883690, 2.8519856296043488, -6.7065553751532712,
      # beta
      1.9870627284547445, 1.4886599950375670, 0.0076781035330373560,
      -0.046397988675990002, -0.36478325818387228, 0.85906851661291888,
      # eta
      -0.012747398476564869, 0.0076781035330373560, 0.016681432220990243,
      -0.037754568939522141, -0.0055461232427009715, -0.056647742370221930,
      # phi
      0.30791044040883690, -0.046397988675990002, -0.037754568939522141,
      0.13344876935609788, -0.038587395248642860, 0.14652892520075174,
      # log_nu
      2.8519856296043488, -0.36478325818387228, -0.0055461232427009715,
      -0.038587395248642860, 0.92894109931730305, -2.3526424176254161,
      # log_xi
      -6.7065553751532712, 0.85906851661291888, -0.056647742370221930,
      0.14652892520075174, -2.3526424176254161, -1.4043333859738446
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

  gradient_hessian <- rss_gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6)
  expect_length(gradient_hessian$H, 6 * 6)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2), NA)
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(c(9 / 10, -2, log(3 / 2)))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 3)
  expect_length(gradient_hessian$H, 3 * 3)

  expect_equal(gradient_hessian$G, true_gradient[c(2, 3, 6)])
  expect_equal(gradient_hessian$H, true_hessian[c(2, 3, 6), c(2, 3, 6)])
})

context("6-parameter logistic - support functions")

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
    0, 1, -1.7617462932355762, -0.48082331303125517, 1.3765334489390729,
    0.0043226032383796241
  )

  true_value <- c(
    0.093212121358460317, 0.92119692151923081, -1.7617462932355762,
    -0.48082331303125517, 1.3765334489390729, 0.0043226032383796241
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 6)
  expect_equal(result, true_value)
})

context("6-parameter logistic - fit")

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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.024883087882351184

  fitted_values <- c(
    rep(0.9202839186685335, 3), rep(0.9197191693822362, 2),
    rep(0.890027420820781, 2), rep(0.5533792934125556, 5),
    rep(0.2627180312034356, 3), rep(0.1541298223060635, 4),
    0.11508521369102623
  )

  residuals <- c(
    0.0077160813314665, -0.0322839186685335, 0.0597160813314665,
    0.0282808306177638, -0.0637191693822362, 0.0069725791792190,
    -0.0070274208207810, -0.0653792934125556, -0.0213792934125556,
    0.0326207065874444, 0.0126207065874444, 0.0456207065874444,
    -0.00371803120343560, 0.00228196879656440, -0.01971803120343560,
    -0.03712982230606350, -0.01112982230606350, 0.02387017769393650,
    0.06487017769393650, -0.02308521369102623
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic6_new(x, y, w, c(0, 1, -1, 0, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
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

  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.024883087882351185

  fitted_values <- c(
    rep(0.9202839184386346, 3), rep(0.9197191691524943, 2),
    rep(0.8900274205992924, 2), rep(0.5533792932846442, 5),
    rep(0.26271803115631854, 3), rep(0.15412982228913039, 4),
    0.11508521368494622
  )

  residuals <- c(
    0.0077160815613654, -0.0322839184386346, 0.0597160815613654,
    0.0282808308475057, -0.0637191691524943, 0.0069725794007076,
    -0.0070274205992924, -0.0653792932846442, -0.0213792932846442,
    0.0326207067153558, 0.0126207067153558, 0.0456207067153558,
    -0.00371803115631854, 0.00228196884368146, -0.01971803115631854,
    -0.03712982228913039, -0.01112982228913039, 0.02387017771086961,
    0.06487017771086961, -0.02308521368494622
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000,
    c(-0.5, 0.9, -2, -1, 0.5, 0), c(0.5, 1.5, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(-0.1, 1.2, -1.3, 0.3, 2, 0.3), 10000,
    c(-0.5, 0.9, -2, -1, 0.5, 0), c(0.5, 1.5, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, -2, 0.1, 2), 10000,
    c(-0.5, 0.9, -2, -1, 0.5, 0), c(0.5, 1.5, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
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
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.027858205619881242

  fitted_values <- c(
    rep(0.91990299490116622, 3), rep(0.91975017005499427, 2),
    rep(0.89071256487960993, 2), rep(0.5492685788947333, 5),
    rep(0.2842905938017210, 3), rep(0.14692068718566469, 4),
    0.07590586225363996
  )

  residuals <- c(
    0.00809700509883378, -0.03190299490116622, 0.06009700509883378,
    0.02824982994500573, -0.06375017005499427, 0.00628743512039007,
    -0.00771256487960993, -0.0612685788947333, -0.0172685788947333,
    0.0367314211052667, 0.0167314211052667, 0.0497314211052667,
    -0.0252905938017210, -0.0192905938017210, -0.0412905938017210,
    -0.02992068718566469, -0.00392068718566469, 0.03107931281433531,
    0.07207931281433531, 0.01609413774636004
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic6_new(
    x, y, w, c(0, 1, -1, 0, 1, 1), 10000,
    c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic6_new(
    x, y, w, c(1, 3, -1, 0, 1, 1), 10000,
    c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.027858205619881242

  fitted_values <- c(
    rep(0.91990299490116622, 3), rep(0.91975017005499427, 2),
    rep(0.89071256487960993, 2), rep(0.5492685788947333, 5),
    rep(0.2842905938017210, 3), rep(0.14692068718566469, 4),
    0.07590586225363996
  )

  residuals <- c(
    0.00809700509883378, -0.03190299490116622, 0.06009700509883378,
    0.02824982994500573, -0.06375017005499427, 0.00628743512039007,
    -0.00771256487960993, -0.0612685788947333, -0.0172685788947333,
    0.0367314211052667, 0.0167314211052667, 0.0497314211052667,
    -0.0252905938017210, -0.0192905938017210, -0.0412905938017210,
    -0.02992068718566469, -0.00392068718566469, 0.03107931281433531,
    0.07207931281433531, 0.01609413774636004
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0, 1, -3, -2, 1, 1), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0, 1, -1.2, -0.3, 2, 1.5), 10000,
    c(0, 1, -3, -2, 1, 1), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, 1, 20, 3), 10000,
    c(0, 1, -3, -2, 1, 1), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

context("6-parameter logistic - weighted fit")

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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.9375852722593340, 3), rep(0.9351351725701264, 1),
    rep(0.8859559720192794, 2), rep(0.5516453155354000, 4),
    rep(0.26479534106274689, 3), rep(0.17495286982311317, 3),
    0.14985594300774279
  )

  residuals <- c(
    -0.0095852722593340, -0.0495852722593340, 0.0424147277406660,
    0.0128648274298736, 0.0110440279807206, -0.0029559720192794,
    -0.0636453155354000, -0.0196453155354000, 0.0143546844646000,
    0.0473546844646000, -0.00579534106274689, 0.00020465893725311,
    -0.02179534106274689, -0.03195286982311317, 0.00304713017688683,
    0.04404713017688683, -0.05785594300774279
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  object <- logistic6_new(x, y, w, c(0, 1, -1, 0, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
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

  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.014141550871844299

  fitted_values <- c(
    rep(0.9375852722593340, 3), rep(0.9351351725701265, 1),
    rep(0.8859559720192794, 2), rep(0.5516453155353999, 4),
    rep(0.26479534106274683, 3), rep(0.17495286982311320, 3),
    0.14985594300774288
  )

  residuals <- c(
    -0.0095852722593340, -0.0495852722593340, 0.0424147277406660,
    0.0128648274298735, 0.0110440279807206, -0.0029559720192794,
    -0.0636453155353999, -0.0196453155353999, 0.0143546844646001,
    0.0473546844646001, -0.00579534106274683, 0.00020465893725317,
    -0.02179534106274683, -0.03195286982311320, 0.00304713017688680,
    0.04404713017688680, -0.05785594300774288
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0, 0.5, -2, -1, 0.5, 0.25), c(0.5, 1, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0.3, 0.6, -1.3, -0.3, 0.7, 0.3), 10000,
    c(0, 0.5, -2, -1, 0.5, 0.25), c(0.5, 1, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, -2, 7, 2), 10000,
    c(0, 0.5, -2, -1, 0.5, 0.25), c(0.5, 1, 0, 1, 5, 1)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
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
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.016942492408652008

  fitted_values <- c(
    rep(0.93741397025414293, 3), rep(0.93679545474108154, 1),
    rep(0.88551343724569236, 2), rep(0.5471760024604304, 4),
    rep(0.2955709497045754, 3), rep(0.15934552408484025, 3),
    0.08588003931268542
  )

  residuals <- c(
    -0.00941397025414293, -0.04941397025414293, 0.04258602974585707,
    0.01120454525891846, 0.01148656275430764, -0.00251343724569236,
    -0.0591760024604304, -0.0151760024604304, 0.0188239975395696,
    0.0518239975395696, -0.0365709497045754, -0.0305709497045754,
    -0.0525709497045754, -0.01634552408484025, 0.01865447591515975,
    0.05965447591515975, 0.00611996068731458
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic6_new(
    x, y, w, c(0, 1, -1, 0, 1, 1), 10000,
    c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic6_new(
    x, y, w, c(1, 3, -1, 0, 1, 1), 10000,
    c(0, 1, rep(-Inf, 4)), c(0, 1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.016942492408652008

  fitted_values <- c(
    rep(0.93741397025414293, 3), rep(0.93679545474108154, 1),
    rep(0.88551343724569236, 2), rep(0.5471760024604304, 4),
    rep(0.2955709497045754, 3), rep(0.15934552408484025, 3),
    0.08588003931268542
  )

  residuals <- c(
    -0.00941397025414293, -0.04941397025414293, 0.04258602974585707,
    0.01120454525891846, 0.01148656275430764, -0.00251343724569236,
    -0.0591760024604304, -0.0151760024604304, 0.0188239975395696,
    0.0518239975395696, -0.0365709497045754, -0.0305709497045754,
    -0.0525709497045754, -0.01634552408484025, 0.01865447591515975,
    0.05965447591515975, 0.00611996068731458
  )

  object <- logistic6_new(
    x, y, w, NULL, 10000, c(0, 1, -3, -2, 0, 0), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic6_new(
    x, y, w, c(0, 1, -1.2, -0.3, 1, 1), 10000,
    c(0, 1, -3, -2, 0, 0), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic6_new(
    x, y, w, c(-1, 2, 1, 1, 20, 3), 10000,
    c(0, 1, -3, -2, 0, 0), c(0, 1, -1, 0, 10, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w)
})

context("6-parameter logistic - general functions")

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
    nu = 1 / 2,
    xi = 3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      4153.2629127765089, 662.45973982797511, 20.966690123749200,
      99.350069363821825, 1969.8255704143403, -1526.0398397171145,
      73336.586414524941,
      # beta
      662.45973982797511, 503.81760756754085, -11.539313577135948,
      -27.254756533043237, -296.44108643671633, 212.74025231384907,
      23986.080835069310,
      # eta
      20.966690123749200, -11.539313577135948, -23.808517105009236,
      -36.920624017648987, -53.963629240443870, -5.1670031045368685,
      188.85343485737145,
      # phi
      99.350069363821825, -27.254756533043237, -36.920624017648987,
      -19.302185112689962, -84.861959254447253, 24.818059120326258,
      1462.0499201047423,
      # nu
      1969.8255704143403, -296.44108643671633, -53.963629240443870,
      -84.861959254447253, 2851.0352393691240, -1063.5701751813548,
      34399.705218858010,
      # xi
      -1526.0398397171145, 212.74025231384907, -5.1670031045368685,
      24.818059120326258, -1063.5701751813548, -685.94762601606976,
      -27016.689384177895,
      # sigma
      73336.586414524941, 23986.080835069310, 188.85343485737145,
      1462.0499201047423, 34399.705218858010, -27016.689384177895,
      1.3696180637620797e+06
    ),
    nrow = 7,
    ncol = 7
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "nu", "xi", "sigma"
  )

  object <- logistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 7 * 7)
  expect_equal(fim, true_value)
})

context("6-parameter logistic - drda fit")

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
      start = c(0, -1, -1, 0, 1, 1)
    ),
    "parameter 'beta' cannot be smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 0, 0, 1, 1)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 0, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, -1, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 1, 0)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 1, -1)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )
})

context("6-parameter logistic - Area under and above the curve")

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

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(nauc(result), 0.53873944885601383)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.48525519463582857)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.52493778668627164)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.47542532439304761
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

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(naac(result), 1 - 0.53873944885601383)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.48525519463582857)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.52493778668627164)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.47542532439304761
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

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(nauc(result), 0.53465218185214561)
  expect_equal(nauc(result, xlim = c(-1, 2)), 0.56806493426957927)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.50308844614993796)
  expect_equal(
    nauc(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 0.61344155711596544
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

  result <- drda(y ~ x, mean_function = "logistic6")

  expect_equal(naac(result), 1 - 0.53465218185214561)
  expect_equal(naac(result, xlim = c(-1, 2)), 1 - 0.56806493426957927)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1 - 0.50308844614993796)
  expect_equal(
    naac(result, xlim = c(-1, 2), ylim = c(0.2, 0.8)), 1 - 0.61344155711596544
  )
})
