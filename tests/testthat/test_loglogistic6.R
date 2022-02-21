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

  start <- c(0, 1, 1, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1, 0, 0.5)
  upper_bound <- c(3, 2, 2, 5, 2, 1)

  object <- loglogistic6_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic6"))
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

  object <- loglogistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic6"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5), log(2), 0))

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

  object <- loglogistic6_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic6"))
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

  object <- loglogistic6_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic6"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, 7)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, 0, 0, 0, 0))
  expect_equal(object$lower_bound, c(0, -1, log(0.5), 0, -Inf, log(0.5)))
  expect_equal(object$upper_bound, c(3, 2, log(2), log(5), log(2), 0))
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
    loglogistic6_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 6"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 0, 1, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, -1, 1, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, -1, 1, 1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, c(0, 1, 1, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 5)),
    "'lower_bound' must be of length 6"
  )

  expect_error(
    loglogistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 5), rep(Inf, 6)),
    "'lower_bound' must be of length 6"
  )

  expect_error(
    loglogistic6_new(x, y, w, NULL, max_iter, rep(-Inf, 6), rep(Inf, 5)),
    "'upper_bound' must be of length 6"
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, 1, 0, rep(Inf, 3))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, 1, -1, rep(Inf, 3))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, 1, Inf, 0, Inf, Inf)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, 1, Inf, -1, Inf, Inf)
    ),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, Inf, Inf, 0, Inf)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, Inf, Inf, Inf, -1, Inf)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, rep(Inf, 4), 0)
    ),
    "'upper_bound[6]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic6_new(
      x, y, w, NULL, max_iter, rep(-Inf, 6), c(1, rep(Inf, 4), -1)
    ),
    "'upper_bound[6]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2, 1 / 2, 3 / 2)

  true_value <- c(
    0.040000000000000000, -0.18500000000000000, -0.32864000000000000,
    -0.35030339083878644, -0.35586566082310935, -0.35787516976007243
  )

  value <- loglogistic6_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic6"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic6_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- c(0, 2, 4, 6, 8, 10)
  theta <- c(4 / 100, -9 / 10, 3, 2, 1 / 2, 3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # delta
      0.0, 0.25000000000000000, 0.40960000000000000, 0.43367043426531826,
      0.43985073424789927, 0.44208352195563603,
      # log_eta
      0.0, 0, -0.061325226393988377, -0.031375031815926065,
      -0.017060715026738104, -0.010218436956830015,
      # log_phi
      0.0, 0.33750000000000000, 0.088473600000000000, 0.028558784695520959,
      0.012306704481547441, 0.0063490718578734962,
      # log_nu
      0.0, -0.19941623125197539, -0.29954735502588114, -0.31656733571333840,
      -0.32103009032854108, -0.32265181674542238,
      # log_xi
      0.0, 0.33750000000000000, 0.70778880000000000, 0.77108718677906589,
      0.78762908681903621, 0.79363398223418703
    ),
    nrow = 6,
    ncol = 6
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
      # (alpha, log_xi)
      rep(0, 6),
      # (delta, alpha)
      rep(0, 6),
      # (delta, delta)
      rep(0, 6),
      # (delta, log_eta)
      0, 0, 0.068139140437764864, 0.034861146462140072, 0.018956350029709005,
      0.011353818840922239,
      # (delta, log_phi)
      0, -0.37500000000000000, -0.098304000000000000, -0.031731982995023288,
      -0.013674116090608268, -0.0070545242865261069,
      # (delta, log_nu)
      0, 0.22157359027997265, 0.33283039447320126, 0.35174148412593155,
      0.35670010036504565, 0.35850201860602487,
      # (delta, log_xi)
      0, -0.37500000000000000, -0.78643200000000000, -0.85676354086562876,
      -0.87514342979892912, -0.88181553581576337,
      # (log_eta, alpha)
      rep(0, 6),
      # (log_eta, delta)
      0, 0, 0.068139140437764864, 0.034861146462140072, 0.018956350029709005,
      0.011353818840922239,
      # (log_eta, log_eta)
      0, 0, 0.050894330124602161, 0.068248772281001278, 0.052789899645785409,
      0.038725729760171624,
      # (log_eta, log_phi)
      0, 0.33750000000000000, -0.073424997680129316, -0.062122709699289740,
      -0.038079863214000328, -0.024061648766309447,
      # (log_eta, log_nu)
      0, 0, -0.047378288043531299, -0.025065043407532354,
      -0.013747111444595107, -0.0082593350036119815,
      # (log_eta, log_xi)
      0, 0, 0.17661665201468653, 0.092977228430122362, 0.050916952618658798,
      0.030573780788387678,
      # (log_phi, alpha)
      rep(0, 6),
      # (log_phi, delta)
      0, -0.37500000000000000, -0.098304000000000000, -0.031731982995023288,
      -0.013674116090608268, -0.0070545242865261069,
      # (log_phi, log_eta)
      0, 0.33750000000000000, -0.073424997680129316, -0.062122709699289740,
      -0.038079863214000328, -0.024061648766309447,
      # (log_phi, log_phi)
      0, 0.25312500000000000, 0.23357030400000000, 0.082541853327298381,
      0.036346225670891405, 0.018895243108937347,
      # (log_phi, log_nu)
      0, 0.21474934687796308, 0.068352421206211473, 0.022815185726640359,
      0.0099164447538326376, 0.0051318133739751617,
      # (log_phi, log_xi)
      0, -0.75937500000000000, -0.25480396800000000, -0.084631520500141378,
      -0.036728817520058683, -0.018996558085392775,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, delta)
      0, 0.22157359027997265, 0.33283039447320126, 0.35174148412593155,
      0.35670010036504565, 0.35850201860602487,
      # (log_nu, log_eta)
      0, 0, -0.047378288043531299, -0.025065043407532354,
      -0.013747111444595107, -0.0082593350036119815,
      # (log_nu, log_phi)
      0, 0.21474934687796308, 0.068352421206211473, 0.022815185726640359,
      0.0099164447538326376, 0.0051318133739751617,
      # (log_nu, log_nu)
      0, -0.0054502500224305037, 0.054963253078937568, 0.059689758583680598,
      0.060667186288486732, 0.060995794340251851,
      # (log_nu, log_xi)
      0, -0.12275065312203692, -0.16096943035030822, -0.15507717215977620,
      -0.15297662257374740, -0.15215731048729182,
      # (log_xi, alpha)
      rep(0, 6),
      # (log_xi, delta)
      0, -0.37500000000000000, -0.78643200000000000, -0.85676354086562876,
      -0.87514342979892912, -0.88181553581576337,
      # (log_xi, log_eta)
      0, 0, 0.17661665201468653, 0.092977228430122362, 0.050916952618658798,
      0.030573780788387678,
      # (log_xi, log_phi)
      0, -0.75937500000000000, -0.25480396800000000, -0.084631520500141378,
      -0.036728817520058683, -0.018996558085392775,
      # (log_xi, log_nu)
      0, -0.12275065312203692, -0.16096943035030822, -0.15507717215977620,
      -0.15297662257374740, -0.15215731048729182,
      # (log_xi, log_xi)
      0, -0.42187500000000000, -1.3306429440000000, -1.5139638667247513,
      -1.5630152344647195, -1.5809357784399098
    ),
    dim = c(6, 6, 6)
  )

  object <- structure(
    list(stats = matrix(x, nrow = 6, ncol = 1)),
    class = "loglogistic6"
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

test_that("Value of the RSS", {
  x <- c(0, 2, 4, 6, 8)
  n <- c(3, 3, 2, 4, 3)
  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)
  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  theta <- c(4 / 100, -9 / 10, log(3), log(2), -log(2), log(3 / 2))

  true_value <- 10.4307168300795

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic6"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, log(2), -log(2), NA)
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(c(-9 / 10, log(3), log(3 / 2)))

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

  theta <- c(4 / 100, -9 / 10, log(3), log(2), -log(2), log(3 / 2))

  true_gradient <- c(
    -11.693490545824474, -3.2117663784095942, 0.23558698987525352,
    -1.3121619912391382, 2.3970268343462383, -5.3437921507242235
  )

  true_hessian <- matrix(
    c(
      # alpha
      15, 4.6234339398049709, -0.29933272513189533, 1.3406024522267262,
      -3.4267030176466653, 7.8753136075733722,
      # beta
      4.6234339398049709, 1.8557305074787272, -0.38893914633028841,
      1.8493400824146138, -4.0310740546018891, 9.1474009761426928,
      # eta
      -0.29933272513189533, -0.38893914633028841, -0.31862636550829760,
      -0.68255388120233939, 0.27730331196932592, -0.90995983340678540,
      # phi
      1.3406024522267262, 1.8493400824146138, -0.68255388120233939,
      -1.1909872915340810, -1.1794454814895983, 3.7258079121909984,
      # log_nu
      -3.4267030176466653, -4.0310740546018891, 0.27730331196932592,
      -1.1794454814895983, 0.69072317533181460, -1.0816053233162648,
      # log_xi
      7.8753136075733722, 9.1474009761426928, -0.90995983340678540,
      3.7258079121909984, -1.0816053233162648, 15.223385786326514
    ),
    nrow = 6,
    ncol = 6
  )

  object <- structure(
    list(stats = cbind(x, n, m, v), m = 5),
    class = "loglogistic6"
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

  known_param <- c(4 / 100, NA, NA, log(2), -log(2), NA)
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gradient_hessian <- rss_gh(c(-9 / 10, log(3), log(3 / 2)))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 3)
  expect_length(gradient_hessian$H, 3 * 3)

  expect_equal(gradient_hessian$G, true_gradient[c(2, 3, 6)])
  expect_equal(gradient_hessian$H, true_hessian[c(2, 3, 6), c(2, 3, 6)])
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

  theta <- c(
    0, 1, 1.5969577078711564, 1.8658702431219855, -0.93859189706413948,
    0.42129720705831475
  )

  true_value <- c(
    0.91962793078907831, -2.4062725028857457, 1.5969577078711564,
    1.8658702431219855, -0.93859189706413948, 0.42129720705831475
  )

  object <- loglogistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 6)
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

  # loglogistic6 model is basically unidentifiable: many parameters are
  # associated with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.024098464516913181

  fitted_values <- c(
    rep(0.91962793078907831, 3), rep(0.91961835703802008, 2),
    rep(0.8915156044609054, 2), rep(0.553185547860874, 5),
    rep(0.261160239126468, 3), rep(0.15910736421272, 4),
    0.10001037110104
  )

  residuals <- c(
    0.00837206921092169, -0.03162793078907831, 0.06037206921092169,
    0.02838164296197992, -0.06361835703802008, 0.0054843955390946,
    -0.0085156044609054, -0.065185547860874, -0.021185547860874,
    0.032814452139126, 0.012814452139126, 0.045814452139126, -0.002160239126468,
    0.003839760873532, -0.018160239126468, -0.04210736421272, -0.01610736421272,
    0.01889263578728, 0.05989263578728, -0.00801037110104
  )

  object <- loglogistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic6_new(x, y, w, c(0, 1, 1, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.024098464516913181

  fitted_values <- c(
    rep(0.91962793078907831, 3), rep(0.91961835703802008, 2),
    rep(0.8915156044609054, 2), rep(0.553185547860874, 5),
    rep(0.261160239126468, 3), rep(0.15910736421272, 4),
    0.10001037110104
  )

  residuals <- c(
    0.00837206921092169, -0.03162793078907831, 0.06037206921092169,
    0.02838164296197992, -0.06361835703802008, 0.0054843955390946,
    -0.0085156044609054, -0.065185547860874, -0.021185547860874,
    0.032814452139126, 0.012814452139126, 0.045814452139126, -0.002160239126468,
    0.003839760873532, -0.018160239126468, -0.04210736421272, -0.01610736421272,
    0.01889263578728, 0.05989263578728, -0.00801037110104
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1, 0.2, 0), c(1, 3, 5, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic6_new(
    x, y, w, c(0, 0, 2, 2, 2, 1.5), 10000,
    c(-1, -3, 1, 1, 0.2, 0), c(1, 3, 5, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic6_new(
    x, y, w, c(-2, -5, 0.5, 20, 0.1, 3), 10000,
    c(-1, -3, 1, 1, 0.2, 0), c(1, 3, 5, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.052758980824378086

  fitted_values <- c(
    rep(1, 3), rep(0.98349288740338228, 2), rep(0.8603417073342989, 2),
    rep(0.556567665608271, 5), rep(0.261860304977413, 3),
    rep(0.154175775776116, 4), 0.115415410988357
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03549288740338228, -0.12749288740338228,
    0.0366582926657011, 0.0226582926657011, -0.068567665608271,
    -0.024567665608271, 0.029432334391729, 0.009432334391729, 0.042432334391729,
    -0.002860304977413, 0.003139695022587, -0.018860304977413,
    -0.037175775776116, -0.011175775776116, 0.023824224223884,
    0.064824224223884, -0.023415410988357
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic6_new(
    x, y, w, c(1, -1, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic6_new(
    x, y, w, c(0, 1, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.052758980824378086

  fitted_values <- c(
    rep(1, 3), rep(0.98349288740338228, 2),
    rep(0.8603417073342989, 2), rep(0.556567665608271, 5),
    rep(0.261860304977413, 3), rep(0.154175775776116, 4),
    0.115415410988357
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.03549288740338228, -0.12749288740338228,
    0.0366582926657011, 0.0226582926657011, -0.068567665608271,
    -0.024567665608271, 0.029432334391729, 0.009432334391729, 0.042432334391729,
    -0.002860304977413, 0.003139695022587, -0.018860304977413,
    -0.037175775776116, -0.011175775776116, 0.023824224223884,
    0.064824224223884, -0.023415410988357
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic6_new(
    x, y, w, c(1, -1, 2, 2, 2, 1), 10000,
    c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic6_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10, 3), 10000,
    c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  # loglogistic6 model is basically unidentifiable: many parameters are
  # associated with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.013857640726518951

  fitted_values <- c(
    rep(0.93735105801194196, 3), rep(0.93691762149637618, 2),
    rep(0.8856728142723709, 2), rep(0.551303815084240, 5),
    rep(0.266468096196843, 3), rep(0.175184022895764, 4),
    0.13295499055321
  )

  residuals <- c(
    -0.00935105801194196, -0.04935105801194196, 0.04264894198805804,
    0.01108237850362382, -0.08091762149637618, 0.0113271857276291,
    -0.0026728142723709, -0.063303815084240, -0.019303815084240,
    0.034696184915760, 0.014696184915760, 0.047696184915760,
    -0.007468096196843, -0.001468096196843, -0.023468096196843,
    -0.058184022895764, -0.032184022895764, 0.002815977104236,
    0.043815977104236, -0.04095499055321
  )

  object <- loglogistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic6_new(x, y, w, c(1, -1, 1, 1, 1, 1), 10000, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.013857640726518951

  fitted_values <- c(
    rep(0.93735105801194196, 3), rep(0.93691762149637618, 2),
    rep(0.8856728142723709, 2), rep(0.551303815084240, 5),
    rep(0.266468096196843, 3), rep(0.175184022895764, 4),
    0.13295499055321
  )

  residuals <- c(
    -0.00935105801194196, -0.04935105801194196, 0.04264894198805804,
    0.01108237850362382, -0.08091762149637618, 0.0113271857276291,
    -0.0026728142723709, -0.063303815084240, -0.019303815084240,
    0.034696184915760, 0.014696184915760, 0.047696184915760,
    -0.007468096196843, -0.001468096196843, -0.023468096196843,
    -0.058184022895764, -0.032184022895764, 0.002815977104236,
    0.043815977104236, -0.04095499055321
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(-1, -3, 1, 1, 0.5, 0), c(1, 3, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic6_new(
    x, y, w, c(0.3, 0.6, 2, 3, 0.7, 0.3), 10000,
    c(-1, -3, 1, 1, 0.5, 0), c(1, 3, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic6_new(
    x, y, w, c(0, 0, 2, 2, 2, 1.5), 10000,
    c(-1, -3, 1, 1, 0.5, 0), c(1, 3, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.034604471317882127

  fitted_values <- c(
    rep(1, 3), rep(0.98849829304962502, 2), rep(0.8739525176987930, 2),
    rep(0.552619881220989, 5), rep(0.264042171991040, 3),
    rep(0.175313804732794, 4), 0.146987373442839
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04049829304962502, -0.13249829304962502,
    0.0230474823012070, 0.0090474823012070, -0.064619881220989,
    -0.020619881220989, 0.033380118779011, 0.013380118779011, 0.046380118779011,
    -0.005042171991040, 0.000957828008960, -0.021042171991040,
    -0.058313804732794, -0.032313804732794, 0.002686195267206,
    0.043686195267206, -0.054987373442839
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic6_new(
    x, y, w, c(1, -1, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic6_new(
    x, y, w, c(0, 1, 1, 1, 1, 1), 10000,
    c(1, -1, rep(-Inf, 4)), c(1, -1, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.034604471317882127

  fitted_values <- c(
    rep(1, 3), rep(0.98849829304962502, 2), rep(0.8739525176987930, 2),
    rep(0.552619881220989, 5), rep(0.264042171991040, 3),
    rep(0.175313804732794, 4), 0.146987373442839
  )

  residuals <- c(
    -9 / 125, -14 / 125, -1 / 50, -0.04049829304962502, -0.13249829304962502,
    0.0230474823012070, 0.0090474823012070, -0.064619881220989,
    -0.020619881220989, 0.033380118779011, 0.013380118779011, 0.046380118779011,
    -0.005042171991040, 0.000957828008960, -0.021042171991040,
    -0.058313804732794, -0.032313804732794, 0.002686195267206,
    0.043686195267206, -0.054987373442839
  )

  object <- loglogistic6_new(
    x, y, w, NULL, 10000, c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic6_new(
    x, y, w, c(1, -1, 2, 2, 2, 1), 10000,
    c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic6_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10, 3), 10000,
    c(1, -1, 1, 1, 0.2, 0), c(1, -1, 10, 10, 5, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic6_fit"))
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
    alpha = 4 / 100, delta = -9 / 10, eta = 3, phi = 2, nu = 1 / 2, xi = 3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5982, 1737.6411673917355, -32.563864073319139, 113.91505984214579,
      -2550.1562213589091, 2034.5404853846844, 199997.69002610248,
      # beta
      1737.6411673917355, 737.62820748106778, -47.042505808076308,
      172.01641282586631, -3330.1732683364612, 2644.1178866953831,
      55062.473302554603,
      # eta
      -32.563864073319139, -47.042505808076308, -27.927845904655816,
      0.64795164899965898, 67.445838938087329, -74.731122594840330,
      -1191.1302018075331,
      # phi
      113.91505984214579, 172.01641282586631, 0.64795164899965898,
      -21.195434262100594, -225.31240521880504, 242.82809406814098,
      4857.8770209102063,
      # nu
      -2550.1562213589091, -3330.1732683364612, 67.445838938087329,
      -225.31240521880504, -3153.3560290660310, -592.32067708656387,
      -80996.807471355637,
      # xi
      2034.5404853846844, 2644.1178866953831, -74.731122594840330,
      242.82809406814098, -592.32067708656387, 4130.2830584929690,
      63915.911509327654,
      # sigma
      199997.69002610248, 55062.473302554603, -1191.1302018075331,
      4857.8770209102063, -80996.807471355637, 63915.911509327654,
      5309646.9785863507
    ),
    nrow = 7,
    ncol = 7
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "xi", "sigma"
  )

  object <- loglogistic6_new(x, y, w, NULL, 10000, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 7 * 7)
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
      y ~ x, mean_function = "loglogistic6",
      lower_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = matrix(-Inf, nrow = 6, ncol = 2),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
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
      y ~ x, mean_function = "loglogistic6",
      upper_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = rep(-Inf, 6),
      upper_bound = matrix(Inf, nrow = 6, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
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
      y ~ x, mean_function = "loglogistic6",
      start = c("a", "b", "c", "d", "e", "f")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, Inf, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(-Inf, 1, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(1, 1, 1, 1, 1, 1, 1)
    ),
    "'start' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, -1, 1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 0, 1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, -1, 1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, 0, 1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, 1, -1, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, 1, 0, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, 1, 1, -1)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic6",
      start = c(0, 1, 1, 1, 1, 0)
    ),
    "parameter 'xi' cannot be negative nor zero"
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

  result <- drda(y ~ x, mean_function = "loglogistic6")

  expect_equal(nauc(result), 0.63385859186507327)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.91962722237822136)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.64213201737475438)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.46158084161668578
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

  result <- drda(y ~ x, mean_function = "loglogistic6")

  expect_equal(naac(result), 1.0 - 0.63385859186507327)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.91962722237822136)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.64213201737475438)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.46158084161668578
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

  result <- drda(y ~ x, mean_function = "loglogistic6")

  expect_equal(nauc(result), 0.44217263239185015)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.14940394886168675)
  expect_equal(nauc(result, ylim = c(0.2, 0.8)), 0.40787801585973162)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 0.62476097129796916
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

  result <- drda(y ~ x, mean_function = "loglogistic6")

  expect_equal(naac(result), 1.0 - 0.44217263239185015)
  expect_equal(naac(result, xlim = c(0, 2)), 1.0 - 0.14940394886168675)
  expect_equal(naac(result, ylim = c(0.2, 0.8)), 1.0 - 0.40787801585973162)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.2, 0.8)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.2, 0.8)), 1.0 - 0.62476097129796916
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.2, 0.8)), 0.0)
})
