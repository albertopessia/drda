test_that("Constructor (decreasing)", {
  x <- ltd$D$x
  y <- ltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- ltd$stats_1

  start <- c(1, 1)

  lower_bound <- c(0.5, 1)
  upper_bound <- c(2, 5)

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, log(start[1]), start[2]))
  expect_equal(object$lower_bound, c(log(lower_bound[1]), lower_bound[2]))
  expect_equal(object$upper_bound, c(log(upper_bound[1]), upper_bound[2]))

  w <- ltd$D$w
  stats <- ltd$stats_2

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, log(start[1]), start[2]))
  expect_equal(object$lower_bound, c(log(lower_bound[1]), lower_bound[2]))
  expect_equal(object$upper_bound, c(log(upper_bound[1]), upper_bound[2]))
})

test_that("Constructor (increasing)", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- ltd$stats_1_i

  start <- c(1, 1)

  lower_bound <- c(0.5, 1)
  upper_bound <- c(2, 5)

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, log(start[1]), start[2]))
  expect_equal(object$lower_bound, c(log(lower_bound[1]), lower_bound[2]))
  expect_equal(object$upper_bound, c(log(upper_bound[1]), upper_bound[2]))

  w <- ltd$D$w
  stats <- ltd$stats_2_i

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, NA_real_, NA_real_))
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- logistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, log(start[1]), start[2]))
  expect_equal(object$lower_bound, c(log(lower_bound[1]), lower_bound[2]))
  expect_equal(object$upper_bound, c(log(upper_bound[1]), upper_bound[2]))
})

test_that("Constructor: errors", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w
  max_iter <- 10000

  expect_error(
    logistic2_new(x, y, w, 1, max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, c(1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, c(0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic2_new(x, y, w, c(-1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, -Inf, Inf),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, -Inf, rep(Inf, 2)),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), Inf),
    "'upper_bound' must be of length 2"
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(0, Inf)),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(-1, Inf)),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value (decreasing)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d
  theta[2] <- -0.5

  m <- length(x)

  true_value <- c(
    0.99993227585038023, 0.93086157965665318, 0.54983399731247791,
    0.45016600268752209, 0.35434369377420455, 0.26894142136999512,
    0.0044962731609411802, 0.000030431556900565341
  )

  value <- logistic2_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic2")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic2_fit")

  value <- fn(object, object$stats[, 1], c(1, -1, theta[3:4]))

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Function value (increasing)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_i
  theta[2] <- 0.5

  m <- length(x)

  true_value <- c(
    0.000067724149619770208, 0.069138420343346818, 0.45016600268752209,
    0.54983399731247791, 0.64565630622579545, 0.73105857863000488,
    0.99550372683905882, 0.99996956844309943
  )

  value <- logistic2_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic2")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic2_fit")

  value <- fn(object, object$stats[, 1], c(0, 1, theta[3:4]))

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # eta
      0.0065010780536955348, 1.6733157785701111, 0.49503314542371989,
      -0.49503314542371989, -1.3727054427399438, -1.9661193324148185,
      -0.24170706118458253, -0.0031647856053746346,
      # phi
      6.7719563059328487e-06, 0.0064358299175773505, 0.024751657271185994,
      0.024751657271185994, 0.022878424045665730, 0.019661193324148185,
      0.00044760566886033802, 3.0430630820909948e-06
    ),
    nrow = m,
    ncol = 2
  )

  G <- logistic2_gradient(x, theta, -1)

  expect_type(G, "double")
  expect_length(G, m * 2)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_hessian <- array(
    c(
      # (eta, eta)
      -0.62401895939807407, -37.490308940194270, -0.098677921754532571,
      0.098677921754532571, 2.3993184519065446, 9.0857747672948409,
      12.934808958989922, 0.32911767061349057,
      # (eta, phi)
      -0.00058230018631366534, -0.079835196748050610, 0.24258267662413332,
      0.24258267662413332, 0.18879559959154822, 0.10575418556853344,
      -0.019477293235452032, -0.00028602866784590791,
      # (phi, eta)
      -0.00058230018631366534, -0.079835196748050610, 0.24258267662413332,
      0.24258267662413332, 0.18879559959154822, 0.10575418556853344,
      -0.019477293235452032, -0.00028602866784590791,
      # (phi, phi)
      -6.7710390559686857e-07, -0.00055459036893778505, -0.00024669480438633143,
      0.00024669480438633143, 0.00066647734775181793, 0.00090857747672948409,
      0.000044358055414917430, 3.0428778717963255e-07
    ),
    dim = c(m, 2, 2)
  )

  H <- logistic2_hessian(x, theta, -1)

  expect_type(H, "double")
  expect_length(H, m * 2 * 2)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # eta
      0.0065010780536955348, 1.6733157785701111, 0.49503314542371989,
      -0.49503314542371989, -1.3727054427399438, -1.9661193324148185,
      -0.24170706118458253, -0.0031647856053746346,
      # phi
      6.7719563059328487e-06, 0.0064358299175773505, 0.024751657271185994,
      0.024751657271185994, 0.022878424045665730, 0.019661193324148185,
      0.00044760566886033802, 3.0430630820909948e-06
    ),
    nrow = m,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (eta, eta)
      -0.62401895939807407, -37.490308940194270, -0.098677921754532571,
      0.098677921754532571, 2.3993184519065446, 9.0857747672948409,
      12.934808958989922, 0.32911767061349057,
      # (eta, phi)
      -0.00058230018631366534, -0.079835196748050610, 0.24258267662413332,
      0.24258267662413332, 0.18879559959154822, 0.10575418556853344,
      -0.019477293235452032, -0.00028602866784590791,
      # (phi, eta)
      -0.00058230018631366534, -0.079835196748050610, 0.24258267662413332,
      0.24258267662413332, 0.18879559959154822, 0.10575418556853344,
      -0.019477293235452032, -0.00028602866784590791,
      # (phi, phi)
      -6.7710390559686857e-07, -0.00055459036893778505, -0.00024669480438633143,
      0.00024669480438633143, 0.00066647734775181793, 0.00090857747672948409,
      0.000044358055414917430, 3.0428778717963255e-07
    ),
    dim = c(m, 2, 2)
  )

  gh <- logistic2_gradient_hessian(x, theta, -1)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 2)
  expect_length(gh$H, m * 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # log_eta
      0.00065010780536955348, 0.16733157785701111, 0.049503314542371989,
      -0.049503314542371989, -0.13727054427399438, -0.19661193324148185,
      -0.024170706118458253, -0.00031647856053746346,
      # phi
      6.7719563059328487e-06, 0.0064358299175773505, 0.024751657271185994,
      0.024751657271185994, 0.022878424045665730, 0.019661193324148185,
      0.00044760566886033802, 3.0430630820909948e-06
    ),
    nrow = m,
    ncol = 2
  )

  G <- logistic2_gradient_2(x, theta, -1)

  expect_type(G, "double")
  expect_length(G, m * 2)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_hessian <- array(
    c(
      # (log_eta, log_eta)
      -0.0055900817886111872, -0.20757151154493158, 0.048516535324826663,
      -0.048516535324826663, -0.11327735975492893, -0.10575418556853344,
      0.10517738347144097, 0.0029746981455974422,
      # (log_eta, phi)
      -0.000058230018631366534, -0.0079835196748050610, 0.024258267662413332,
      0.024258267662413332, 0.018879559959154822, 0.010575418556853344,
      -0.0019477293235452032, -0.000028602866784590791,
      # (phi, log_eta)
      -0.000058230018631366534, -0.0079835196748050610, 0.024258267662413332,
      0.024258267662413332, 0.018879559959154822, 0.010575418556853344,
      -0.0019477293235452032, -0.000028602866784590791,
      # (phi, phi)
      -6.7710390559686857e-07, -0.00055459036893778505, -0.00024669480438633143,
      0.00024669480438633143, 0.00066647734775181793, 0.00090857747672948409,
      0.000044358055414917430, 3.0428778717963255e-07
    ),
    dim = c(m, 2, 2)
  )

  H <- logistic2_hessian_2(x, theta, -1)

  expect_type(H, "double")
  expect_length(H, m * 2 * 2)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # log_eta
      0.00065010780536955348, 0.16733157785701111, 0.049503314542371989,
      -0.049503314542371989, -0.13727054427399438, -0.19661193324148185,
      -0.024170706118458253, -0.00031647856053746346,
      # phi
      6.7719563059328487e-06, 0.0064358299175773505, 0.024751657271185994,
      0.024751657271185994, 0.022878424045665730, 0.019661193324148185,
      0.00044760566886033802, 3.0430630820909948e-06
    ),
    nrow = m,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (log_eta, log_eta)
      -0.0055900817886111872, -0.20757151154493158, 0.048516535324826663,
      -0.048516535324826663, -0.11327735975492893, -0.10575418556853344,
      0.10517738347144097, 0.0029746981455974422,
      # (log_eta, phi)
      -0.000058230018631366534, -0.0079835196748050610, 0.024258267662413332,
      0.024258267662413332, 0.018879559959154822, 0.010575418556853344,
      -0.0019477293235452032, -0.000028602866784590791,
      # (phi, log_eta)
      -0.000058230018631366534, -0.0079835196748050610, 0.024258267662413332,
      0.024258267662413332, 0.018879559959154822, 0.010575418556853344,
      -0.0019477293235452032, -0.000028602866784590791,
      # (phi, phi)
      -6.7710390559686857e-07, -0.00055459036893778505, -0.00024669480438633143,
      0.00024669480438633143, 0.00066647734775181793, 0.00090857747672948409,
      0.000044358055414917430, 3.0428778717963255e-07
    ),
    dim = c(m, 2, 2)
  )

  gh <- logistic2_gradient_hessian_2(x, theta, -1)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 2)
  expect_length(gh$H, m * 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(
    list(stats = ltd$stats_1, start = c(1, -1, NA_real_, NA_real_)),
    class = "logistic2"
  )

  gh <- gradient_hessian(object, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 2)
  expect_length(gh$H, m * 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Value of the RSS", {
  theta <- c(log(ltd$theta_2_d[3]), ltd$theta_2_d[4])

  true_value <- 0.16502893044345798

  object <- structure(
    list(
      stats = ltd$stats_1,
      m = nrow(ltd$stats_1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "logistic2"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(NA, theta[2])
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta[1])

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  theta <- c(log(ltd$theta_2_d[3]), ltd$theta_2_d[4])

  true_gradient <- c(0.025653880019264530, -0.0036804593089360313)

  true_hessian <- matrix(
    c(
      # log_eta
      0.097618662721905987, -0.021398439978182300,
      # phi
      -0.021398439978182300, 0.0057166893653613790
    ),
    nrow = 2,
    ncol = 2
  )

  object <- structure(
    list(
      stats = ltd$stats_1,
      m = nrow(ltd$stats_1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "logistic2"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta)

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 2)
  expect_length(gh$H, 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  known_param <- c(NA, theta[2])
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta[1])

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 1)
  expect_length(gh$H, 1)

  expect_equal(gh$G, true_gradient[1])
  expect_equal(gh$H, true_hessian[1, 1, drop = FALSE])
})

test_that("mle_asy", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(-2.5338955025359479, -3.9372918905117386)

  true_value <- c(-2.5338955025359479, -3.9372918905117386)

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.20727384125777396

  theta <- c(
    alpha = 1, delta = -1, eta = exp(-2.5338955025359479),
    phi = -3.9372918905117386
  )

  fitted_values <- rep(
    c(
      0.99951092709497446, 0.88775959713220388, 0.54082751285442484,
      0.46164480559169129, 0.38435118260736283, 0.31248712794063619,
      0.013655180945348409, 0.00026188252852923363
    ),
    k
  )

  residuals <- c(
    -0.14681092709497446, -0.24241092709497446, -0.062510927094974457,
    -0.060259597132203877, -0.10995959713220388, -0.0031595971322038767,
    0.015272487145575160, 0.16397248714557516, 0.083655194408308708,
    0.024355194408308708, 0.073455194408308708, -0.021744805591691292,
    -0.029051182607362831, -0.071151182607362831, -0.034451182607362831,
    -0.17038712794063619, 0.0031448190546515905, 0.12494481905465159,
    0.12243811747147077
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic2_new(x, y, w, c(1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.20727384125777396

  theta <- c(
    alpha = 1, delta = -1, eta = exp(-2.5338955025359479),
    phi = -3.9372918905117386
  )

  fitted_values <- rep(
    c(
      0.99951092709497446, 0.88775959713220388, 0.54082751285442484,
      0.46164480559169129, 0.38435118260736283, 0.31248712794063619,
      0.013655180945348409, 0.00026188252852923363
    ),
    k
  )

  residuals <- c(
    -0.14681092709497446, -0.24241092709497446, -0.062510927094974457,
    -0.060259597132203877, -0.10995959713220388, -0.0031595971322038767,
    0.015272487145575160, 0.16397248714557516, 0.083655194408308708,
    0.024355194408308708, 0.073455194408308708, -0.021744805591691292,
    -0.029051182607362831, -0.071151182607362831, -0.034451182607362831,
    -0.17038712794063619, 0.0031448190546515905, 0.12494481905465159,
    0.12243811747147077
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(
    x, y, w, c(1, 0), max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(
    x, y, w, c(7, -8), max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.51514867830805084

  theta <- c(
    alpha = 1, delta = -1, eta = 0.5, phi = -1.8048595456927261
  )

  fitted_values <- rep(
    c(
      1, 0.99999924577188988, 0.89066679334583969, 0.52437322386786472,
      0.12983371911440877, 0.019793109833897495, 5.6327137560357419e-12,
      7.8226812451675656e-23
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.17249924577188988, -0.22219924577188988,
    -0.11539924577188988, -0.33456679334583969, -0.18586679334583969,
    0.020926776132135283, -0.038373223867864717, 0.010726776132135283,
    -0.084473223867864717, 0.22546628088559123, 0.18336628088559123,
    0.22006628088559123, 0.12230689016610250, 0.016799999994367286,
    0.13859999999436729, 0.1227
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic2_new(
    x, y, w, c(0.5, 0), max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic2_new(
    x, y, w, c(1, 0), max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: equalities and inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.51514867830805084

  theta <- c(
    alpha = 1, delta = -1, eta = 0.5, phi = -1.8048595456927261
  )

  fitted_values <- rep(
    c(
      1, 0.99999924577188988, 0.89066679334583969, 0.52437322386786472,
      0.12983371911440877, 0.019793109833897495, 5.6327137560357419e-12,
      7.8226812451675656e-23
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.17249924577188988, -0.22219924577188988,
    -0.11539924577188988, -0.33456679334583969, -0.18586679334583969,
    0.020926776132135283, -0.038373223867864717, 0.010726776132135283,
    -0.084473223867864717, 0.22546628088559123, 0.18336628088559123,
    0.22006628088559123, 0.12230689016610250, 0.016799999994367286,
    0.13859999999436729, 0.1227
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(
    x, y, w, c(0.5, 0), max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(
    x, y, w, c(8, -8), max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit (weighted)", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.14802210304519822

  theta <- c(
    alpha = 1, delta = -1, eta = exp(-1.7118497175241647),
    phi = -2.0657104870916123
  )

  fitted_values <- rep(
    c(
      0.99999997903117914, 0.99358716675353491, 0.67045935646748149,
      0.49703433073752424, 0.32431997533986986, 0.18906221049072057,
      0.000082758300470677611, 9.9461424961556887e-09
    ),
    k
  )

  residuals <- c(
    -0.14729997903117914, -0.24289997903117914, -0.062999979031179143,
    -0.16608716675353491, -0.21578716675353491, -0.10898716675353491,
    -0.11435935646748149, 0.034340643532518513, 0.048265669262475756,
    -0.011034330737524244, 0.038065669262475756, -0.057134330737524244,
    0.030980024660130145, -0.011119975339869855, 0.025580024660130145,
    -0.046962210490720567, 0.016717241699529322, 0.13851724169952932,
    0.12269999005385750
  )

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic2_new(x, y, w, c(1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.14802210304519822

  theta <- c(
    alpha = 1, delta = -1, eta = exp(-1.7118497175241647),
    phi = -2.0657104870916123
  )

  fitted_values <- rep(
    c(
      0.99999997903117914, 0.99358716675353491, 0.67045935646748149,
      0.49703433073752424, 0.32431997533986986, 0.18906221049072057,
      0.000082758300470677611, 9.9461424961556887e-09
    ),
    k
  )

  residuals <- c(
    -0.14729997903117914, -0.24289997903117914, -0.062999979031179143,
    -0.16608716675353491, -0.21578716675353491, -0.10898716675353491,
    -0.11435935646748149, 0.034340643532518513, 0.048265669262475756,
    -0.011034330737524244, 0.038065669262475756, -0.057134330737524244,
    0.030980024660130145, -0.011119975339869855, 0.025580024660130145,
    -0.046962210490720567, 0.016717241699529322, 0.13851724169952932,
    0.12269999005385750
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(
    x, y, w, c(1, 0), max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(
    x, y, w, c(7, -8), max_iter,
    c(0.01, -5), c(2, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.33123657223639616

  theta <- c(
    alpha = 1, delta = -1, eta = 0.5, phi = -1.6708752737442977
  )

  fitted_values <- rep(
    c(
      1, 0.99999929464390932, 0.89702175505783073, 0.54104799827715088,
      0.13759176925344944, 0.021135530999337022, 6.0229879840594822e-12,
      8.3646919022474704e-23
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.17249929464390932, -0.22219929464390932,
    -0.11539929464390932, -0.34092175505783073, -0.19222175505783073,
    0.0042520017228491238, -0.055047998277150876, -0.0059479982771508762,
    -0.10114799827715088, 0.21770823074655056, 0.17560823074655056,
    0.21230823074655056, 0.12096446900066298, 0.016799999993977012,
    0.13859999999397701, 0.1227
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic2_new(
    x, y, w, c(0.5, 0), max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic2_new(
    x, y, w, c(1, 0), max_iter,
    c(0.5, -Inf), c(0.5, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fit_constrained (weighted): equalities and inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.33123657223639616

  theta <- c(
    alpha = 1, delta = -1, eta = 0.5, phi = -1.6708752737442977
  )

  fitted_values <- rep(
    c(
      1, 0.99999929464390932, 0.89702175505783073, 0.54104799827715088,
      0.13759176925344944, 0.021135530999337022, 6.0229879840594822e-12,
      8.3646919022474704e-23
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.17249929464390932, -0.22219929464390932,
    -0.11539929464390932, -0.34092175505783073, -0.19222175505783073,
    0.0042520017228491238, -0.055047998277150876, -0.0059479982771508762,
    -0.10114799827715088, 0.21770823074655056, 0.17560823074655056,
    0.21230823074655056, 0.12096446900066298, 0.016799999993977012,
    0.13859999999397701, 0.1227
  )

  object <- logistic2_new(
    x, y, w, NULL, max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic2_new(
    x, y, w, c(0.5, 0), max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic2_new(
    x, y, w, c(8, -8), max_iter,
    c(0.5, -5), c(0.5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic2_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  max_iter <- 10000

  theta <- ltd$theta_2_d
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- ltd$sigma

  true_value <- matrix(c(
      # eta
      3671.8976437329046, -79.729113965633577, -2000.7114186570950,
      # phi
      -79.729113965633577, 1.9250613126586459, 34.261360313645943,
      # sigma
      -2000.7114186570950, 34.261360313645943, 65475.726403465305
    ),
    nrow = 3,
    ncol = 3
  )

  rownames(true_value) <- colnames(true_value) <- c("eta", "phi", "sigma")

  object <- logistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 3 * 3)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- ltd$D$x
  y <- ltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c("c", "d")
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
      upper_bound = c(-1, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      lower_bound = c(Inf, -Inf),
      upper_bound = c(Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
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

test_that("drda: 'upper_bound' argument errors", {
  x <- ltd$D$x
  y <- ltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      upper_bound = c("c", "d")
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
      upper_bound = c(-Inf, Inf)
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
  x <- ltd$D$x
  y <- ltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c("c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(1, Inf)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(-Inf, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = rep(1, 3)
    ),
    "'start' must be of length 2"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(-1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic2",
      start = c(0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic2")

  expect_equal(nauc(result), 0.42629872039669964)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.40878964077573577)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.39671447564542112)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 0.032549163262074700
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic2")

  expect_equal(naac(result), 1 - 0.42629872039669964)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.40878964077573577)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.39671447564542112)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 1 - 0.032549163262074700
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic2")

  expect_equal(nauc(result), 0.64408118428934661)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.64801389875761971)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.83138849342085760)
  expect_equal(nauc(result, xlim = c(-40, -25), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 0.76263081184979702
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic2")

  expect_equal(naac(result), 1 - 0.64408118428934661)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.64801389875761971)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.83138849342085760)
  expect_equal(naac(result, xlim = c(-40, -25), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 1 - 0.76263081184979702
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})
