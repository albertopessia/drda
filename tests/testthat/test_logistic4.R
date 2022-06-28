test_that("Constructor", {
  x <- ltd$D$x
  y <- ltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- ltd$stats_1

  start <- c(0, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1)
  upper_bound <- c(3, 2, 2, 5)

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_false(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_null(object$start)
  expect_null(object$lower_bound)
  expect_null(object$upper_bound)

  object <- logistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  i <- c(1, 2)

  expect_true(inherits(object, "logistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(start[i], log(start[3]), start[4]))
  expect_equal(
    object$lower_bound, c(lower_bound[i], log(lower_bound[3]), lower_bound[4])
  )
  expect_equal(
    object$upper_bound, c(upper_bound[i], log(upper_bound[3]), upper_bound[4])
  )

  w <- ltd$D$w
  stats <- ltd$stats_2

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic4"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
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
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(start[i], log(start[3]), start[4]))
  expect_equal(
    object$lower_bound, c(lower_bound[i], log(lower_bound[3]), lower_bound[4])
  )
  expect_equal(
    object$upper_bound, c(upper_bound[i], log(upper_bound[3]), upper_bound[4])
  )
})

test_that("Constructor: errors", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w
  max_iter <- 10000

  expect_error(
    logistic4_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic4_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_value <- c(
    0.89995259309526616, 0.85160310575965723, 0.58488379811873454,
    0.51511620188126546, 0.44804058564194318, 0.38825899495899658,
    0.20314739121265883, 0.20002130208983040
  )

  value <- logistic4_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic4")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic4_fit")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.000067724149619770208, 0.069138420343346818, 0.45016600268752209,
      0.54983399731247791, 0.64565630622579545, 0.73105857863000488,
      0.99550372683905882, 0.99996956844309943,
      # eta
      0.0045507546375868744, 1.1713210449990778, 0.34652320179660392,
      -0.34652320179660392, -0.96089380991796066, -1.3762835326903730,
      -0.16919494282920777, -0.0022153499237622442,
      # phi
      4.7403694141529941e-06, 0.0045050809423041453, 0.017326160089830196,
      0.017326160089830196, 0.016014896831966011, 0.013762835326903730,
      0.00031332396820223662, 2.1301441574636964e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- logistic4_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, m),
      # (alpha, delta)
      rep(0, m),
      # (alpha, eta)
      rep(0, m),
      # (alpha, phi)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, eta)
      -0.0065010780536955348, -1.6733157785701111, -0.49503314542371989,
      0.49503314542371989, 1.3727054427399438, 1.9661193324148185,
      0.24170706118458253, 0.0031647856053746346,
      # (delta, phi)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      -0.0065010780536955348, -1.6733157785701111, -0.49503314542371989,
      0.49503314542371989, 1.3727054427399438, 1.9661193324148185,
      0.24170706118458253, 0.0031647856053746346,
      # (eta, eta)
      -0.43681327157865185, -26.243216258135989, -0.069074545228172799,
      0.069074545228172799, 1.6795229163345812, 6.3600423371063887,
      9.0543662712929457, 0.23038236942944340,
      # (eta, phi)
      -0.00040761013041956573, -0.055884637723635427, 0.16980787363689332,
      0.16980787363689332, 0.13215691971408376, 0.074027929897973410,
      -0.013634105264816422, -0.00020022006749213553,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (phi, eta)
      -0.00040761013041956573, -0.055884637723635427, 0.16980787363689332,
      0.16980787363689332, 0.13215691971408376, 0.074027929897973410,
      -0.013634105264816422, -0.00020022006749213553,
      # (phi, phi)
      -4.7397273391780800e-07, -0.00038821325825644954, -0.00017268636307043200,
      0.00017268636307043200, 0.00046653414342627255, 0.00063600423371063887,
      0.000031050638790442201, 2.1300145102574279e-07
    ),
    dim = c(m, 4, 4)
  )

  H <- logistic4_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.000067724149619770208, 0.069138420343346818, 0.45016600268752209,
      0.54983399731247791, 0.64565630622579545, 0.73105857863000488,
      0.99550372683905882, 0.99996956844309943,
      # eta
      0.0045507546375868744, 1.1713210449990778, 0.34652320179660392,
      -0.34652320179660392, -0.96089380991796066, -1.3762835326903730,
      -0.16919494282920777, -0.0022153499237622442,
      # phi
      4.7403694141529941e-06, 0.0045050809423041453, 0.017326160089830196,
      0.017326160089830196, 0.016014896831966011, 0.013762835326903730,
      0.00031332396820223662, 2.1301441574636964e-06
    ),
    nrow = m,
    ncol = 4
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, m),
      # (alpha, delta)
      rep(0, m),
      # (alpha, eta)
      rep(0, m),
      # (alpha, phi)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, eta)
      -0.0065010780536955348, -1.6733157785701111, -0.49503314542371989,
      0.49503314542371989, 1.3727054427399438, 1.9661193324148185,
      0.24170706118458253, 0.0031647856053746346,
      # (delta, phi)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      -0.0065010780536955348, -1.6733157785701111, -0.49503314542371989,
      0.49503314542371989, 1.3727054427399438, 1.9661193324148185,
      0.24170706118458253, 0.0031647856053746346,
      # (eta, eta)
      -0.43681327157865185, -26.243216258135989, -0.069074545228172799,
      0.069074545228172799, 1.6795229163345812, 6.3600423371063887,
      9.0543662712929457, 0.23038236942944340,
      # (eta, phi)
      -0.00040761013041956573, -0.055884637723635427, 0.16980787363689332,
      0.16980787363689332, 0.13215691971408376, 0.074027929897973410,
      -0.013634105264816422, -0.00020022006749213553,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (phi, eta)
      -0.00040761013041956573, -0.055884637723635427, 0.16980787363689332,
      0.16980787363689332, 0.13215691971408376, 0.074027929897973410,
      -0.013634105264816422, -0.00020022006749213553,
      # (phi, phi)
      -4.7397273391780800e-07, -0.00038821325825644954, -0.00017268636307043200,
      0.00017268636307043200, 0.00046653414342627255, 0.00063600423371063887,
      0.000031050638790442201, 2.1300145102574279e-07
    ),
    dim = c(m, 4, 4)
  )

  gh <- logistic4_gradient_hessian(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.000067724149619770208, 0.069138420343346818, 0.45016600268752209,
      0.54983399731247791, 0.64565630622579545, 0.73105857863000488,
      0.99550372683905882, 0.99996956844309943,
      # log_eta
      0.00045507546375868744, 0.11713210449990778, 0.034652320179660392,
      -0.034652320179660392, -0.096089380991796066, -0.13762835326903730,
      -0.016919494282920777, -0.00022153499237622442,
      # phi
      4.7403694141529941e-06, 0.0045050809423041453, 0.017326160089830196,
      0.017326160089830196, 0.016014896831966011, 0.013762835326903730,
      0.00031332396820223662, 2.1301441574636964e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- logistic4_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, m),
      # (alpha, delta)
      rep(0, m),
      # (alpha, log_eta)
      rep(0, m),
      # (alpha, phi)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      -0.00065010780536955348, -0.16733157785701111, -0.049503314542371989,
      0.049503314542371989, 0.13727054427399438, 0.19661193324148185,
      0.024170706118458253, 0.00031647856053746346,
      # (delta, phi)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      -0.00065010780536955348, -0.16733157785701111, -0.049503314542371989,
      0.049503314542371989, 0.13727054427399438, 0.19661193324148185,
      0.024170706118458253, 0.00031647856053746346,
      # (log_eta, log_eta)
      -0.0039130572520278311, -0.14530005808145211, 0.033961574727378664,
      -0.033961574727378664, -0.079294151828450254, -0.074027929897973410,
      0.073624168430008680, 0.0020822887019182096,
      # (log_eta, phi)
      -0.000040761013041956573, -0.0055884637723635427, 0.016980787363689332,
      0.016980787363689332, 0.013215691971408376, 0.0074027929897973410,
      -0.0013634105264816422, -0.000020022006749213553,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (phi, log_eta)
      -0.000040761013041956573, -0.0055884637723635427, 0.016980787363689332,
      0.016980787363689332, 0.013215691971408376, 0.0074027929897973410,
      -0.0013634105264816422, -0.000020022006749213553,
      # (phi, phi)
      -4.7397273391780800e-07, -0.00038821325825644954, -0.000172686363070432,
      0.000172686363070432, 0.00046653414342627255, 0.00063600423371063887,
      0.000031050638790442201, 2.1300145102574279e-07
    ),
    dim = c(m, 4, 4)
  )

  H <- logistic4_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.000067724149619770208, 0.069138420343346818, 0.45016600268752209,
      0.54983399731247791, 0.64565630622579545, 0.73105857863000488,
      0.99550372683905882, 0.99996956844309943,
      # log_eta
      0.00045507546375868744, 0.11713210449990778, 0.034652320179660392,
      -0.034652320179660392, -0.096089380991796066, -0.13762835326903730,
      -0.016919494282920777, -0.00022153499237622442,
      # phi
      4.7403694141529941e-06, 0.0045050809423041453, 0.017326160089830196,
      0.017326160089830196, 0.016014896831966011, 0.013762835326903730,
      0.00031332396820223662, 2.1301441574636964e-06
    ),
    nrow = m,
    ncol = 4
  )

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, m),
      # (alpha, delta)
      rep(0, m),
      # (alpha, log_eta)
      rep(0, m),
      # (alpha, phi)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      -0.00065010780536955348, -0.16733157785701111, -0.049503314542371989,
      0.049503314542371989, 0.13727054427399438, 0.19661193324148185,
      0.024170706118458253, 0.00031647856053746346,
      # (delta, phi)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      -0.00065010780536955348, -0.16733157785701111, -0.049503314542371989,
      0.049503314542371989, 0.13727054427399438, 0.19661193324148185,
      0.024170706118458253, 0.00031647856053746346,
      # (log_eta, log_eta)
      -0.0039130572520278311, -0.14530005808145211, 0.033961574727378664,
      -0.033961574727378664, -0.079294151828450254, -0.074027929897973410,
      0.073624168430008680, 0.0020822887019182096,
      # (log_eta, phi)
      -0.000040761013041956573, -0.0055884637723635427, 0.016980787363689332,
      0.016980787363689332, 0.013215691971408376, 0.0074027929897973410,
      -0.0013634105264816422, -0.000020022006749213553,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -6.7719563059328487e-06, -0.0064358299175773505, -0.024751657271185994,
      -0.024751657271185994, -0.022878424045665730, -0.019661193324148185,
      -0.00044760566886033802, -3.0430630820909948e-06,
      # (phi, log_eta)
      -0.000040761013041956573, -0.0055884637723635427, 0.016980787363689332,
      0.016980787363689332, 0.013215691971408376, 0.0074027929897973410,
      -0.0013634105264816422, -0.000020022006749213553,
      # (phi, phi)
      -4.7397273391780800e-07, -0.00038821325825644954, -0.000172686363070432,
      0.000172686363070432, 0.00046653414342627255, 0.00063600423371063887,
      0.000031050638790442201, 2.1300145102574279e-07
    ),
    dim = c(m, 4, 4)
  )

  gh <- logistic4_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(list(stats = ltd$stats_1), class = "logistic4")

  gh <- gradient_hessian(object, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Value of the RSS", {
  theta <- ltd$theta_4
  theta[3] <- log(theta[3])

  true_value <- 0.14751113315794743

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "logistic4"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(theta[1], NA, NA, theta[4])
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta[2:3])

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  theta <- ltd$theta_4
  theta[3] <- log(theta[3])

  true_gradient <- c(
    1.0810963357272753, 0.71060162681045976, -0.066801269060968521,
    0.0083356628669064596
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 8.9662909475324639, -0.17650012027096466, 0.17992272837749977,
      # delta
      8.9662909475324639, 6.3959662564253243, -0.24593293929902169,
      0.084439434259839158,
      # log_eta
      -0.17650012027096466, -0.24593293929902169, 0.055196449871953286,
      -0.0013524912389573294,
      # phi
      0.17992272837749977, 0.084439434259839158, -0.0013524912389573294,
      0.0031372500580327936
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "logistic4"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta)

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 4)
  expect_length(gh$H, 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  known_param <- c(theta[1], NA, NA, theta[4])
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta[2:3])

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 2)
  expect_length(gh$H, 2 * 2)

  expect_equal(gh$G, true_gradient[2:3])
  expect_equal(gh$H, true_hessian[2:3, 2:3])
})

test_that("mle_asy", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(0, 1, -1.4610997576603913, -1.2466796608100031)

  true_value <- c(
    0.83800231350869869, -0.75262962425256809, -1.4610997576603913,
    -1.2466796608100031
  )

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.055828186483232251

  theta <- c(
    alpha = 0.83800231350869869, delta = -0.75262962425256809,
    eta = exp(-1.4610997576603913), phi = -1.2466796608100031
  )

  fitted_values <- rep(
    c(
      0.83800231342409788, 0.83704910194951189, 0.6504183594995755,
      0.4944857116477287, 0.3263130606134854, 0.2034985966434521,
      0.0853778602508113, 0.0853726893035737
    ),
    k
  )

  residuals <- c(
    0.01469768657590212, -0.08090231342409788, 0.09899768657590212,
    -0.00954910194951189, -0.05924910194951189, 0.04755089805048811,
    -0.0943183594995755, 0.0543816405004245, 0.0508142883522713,
    -0.0084857116477287, 0.0406142883522713, -0.0545857116477287,
    0.0289869393865146, -0.0131130606134854, 0.0235869393865146,
    -0.0613985966434521, -0.0685778602508113, 0.0532221397491887,
    0.0373273106964263
  )

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic4_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.055828186483232251

  theta <- c(
    alpha = 0.83800231350869869, delta = -0.75262962425256809,
    eta = exp(-1.4610997576603913), phi = -1.2466796608100031
  )

  fitted_values <- rep(
    c(
      0.83800231342409788, 0.83704910194951189, 0.6504183594995755,
      0.4944857116477287, 0.3263130606134854, 0.2034985966434521,
      0.0853778602508113, 0.0853726893035737
    ),
    k
  )

  residuals <- c(
    0.01469768657590212, -0.08090231342409788, 0.09899768657590212,
    -0.00954910194951189, -0.05924910194951189, 0.04755089805048811,
    -0.0943183594995755, 0.0543816405004245, 0.0508142883522713,
    -0.0084857116477287, 0.0406142883522713, -0.0545857116477287,
    0.0289869393865146, -0.0131130606134854, 0.0235869393865146,
    -0.0613985966434521, -0.0685778602508113, 0.0532221397491887,
    0.0373273106964263
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic4_new(
    x, y, w, c(0.6, -0.6, 2, 2), max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic4_new(
    x, y, w, c(-2, 2, 7, -5), max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.17219084390371050

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-1.6104963245593389),
    phi = 1.5491035993317962
  )

  fitted_values <- rep(
    c(
      0.799999998609630709462304, 0.798355526885925310, 0.63691800939495318,
      0.50317603519210181, 0.32974482203967835, 0.16214046619972323,
      -0.09994373014167510483, -0.099999997418023120612338
    ),
    k
  )

  residuals <- c(
    0.052700001390369290537696, -0.042899998609630709462304,
    0.137000001390369290537696, 0.029144473114074690, -0.020555526885925310,
    0.086244473114074690, -0.08081800939495318, 0.06788199060504682,
    0.04212396480789819, -0.01717603519210181, 0.03192396480789819,
    -0.06327603519210181, 0.02555517796032165, -0.01654482203967835,
    0.02015517796032165, -0.02004046619972323, 0.11674373014167510483,
    0.23854373014167510483, 0.222699997418023120612338
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

  # initial values with same equalities
  object <- logistic4_new(
    x, y, w, c(0.8, -0.9, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

  # initial values with different equalities
  object <- logistic4_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fit_constrained: equalities and inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.17219084390371050

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-1.6104963245593389),
    phi = 1.5491035993317962
  )

  fitted_values <- rep(
    c(
      0.799999998609630709462304, 0.798355526885925310, 0.63691800939495318,
      0.50317603519210181, 0.32974482203967835, 0.16214046619972323,
      -0.09994373014167510483, -0.099999997418023120612338
    ),
    k
  )

  residuals <- c(
    0.052700001390369290537696, -0.042899998609630709462304,
    0.137000001390369290537696, 0.029144473114074690, -0.020555526885925310,
    0.086244473114074690, -0.08081800939495318, 0.06788199060504682,
    0.04212396480789819, -0.01717603519210181, 0.03192396480789819,
    -0.06327603519210181, 0.02555517796032165, -0.01654482203967835,
    0.02015517796032165, -0.02004046619972323, 0.11674373014167510483,
    0.23854373014167510483, 0.222699997418023120612338
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0.8, -0.9, 0.5, 2), max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0, 1, 8, -5), max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fit (weighted)", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.026367978674860599

  theta <- c(
    alpha = 0.85067612224927024, delta = -0.76270152717976731,
    eta = exp(-1.2827721627260799), phi = -1.2499066142769304
  )

  fitted_values <- rep(
    c(
      0.85067612224829011, 0.85041295034298258, 0.6895096518932521,
      0.5088390499570532, 0.3082624505000658, 0.1780800518792573,
      0.0879751091798297, 0.0879745950699930
    ),
    k
  )

  residuals <- c(
    0.00202387775170989, -0.09357612224829011, 0.08632387775170989,
    -0.02291295034298258, -0.07261295034298258, 0.03418704965701742,
    -0.1334096518932521, 0.0152903481067479, 0.0364609500429468,
    -0.0228390499570532, 0.0262609500429468, -0.0689390499570532,
    0.0470375494999342, 0.0049375494999342, 0.0416375494999342,
    -0.0359800518792573, -0.0711751091798297, 0.0506248908201703,
    0.0347254049300070
  )

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic4_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.026367978674860599

  theta <- c(
    alpha = 0.85067612224927024, delta = -0.76270152717976731,
    eta = exp(-1.2827721627260799), phi = -1.2499066142769304
  )

  fitted_values <- rep(
    c(
      0.85067612224829011, 0.85041295034298258, 0.6895096518932521,
      0.5088390499570532, 0.3082624505000658, 0.1780800518792573,
      0.0879751091798297, 0.0879745950699930
    ),
    k
  )

  residuals <- c(
    0.00202387775170989, -0.09357612224829011, 0.08632387775170989,
    -0.02291295034298258, -0.07261295034298258, 0.03418704965701742,
    -0.1334096518932521, 0.0152903481067479, 0.0364609500429468,
    -0.0228390499570532, 0.0262609500429468, -0.0689390499570532,
    0.0470375494999342, 0.0049375494999342, 0.0416375494999342,
    -0.0359800518792573, -0.0711751091798297, 0.0506248908201703,
    0.0347254049300070
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic4_new(
    x, y, w, c(0.6, -0.6, 2, 2), max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic4_new(
    x, y, w, c(2, -2, 7, -5), max_iter,
    c(0.5, -1, 0.05, -3), c(1, -0.5, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.056056044513102930

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-1.5066835513852191),
    phi = -0.80379896247074719
  )

  fitted_values <- rep(
    c(
      0.89999999974546595, 0.89860946205304318, 0.68384013229598331,
      0.50930737309648215, 0.31450764473712351, 0.16310916420690342,
      0.00001158589661434567, 1.782381258432753e-10
    ),
    k
  )

  residuals <- c(
    -0.047299999745465954, -0.14289999974546595, 0.037000000254534046,
    -0.071109462053043175, -0.12080946205304318, -0.014009462053043175,
    -0.12774013229598331, 0.02095986770401669, 0.03599262690351785,
    -0.02330737309648215, 0.02579262690351785, -0.06940737309648215,
    0.04079235526287649, -0.00130764473712351, 0.03539235526287649,
    -0.02100916420690342, 0.016788414103385654, 0.13858841410338565,
    0.12269999982176187
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

  # initial values with same equalities
  object <- logistic4_new(
    x, y, w, c(0.9, -0.9, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

  # initial values with different equalities
  object <- logistic4_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fit_constrained (weighted): equalities and inequalities", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.056056044513102930

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-1.5066835513852191),
    phi = -0.80379896247074719
  )

  fitted_values <- rep(
    c(
      0.89999999974546595, 0.89860946205304318, 0.68384013229598331,
      0.50930737309648215, 0.31450764473712351, 0.16310916420690342,
      0.00001158589661434567, 1.782381258432753e-10
    ),
    k
  )

  residuals <- c(
    -0.047299999745465954, -0.14289999974546595, 0.037000000254534046,
    -0.071109462053043175, -0.12080946205304318, -0.014009462053043175,
    -0.12774013229598331, 0.02095986770401669, 0.03599262690351785,
    -0.02330737309648215, 0.02579262690351785, -0.06940737309648215,
    0.04079235526287649, -0.00130764473712351, 0.03539235526287649,
    -0.02100916420690342, 0.016788414103385654, 0.13858841410338565,
    0.12269999982176187
  )

  object <- logistic4_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 0.05, -3), c(0.9, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0.9, -0.9, 0.5, 2), max_iter,
    c(0.9, -0.9, 0.05, -3), c(0.9, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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
  object <- logistic4_new(
    x, y, w, c(0, 1, 8, -5), max_iter,
    c(0.9, -0.9, 0.05, -3), c(0.9, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic4_fit"))
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

test_that("fisher_info", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  max_iter <- 10000

  theta <- ltd$theta_4
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- ltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.96, 2988.8313875883050, -372.86733266468586, 61.788189997723084,
      -13304.137147527459,
      # delta
      2988.8313875883050, 2106.1429590494038, -837.13692507712205,
      28.477674031988596, -10175.765394964461,
      # eta
      -372.86733266468586, -837.13692507712205, 5131.4544641662411,
      -6.6339556209875097, 12518.932610111821,
      # phi
      61.788189997723084, 28.477674031988596, -6.6339556209875097,
      1.0877076030176040, -123.54365830348264,
      # sigma
      -13304.137147527459, -10175.765394964461, 12518.932610111821,
      -123.54365830348264, 82063.580956683193
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- logistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- ltd$D$x
  y <- ltd$D$y

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
  x <- ltd$D$x
  y <- ltd$D$y

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
  x <- ltd$D$x
  y <- ltd$D$y

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
      start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = rep(1, 5)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic4")

  expect_equal(nauc(result), 0.42736069177317891)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.40547947178341702)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.40516512132073220)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 0.019739810333290101
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic4")

  expect_equal(naac(result), 1 - 0.42736069177317891)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.40547947178341702)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.40516512132073220)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 1 - 0.019739810333290101
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic4")

  expect_equal(nauc(result), 0.62568927373674172)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.66856228471375872)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.69696125124053229)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 0.59489970622431733
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic4")

  expect_equal(naac(result), 1 - 0.62568927373674172)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.66856228471375872)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.69696125124053229)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 1 - 0.59489970622431733
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})
