test_that("Constructor", {
  x <- ltd$D$x
  y <- ltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- ltd$stats_1

  start <- c(0, 1, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1, 0)
  upper_bound <- c(3, 2, 2, 5, 2)

  i <- c(1, 2, 4)

  s <- start
  s[-i] <- log(s[-i])

  lb <- lower_bound
  lb[-i] <- log(lb[-i])

  ub <- upper_bound
  ub[-i] <- log(ub[-i])

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic5"))
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

  object <- logistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, s)
  expect_equal(object$lower_bound, lb)
  expect_equal(object$upper_bound, ub)

  w <- ltd$D$w
  stats <- ltd$stats_2

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "logistic5"))
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

  object <- logistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "logistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, s)
  expect_equal(object$lower_bound, lb)
  expect_equal(object$upper_bound, ub)
})

test_that("Constructor: errors", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w
  max_iter <- 10000

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 5"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 0, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, -1, 1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    logistic5_new(x, y, w, c(0, 1, 1, 1, -1), max_iter, NULL, NULL),
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
    logistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, 0, rep(Inf, 2))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, 1, -1, rep(Inf, 2))
    ),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, Inf, Inf, Inf, 0)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    logistic5_new(
      x, y, w, NULL, max_iter, rep(-Inf, 5), c(1, Inf, Inf, Inf, -1)
    ),
    "'upper_bound[5]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

  m <- length(x)

  true_value <- c(
    0.89592655200971904, 0.76754077923524783, 0.52273910133119623,
    0.46897250433869709, 0.41668052654570241, 0.36868345260545831,
    0.20314034715510961, 0.20002130176571238
  )

  value <- logistic5_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic5")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "logistic5_fit")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0058192114146870848, 0.18922745823536025, 0.53894414095543395,
      0.61575356523043272, 0.69045639064899656, 0.75902363913505955,
      0.99551378977841484, 0.99996956890612517,
      # eta
      0.19551888241005318, 1.6603112415134353, 0.26768140111065558,
      -0.26760236333036131, -0.75871970291177206, -1.1260805664253117,
      -0.16843930399577669, -0.0022152825102922767,
      # phi
      0.00020366550251047206, 0.0063858124673593664, 0.013384070055532779,
      0.013380118166518065, 0.012645328381862868, 0.011260805664253117,
      0.00031192463702921610, 2.1300793368194968e-06,
      # nu
      -0.0094638568971289273, -0.078330357608578664, -0.049680307215175553,
      -0.037603846999036000, -0.026284720101306670, -0.016943897050071443,
      -7.0230492788902968e-06, -3.2411144051098717e-10
    ),
    nrow = m,
    ncol = 5
  )

  G <- logistic5_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 5)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

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
      # (alpha, nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, eta)
      -0.27931268915721883, -2.3718732021620504, -0.38240200158665083,
      0.38228909047194473, 1.0838852898739601, 1.6086865234647310,
      0.24062757713682385, 0.0031646893004175382,
      # (delta, phi)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (delta, nu)
      0.013519795567327039, 0.11190051086939809, 0.070971867450250790,
      0.053719781427194286, 0.037549600144723814, 0.024205567214387776,
      0.000010032927541271853, 4.6301634358712452e-10,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      -0.27931268915721883, -2.3718732021620504, -0.38240200158665083,
      0.38228909047194473, 1.0838852898739601, 1.6086865234647310,
      0.24062757713682385, 0.0031646893004175382,
      # (eta, eta)
      -9.3839529461948296, -19.265464315212061, -0.034428549758698295,
      0.036783914201004759, 0.97918107400545410, 4.1009075867378045,
      8.9735810390337902, 0.23036834838777607,
      # (eta, phi)
      -0.0077382959605148935, -0.010239815000298878, 0.13211927306739288,
      0.13196198595513042, 0.11013359925187111, 0.071598980775153122,
      -0.013498496294585228, -0.00020020723392774356,
      # (eta, nu)
      0.35649361742602426, 0.18140232815072252, -0.059715099250434786,
      0.059723879246372504, 0.15724560477882744, 0.20275210992155360,
      0.00075226156074727290, 6.7411418553929792e-08,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (phi, eta)
      -0.0077382959605148935, -0.010239815000298878, 0.13211927306739288,
      0.13196198595513042, 0.11013359925187111, 0.071598980775153122,
      -0.013498496294585228, -0.00020020723392774356,
      # (phi, phi)
      -0.000010182240610020431, -0.00028499207566881747,
      -0.000086071374396745737, 0.000091959785502511897,
      0.00027199474277929280, 0.00041009075867378045, 0.000030773597527550721,
      2.1298848778455627e-07,
      # (phi, nu)
      0.00037134751815210861, 0.00069770126211816355, -0.0029857549625217393,
      -0.0029861939623186252, -0.0026207600796471241, -0.0020275210992155360,
      -1.3930769643468017e-06, -6.4818671686470954e-11,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0.013519795567327039, 0.11190051086939809, 0.070971867450250790,
      0.053719781427194286, 0.037549600144723814, 0.024205567214387776,
      0.000010032927541271853, 4.6301634358712452e-10,
      # (nu, eta)
      0.35649361742602426, 0.18140232815072252, -0.059715099250434786,
      0.059723879246372504, 0.15724560477882744, 0.20275210992155360,
      0.00075226156074727290, 6.7411418553929792e-08,
      # (nu, phi)
      0.00037134751815210861, 0.00069770126211816355, -0.0029857549625217393,
      -0.0029861939623186252, -0.0026207600796471241, -0.0020275210992155360,
      -1.3930769643468017e-06, -6.4818671686470954e-11,
      # (nu, nu)
      -0.013032703617114021, 0.016616458649290703, 0.019396756685878632,
      0.013555665569860318, 0.0083129571533303636, 0.0044703859385760703,
      4.1875576766906138e-08, 1.3150470758856414e-14
    ),
    dim = c(m, 5, 5)
  )

  H <- logistic5_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 5 * 5)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0058192114146870848, 0.18922745823536025, 0.53894414095543395,
      0.61575356523043272, 0.69045639064899656, 0.75902363913505955,
      0.99551378977841484, 0.99996956890612517,
      # eta
      0.19551888241005318, 1.6603112415134353, 0.26768140111065558,
      -0.26760236333036131, -0.75871970291177206, -1.1260805664253117,
      -0.16843930399577669, -0.0022152825102922767,
      # phi
      0.00020366550251047206, 0.0063858124673593664, 0.013384070055532779,
      0.013380118166518065, 0.012645328381862868, 0.011260805664253117,
      0.00031192463702921610, 2.1300793368194968e-06,
      # nu
      -0.0094638568971289273, -0.078330357608578664, -0.049680307215175553,
      -0.037603846999036000, -0.026284720101306670, -0.016943897050071443,
      -7.0230492788902968e-06, -3.2411144051098717e-10
    ),
    nrow = m,
    ncol = 5
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
      # (alpha, nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, eta)
      -0.27931268915721883, -2.3718732021620504, -0.38240200158665083,
      0.38228909047194473, 1.0838852898739601, 1.6086865234647310,
      0.24062757713682385, 0.0031646893004175382,
      # (delta, phi)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (delta, nu)
      0.013519795567327039, 0.11190051086939809, 0.070971867450250790,
      0.053719781427194286, 0.037549600144723814, 0.024205567214387776,
      0.000010032927541271853, 4.6301634358712452e-10,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      -0.27931268915721883, -2.3718732021620504, -0.38240200158665083,
      0.38228909047194473, 1.0838852898739601, 1.6086865234647310,
      0.24062757713682385, 0.0031646893004175382,
      # (eta, eta)
      -9.3839529461948296, -19.265464315212061, -0.034428549758698295,
      0.036783914201004759, 0.97918107400545410, 4.1009075867378045,
      8.9735810390337902, 0.23036834838777607,
      # (eta, phi)
      -0.0077382959605148935, -0.010239815000298878, 0.13211927306739288,
      0.13196198595513042, 0.11013359925187111, 0.071598980775153122,
      -0.013498496294585228, -0.00020020723392774356,
      # (eta, nu)
      0.35649361742602426, 0.18140232815072252, -0.059715099250434786,
      0.059723879246372504, 0.15724560477882744, 0.20275210992155360,
      0.00075226156074727290, 6.7411418553929792e-08,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (phi, eta)
      -0.0077382959605148935, -0.010239815000298878, 0.13211927306739288,
      0.13196198595513042, 0.11013359925187111, 0.071598980775153122,
      -0.013498496294585228, -0.00020020723392774356,
      # (phi, phi)
      -0.000010182240610020431, -0.00028499207566881747,
      -0.000086071374396745737, 0.000091959785502511897,
      0.00027199474277929280, 0.00041009075867378045, 0.000030773597527550721,
      2.1298848778455627e-07,
      # (phi, nu)
      0.00037134751815210861, 0.00069770126211816355, -0.0029857549625217393,
      -0.0029861939623186252, -0.0026207600796471241, -0.0020275210992155360,
      -1.3930769643468017e-06, -6.4818671686470954e-11,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0.013519795567327039, 0.11190051086939809, 0.070971867450250790,
      0.053719781427194286, 0.037549600144723814, 0.024205567214387776,
      0.000010032927541271853, 4.6301634358712452e-10,
      # (nu, eta)
      0.35649361742602426, 0.18140232815072252, -0.059715099250434786,
      0.059723879246372504, 0.15724560477882744, 0.20275210992155360,
      0.00075226156074727290, 6.7411418553929792e-08,
      # (nu, phi)
      0.00037134751815210861, 0.00069770126211816355, -0.0029857549625217393,
      -0.0029861939623186252, -0.0026207600796471241, -0.0020275210992155360,
      -1.3930769643468017e-06, -6.4818671686470954e-11,
      # (nu, nu)
      -0.013032703617114021, 0.016616458649290703, 0.019396756685878632,
      0.013555665569860318, 0.0083129571533303636, 0.0044703859385760703,
      4.1875576766906138e-08, 1.3150470758856414e-14
    ),
    dim = c(m, 5, 5)
  )

  gh <- logistic5_gradient_hessian(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 5)
  expect_length(gh$H, m * 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0058192114146870848, 0.18922745823536025, 0.53894414095543395,
      0.61575356523043272, 0.69045639064899656, 0.75902363913505955,
      0.99551378977841484, 0.99996956890612517,
      # log_eta
      0.019551888241005318, 0.16603112415134353, 0.026768140111065558,
      -0.026760236333036131, -0.075871970291177206, -0.11260805664253117,
      -0.016843930399577669, -0.00022152825102922767,
      # phi
      0.00020366550251047206, 0.0063858124673593664, 0.013384070055532779,
      0.013380118166518065, 0.012645328381862868, 0.011260805664253117,
      0.00031192463702921610, 2.1300793368194968e-06,
      # log_nu
      -0.018927713794257855, -0.15666071521715733, -0.099360614430351106,
      -0.075207693998072001, -0.052569440202613339, -0.033887794100142886,
      -0.000014046098557780594, -6.4822288102197433e-10
    ),
    nrow = m,
    ncol = 5
  )

  G <- logistic5_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 5)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

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
      # (alpha, log_nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      -0.027931268915721883, -0.23718732021620504, -0.038240200158665083,
      0.038228909047194473, 0.10838852898739601, 0.16086865234647310,
      0.024062757713682385, 0.00031646893004175382,
      # (delta, phi)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (delta, log_nu)
      0.027039591134654078, 0.22380102173879618, 0.14194373490050158,
      0.10743956285438857, 0.075099200289447628, 0.048411134428775551,
      0.000020065855082543705, 9.2603268717424905e-10,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      -0.027931268915721883, -0.23718732021620504, -0.038240200158665083,
      0.038228909047194473, 0.10838852898739601, 0.16086865234647310,
      0.024062757713682385, 0.00031646893004175382,
      # (log_eta, log_eta)
      -0.074287641220942978, -0.026623519000777082, 0.026423854613478575,
      -0.026392397191026083, -0.066080159551122665, -0.071598980775153122,
      0.072891879990760233, 0.0020821552328485330,
      # (log_eta, phi)
      -0.00077382959605148935, -0.0010239815000298878, 0.013211927306739288,
      0.013196198595513042, 0.011013359925187111, 0.0071598980775153122,
      -0.0013498496294585228, -0.000020020723392774356,
      # (log_eta, log_nu)
      0.071298723485204852, 0.036280465630144505, -0.011943019850086957,
      0.011944775849274501, 0.031449120955765489, 0.040550421984310720,
      0.00015045231214945458, 1.3482283710785958e-08,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (phi, log_eta)
      -0.00077382959605148935, -0.0010239815000298878, 0.013211927306739288,
      0.013196198595513042, 0.011013359925187111, 0.0071598980775153122,
      -0.0013498496294585228, -0.000020020723392774356,
      # (phi, phi)
      -0.000010182240610020431, -0.00028499207566881747,
      -0.000086071374396745737, 0.000091959785502511897, 0.00027199474277929280,
      0.00041009075867378045, 0.000030773597527550721, 2.1298848778455627e-07,
      # (phi, log_nu)
      0.00074269503630421721, 0.0013954025242363271, -0.0059715099250434786,
      -0.0059723879246372504, -0.0052415201592942481, -0.0040550421984310720,
      -2.7861539286936033e-06, -1.2963734337294191e-10,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0.027039591134654078, 0.22380102173879618, 0.14194373490050158,
      0.10743956285438857, 0.075099200289447628, 0.048411134428775551,
      0.000020065855082543705, 9.2603268717424905e-10,
      # (log_nu, log_eta)
      0.071298723485204852, 0.036280465630144505, -0.011943019850086957,
      0.011944775849274501, 0.031449120955765489, 0.040550421984310720,
      0.00015045231214945458, 1.3482283710785958e-08,
      # (log_nu, phi)
      0.00074269503630421721, 0.0013954025242363271, -0.0059715099250434786,
      -0.0059723879246372504, -0.0052415201592942481, -0.0040550421984310720,
      -2.7861539286936033e-06, -1.2963734337294191e-10,
      # (log_nu, log_nu)
      -0.071058528262713939, -0.090194880619994517, -0.021773587686836579,
      -0.020985031718630730, -0.019317611589291885, -0.016006250345838605,
      -0.000013878596250712969, -6.4817027913893891e-10
    ),
    dim = c(m, 5, 5)
  )

  H <- logistic5_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 5 * 5)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0058192114146870848, 0.18922745823536025, 0.53894414095543395,
      0.61575356523043272, 0.69045639064899656, 0.75902363913505955,
      0.99551378977841484, 0.99996956890612517,
      # log_eta
      0.019551888241005318, 0.16603112415134353, 0.026768140111065558,
      -0.026760236333036131, -0.075871970291177206, -0.11260805664253117,
      -0.016843930399577669, -0.00022152825102922767,
      # phi
      0.00020366550251047206, 0.0063858124673593664, 0.013384070055532779,
      0.013380118166518065, 0.012645328381862868, 0.011260805664253117,
      0.00031192463702921610, 2.1300793368194968e-06,
      # log_nu
      -0.018927713794257855, -0.15666071521715733, -0.099360614430351106,
      -0.075207693998072001, -0.052569440202613339, -0.033887794100142886,
      -0.000014046098557780594, -6.4822288102197433e-10
    ),
    nrow = m,
    ncol = 5
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
      # (alpha, log_nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      -0.027931268915721883, -0.23718732021620504, -0.038240200158665083,
      0.038228909047194473, 0.10838852898739601, 0.16086865234647310,
      0.024062757713682385, 0.00031646893004175382,
      # (delta, phi)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (delta, log_nu)
      0.027039591134654078, 0.22380102173879618, 0.14194373490050158,
      0.10743956285438857, 0.075099200289447628, 0.048411134428775551,
      0.000020065855082543705, 9.2603268717424905e-10,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      -0.027931268915721883, -0.23718732021620504, -0.038240200158665083,
      0.038228909047194473, 0.10838852898739601, 0.16086865234647310,
      0.024062757713682385, 0.00031646893004175382,
      # (log_eta, log_eta)
      -0.074287641220942978, -0.026623519000777082, 0.026423854613478575,
      -0.026392397191026083, -0.066080159551122665, -0.071598980775153122,
      0.072891879990760233, 0.0020821552328485330,
      # (log_eta, phi)
      -0.00077382959605148935, -0.0010239815000298878, 0.013211927306739288,
      0.013196198595513042, 0.011013359925187111, 0.0071598980775153122,
      -0.0013498496294585228, -0.000020020723392774356,
      # (log_eta, log_nu)
      0.071298723485204852, 0.036280465630144505, -0.011943019850086957,
      0.011944775849274501, 0.031449120955765489, 0.040550421984310720,
      0.00015045231214945458, 1.3482283710785958e-08,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      -0.00029095071787210295, -0.0091225892390848092, -0.019120100079332541,
      -0.019114454523597236, -0.018064754831232668, -0.016086865234647310,
      -0.00044560662432745157, -3.0429704811707098e-06,
      # (phi, log_eta)
      -0.00077382959605148935, -0.0010239815000298878, 0.013211927306739288,
      0.013196198595513042, 0.011013359925187111, 0.0071598980775153122,
      -0.0013498496294585228, -0.000020020723392774356,
      # (phi, phi)
      -0.000010182240610020431, -0.00028499207566881747,
      -0.000086071374396745737, 0.000091959785502511897, 0.00027199474277929280,
      0.00041009075867378045, 0.000030773597527550721, 2.1298848778455627e-07,
      # (phi, log_nu)
      0.00074269503630421721, 0.0013954025242363271, -0.0059715099250434786,
      -0.0059723879246372504, -0.0052415201592942481, -0.0040550421984310720,
      -2.7861539286936033e-06, -1.2963734337294191e-10,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0.027039591134654078, 0.22380102173879618, 0.14194373490050158,
      0.10743956285438857, 0.075099200289447628, 0.048411134428775551,
      0.000020065855082543705, 9.2603268717424905e-10,
      # (log_nu, log_eta)
      0.071298723485204852, 0.036280465630144505, -0.011943019850086957,
      0.011944775849274501, 0.031449120955765489, 0.040550421984310720,
      0.00015045231214945458, 1.3482283710785958e-08,
      # (log_nu, phi)
      0.00074269503630421721, 0.0013954025242363271, -0.0059715099250434786,
      -0.0059723879246372504, -0.0052415201592942481, -0.0040550421984310720,
      -2.7861539286936033e-06, -1.2963734337294191e-10,
      # (log_nu, log_nu)
      -0.071058528262713939, -0.090194880619994517, -0.021773587686836579,
      -0.020985031718630730, -0.019317611589291885, -0.016006250345838605,
      -0.000013878596250712969, -6.4817027913893891e-10
    ),
    dim = c(m, 5, 5)
  )

  gh <- logistic5_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 5)
  expect_length(gh$H, m * 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(list(stats = ltd$stats_1), class = "logistic5")

  gh <- gradient_hessian(object, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 5)
  expect_length(gh$H, m * 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Value of the RSS", {
  theta <- ltd$theta_5
  theta[c(3, 5)] <- log(theta[c(3, 5)])

  true_value <- 0.15244617497377583

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "logistic5"
  )

  rss_fn <- rss(object)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(theta[1], NA, NA, theta[4], theta[5])
  rss_fn <- rss_fixed(object, known_param)

  expect_type(rss_fn, "closure")

  value <- rss_fn(theta[c(2, 3)])

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  theta <- ltd$theta_5
  theta[c(3, 5)] <- log(theta[c(3, 5)])

  true_gradient <- c(
    0.39429724107057858, 0.42797633159869274, -0.077947312409148747,
    -0.00023621212406917749, 0.038023700840970313
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 9.9474325113277449, 0.12911101550078578, 0.14987981684998431,
      -1.2179414994405571,
      # delta
      9.9474325113277449, 7.1933955859776619, -0.10749949187909019,
      0.086711610373157204, -0.57056111883848218,
      # log_eta
      0.12911101550078578, -0.10749949187909019, 0.097841672037622319,
      -0.0023309653544513332, -0.039845880993324202,
      # phi
      0.14987981684998431, 0.086711610373157204, -0.0023309653544513332,
      0.0020256988639591056, -0.012298547815014724,
      # log_nu
      -1.2179414994405571, -0.57056111883848218, -0.039845880993324202,
      -0.012298547815014724, 0.13270723470839286
    ),
    nrow = 5,
    ncol = 5
  )

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "logistic5"
  )

  rss_gh <- rss_gradient_hessian(object)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta)

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 5)
  expect_length(gh$H, 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  known_param <- c(theta[1], NA, NA, theta[4], theta[5])
  rss_gh <- rss_gradient_hessian_fixed(object, known_param)

  expect_type(rss_gh, "closure")

  gh <- rss_gh(theta[c(2, 3)])

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 2)
  expect_length(gh$H, 2 * 2)

  expect_equal(gh$G, true_gradient[c(2, 3)])
  expect_equal(gh$H, true_hessian[c(2, 3), c(2, 3)])
})

test_that("mle_asy", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(
    0, 1, -0.13948206931816815, 3.2950883977695700, 2.1396800314205292
  )

  true_value <- c(
    0.84948649909664558, -0.75695416061024096, -0.13948206931816815,
    3.2950883977695700, 2.1396800314205292
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.050168602337768019

  theta <- c(
    alpha = 0.84948649909664558, delta = -0.75695416061024096,
    eta = exp(-0.13948206931816815), phi = 3.29508839776957,
    nu = exp(2.1396800314205292)
  )

  fitted_values <- rep(
    c(
      0.84947145597244308, 0.83001285437491464, 0.62226246927410375,
      0.50732524066346644, 0.33637480874558723, 0.14349972605987800,
      0.092532338486404620, 0.092532338486404619
    ),
    k
  )

  residuals <- c(
    0.0032285440275569165, -0.092371455972443084, 0.087528544027556916,
    -0.0025128543749146426, -0.052212854374914643, 0.054587145625085357,
    -0.066162469274103751, 0.082537530725896249, 0.037974759336533557,
    -0.021325240663466443, 0.027774759336533557, -0.067425240663466443,
    0.018925191254412770, -0.023174808745587230, 0.013525191254412770,
    -0.0013997260598779998, -0.075732338486404620, 0.046067661513595380,
    0.030167661513595381
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic5_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.050168602337768019

  theta <- c(
    alpha = 0.84948649909664558, delta = -0.75695416061024096,
    eta = exp(-0.13948206931816815), phi = 3.29508839776957,
    nu = exp(2.1396800314205292)
  )

  fitted_values <- rep(
    c(
      0.84947145597244308, 0.83001285437491464, 0.62226246927410375,
      0.50732524066346644, 0.33637480874558723, 0.14349972605987800,
      0.092532338486404620, 0.092532338486404619
    ),
    k
  )

  residuals <- c(
    0.0032285440275569165, -0.092371455972443084, 0.087528544027556916,
    -0.0025128543749146426, -0.052212854374914643, 0.054587145625085357,
    -0.066162469274103751, 0.082537530725896249, 0.037974759336533557,
    -0.021325240663466443, 0.027774759336533557, -0.067425240663466443,
    0.018925191254412770, -0.023174808745587230, 0.013525191254412770,
    -0.0013997260598779998, -0.075732338486404620, 0.046067661513595380,
    0.030167661513595381
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -5, 3), c(1, -0.5, 5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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
    x, y, w, c(0.6, -0.6, 2, 2, 8), max_iter,
    c(0.5, -1, 0.05, -5, 3), c(1, -0.5, 5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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
    x, y, w, c(-2, 2, 7, -3, 4), max_iter,
    c(0.5, -1, 0.05, -5, 3), c(1, -0.5, 5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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
  x <- ltd$D$x
  y <- ltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.17215094563016538

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-1.5326708776200168),
    phi = 1.9358453966983796, nu = exp(0.18317813414869142)
  )

  fitted_values <- rep(
    c(
      0.79999999153278315, 0.79752371287734322, 0.63490156174551599,
      0.50449794753818087, 0.33044038292270020, 0.15762370921429484,
      -0.099972054046131580, -0.099999999428701282
    ),
    k
  )

  residuals <- c(
    0.052700008467216853, -0.042899991532783147, 0.13700000846721685,
    0.029976287122656784, -0.019723712877343216, 0.087076287122656784,
    -0.078801561745515992, 0.069898438254484008, 0.040802052461819132,
    -0.018497947538180868, 0.030602052461819132, -0.064597947538180868,
    0.024859617077299801, -0.017240382922700199, 0.019459617077299801,
    -0.015523709214294845, 0.11677205404613158, 0.23857205404613158,
    0.22269999942870128
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 3)), c(0.8, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic5_new(
    x, y, w, c(0.8, -0.9, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 3)), c(0.8, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 3)), c(0.8, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.17215094563016538

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-1.5326708776200168),
    phi = 1.9358453966983796, nu = exp(0.18317813414869142)
  )

  fitted_values <- rep(
    c(
      0.79999999153278315, 0.79752371287734322, 0.63490156174551599,
      0.50449794753818087, 0.33044038292270020, 0.15762370921429484,
      -0.099972054046131580, -0.099999999428701282
    ),
    k
  )

  residuals <- c(
    0.052700008467216853, -0.042899991532783147, 0.13700000846721685,
    0.029976287122656784, -0.019723712877343216, 0.087076287122656784,
    -0.078801561745515992, 0.069898438254484008, 0.040802052461819132,
    -0.018497947538180868, 0.030602052461819132, -0.064597947538180868,
    0.024859617077299801, -0.017240382922700199, 0.019459617077299801,
    -0.015523709214294845, 0.11677205404613158, 0.23857205404613158,
    0.22269999942870128
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 0.05, -3, 0.5), c(0.8, -0.9, 2, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0.8, -0.9, 0.5, 2, 1.5), max_iter,
    c(0.8, -0.9, 0.05, -3, 0.5), c(0.8, -0.9, 2, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, 8, -8, 4), max_iter,
    c(0.8, -0.9, 0.05, -3, 0.5), c(0.8, -0.9, 2, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.026367789609414469

  theta <- c(
    alpha = 0.85062856577509035, delta = -0.76270201072171628,
    eta = exp(-1.2882638383496943), phi = -1.2804646174528685,
    nu = exp(-0.019767200884507061)
  )

  fitted_values <- rep(
    c(
      0.85062856577440989, 0.85038705397619911, 0.68988971453996974,
      0.50866177107419315, 0.30821950267475678, 0.17833789292363378,
      0.087927106108376473, 0.087926555053940805
    ),
    k
  )

  residuals <- c(
    0.0020714342255901108, -0.093528565774409889, 0.086371434225590111,
    -0.022887053976199108, -0.072587053976199108, 0.034212946023800892,
    -0.13378971453996974, 0.014910285460030258, 0.036638228925806847,
    -0.022661771074193153, 0.026438228925806847, -0.068761771074193153,
    0.047080497325243221, 0.0049804973252432212, 0.041680497325243221,
    -0.036237892923633776, -0.071127106108376473, 0.050672893891623527,
    0.034773444946059195
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- logistic5_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  rss_value <- 0.026367789609414469

  theta <- c(
    alpha = 0.85062856577509035, delta = -0.76270201072171628,
    eta = exp(-1.2882638383496943), phi = -1.2804646174528685,
    nu = exp(-0.019767200884507061)
  )

  fitted_values <- rep(
    c(
      0.85062856577440989, 0.85038705397619911, 0.68988971453996974,
      0.50866177107419315, 0.30821950267475678, 0.17833789292363378,
      0.087927106108376473, 0.087926555053940805
    ),
    k
  )

  residuals <- c(
    0.0020714342255901108, -0.093528565774409889, 0.086371434225590111,
    -0.022887053976199108, -0.072587053976199108, 0.034212946023800892,
    -0.13378971453996974, 0.014910285460030258, 0.036638228925806847,
    -0.022661771074193153, 0.026438228925806847, -0.068761771074193153,
    0.047080497325243221, 0.0049804973252432212, 0.041680497325243221,
    -0.036237892923633776, -0.071127106108376473, 0.050672893891623527,
    0.034773444946059195
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -3, 0.5), c(1, -0.5, 2, 3, 4)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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
    x, y, w, c(0.6, -0.6, 0.5, 2, 3), max_iter,
    c(0.5, -1, 0.05, -3, 0.5), c(1, -0.5, 2, 3, 4)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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
    x, y, w, c(2, -2, 7, -5, 7), max_iter,
    c(0.5, -1, 0.05, -3, 0.5), c(1, -0.5, 2, 3, 4)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
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

test_that("fit_constrained (weighted): equalities", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.052161993061533667

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-0.70834912297949987),
    phi = 3.3026845129113375, nu = exp(1.6191105670237828)
  )

  fitted_values <- rep(
    c(
      0.89997253346668808, 0.87463833970460955, 0.63654713201903529,
      0.51177369572845392, 0.33604579265372012, 0.13931961659096753,
      9.2690269562236922e-11, 1.8770544051687683e-21
    ),
    k
  )

  residuals <- c(
    -0.047272533466688080, -0.14287253346668808, 0.037027466533311920,
    -0.047138339704609554, -0.096838339704609554, 0.0099616602953904463,
    -0.080447132019035287, 0.068252867980964713, 0.033526304271546084,
    -0.025773695728453916, 0.023326304271546084, -0.071873695728453916,
    0.019254207346279876, -0.022845792653720124, 0.013854207346279876,
    0.0027803834090324655, 0.016799999907309730, 0.13859999990730973,
    0.1227
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- logistic5_new(
    x, y, w, c(0.9, -0.9, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- logistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  rss_value <- 0.052161993061533667

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-0.70834912297949987),
    phi = 3.3026845129113375, nu = exp(1.6191105670237828)
  )

  fitted_values <- rep(
    c(
      0.89997253346668808, 0.87463833970460955, 0.63654713201903529,
      0.51177369572845392, 0.33604579265372012, 0.13931961659096753,
      9.2690269562236922e-11, 1.8770544051687683e-21
    ),
    k
  )

  residuals <- c(
    -0.047272533466688080, -0.14287253346668808, 0.037027466533311920,
    -0.047138339704609554, -0.096838339704609554, 0.0099616602953904463,
    -0.080447132019035287, 0.068252867980964713, 0.033526304271546084,
    -0.025773695728453916, 0.023326304271546084, -0.071873695728453916,
    0.019254207346279876, -0.022845792653720124, 0.013854207346279876,
    0.0027803834090324655, 0.016799999907309730, 0.13859999990730973,
    0.1227
  )

  object <- logistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 0.05, -5, 2), c(0.9, -0.9, 3, 5, 7)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- logistic5_new(
    x, y, w, c(0.9, -0.9, 0.5, 2, 3), max_iter,
    c(0.9, -0.9, 0.05, -5, 2), c(0.9, -0.9, 3, 5, 7)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- logistic5_new(
    x, y, w, c(0, 1, 8, -8, 0.5), max_iter,
    c(0.9, -0.9, 0.05, -5, 2), c(0.9, -0.9, 3, 5, 7)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "logistic5_fit"))
  expect_true(inherits(result, "logistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  max_iter <- 10000

  theta <- ltd$theta_5
  names(theta) <- c("alpha", "delta", "eta", "phi", "nu")

  sigma <- ltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.96, 3355.2518181039595, 809.56316763805272, 52.581615187017823,
      -226.81180476705522, -3044.3650930891329,
      # delta
      3355.2518181039595, 2379.6377719606501, -189.72541160361480,
      29.695154858733559, -104.20785201427244, -6183.1689982196366,
      # eta
      809.56316763805272, -189.72541160361480, 8275.1143459227209,
      -3.7929172793373654, -122.28116297367687, 15645.875128100754,
      # phi
      52.581615187017823, 29.695154858733559, -3.7929172793373654,
      0.72389494640083248, -2.2876753231090397, -1.8946819202146866,
      # nu
      -226.81180476705522, -104.20785201427244, -122.28116297367687,
      -2.2876753231090397, 9.9150478445799180, -380.37357037713062,
      # sigma
      -3044.3650930891329, -6183.1689982196366, 15645.875128100754,
      -1.8946819202146866, -380.37357037713062, 88937.679628153
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "sigma"
  )

  object <- logistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- ltd$D$x
  y <- ltd$D$y

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
  x <- ltd$D$x
  y <- ltd$D$y

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
  x <- ltd$D$x
  y <- ltd$D$y

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
      start = c(0, Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(-Inf, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(1, 1, 1, 1, 1, 1)
    ),
    "'start' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, -1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, 0, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, 1, 1, -1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, 1, 1, 0)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic5")

  expect_equal(nauc(result), 0.42753516878441386)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.40530818958835738)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.40513264312426204)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 0.019678613960950707
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic5")

  expect_equal(naac(result), 1 - 0.42753516878441386)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.40530818958835738)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.40513264312426204)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 1 - 0.019678613960950707
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic5")

  expect_equal(nauc(result), 0.65302073661956755)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.65733232476226312)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.74700167012689969)
  expect_equal(nauc(result, xlim = c(-25, -18), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 0.62649140813097456
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "logistic5")

  expect_equal(naac(result), 1 - 0.65302073661956755)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.65733232476226312)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.74700167012689969)
  expect_equal(naac(result, xlim = c(-25, -18), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 1 - 0.62649140813097456
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
