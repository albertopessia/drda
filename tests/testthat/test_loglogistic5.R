test_that("Constructor", {
  x <- lltd$D$x
  y <- lltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- lltd$stats_1

  start <- c(0, 1, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1, 0)
  upper_bound <- c(3, 2, 2, 5, 2)

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic5"))
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

  object <- loglogistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  i <- c(1, 2)

  expect_true(inherits(object, "loglogistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(start[i], log(start[-i])))
  expect_equal(object$lower_bound, c(lower_bound[i], log(lower_bound[-i])))
  expect_equal(object$upper_bound, c(upper_bound[i], log(upper_bound[-i])))

  w <- lltd$D$w
  stats <- lltd$stats_2

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic5"))
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

  object <- loglogistic5_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic5"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(start[i], log(start[-i])))
  expect_equal(object$lower_bound, c(lower_bound[i], log(lower_bound[-i])))
  expect_equal(object$upper_bound, c(upper_bound[i], log(upper_bound[-i])))
})

test_that("Constructor: errors", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w
  max_iter <- 10000

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 5"
  )

  expect_error(
    loglogistic5_new(x, y, w, c(0, 1, 1, 1, 1, 1), max_iter, NULL, NULL),
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
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_value <- c(
    0.9, 0.18445824720006730, 0.11221445773227051, 0.10367818040256333,
    0.10155793721368916, 0.10079880199650629, 0.10000079999880000, 0.1
  )

  value <- loglogistic5_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = m, ncol = 1)),
    class = "loglogistic5"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(
    list(stats = matrix(x, nrow = m, ncol = 1)),
    class = "loglogistic5_fit"
  )

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.89442719099991588, 0.98473192783466186, 0.99540227449679584,
      0.99805257848288855, 0.99900149750436714, 0.99999900000150000,
      0.99999999900000000,
      # eta
      0, -0.049597574852619471, -0.016547009924508010, -0.0065450328477826702,
      -0.0032301820894599116, -0.0018365556535187897, -3.6841250964096577e-06,
      -5.5262042066070970e-09,
      # phi
      0, 0.21466252583997981, 0.035808433739442250, 0.010958557150423440,
      0.0046601676816321644, 0.0023928179580942925, 2.3999928000180000e-06,
      2.3999999928000000e-09,
      # nu
      0, -0.0041400443183462087, -0.000092294662409466473,
      -8.4296992015322137e-06, -1.5150102795175839e-06, -3.9853738990374935e-07,
      -3.9999853333739999e-13, -3.9999999853333334e-19
    ),
    nrow = m,
    ncol = 5
  )

  G <- loglogistic5_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 5)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

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
      0, 0.061996968565774338, 0.020683762405635012, 0.0081812910597283377,
      0.0040377276118248895, 0.0022956945668984871, 4.6051563705120721e-06,
      6.9077552582588713e-09,
      # (delta, phi)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (delta, nu)
      0, 0.0051750553979327608, 0.00011536832801183309, 0.000010537124001915267,
      1.8937628493969799e-06, 4.9817173737968668e-07, 4.9999816667174999e-13,
      4.9999999816666667e-19,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, 0.061996968565774338, 0.020683762405635012, 0.0081812910597283377,
      0.0040377276118248895, 0.0022956945668984871, 4.6051563705120721e-06,
      6.9077552582588713e-09,
      # (eta, eta)
      0, 0.024064893420192822, 0.021896343526661215, 0.011565742133055052,
      0.0066777706907990280, 0.0042161645155449665, 0.000016965972157468888,
      3.8173666166402468e-08,
      # (eta, phi)
      0, -0.032600731910507618, -0.035448474749458854, -0.015712038381754062,
      -0.0080805974325987036, -0.0046955649656609970, -0.000010252344532163419,
      -1.5778612572485453e-08,
      # (eta, nu)
      0, 0.0046727914805964575, 0.00024877367203709549, 0.000029953802084475141,
      6.2782718174097486e-06, 1.8319740391874217e-06, 3.6841158861165653e-12,
      5.5262041927915865e-18,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (phi, eta)
      0, -0.032600731910507618, -0.035448474749458854, -0.015712038381754062,
      -0.0080805974325987036, -0.0046955649656609970, -0.000010252344532163419,
      -1.5778612572485453e-08,
      # (phi, phi)
      0, 0.23612877842397779, 0.066733899241687829, 0.021464696803811051,
      0.0092387370964653221, 0.0047641435393194747, 4.7999640001439995e-06,
      4.7999999640000001e-09,
      # (phi, nu)
      0, -0.020224239288494118, -0.00053835681442748258,
      -0.000050152605746875430, -9.0576316163299838e-06,
      -2.3868486486273266e-06, -2.3999868000487998e-12, -2.3999999868000000e-18,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0, 0.0051750553979327608, 0.00011536832801183309, 0.000010537124001915267,
      1.8937628493969799e-06, 4.9817173737968668e-07, 4.9999816667174999e-13,
      4.9999999816666667e-19,
      # (nu, eta)
      0, 0.0046727914805964575, 0.00024877367203709549, 0.000029953802084475141,
      6.2782718174097486e-06, 1.8319740391874217e-06, 3.6841158861165653e-12,
      5.5262041927915865e-18,
      # (nu, phi)
      0, -0.020224239288494118, -0.00053835681442748258,
      -0.000050152605746875430, -9.0576316163299838e-06,
      -2.3868486486273266e-06, -2.3999868000487998e-12, -2.3999999868000000e-18,
      # (nu, nu)
      0, 0.00053838172231535815, 1.8585117094458687e-06, 5.1508169508281362e-08,
      3.9283900527873613e-09, 5.3021210764429046e-10, 5.3333020001214663e-19,
      5.3333333020000001e-28
    ),
    dim = c(m, 5, 5)
  )

  H <- loglogistic5_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 5 * 5)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.89442719099991588, 0.98473192783466186, 0.99540227449679584,
      0.99805257848288855, 0.99900149750436714, 0.99999900000150000,
      0.99999999900000000,
      # eta
      0, -0.049597574852619471, -0.016547009924508010, -0.0065450328477826702,
      -0.0032301820894599116, -0.0018365556535187897, -3.6841250964096577e-06,
      -5.5262042066070970e-09,
      # phi
      0, 0.21466252583997981, 0.035808433739442250, 0.010958557150423440,
      0.0046601676816321644, 0.0023928179580942925, 2.3999928000180000e-06,
      2.3999999928000000e-09,
      # nu
      0, -0.0041400443183462087, -0.000092294662409466473,
      -8.4296992015322137e-06, -1.5150102795175839e-06, -3.9853738990374935e-07,
      -3.9999853333739999e-13, -3.9999999853333334e-19
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
      0, 0.061996968565774338, 0.020683762405635012, 0.0081812910597283377,
      0.0040377276118248895, 0.0022956945668984871, 4.6051563705120721e-06,
      6.9077552582588713e-09,
      # (delta, phi)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (delta, nu)
      0, 0.0051750553979327608, 0.00011536832801183309, 0.000010537124001915267,
      1.8937628493969799e-06, 4.9817173737968668e-07, 4.9999816667174999e-13,
      4.9999999816666667e-19,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, 0.061996968565774338, 0.020683762405635012, 0.0081812910597283377,
      0.0040377276118248895, 0.0022956945668984871, 4.6051563705120721e-06,
      6.9077552582588713e-09,
      # (eta, eta)
      0, 0.024064893420192822, 0.021896343526661215, 0.011565742133055052,
      0.0066777706907990280, 0.0042161645155449665, 0.000016965972157468888,
      3.8173666166402468e-08,
      # (eta, phi)
      0, -0.032600731910507618, -0.035448474749458854, -0.015712038381754062,
      -0.0080805974325987036, -0.0046955649656609970, -0.000010252344532163419,
      -1.5778612572485453e-08,
      # (eta, nu)
      0, 0.0046727914805964575, 0.00024877367203709549, 0.000029953802084475141,
      6.2782718174097486e-06, 1.8319740391874217e-06, 3.6841158861165653e-12,
      5.5262041927915865e-18,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (phi, eta)
      0, -0.032600731910507618, -0.035448474749458854, -0.015712038381754062,
      -0.0080805974325987036, -0.0046955649656609970, -0.000010252344532163419,
      -1.5778612572485453e-08,
      # (phi, phi)
      0, 0.23612877842397779, 0.066733899241687829, 0.021464696803811051,
      0.0092387370964653221, 0.0047641435393194747, 4.7999640001439995e-06,
      4.7999999640000001e-09,
      # (phi, nu)
      0, -0.020224239288494118, -0.00053835681442748258,
      -0.000050152605746875430, -9.0576316163299838e-06,
      -2.3868486486273266e-06, -2.3999868000487998e-12, -2.3999999868000000e-18,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0, 0.0051750553979327608, 0.00011536832801183309, 0.000010537124001915267,
      1.8937628493969799e-06, 4.9817173737968668e-07, 4.9999816667174999e-13,
      4.9999999816666667e-19,
      # (nu, eta)
      0, 0.0046727914805964575, 0.00024877367203709549, 0.000029953802084475141,
      6.2782718174097486e-06, 1.8319740391874217e-06, 3.6841158861165653e-12,
      5.5262041927915865e-18,
      # (nu, phi)
      0, -0.020224239288494118, -0.00053835681442748258,
      -0.000050152605746875430, -9.0576316163299838e-06,
      -2.3868486486273266e-06, -2.3999868000487998e-12, -2.3999999868000000e-18,
      # (nu, nu)
      0, 0.00053838172231535815, 1.8585117094458687e-06, 5.1508169508281362e-08,
      3.9283900527873613e-09, 5.3021210764429046e-10, 5.3333020001214663e-19,
      5.3333333020000001e-28
    ),
    dim = c(m, 5, 5)
  )

  gh <- loglogistic5_gradient_hessian(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 5)
  expect_length(gh$H, m * 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.89442719099991588, 0.98473192783466186, 0.99540227449679584,
      0.99805257848288855, 0.99900149750436714, 0.9999990000015, 0.999999999,
      # log_eta
      0, -0.14879272455785841, -0.049641029773524029, -0.019635098543348011,
      -0.0096905462683797348, -0.0055096669605563691, -0.000011052375289228973,
      -1.6578612619821291e-08,
      # log_phi
      0, 0.21466252583997981, 0.035808433739442250, 0.010958557150423440,
      0.0046601676816321644, 0.0023928179580942925, 2.399992800018e-06,
      2.3999999928000000e-09,
      # log_nu
      0, -0.0082800886366924173, -0.00018458932481893295,
      -0.000016859398403064427, -3.0300205590351678e-06,
      -7.9707477980749869e-07, -7.9999706667479998e-13, -7.9999999706666667e-19
    ),
    nrow = m,
    ncol = 5
  )

  G <- loglogistic5_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 5)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_hessian <- array(
    c(
      # (alpha, alpha)
      rep(0, m),
      # (alpha, delta)
      rep(0, m),
      # (alpha, log_eta)
      rep(0, m),
      # (alpha, log_phi)
      rep(0, m),
      # (alpha, log_nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      0, 0.18599090569732301, 0.062051287216905036, 0.024543873179185013,
      0.012113182835474669, 0.0068870837006954614, 0.000013815469111536216,
      2.0723265774776614e-08,
      # (delta, log_phi)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (delta, log_nu)
      0, 0.010350110795865522, 0.00023073665602366618, 0.000021074248003830534,
      3.7875256987939597e-06, 9.9634347475937337e-07, 9.9999633334349997e-13,
      9.9999999633333334e-19,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, 0.18599090569732301, 0.062051287216905036, 0.024543873179185013,
      0.012113182835474669, 0.0068870837006954614, 0.000013815469111536216,
      2.0723265774776614e-08,
      # (log_eta, log_eta)
      0, 0.067791316223876984, 0.14742606196642691, 0.084456580654147461,
      0.050409389948811517, 0.032435813679348329, 0.00014164137412799102,
      3.2698438287780092e-07,
      # (log_eta, log_phi)
      0, -0.097802195731522854, -0.10634542424837656, -0.047136115145262187,
      -0.024241792297796111, -0.014086694896982991, -0.000030757033596490258,
      -4.7335837717456360e-08,
      # (log_eta, log_nu)
      0, 0.028036748883578745, 0.0014926420322225729, 0.00017972281250685084,
      0.000037669630904458492, 0.000010991844235124530, 2.2104695316699392e-11,
      3.3157225156749519e-17,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.999999991e-09,
      # (log_phi, log_eta)
      0, -0.097802195731522854, -0.10634542424837656, -0.047136115145262187,
      -0.024241792297796111, -0.014086694896982991, -0.000030757033596490258,
      -4.733583771745636e-08,
      # (log_phi, log_phi)
      0, 0.45079130426395760, 0.10254233298113008, 0.032423253954234491,
      0.013898904778097487, 0.0071569614974137672, 7.1999568001619995e-06,
      7.1999999568000002e-09,
      # (log_phi, log_nu)
      0, -0.040448478576988237, -0.0010767136288549652, -0.00010030521149375086,
      -0.000018115263232659968, -4.7736972972546533e-06,
      -4.7999736000975997e-12, -4.7999999736000001e-18,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0, 0.010350110795865522, 0.00023073665602366618, 0.000021074248003830534,
      3.7875256987939597e-06, 9.9634347475937337e-07, 9.9999633334349997e-13,
      9.9999999633333334e-19,
      # (log_nu, log_eta)
      0, 0.028036748883578745, 0.0014926420322225729, 0.00017972281250685084,
      0.000037669630904458492, 0.000010991844235124530, 2.2104695316699392e-11,
      3.3157225156749519e-17,
      # (log_nu, log_phi)
      0, -0.040448478576988237, -0.0010767136288549652, -0.00010030521149375086,
      -0.000018115263232659968, -4.7736972972546533e-06,
      -4.7999736000975997e-12, -4.7999999736000001e-18,
      # (log_nu, log_nu)
      0, -0.0061265617474309847, -0.00017715527798114947,
      -0.000016653365725031302, -3.0143069988240183e-06,
      -7.9495393137692153e-07, -7.9999493335399993e-13, -7.9999999493333335e-19
    ),
    dim = c(m, 5, 5)
  )

  H <- loglogistic5_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 5 * 5)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_5

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.89442719099991588, 0.98473192783466186, 0.99540227449679584,
      0.99805257848288855, 0.99900149750436714, 0.9999990000015, 0.999999999,
      # log_eta
      0, -0.14879272455785841, -0.049641029773524029, -0.019635098543348011,
      -0.0096905462683797348, -0.0055096669605563691, -0.000011052375289228973,
      -1.6578612619821291e-08,
      # log_phi
      0, 0.21466252583997981, 0.035808433739442250, 0.010958557150423440,
      0.0046601676816321644, 0.0023928179580942925, 2.399992800018e-06,
      2.3999999928000000e-09,
      # log_nu
      0, -0.0082800886366924173, -0.00018458932481893295,
      -0.000016859398403064427, -3.0300205590351678e-06,
      -7.9707477980749869e-07, -7.9999706667479998e-13, -7.9999999706666667e-19
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
      # (alpha, log_phi)
      rep(0, m),
      # (alpha, log_nu)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      0, 0.18599090569732301, 0.062051287216905036, 0.024543873179185013,
      0.012113182835474669, 0.0068870837006954614, 0.000013815469111536216,
      2.0723265774776614e-08,
      # (delta, log_phi)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.9999999910000000e-09,
      # (delta, log_nu)
      0, 0.010350110795865522, 0.00023073665602366618, 0.000021074248003830534,
      3.7875256987939597e-06, 9.9634347475937337e-07, 9.9999633334349997e-13,
      9.9999999633333334e-19,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, 0.18599090569732301, 0.062051287216905036, 0.024543873179185013,
      0.012113182835474669, 0.0068870837006954614, 0.000013815469111536216,
      2.0723265774776614e-08,
      # (log_eta, log_eta)
      0, 0.067791316223876984, 0.14742606196642691, 0.084456580654147461,
      0.050409389948811517, 0.032435813679348329, 0.00014164137412799102,
      3.2698438287780092e-07,
      # (log_eta, log_phi)
      0, -0.097802195731522854, -0.10634542424837656, -0.047136115145262187,
      -0.024241792297796111, -0.014086694896982991, -0.000030757033596490258,
      -4.7335837717456360e-08,
      # (log_eta, log_nu)
      0, 0.028036748883578745, 0.0014926420322225729, 0.00017972281250685084,
      0.000037669630904458492, 0.000010991844235124530, 2.2104695316699392e-11,
      3.3157225156749519e-17,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.26832815729997476, -0.044760542174302812, -0.013698196438029301,
      -0.0058252096020402056, -0.0029910224476178657, -2.9999910000224999e-06,
      -2.999999991e-09,
      # (log_phi, log_eta)
      0, -0.097802195731522854, -0.10634542424837656, -0.047136115145262187,
      -0.024241792297796111, -0.014086694896982991, -0.000030757033596490258,
      -4.733583771745636e-08,
      # (log_phi, log_phi)
      0, 0.45079130426395760, 0.10254233298113008, 0.032423253954234491,
      0.013898904778097487, 0.0071569614974137672, 7.1999568001619995e-06,
      7.1999999568000002e-09,
      # (log_phi, log_nu)
      0, -0.040448478576988237, -0.0010767136288549652, -0.00010030521149375086,
      -0.000018115263232659968, -4.7736972972546533e-06,
      -4.7999736000975997e-12, -4.7999999736000001e-18,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0, 0.010350110795865522, 0.00023073665602366618, 0.000021074248003830534,
      3.7875256987939597e-06, 9.9634347475937337e-07, 9.9999633334349997e-13,
      9.9999999633333334e-19,
      # (log_nu, log_eta)
      0, 0.028036748883578745, 0.0014926420322225729, 0.00017972281250685084,
      0.000037669630904458492, 0.000010991844235124530, 2.2104695316699392e-11,
      3.3157225156749519e-17,
      # (log_nu, log_phi)
      0, -0.040448478576988237, -0.0010767136288549652, -0.00010030521149375086,
      -0.000018115263232659968, -4.7736972972546533e-06,
      -4.7999736000975997e-12, -4.7999999736000001e-18,
      # (log_nu, log_nu)
      0, -0.0061265617474309847, -0.00017715527798114947,
      -0.000016653365725031302, -3.0143069988240183e-06,
      -7.9495393137692153e-07, -7.9999493335399993e-13, -7.9999999493333335e-19
    ),
    dim = c(m, 5, 5)
  )

  gh <- loglogistic5_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 5)
  expect_length(gh$H, m * 5 * 5)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(
    list(stats = matrix(x, nrow = m, ncol = 1)),
    class = "loglogistic5"
  )

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
  theta <- lltd$theta_5
  theta[3:5] <- log(theta[3:5])

  true_value <- 2.6013078656623450

  object <- structure(
    list(stats = lltd$stats_1, m = 5),
    class = "loglogistic5"
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

  value <- rss_fn(theta[2:3])

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  theta <- lltd$theta_5
  theta[3:5] <- log(theta[3:5])

  true_gradient <- c(
    -5.1445094078898300, -5.0686912803100457, 0.37798679446927383,
    -0.47368011281965922, 0.016254951021867376
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 15.627511758612288, -0.65880405448890199, 0.77581679698910851,
      -0.025286769291384284,
      # delta
      15.627511758612288, 15.289023529960513, -1.0822230809079753,
      1.2986047560701123, -0.042977036216411158,
      # log_eta
      -0.65880405448890199, -1.0822230809079753, -0.38261378223345501,
      0.29213889647528679, -0.052438543155035400,
      # log_phi
      0.77581679698910851, 1.2986047560701123, 0.29213889647528679,
      -0.89971490664874683, 0.074271976787634170,
      # log_nu
      -0.025286769291384284, -0.042977036216411158, -0.052438543155035400,
      0.074271976787634170, 0.012282296322732956
    ),
    nrow = 5,
    ncol = 5
  )

  object <- structure(
    list(stats = lltd$stats_1, m = 5),
    class = "loglogistic5"
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

  gh <- rss_gh(theta[2:3])

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 2)
  expect_length(gh$H, 2 * 2)

  expect_equal(gh$G, true_gradient[2:3])
  expect_equal(gh$H, true_hessian[2:3, 2:3])
})

test_that("mle_asy", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(
    0, 1, 3.3501562081181829, 2.2334446013195606, 3.0041670128832310
  )

  true_value <- c(
    0.86567490767345473, -0.77298635954637814, 3.3501562081181829,
    2.2334446013195606, 3.0041670128832310
  )

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 5)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.86567490767345473, delta = -0.77298635954637814,
    eta = 28.50718634991821, phi = 9.331955646606041, nu = 20.169408239379021
  )

  rss_value <- 0.058276476180507353

  fitted_values <- rep(
    c(
      0.86567490767345473, 0.7901639477306724, 0.6645425346912730,
      0.508922320859151, 0.329951365642644, 0.142289339659307,
      0.092688548127077, 0.092688548127077
    ),
    k
  )

  residuals <- c(
    -0.01297490767345473, -0.10857490767345473, 0.07132509232654527,
    0.0373360522693276, -0.0123639477306724, 0.0944360522693276,
    -0.1084425346912730, 0.0402574653087270, 0.036377679140849,
    -0.022922320859151, 0.026177679140849, -0.069022320859151,
    0.025348634357356, -0.016751365642644, 0.019948634357356,
    -0.000189339659307, -0.075888548127077, 0.045911451872923,
    0.030011451872923
  )

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic5_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fit_constrained: inequalities", {
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.058276476180507351

  fitted_values <- rep(
    c(
      0.86567490772243801, 0.7901639477115967, 0.6645425346486514,
      0.508922320838699, 0.329951365685670, 0.142289339645596,
      0.092688548125015, 0.092688548125015
    ),
    k
  )

  residuals <- c(
    -0.01297490772243801, -0.10857490772243801, 0.07132509227756199,
    0.0373360522884033, -0.0123639477115967, 0.0944360522884033,
    -0.1084425346486514, 0.0402574653513486, 0.036377679161301,
    -0.022922320838699, 0.026177679161301, -0.069022320838699,
    0.025348634314330, -0.016751365685670, 0.019948634314330,
    -0.000189339645596, -0.075888548125015, 0.045911451874985, 0.030011451874985
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 25, 8, 15, 20), c(1, -0.5, 30, 12, 30, 40)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.7, -0.6, 29, 11, 16, 38), max_iter,
    c(0.5, -1, 25, 8, 15, 20), c(1, -0.5, 30, 12, 30, 40)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, 2, 1, 1, 1, 1), max_iter,
    c(0.5, -1, 25, 8, 15, 20), c(1, -0.5, 30, 12, 30, 40)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.074747066165150929

  fitted_values <- rep(
    c(
      0.8, 0.77225113247950056, 0.6741844832087593, 0.5057443762603047,
      0.317532596086360, 0.190181348756377, 0.088632229933508, 0.088631934629668
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.1370, 0.05524886752049944, 0.00554886752049944,
    0.11234886752049944, -0.1180844832087593, 0.0306155167912407,
    0.0395556237396953, -0.0197443762603047, 0.0293556237396953,
    -0.0658443762603047, 0.037767403913640, -0.004332596086360,
    0.032367403913640, -0.048081348756377, -0.071832229933508,
    0.049967770066492, 0.034068065370332
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 4)), c(0.8, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(0.8, -0.9, 1, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 4)), c(0.8, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 4)), c(0.8, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.074747066165150929

  fitted_values <- rep(
    c(
      0.8, 0.77225113247950056, 0.6741844832087593, 0.5057443762603047,
      0.317532596086360, 0.190181348756377, 0.088632229933508, 0.088631934629668
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.1370, 0.05524886752049944, 0.00554886752049944,
    0.11234886752049944, -0.1180844832087593, 0.0306155167912407,
    0.0395556237396953, -0.0197443762603047, 0.0293556237396953,
    -0.0658443762603047, 0.037767403913640, -0.004332596086360,
    0.032367403913640, -0.048081348756377, -0.071832229933508,
    0.049967770066492, 0.034068065370332
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 4, 5, 2, 1), c(0.8, -0.9, 8, 10, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.8, -0.9, 7, 9, 2.1, 1.2), max_iter,
    c(0.8, -0.9, 4, 5, 2, 1), c(0.8, -0.9, 8, 10, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10, 5), max_iter,
    c(0.8, -0.9, 4, 5, 2, 1), c(0.8, -0.9, 8, 10, 3, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  # loglogistic5 model is basically unidentifiable: many parameters are
  # associated with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022085036915158753

  fitted_values <- rep(
    c(
      0.90681779366690136, 0.82534227701819370, 0.68768283546185108,
      0.51593867110192256, 0.31929217045238740, 0.14655690461413620,
      0.093029886494351149, 0.093029886494351068
    ),
    k
  )

  residuals <- c(
    -0.054117793666901360, -0.14971779366690136, 0.030182206333098640,
    0.0021577229818062970, -0.047542277018193703, 0.059257722981806297,
    -0.13158283546185108, 0.017117164538148918, 0.029361328898077444,
    -0.029938671101922556, 0.019161328898077444, -0.076038671101922556,
    0.036007829547612596, -0.0060921704523874042, 0.030607829547612596,
    -0.0044569046141361993, -0.076229886494351149, 0.045570113505648851,
    0.029670113505648932
  )

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic5_new(
    x, y, w, c(1, -1, 1, 1, 1, 1), max_iter, NULL, NULL
  )

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022085036915158753

  fitted_values <- rep(
    c(
      0.90681779366690136, 0.82534227701819370, 0.68768283546185108,
      0.51593867110192256, 0.31929217045238740, 0.14655690461413620,
      0.093029886494351149, 0.093029886494351068
    ),
    k
  )

  residuals <- c(
    -0.054117793666901360, -0.14971779366690136, 0.030182206333098640,
    0.0021577229818062970, -0.047542277018193703, 0.059257722981806297,
    -0.13158283546185108, 0.017117164538148918, 0.029361328898077444,
    -0.029938671101922556, 0.019161328898077444, -0.076038671101922556,
    0.036007829547612596, -0.0060921704523874042, 0.030607829547612596,
    -0.0044569046141361993, -0.076229886494351149, 0.045570113505648851,
    0.029670113505648932
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 10, 8, 5, 0), c(1, -0.5, 20, 10, 15, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.7, -0.6, 18, 8.5, 7, 1.8), max_iter,
    c(0.5, -1, 10, 8, 5, 0), c(1, -0.5, 20, 10, 15, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, -5, 0.5, 20, 0.1, 3), 10000,
    c(0.5, -1, 10, 8, 5, 0), c(1, -0.5, 20, 10, 15, 2)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022165226781174076

  fitted_values <- rep(
    c(
      0.9, 0.82394807262993030, 0.68895716760665773, 0.51663843729538372,
      0.31768368751798345, 0.14888199437111868, 0.092568060779228249,
      0.092568060779225243
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, 0.0035519273700696982, -0.046148072629930302,
    0.060651927370069698, -0.13285716760665773, 0.015842832393342270,
    0.028661562704616279, -0.030638437295383721, 0.018461562704616279,
    -0.076738437295383721, 0.037616312482016546, -0.0044836875179834541,
    0.032216312482016546, -0.0067819943711186764, -0.075768060779228249,
    0.046031939220771751, 0.030131939220774757
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 4)), c(0.9, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(0.9, -0.9, 1, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 4)), c(0.9, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 4)), c(0.9, -0.9, rep(Inf, 4))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss_value <- 0.022165226781174076

  fitted_values <- rep(
    c(
      0.9, 0.82394807262993030, 0.68895716760665773, 0.51663843729538372,
      0.31768368751798345, 0.14888199437111868, 0.092568060779228249,
      0.092568060779225243
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, 0.0035519273700696982, -0.046148072629930302,
    0.060651927370069698, -0.13285716760665773, 0.015842832393342270,
    0.028661562704616279, -0.030638437295383721, 0.018461562704616279,
    -0.076738437295383721, 0.037616312482016546, -0.0044836875179834541,
    0.032216312482016546, -0.0067819943711186764, -0.075768060779228249,
    0.046031939220771751, 0.030131939220774757
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 5, 5, 7, 0), c(0.9, -0.9, 15, 10, 12, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.9, -0.9, 6, 7, 11, 0.5), max_iter,
    c(0.9, -0.9, 5, 5, 7, 0), c(0.9, -0.9, 15, 10, 12, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 5, 6), max_iter,
    c(0.9, -0.9, 5, 5, 7, 0), c(0.9, -0.9, 15, 10, 12, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  max_iter <- 10000

  theta <- lltd$theta_6
  names(theta) <- c("alpha", "delta", "eta", "phi", "nu", "xi")

  sigma <- lltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.9600000000000, 4033.6259327887426, -78.699146035632538,
      198.77058917560704, -293.58875257662276, 269.87863134052620,
      89181.461311549703,
      # delta
      4033.6259327887426, 2878.4597528145919, -88.705348162354274,
      244.00408263047195, -320.80516117638226, 282.55798267636782,
      59530.094792948388,
      # eta
      -78.699146035632538, -88.705348162354274, -24.139034413674205,
      -33.757408298046454, 4.9250747626882528, -15.028191756557660,
      -1242.7147074919567,
      # phi
      198.77058917560704, 244.00408263047195, -33.757408298046454,
      21.182660748584033, -6.5834237992925952, 32.161335478721606,
      4533.2517221704692,
      # nu
      -293.58875257662276, -320.80516117638226, 4.9250747626882528,
      -6.5834237992925952, -29.341061059068340, 6.4659688730596288,
      -4510.5256202602449,
      # xi
      269.87863134052620, 282.55798267636782, -15.028191756557660,
      32.161335478721606, 6.4659688730596288, 41.574031108086305,
      3449.7573253555426,
      # sigma
      89181.461311549703, 59530.094792948388, -1242.7147074919567,
      4533.2517221704692, -4510.5256202602449, 3449.7573253555426,
      1.3580467956656330e+06
    ),
    nrow = 7,
    ncol = 7
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "xi", "sigma"
  )

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 7 * 7)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = matrix(-Inf, nrow = 6, ncol = 2),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      upper_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = matrix(Inf, nrow = 6, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' must be of length 6"
  )
})

test_that("drda: 'start' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c("a", "b", "c", "d", "e", "f")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, Inf, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(-Inf, 1, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(1, 1, 1, 1, 1, 1, 1)
    ),
    "'start' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, -1, 1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 0, 1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, -1, 1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 0, 1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, -1, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, 0, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, 1, -1)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, 1, 0)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.57634025117613155)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.87325260264665631)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.61191350911591019)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.41629510891428020
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(naac(result), 1 - 0.57634025117613155)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.87325260264665631)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.61191350911591019)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.41629510891428020
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.48956965931169600)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.17502499958193786)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.49142869302873806)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.77685613973089947
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, mean_function = "loglogistic5")

  expect_equal(naac(result), 1 - 0.48956965931169600)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.17502499958193786)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.49142869302873806)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.77685613973089947
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
