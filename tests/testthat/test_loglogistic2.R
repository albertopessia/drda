test_that("Constructor (decreasing)", {
  x <- lltd$D$x
  y <- lltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- lltd$stats_1

  start <- c(1, 1)

  lower_bound <- c(0.5, 1)
  upper_bound <- c(2, 5)

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
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

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, log(start)))
  expect_equal(object$lower_bound, log(lower_bound))
  expect_equal(object$upper_bound, log(upper_bound))

  w <- lltd$D$w
  stats <- lltd$stats_2

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
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

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(1, -1, log(start)))
  expect_equal(object$lower_bound, log(lower_bound))
  expect_equal(object$upper_bound, log(upper_bound))
})

test_that("Constructor (increasing)", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- lltd$stats_1_i

  start <- c(1, 1)

  lower_bound <- c(0.5, 1)
  upper_bound <- c(2, 5)

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
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

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, log(start)))
  expect_equal(object$lower_bound, log(lower_bound))
  expect_equal(object$upper_bound, log(upper_bound))

  w <- lltd$D$w
  stats <- lltd$stats_2_i

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic2"))
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

  object <- loglogistic2_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic2"))
  expect_equal(object$x, x)
  expect_equal(object$y, y)
  expect_equal(object$w, w)
  expect_equal(object$n, n)
  expect_equal(object$m, m)
  expect_equal(object$stats, stats)
  expect_true(object$constrained)
  expect_equal(object$max_iter, max_iter)
  expect_equal(object$start, c(0, 1, log(start)))
  expect_equal(object$lower_bound, log(lower_bound))
  expect_equal(object$upper_bound, log(upper_bound))
})

test_that("Constructor: errors", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w
  max_iter <- 10000

  expect_error(
    loglogistic2_new(x, y, w, 1, max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(-1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, c(1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, -Inf, Inf),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, -Inf, rep(Inf, 2)),
    "'lower_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), Inf),
    "'upper_bound' must be of length 2"
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(0, Inf)),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(-1, Inf)),
    "'upper_bound[1]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(Inf, 0)),
    "'upper_bound[2]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic2_new(x, y, w, NULL, max_iter, rep(-Inf, 2), c(Inf, -1)),
    "'upper_bound[2]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value (decreasing)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d
  theta[2] <- -0.5

  m <- length(x)

  true_value <- c(
    1, 0.86206896551724138, 0.60975609756097561, 0.40983606557377049,
    0.28089887640449438, 0.20000000000000000, 0.0024937655860349127,
    0.000024999375015624609
  )

  value <- loglogistic2_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic2")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic2_fit")

  value <- fn(object, object$stats[, 1], c(1, -1, theta[3:4]))

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Function value (increasing)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_i
  theta[2] <- 0.5

  m <- length(x)

  true_value <- c(
    0, 0.13793103448275862, 0.39024390243902439, 0.59016393442622951,
    0.71910112359550562, 0.80000000000000000, 0.99750623441396509,
    0.99997500062498438
  )

  value <- loglogistic2_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic2")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic2_fit")

  value <- fn(object, object$stats[, 1], c(0, 1, theta[3:4]))

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # eta
      0, 0.10895252459859157, 0.053097811139609698, -0.044098199708293245,
      -0.094938240978812888, -0.11090354888959125, -0.0074520239887910921,
      -0.00013245131151534308,
      # phi
      0, 0.047562425683709869, 0.095181439619274242, 0.096748185971513034,
      0.080797879055674789, 0.064, 0.00099501868769472827,
      9.9995000187493750e-06
    ),
    nrow = m,
    ncol = 2
  )

  G <- loglogistic2_gradient(x, theta, -1)

  expect_type(G, "double")
  expect_length(G, m * 2)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_hessian <- array(
    c(
      # (eta, eta)
      0, -0.072292274433916804, -0.0026008757878617911, 0.0014498455188343302,
      0.019553161738827894, 0.046123489336147337, 0.022212925780153336,
      0.00070173399644679443,
      # (eta, phi)
      0, -0.0077774494556681395, 0.042928472977866513, 0.045193239236305857,
      0.023758079311326371, 0.0053831482664981001, -0.0024684333708934296,
      -0.000047978125636756238,
      # (phi, eta)
      0, -0.0077774494556681395, 0.042928472977866513, 0.045193239236305857,
      0.023758079311326371, 0.0053831482664981001, -0.0024684333708934296,
      -0.000047978125636756238,
      # (phi, phi)
      0, -0.023289187748575177, -0.027393682622132587, -0.01237107951766888,
      -0.0019972509429492644, 0.00256, 0.00019701866285027787,
      1.9997000187491250e-06
    ),
    dim = c(m, 2, 2)
  )

  H <- loglogistic2_hessian(x, theta, -1)

  expect_type(H, "double")
  expect_length(H, m * 2 * 2)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # eta
      0, 0.10895252459859157, 0.053097811139609698, -0.044098199708293245,
      -0.094938240978812888, -0.11090354888959125, -0.0074520239887910921,
      -0.00013245131151534308,
      # phi
      0, 0.047562425683709869, 0.095181439619274242, 0.096748185971513034,
      0.080797879055674789, 0.064, 0.00099501868769472827,
      9.9995000187493750e-06
    ),
    nrow = m,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (eta, eta)
      0, -0.072292274433916804, -0.0026008757878617911, 0.0014498455188343302,
      0.019553161738827894, 0.046123489336147337, 0.022212925780153336,
      0.00070173399644679443,
      # (eta, phi)
      0, -0.0077774494556681395, 0.042928472977866513, 0.045193239236305857,
      0.023758079311326371, 0.0053831482664981001, -0.0024684333708934296,
      -0.000047978125636756238,
      # (phi, eta)
      0, -0.0077774494556681395, 0.042928472977866513, 0.045193239236305857,
      0.023758079311326371, 0.0053831482664981001, -0.0024684333708934296,
      -0.000047978125636756238,
      # (phi, phi)
      0, -0.023289187748575177, -0.027393682622132587, -0.01237107951766888,
      -0.0019972509429492644, 0.00256, 0.00019701866285027787,
      1.9997000187491250e-06
    ),
    dim = c(m, 2, 2)
  )

  gh <- loglogistic2_gradient_hessian(x, theta, -1)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 2)
  expect_length(gh$H, m * 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # log_eta
      0, 0.21790504919718313, 0.10619562227921940, -0.088196399416586489,
      -0.18987648195762578, -0.22180709777918250, -0.014904047977582184,
      -0.00026490262303068616,
      # log_phi
      0, 0.23781212841854935, 0.47590719809637121, 0.48374092985756517,
      0.40398939527837394, 0.32000000000000000, 0.0049750934384736413,
      0.000049997500093746875
    ),
    nrow = m,
    ncol = 2
  )

  G <- loglogistic2_gradient_2(x, theta, -1)

  expect_type(G, "double")
  expect_length(G, m * 2)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_hessian <- array(
    c(
      # (log_eta, log_eta)
      0, -0.071264048538484085, 0.095792119127772232, -0.082397017341249169,
      -0.11166383500231420, -0.037313140434593152, 0.073947655143031158,
      0.0025420333627564916,
      # (log_eta, log_phi)
      0, -0.077774494556681395, 0.42928472977866513, 0.45193239236305857,
      0.23758079311326371, 0.053831482664981001, -0.024684333708934296,
      -0.00047978125636756238,
      # (log_phi, log_eta)
      0, -0.077774494556681395, 0.42928472977866513, 0.45193239236305857,
      0.23758079311326371, 0.053831482664981001, -0.024684333708934296,
      -0.00047978125636756238,
      # (log_phi, log_phi)
      0, -0.34441756529583009, -0.20893486745694346, 0.17446394191584318,
      0.35405812170464233, 0.384, 0.009900560009730588, 0.000099990000562475001
    ),
    dim = c(m, 2, 2)
  )

  H <- loglogistic2_hessian_2(x, theta, -1)

  expect_type(H, "double")
  expect_length(H, m * 2 * 2)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_2_d[3:4]

  m <- length(x)

  true_gradient <- matrix(
    c(
      # log_eta
      0, 0.21790504919718313, 0.10619562227921940, -0.088196399416586489,
      -0.18987648195762578, -0.22180709777918250, -0.014904047977582184,
      -0.00026490262303068616,
      # log_phi
      0, 0.23781212841854935, 0.47590719809637121, 0.48374092985756517,
      0.40398939527837394, 0.32000000000000000, 0.0049750934384736413,
      0.000049997500093746875
    ),
    nrow = m,
    ncol = 2
  )

  true_hessian <- array(
    c(
      # (log_eta, log_eta)
      0, -0.071264048538484085, 0.095792119127772232, -0.082397017341249169,
      -0.11166383500231420, -0.037313140434593152, 0.073947655143031158,
      0.0025420333627564916,
      # (log_eta, log_phi)
      0, -0.077774494556681395, 0.42928472977866513, 0.45193239236305857,
      0.23758079311326371, 0.053831482664981001, -0.024684333708934296,
      -0.00047978125636756238,
      # (log_phi, log_eta)
      0, -0.077774494556681395, 0.42928472977866513, 0.45193239236305857,
      0.23758079311326371, 0.053831482664981001, -0.024684333708934296,
      -0.00047978125636756238,
      # (log_phi, log_phi)
      0, -0.34441756529583009, -0.20893486745694346, 0.17446394191584318,
      0.35405812170464233, 0.384, 0.009900560009730588, 0.000099990000562475001
    ),
    dim = c(m, 2, 2)
  )

  gh <- loglogistic2_gradient_hessian_2(x, theta, -1)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 2)
  expect_length(gh$H, m * 2 * 2)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(
    list(stats = lltd$stats_1, start = c(1, -1, NA_real_, NA_real_)),
    class = "loglogistic2"
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
  theta <- log(lltd$theta_2_d[3:4])

  true_value <- 0.14607976935420221

  object <- structure(
    list(
      stats = lltd$stats_1,
      m = nrow(lltd$stats_1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
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
  theta <- log(lltd$theta_2_d[3:4])

  true_gradient <- c(0.071748288653048041, -0.22751406466277087)

  true_hessian <- matrix(
    c(
      # log_eta
      0.37935222818782687, -0.44131953032097426,
      # log_phi
      -0.44131953032097426, 2.0207120814395956
    ),
    nrow = 2,
    ncol = 2
  )

  object <- structure(
    list(
      stats = lltd$stats_1,
      m = nrow(lltd$stats_1),
      start = c(1, -1, NA_real_, NA_real_)
    ),
    class = "loglogistic2"
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(0.54763114884776605, 1.7084450404687743)

  true_value <- c(0.54763114884776605, 1.7084450404687743)

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 2)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.16419461331342958

  theta <- c(
    alpha = 1, delta = -1, eta = exp(0.54763114884776605),
    phi = exp(1.7084450404687743)
  )

  fitted_values <- rep(
    c(
      1, 0.85265828525963440, 0.6357691163194841, 0.4640463320849631,
      0.3449054253807791, 0.2635965354058204, 0.0066340855843002,
      0.000124584968709
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.02515828525963440, -0.07485828525963440,
    0.03194171474036560, -0.0796691163194841, 0.0690308836805159,
    0.0812536679150369, 0.0219536679150369, 0.0710536679150369,
    -0.0241463320849631, 0.0103945746192209, -0.0317054253807791,
    0.0049945746192209, -0.1214965354058204, 0.0101659144156998,
    0.1319659144156998, 0.122575415031291
  )

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic2_new(x, y, w, c(1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.16419461331342958

  theta <- c(
    alpha = 1, delta = -1, eta = exp(0.54763114884776605),
    phi = exp(1.7084450404687743)
  )

  fitted_values <- rep(
    c(
      1, 0.85265828525963440, 0.6357691163194841, 0.4640463320849631,
      0.3449054253807791, 0.2635965354058204, 0.0066340855843002,
      0.000124584968709
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.02515828525963440, -0.07485828525963440,
    0.03194171474036560, -0.0796691163194841, 0.0690308836805159,
    0.0812536679150369, 0.0219536679150369, 0.0710536679150369,
    -0.0241463320849631, 0.0103945746192209, -0.0317054253807791,
    0.0049945746192209, -0.1214965354058204, 0.0101659144156998,
    0.1319659144156998, 0.122575415031291
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(
    x, y, w, c(4, 8), max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(
    x, y, w, c(7, 1), max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.16968996645062665

  theta <- c(
    alpha = 1, delta = -1, eta = 2, phi = exp(1.7176599758908608)
  )

  fitted_values <- rep(
    c(
      1, 0.88584912304501616, 0.65987365414469592, 0.46301791484323224,
      0.32660884016084947, 0.23688205621120481, 0.003094528439196051,
      0.00003104037917471891
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05834912304501616, -0.10804912304501616,
    -0.00124912304501616, -0.10377365414469592, 0.04492634585530408,
    0.08228208515676776, 0.02298208515676776, 0.07208208515676776,
    -0.02311791484323224, 0.02869115983915053, -0.01340884016084947,
    0.02329115983915053, -0.09478205621120481, 0.013705471560803949,
    0.135505471560803949, 0.12266895962082528109
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic2_new(
    x, y, w, c(2, 1), max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic2_new(
    x, y, w, c(1, 1), max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.16968996645062665

  theta <- c(
    alpha = 1, delta = -1, eta = 2, phi = exp(1.7176599758908608)
  )

  fitted_values <- rep(
    c(
      1, 0.88584912304501616, 0.65987365414469592, 0.46301791484323224,
      0.32660884016084947, 0.23688205621120481, 0.003094528439196051,
      0.00003104037917471891
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05834912304501616, -0.10804912304501616,
    -0.00124912304501616, -0.10377365414469592, 0.04492634585530408,
    0.08228208515676776, 0.02298208515676776, 0.07208208515676776,
    -0.02311791484323224, 0.02869115983915053, -0.01340884016084947,
    0.02329115983915053, -0.09478205621120481, 0.013705471560803949,
    0.135505471560803949, 0.12266895962082528109
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(
    x, y, w, c(2, 7), max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(
    x, y, w, c(8, 1), max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.080570305967227270

  theta <- c(
    alpha = 1, delta = -1, eta = exp(0.66309884142824792),
    phi = exp(1.7153886937325806)
  )

  fitted_values <- rep(
    c(
      1, 0.87910304144164237, 0.6544607502962267, 0.4630126427820133,
      0.3303609999355391, 0.2423884549177716, 0.003653243406751,
      0.000042019666396
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05160304144164237, -0.10130304144164237,
    0.00549695855835763, -0.0983607502962267, 0.0503392497037733,
    0.0822873572179867, 0.0229873572179867, 0.0720873572179867,
    -0.0231126427820133, 0.0249390000644609, -0.0171609999355391,
    0.0195390000644609, -0.1002884549177716, 0.013146756593249,
    0.134946756593249, 0.122657980333604
  )

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic2_new(x, y, w, c(1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.080570305967227270

  theta <- c(
    alpha = 1, delta = -1, eta = exp(0.66309884142824792),
    phi = exp(1.7153886937325806)
  )

  fitted_values <- rep(
    c(
      1, 0.87910304144164237, 0.6544607502962267, 0.4630126427820133,
      0.3303609999355391, 0.2423884549177716, 0.003653243406751,
      0.000042019666396
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05160304144164237, -0.10130304144164237,
    0.00549695855835763, -0.0983607502962267, 0.0503392497037733,
    0.0822873572179867, 0.0229873572179867, 0.0720873572179867,
    -0.0231126427820133, 0.0249390000644609, -0.0171609999355391,
    0.0195390000644609, -0.1002884549177716, 0.013146756593249,
    0.134946756593249, 0.122657980333604
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(
    x, y, w, c(4, 8), max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(
    x, y, w, c(7, 1), max_iter,
    c(1, 3), c(5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 2)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.080806739950738768

  theta <- c(
    alpha = 1, delta = -1, eta = 2, phi = exp(1.7179698864941254)
  )

  fitted_values <- rep(
    c(
      1, 0.88591178465897222, 0.66001275328511998, 0.46317202595877424,
      0.32674517549914561, 0.23699411891072304, 0.003096441147111687,
      0.00003105962402598627
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05841178465897222, -0.10811178465897222,
    -0.00131178465897222, -0.10391275328511998, 0.04478724671488002,
    0.08212797404122576, 0.02282797404122576, 0.07192797404122576,
    -0.02327202595877424, 0.02855482450085439, -0.01354517549914561,
    0.02315482450085439, -0.09489411891072304, 0.013703558852888313,
    0.13550355885288831, 0.12266894037597401
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic2_new(
    x, y, w, c(2, 1), max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic2_new(
    x, y, w, c(1, 1), max_iter,
    c(2, -Inf), c(2, Inf)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = FALSE, phi = TRUE)

  rss_value <- 0.080806739950738768

  theta <- c(
    alpha = 1, delta = -1, eta = 2, phi = exp(1.7179698864941254)
  )

  fitted_values <- rep(
    c(
      1, 0.88591178465897222, 0.66001275328511998, 0.46317202595877424,
      0.32674517549914561, 0.23699411891072304, 0.003096441147111687,
      0.00003105962402598627
    ),
    k
  )

  residuals <- c(
    -0.1473, -0.2429, -0.063, -0.05841178465897222, -0.10811178465897222,
    -0.00131178465897222, -0.10391275328511998, 0.04478724671488002,
    0.08212797404122576, 0.02282797404122576, 0.07192797404122576,
    -0.02327202595877424, 0.02855482450085439, -0.01354517549914561,
    0.02315482450085439, -0.09489411891072304, 0.013703558852888313,
    0.13550355885288831, 0.12266894037597401
  )

  object <- loglogistic2_new(
    x, y, w, NULL, max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic2_new(
    x, y, w, c(2, 7), max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic2_new(
    x, y, w, c(8, 1), max_iter,
    c(2, 3), c(2, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic2_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 1)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  max_iter <- 10000

  theta <- lltd$theta_2_d
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- lltd$sigma

  true_value <- matrix(c(
      # eta
      33.710156780065974, -14.980681134843766, -303.36552779699696,
      # phi
      -14.980681134843766, 30.964961924076123, 628.42233185708605,
      # sigma
      -303.36552779699696, 628.42233185708605, 42740.806385880433
    ),
    nrow = 3,
    ncol = 3
  )

  rownames(true_value) <- colnames(true_value) <- c("eta", "phi", "sigma")

  object <- loglogistic2_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 3 * 3)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c("c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = matrix(-Inf, nrow = 2, ncol = 2),
      upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep(Inf, 2)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c( 0, -Inf),
      upper_bound = c(-1, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c(Inf, -Inf),
      upper_bound = c(Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep(Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      upper_bound = c("c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 2),
      upper_bound = matrix(Inf, nrow = 2, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = c(-Inf, -Inf),
      upper_bound = c(-Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      lower_bound = rep(-Inf, 3),
      upper_bound = rep(Inf, 3)
    ),
    "'lower_bound' must be of length 2"
  )
})

test_that("drda: 'start' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c("c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(1, Inf)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(-Inf, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = rep(1, 3)
    ),
    "'start' must be of length 2"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(-1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(1, -1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic2",
      start = c(1, 0)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic2")

  expect_equal(nauc(result), 0.0089638310449171108)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.95676585686209365)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.0057288924161456609)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.32521005090585268
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic2")

  expect_equal(naac(result), 1 - 0.0089638310449171108)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.95676585686209365)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.0057288924161456609)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.32521005090585268
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic2")

  expect_equal(nauc(result), 0.99027624631427548)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.099216635018969768)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.99536629407136298)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.84693307096743735
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic2")

  expect_equal(naac(result), 1 - 0.99027624631427548)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.099216635018969768)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.99536629407136298)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.84693307096743735
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
