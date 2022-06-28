test_that("Constructor", {
  x <- lltd$D$x
  y <- lltd$D$y

  m <- length(unique(x))
  n <- length(y)

  w <- rep(1, n)

  max_iter <- 10000

  stats <- lltd$stats_1

  start <- c(0, 1, 1, 1)

  lower_bound <- c(0, -1, 0.5, 1)
  upper_bound <- c(3, 2, 2, 5)

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic4"))
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

  object <- loglogistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  i <- c(1, 2)

  expect_true(inherits(object, "loglogistic4"))
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

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loglogistic4"))
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

  object <- loglogistic4_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loglogistic4"))
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
    loglogistic4_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, c(0, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, 0)),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loglogistic4_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, -1)),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

  m <- length(x)

  true_value <- c(
    1, 0.88275862068965517, 0.66829268292682927, 0.49836065573770492,
    0.38876404494382022, 0.32, 0.15211970074812968, 0.15002124946876328
  )

  value <- loglogistic4_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic4")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic4_fit")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.13793103448275862, 0.39024390243902439, 0.59016393442622951,
      0.71910112359550562, 0.8, 0.99750623441396509, 0.99997500062498438,
      # eta
      0, 0.092609645908802831, 0.045133139468668243, -0.037483469752049258,
      -0.080697504831990955, -0.094268016556152562, -0.0063342203904724283,
      -0.00011258361478804162,
      # phi
      0, 0.040428061831153389, 0.080904223676383105, 0.082235958075786079,
      0.068678197197323570, 0.0544, 0.00084576588454051903,
      8.4995750159369688e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- loglogistic4_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

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
      0, -0.10895252459859157, -0.053097811139609698, 0.044098199708293245,
      0.094938240978812888, 0.11090354888959125, 0.0074520239887910921,
      0.00013245131151534308,
      # (delta, phi)
      0, -0.047562425683709869, -0.095181439619274242, -0.096748185971513034,
      -0.080797879055674789, -0.064, -0.00099501868769472827,
      -9.9995000187493750e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.10895252459859157, -0.053097811139609698, 0.044098199708293245,
      0.094938240978812888, 0.11090354888959125, 0.0074520239887910921,
      0.00013245131151534308,
      # (eta, eta)
      0, -0.061448433268829283, -0.0022107444196825224, 0.0012323686910091807,
      0.016620187478003710, 0.039204965935725236, 0.018880986913130335,
      0.00059647389697977526,
      # (eta, phi)
      0, -0.0066108320373179186, 0.036489202031186536, 0.038414253350859978,
      0.020194367414627416, 0.0045756760265233851, -0.0020981683652594152,
      -0.000040781406791242802,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.047562425683709869, -0.095181439619274242, -0.096748185971513034,
      -0.080797879055674789, -0.064, -0.00099501868769472827,
      -9.9995000187493750e-06,
      # (phi, eta)
      0, -0.0066108320373179186, 0.036489202031186536, 0.038414253350859978,
      0.020194367414627416, 0.0045756760265233851, -0.0020981683652594152,
      -0.000040781406791242802,
      # (phi, phi)
      0, -0.019795809586288901, -0.023284630228812699, -0.010515417590018548,
      -0.0016976633015068748, 0.002176, 0.00016746586342273619,
      1.6997450159367563e-06
    ),
    dim = c(m, 4, 4)
  )

  H <- loglogistic4_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.13793103448275862, 0.39024390243902439, 0.59016393442622951,
      0.71910112359550562, 0.8, 0.99750623441396509, 0.99997500062498438,
      # eta
      0, 0.092609645908802831, 0.045133139468668243, -0.037483469752049258,
      -0.080697504831990955, -0.094268016556152562, -0.0063342203904724283,
      -0.00011258361478804162,
      # phi
      0, 0.040428061831153389, 0.080904223676383105, 0.082235958075786079,
      0.068678197197323570, 0.0544, 0.00084576588454051903,
      8.4995750159369688e-06
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
      0, -0.10895252459859157, -0.053097811139609698, 0.044098199708293245,
      0.094938240978812888, 0.11090354888959125, 0.0074520239887910921,
      0.00013245131151534308,
      # (delta, phi)
      0, -0.047562425683709869, -0.095181439619274242, -0.096748185971513034,
      -0.080797879055674789, -0.064, -0.00099501868769472827,
      -9.9995000187493750e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.10895252459859157, -0.053097811139609698, 0.044098199708293245,
      0.094938240978812888, 0.11090354888959125, 0.0074520239887910921,
      0.00013245131151534308,
      # (eta, eta)
      0, -0.061448433268829283, -0.0022107444196825224, 0.0012323686910091807,
      0.016620187478003710, 0.039204965935725236, 0.018880986913130335,
      0.00059647389697977526,
      # (eta, phi)
      0, -0.0066108320373179186, 0.036489202031186536, 0.038414253350859978,
      0.020194367414627416, 0.0045756760265233851, -0.0020981683652594152,
      -0.000040781406791242802,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.047562425683709869, -0.095181439619274242, -0.096748185971513034,
      -0.080797879055674789, -0.064, -0.00099501868769472827,
      -9.9995000187493750e-06,
      # (phi, eta)
      0, -0.0066108320373179186, 0.036489202031186536, 0.038414253350859978,
      0.020194367414627416, 0.0045756760265233851, -0.0020981683652594152,
      -0.000040781406791242802,
      # (phi, phi)
      0, -0.019795809586288901, -0.023284630228812699, -0.010515417590018548,
      -0.0016976633015068748, 0.002176, 0.00016746586342273619,
      1.6997450159367563e-06
    ),
    dim = c(m, 4, 4)
  )

  gh <- loglogistic4_gradient_hessian(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)
})

test_that("Gradient (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.13793103448275862, 0.39024390243902439, 0.59016393442622951,
      0.71910112359550562, 0.8, 0.99750623441396509, 0.99997500062498438,
      # log_eta
      0, 0.18521929181760566, 0.090266278937336487, -0.074966939504098516,
      -0.16139500966398191, -0.18853603311230512, -0.012668440780944857,
      -0.00022516722957608324,
      # log_phi
      0, 0.20214030915576694, 0.40452111838191553, 0.41117979037893040,
      0.34339098598661785, 0.272, 0.0042288294227025951, 0.000042497875079684844
    ),
    nrow = m,
    ncol = 4
  )

  G <- loglogistic4_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

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
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      0, -0.21790504919718313, -0.10619562227921940, 0.088196399416586489,
      0.18987648195762578, 0.22180709777918250, 0.014904047977582184,
      0.00026490262303068616,
      # (delta, log_phi)
      0, -0.23781212841854935, -0.47590719809637121, -0.48374092985756517,
      -0.40398939527837394, -0.32, -0.0049750934384736413,
      -0.000049997500093746875,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.21790504919718313, -0.10619562227921940, 0.088196399416586489,
      0.18987648195762578, 0.22180709777918250, 0.014904047977582184,
      0.00026490262303068616,
      # (log_eta, log_eta)
      0, -0.060574441257711472, 0.081423301258606397, -0.070037464740061793,
      -0.094914259751967070, -0.031716169369404179, 0.062855506871576484,
      0.0021607283583430178,
      # (log_eta, log_phi)
      0, -0.066108320373179186, 0.36489202031186536, 0.38414253350859978,
      0.20194367414627416, 0.045756760265233851, -0.020981683652594152,
      -0.00040781406791242802,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.23781212841854935, -0.47590719809637121, -0.48374092985756517,
      -0.40398939527837394, -0.32, -0.0049750934384736413,
      -0.000049997500093746875,
      # (log_phi, log_eta)
      0, -0.066108320373179186, 0.36489202031186536, 0.38414253350859978,
      0.20194367414627416, 0.045756760265233851, -0.020981683652594152,
      -0.00040781406791242802,
      # (log_phi, log_phi)
      0, -0.29275493050145557, -0.17759463733840194, 0.14829435062846670,
      0.30094940344894598, 0.3264, 0.0084154760082709998,
      0.000084991500478103751
    ),
    dim = c(m, 4, 4)
  )

  H <- loglogistic4_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_4

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.13793103448275862, 0.39024390243902439, 0.59016393442622951,
      0.71910112359550562, 0.8, 0.99750623441396509, 0.99997500062498438,
      # log_eta
      0, 0.18521929181760566, 0.090266278937336487, -0.074966939504098516,
      -0.16139500966398191, -0.18853603311230512, -0.012668440780944857,
      -0.00022516722957608324,
      # log_phi
      0, 0.20214030915576694, 0.40452111838191553, 0.41117979037893040,
      0.34339098598661785, 0.272, 0.0042288294227025951, 0.000042497875079684844
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
      # (alpha, log_phi)
      rep(0, m),
      # (delta, alpha)
      rep(0, m),
      # (delta, delta)
      rep(0, m),
      # (delta, log_eta)
      0, -0.21790504919718313, -0.10619562227921940, 0.088196399416586489,
      0.18987648195762578, 0.22180709777918250, 0.014904047977582184,
      0.00026490262303068616,
      # (delta, log_phi)
      0, -0.23781212841854935, -0.47590719809637121, -0.48374092985756517,
      -0.40398939527837394, -0.32, -0.0049750934384736413,
      -0.000049997500093746875,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.21790504919718313, -0.10619562227921940, 0.088196399416586489,
      0.18987648195762578, 0.22180709777918250, 0.014904047977582184,
      0.00026490262303068616,
      # (log_eta, log_eta)
      0, -0.060574441257711472, 0.081423301258606397, -0.070037464740061793,
      -0.094914259751967070, -0.031716169369404179, 0.062855506871576484,
      0.0021607283583430178,
      # (log_eta, log_phi)
      0, -0.066108320373179186, 0.36489202031186536, 0.38414253350859978,
      0.20194367414627416, 0.045756760265233851, -0.020981683652594152,
      -0.00040781406791242802,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.23781212841854935, -0.47590719809637121, -0.48374092985756517,
      -0.40398939527837394, -0.32, -0.0049750934384736413,
      -0.000049997500093746875,
      # (log_phi, log_eta)
      0, -0.066108320373179186, 0.36489202031186536, 0.38414253350859978,
      0.20194367414627416, 0.045756760265233851, -0.020981683652594152,
      -0.00040781406791242802,
      # (log_phi, log_phi)
      0, -0.29275493050145557, -0.17759463733840194, 0.14829435062846670,
      0.30094940344894598, 0.3264, 0.0084154760082709998,
      0.000084991500478103751
    ),
    dim = c(m, 4, 4)
  )

  gh <- loglogistic4_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic4")

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
  theta <- lltd$theta_4
  theta[3:4] <- log(theta[3:4])

  true_value <- 0.13049198890475033

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
    class = "loglogistic4"
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
  theta <- lltd$theta_4
  theta[3:4] <- log(theta[3:4])

  true_gradient <- c(
    1.1763566356699270, 0.46825093563955542, -0.022171348272167373,
    0.15714798395105493
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 9.5072274862706741, -0.26196043558462076, 4.3708554404271919,
      # delta
      9.5072274862706741, 6.9361360507252671, -0.52829868503799506,
      2.1520196116438464,
      # log_eta
      -0.26196043558462076, -0.52829868503799506, 0.24291119321302682,
      -0.10852328156447682,
      # log_phi
      4.3708554404271919, 2.1520196116438464, -0.10852328156447682,
      1.5960201200846326
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
    class = "loglogistic4"
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(0, 1, 1.0597641935779216, 1.8170414285831683)

  true_value <- c(
    0.84803242222048425, -0.76457472791925385, 1.0597641935779216,
    1.8170414285831683
  )

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- mle_asy(object, theta)

  expect_type(result, "double")
  expect_length(result, 4)
  expect_equal(result, true_value)
})

test_that("fit", {
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.065671734266434990

  theta <- c(
    alpha = 0.84803242222048425, delta = -0.76457472791925385,
    eta = exp(1.0597641935779216), phi = exp(1.8170414285831683)
  )

  fitted_values <- rep(
    c(
      0.84803242222048425, 0.8193061170955434, 0.6768336235925299,
      0.4796839372400049, 0.3275506915281221, 0.234565997851723,
      0.083702650178120, 0.083458013114508
    ),
    k
  )

  residuals <- c(
    0.00466757777951575, -0.09093242222048425, 0.08896757777951575,
    0.0081938829044566, -0.0415061170955434, 0.0652938829044566,
    -0.1207336235925299, 0.0279663764074701, 0.0656160627599951,
    0.0063160627599951, 0.0554160627599951, -0.0397839372400049,
    0.0277493084718779, -0.0143506915281221, 0.0223493084718779,
    -0.092465997851723, -0.066902650178120, 0.054897349821880,
    0.039241986885492
  )

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- loglogistic4_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.065671734266434990

  theta <- c(
    alpha = 0.84803242222048425, delta = -0.76457472791925385,
    eta = exp(1.0597641935779216), phi = exp(1.8170414285831683)
  )

  fitted_values <- rep(
    c(
      0.84803242222048425, 0.8193061170955434, 0.6768336235925299,
      0.4796839372400049, 0.3275506915281221, 0.234565997851723,
      0.083702650178120, 0.083458013114508
    ),
    k
  )

  residuals <- c(
    0.00466757777951575, -0.09093242222048425, 0.08896757777951575,
    0.0081938829044566, -0.0415061170955434, 0.0652938829044566,
    -0.1207336235925299, 0.0279663764074701, 0.0656160627599951,
    0.0063160627599951, 0.0554160627599951, -0.0397839372400049,
    0.0277493084718779, -0.0143506915281221, 0.0223493084718779,
    -0.092465997851723, -0.066902650178120, 0.054897349821880,
    0.039241986885492
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.6, -0.6, 4, 8), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(-2, 2, 7, 1), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.18260381517441647

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(1.0311116743448195),
    phi = exp(2.0344389197855765)
  )

  fitted_values <- rep(
    c(
      0.8, 0.77954595995735380, 0.6742409051933649, 0.4974687113710638,
      0.3216436306778372, 0.1883607158194651, -0.099334450265717,
      -0.099998954510080
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.04795404004264620, -0.00174595995735380,
    0.10505404004264620, -0.1181409051933649, 0.0305590948066351,
    0.0478312886289362, -0.0114687113710638, 0.0376312886289362,
    -0.0575687113710638, 0.0336563693221628, -0.0084436306778372,
    0.0282563693221628, -0.0462607158194651, 0.116134450265717,
    0.237934450265717, 0.222698954510080
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.8, -0.9, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.18260381517441647

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(1.0311116743448195),
    phi = exp(2.0344389197855765)
  )

  fitted_values <- rep(
    c(
      0.8, 0.77954595995735380, 0.6742409051933649, 0.4974687113710638,
      0.3216436306778372, 0.1883607158194651, -0.099334450265717,
      -0.099998954510080
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.04795404004264620, -0.00174595995735380,
    0.10505404004264620, -0.1181409051933649, 0.0305590948066351,
    0.0478312886289362, -0.0114687113710638, 0.0376312886289362,
    -0.0575687113710638, 0.0336563693221628, -0.0084436306778372,
    0.0282563693221628, -0.0462607158194651, 0.116134450265717,
    0.237934450265717, 0.222698954510080
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.8, -0.9, 3, 7), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 8, 1), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.028149462100107135

  theta <- c(
    alpha = 0.85946644015940627, delta = -0.77755335460964313,
    eta = exp(1.2410728497889232), phi = exp(1.8216686605281363)
  )

  fitted_values <- rep(
    c(
      0.85946644015940627, 0.8440987075059350, 0.7183269332966290,
      0.4907843357327152, 0.307990679902764, 0.205756559537494,
      0.081964238053554, 0.081913103315067
    ),
    k
  )

  residuals <- c(
    -0.00676644015940627, -0.10236644015940627, 0.07753355984059373,
    -0.0165987075059350, -0.0662987075059350, 0.0405012924940650,
    -0.1622269332966290, -0.0135269332966290, 0.0545156642672848,
    -0.0047843357327152, 0.0443156642672848, -0.0508843357327152,
    0.047309320097236, 0.005209320097236, 0.041909320097236,
    -0.063656559537494, -0.065164238053554, 0.056635761946446,
    0.040786896684933
  )

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-6)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-6)
  expect_equal(result$weights, w)

  object <- loglogistic4_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.028149462100107135

  theta <- c(
    alpha = 0.85946644015940627, delta = -0.77755335460964313,
    eta = exp(1.2410728497889232), phi = exp(1.8216686605281363)
  )

  fitted_values <- rep(
    c(
      0.85946644015940627, 0.8440987075059350, 0.7183269332966290,
      0.4907843357327152, 0.307990679902764, 0.205756559537494,
      0.081964238053554, 0.081913103315067
    ),
    k
  )

  residuals <- c(
    -0.00676644015940627, -0.10236644015940627, 0.07753355984059373,
    -0.0165987075059350, -0.0662987075059350, 0.0405012924940650,
    -0.1622269332966290, -0.0135269332966290, 0.0545156642672848,
    -0.0047843357327152, 0.0443156642672848, -0.0508843357327152,
    0.047309320097236, 0.005209320097236, 0.041909320097236,
    -0.063656559537494, -0.065164238053554, 0.056635761946446,
    0.040786896684933
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.6, -0.6, 4, 8), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(-2, 2, 7, 1), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.047616281998049843

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(0.99935582830934212),
    phi = exp(1.8494626509107548)
  )

  fitted_values <- rep(
    c(
      0.9, 0.86270334460885678, 0.7008462168978897, 0.4851972721174256,
      0.3138327685225718, 0.2034181924346560, 0.000504528021877,
      9.69609809e-07
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.03520334460885678, -0.08490334460885678,
    0.02189665539114322, -0.1447462168978897, 0.0039537831021103,
    0.0601027278825744, 0.0008027278825744, 0.0499027278825744,
    -0.0452972721174256, 0.0414672314774282, -0.0006327685225718,
    0.0360672314774282, -0.0613181924346560, 0.016295471978123,
    0.138095471978123, 0.122699030390191
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.9, -0.9, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.047616281998049843

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(0.99935582830934212),
    phi = exp(1.8494626509107548)
  )

  fitted_values <- rep(
    c(
      0.9, 0.86270334460885678, 0.7008462168978897, 0.4851972721174256,
      0.3138327685225718, 0.2034181924346560, 0.000504528021877,
      9.69609809e-07
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.03520334460885678, -0.08490334460885678,
    0.02189665539114322, -0.1447462168978897, 0.0039537831021103,
    0.0601027278825744, 0.0008027278825744, 0.0499027278825744,
    -0.0452972721174256, 0.0414672314774282, -0.0006327685225718,
    0.0360672314774282, -0.0613181924346560, 0.016295471978123,
    0.138095471978123, 0.122699030390191
  )

  object <- loglogistic4_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 1, 3), c(0.9, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0.9, -0.9, 3, 7), max_iter,
    c(0.9, -0.9, 1, 3), c(0.9, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  object <- loglogistic4_new(
    x, y, w, c(0, 1, 8, 1), max_iter,
    c(0.9, -0.9, 1, 3), c(0.9, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic4_fit"))
  expect_true(inherits(result, "loglogistic"))
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
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  max_iter <- 10000

  theta <- lltd$theta_4
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- lltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.96, 3213.3730235910022, -25.345500911144129, 305.59395175952439,
      -12564.573197905927,
      # delta
      3213.3730235910022, 2314.5672831218128, -91.521827408072732,
      147.85088411543430, -6955.2079876339471,
      # eta
      -25.345500911144129, -91.521827408072732, 27.900844069407526,
      -4.2726065359333406, 353.68532278262947,
      # phi
      305.59395175952439, 147.85088411543430, -4.2726065359333406,
      20.082539137142281, -428.57767745995749,
      # sigma
      -12564.573197905927, -6955.2079876339471, 353.68532278262947,
      -428.57767745995749, 42844.361630297012
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- loglogistic4_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be of length 4"
  )
})

test_that("drda: 'start' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = rep(1, 5)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(0, 1, 1, -1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic4",
      start = c(0, 1, 1, 0)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic4")

  expect_equal(nauc(result), 0.087450408577304738)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.85599012126905068)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.0059897567503640588)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.36063326947370824
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic4")

  expect_equal(naac(result), 1 - 0.087450408577304738)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.85599012126905068)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.0059897567503640588)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.36063326947370824
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic4")

  expect_equal(nauc(result), 0.88198083382092461)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.15122707630633397)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.99526800138759278)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.84700641564881314
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic4")

  expect_equal(naac(result), 1 - 0.88198083382092461)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.15122707630633397)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.99526800138759278)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.84700641564881314
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
