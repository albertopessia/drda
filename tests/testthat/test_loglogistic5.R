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
    1.0, 0.76865930207047762, 0.58148893067026871, 0.45005285666246408,
    0.3631216481244481, 0.30597790621143287, 0.15211706430851203,
    0.1500212492031582
  )

  value <- loglogistic5_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic5")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic5_fit")

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
      0, 0.27216552697590868, 0.49236596391733093, 0.64699663922063049,
      0.74926864926535517, 0.81649658092772603, 0.99750933610763290,
      0.99997500093746094,
      # eta
      0, 0.098136730286166614, 0.035374259952478671, -0.029147447478978971,
      -0.065643670344200000, -0.080176576259309206, -0.0063184832702654069,
      -0.00011258080037357414,
      # phi
      0, 0.042840869986948588, 0.063410768080262317, 0.063947342248550688,
      0.055866522094346658, 0.046268139585904475, 0.00084366461262834624,
      8.4993625398414259e-06,
      # nu
      0, -0.096975924597482464, -0.069000993712605771, -0.039793214116310766,
      -0.022086761931019447, -0.012515261339478399, -2.6320687803095987e-06,
      -2.6560065272938719e-10
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
      0, -0.11545497680725484, -0.041616776414680789, 0.034291114681151731,
      0.077227847463764706, 0.094325383834481419, 0.0074335097297240082,
      0.00013244800043949899,
      # (delta, phi)
      0, -0.050401023514057163, -0.074600903623838020, -0.075232167351236104,
      -0.065725320110996068, -0.054433105395181736, -0.00099254660309217204,
      -9.9992500468722658e-06,
      # (delta, nu)
      0, 0.11408932305586172, 0.081177639661889143, 0.046815546019189137,
      0.025984425801199350, 0.014723836869974587, 3.0965515062465866e-06,
      3.1247135615222022e-10,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.11545497680725484, -0.041616776414680789, 0.034291114681151731,
      0.077227847463764706, 0.094325383834481419, 0.0074335097297240082,
      0.00013244800043949899,
      # (eta, eta)
      0, -0.034969579717974287, -0.0010763915442147436, 0.00067972427918967692,
      0.010554892707478964, 0.027787083890544811, 0.018787226907476051,
      0.00059644407533517483,
      # (eta, phi)
      0, 0.0061547213934039318, 0.029775878951814140, 0.030482406369536885,
      0.018950443000072276, 0.0070987545410903964, -0.0020866998577016452,
      -0.000040779261624360493,
      # (eta, nu)
      0, -0.0042956432234104312, -0.0075670958981963486, 0.0063640368405401223,
      0.012119037736476822, 0.011916943024646478, 0.000015698005713616964,
      2.8143441112200194e-09,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.050401023514057163, -0.074600903623838020, -0.075232167351236104,
      -0.065725320110996068, -0.054433105395181736, -0.00099254660309217204,
      -9.9992500468722658e-06,
      # (phi, eta)
      0, 0.0061547213934039318, 0.029775878951814140, 0.030482406369536885,
      0.018950443000072276, 0.0070987545410903964, -0.0020866998577016452,
      -0.000040779261624360493,
      # (phi, phi)
      0, -0.015232309328692831, -0.016140922784066772, -0.0095177439625749862,
      -0.0035284119217482100, 0, 0.00016621452069692792, 1.6996175398404963e-06,
      # (phi, nu)
      0, -0.0018752315499794457, -0.013564534316371215, -0.013962225756403018,
      -0.010313995026740982, -0.0068770058416855193, -2.0960492167070209e-06,
      -2.1247078395031084e-10,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0, 0.11408932305586172, 0.081177639661889143, 0.046815546019189137,
      0.025984425801199350, 0.014723836869974587, 3.0965515062465866e-06,
      3.1247135615222022e-10,
      # (nu, eta)
      0, -0.0042956432234104312, -0.0075670958981963486, 0.0063640368405401223,
      0.012119037736476822, 0.011916943024646478, 0.000015698005713616964,
      2.8143441112200194e-09,
      # (nu, phi)
      0, -0.0018752315499794457, -0.013564534316371215, -0.013962225756403018,
      -0.010313995026740982, -0.0068770058416855193, -2.0960492167070209e-06,
      -2.1247078395031084e-10,
      # (nu, nu)
      0, 0.031532325488433179, 0.027600584936730370, 0.013677167095710211,
      0.0060065113909322852, 0.0026503785893681148, 8.7253713503965157e-09,
      8.8528663369608776e-15
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
      0, 0.27216552697590868, 0.49236596391733093, 0.64699663922063049,
      0.74926864926535517, 0.81649658092772603, 0.99750933610763290,
      0.99997500093746094,
      # eta
      0, 0.098136730286166614, 0.035374259952478671, -0.029147447478978971,
      -0.065643670344200000, -0.080176576259309206, -0.0063184832702654069,
      -0.00011258080037357414,
      # phi
      0, 0.042840869986948588, 0.063410768080262317, 0.063947342248550688,
      0.055866522094346658, 0.046268139585904475, 0.00084366461262834624,
      8.4993625398414259e-06,
      # nu
      0, -0.096975924597482464, -0.069000993712605771, -0.039793214116310766,
      -0.022086761931019447, -0.012515261339478399, -2.6320687803095987e-06,
      -2.6560065272938719e-10
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
      0, -0.11545497680725484, -0.041616776414680789, 0.034291114681151731,
      0.077227847463764706, 0.094325383834481419, 0.0074335097297240082,
      0.00013244800043949899,
      # (delta, phi)
      0, -0.050401023514057163, -0.074600903623838020, -0.075232167351236104,
      -0.065725320110996068, -0.054433105395181736, -0.00099254660309217204,
      -9.9992500468722658e-06,
      # (delta, nu)
      0, 0.11408932305586172, 0.081177639661889143, 0.046815546019189137,
      0.025984425801199350, 0.014723836869974587, 3.0965515062465866e-06,
      3.1247135615222022e-10,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.11545497680725484, -0.041616776414680789, 0.034291114681151731,
      0.077227847463764706, 0.094325383834481419, 0.0074335097297240082,
      0.00013244800043949899,
      # (eta, eta)
      0, -0.034969579717974287, -0.0010763915442147436, 0.00067972427918967692,
      0.010554892707478964, 0.027787083890544811, 0.018787226907476051,
      0.00059644407533517483,
      # (eta, phi)
      0, 0.0061547213934039318, 0.029775878951814140, 0.030482406369536885,
      0.018950443000072276, 0.0070987545410903964, -0.0020866998577016452,
      -0.000040779261624360493,
      # (eta, nu)
      0, -0.0042956432234104312, -0.0075670958981963486, 0.0063640368405401223,
      0.012119037736476822, 0.011916943024646478, 0.000015698005713616964,
      2.8143441112200194e-09,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.050401023514057163, -0.074600903623838020, -0.075232167351236104,
      -0.065725320110996068, -0.054433105395181736, -0.00099254660309217204,
      -9.9992500468722658e-06,
      # (phi, eta)
      0, 0.0061547213934039318, 0.029775878951814140, 0.030482406369536885,
      0.018950443000072276, 0.0070987545410903964, -0.0020866998577016452,
      -0.000040779261624360493,
      # (phi, phi)
      0, -0.015232309328692831, -0.016140922784066772, -0.0095177439625749862,
      -0.0035284119217482100, 0, 0.00016621452069692792, 1.6996175398404963e-06,
      # (phi, nu)
      0, -0.0018752315499794457, -0.013564534316371215, -0.013962225756403018,
      -0.010313995026740982, -0.0068770058416855193, -2.0960492167070209e-06,
      -2.1247078395031084e-10,
      # (nu, alpha)
      rep(0, m),
      # (nu, delta)
      0, 0.11408932305586172, 0.081177639661889143, 0.046815546019189137,
      0.025984425801199350, 0.014723836869974587, 3.0965515062465866e-06,
      3.1247135615222022e-10,
      # (nu, eta)
      0, -0.0042956432234104312, -0.0075670958981963486, 0.0063640368405401223,
      0.012119037736476822, 0.011916943024646478, 0.000015698005713616964,
      2.8143441112200194e-09,
      # (nu, phi)
      0, -0.0018752315499794457, -0.013564534316371215, -0.013962225756403018,
      -0.010313995026740982, -0.0068770058416855193, -2.0960492167070209e-06,
      -2.1247078395031084e-10,
      # (nu, nu)
      0, 0.031532325488433179, 0.027600584936730370, 0.013677167095710211,
      0.0060065113909322852, 0.0026503785893681148, 8.7253713503965157e-09,
      8.8528663369608776e-15
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
      0, 0.27216552697590868, 0.49236596391733093, 0.64699663922063049,
      0.74926864926535517, 0.81649658092772603, 0.99750933610763290,
      0.99997500093746094,
      # log_eta
      0, 0.19627346057233323, 0.070748519904957342, -0.058294894957957942,
      -0.13128734068840000, -0.16035315251861841, -0.012636966540530814,
      -0.00022516160074714828,
      # log_phi
      0, 0.21420434993474294, 0.31705384040131158, 0.31973671124275344,
      0.27933261047173329, 0.23134069792952238, 0.0042183230631417312,
      0.000042496812699207130,
      # log_nu
      0, -0.19395184919496493, -0.13800198742521154, -0.079586428232621532,
      -0.044173523862038895, -0.025030522678956797, -5.2641375606191973e-06,
      -5.3120130545877438e-10
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
      0, -0.23090995361450968, -0.083233552829361579, 0.068582229362303462,
      0.15445569492752941, 0.18865076766896284, 0.014867019459448016,
      0.00026489600087899798,
      # (delta, log_phi)
      0, -0.25200511757028581, -0.37300451811919010, -0.37616083675618052,
      -0.32862660055498034, -0.27216552697590868, -0.0049627330154608602,
      -0.000049996250234361329,
      # (delta, log_nu)
      0, 0.22817864611172344, 0.16235527932377829, 0.093631092038378273,
      0.051968851602398700, 0.029447673739949173, 6.1931030124931733e-06,
      6.2494271230444044e-10,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.23090995361450968, -0.083233552829361579, 0.068582229362303462,
      0.15445569492752941, 0.18865076766896284, 0.014867019459448016,
      0.00026489600087899798,
      # (log_eta, log_eta)
      0, 0.056395141700436081, 0.066442953728098367, -0.055575997841199235,
      -0.089067769858484146, -0.049204816956439167, 0.062511941089373389,
      0.0021606147005935510,
      # (log_eta, log_phi)
      0, 0.061547213934039318, 0.29775878951814140, 0.30482406369536885,
      0.18950443000072276, 0.070987545410903964, -0.020866998577016452,
      -0.00040779261624360493,
      # (log_eta, log_nu)
      0, -0.017182572893641725, -0.030268383592785394, 0.025456147362160489,
      0.048476150945907289, 0.047667772098585913, 0.000062792022854467856,
      1.1257376444880078e-08,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.25200511757028581, -0.37300451811919010, -0.37616083675618052,
      -0.32862660055498034, -0.27216552697590868, -0.0049627330154608602,
      -0.000049996250234361329,
      # (log_phi, log_eta)
      0, 0.061547213934039318, 0.29775878951814140, 0.30482406369536885,
      0.18950443000072276, 0.070987545410903964, -0.020866998577016452,
      -0.00040779261624360493,
      # (log_phi, log_phi)
      0, -0.16660338328257784, -0.086469229200357705, 0.081793112178378787,
      0.19112231242802804, 0.23134069792952238, 0.0083736860805649291,
      0.000084987251195219538,
      # (log_phi, log_nu)
      0, -0.018752315499794457, -0.13564534316371215, -0.13962225756403018,
      -0.10313995026740982, -0.068770058416855193, -0.000020960492167070209,
      -2.1247078395031084e-09,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0, 0.22817864611172344, 0.16235527932377829, 0.093631092038378273,
      0.051968851602398700, 0.029447673739949173, 6.1931030124931733e-06,
      6.2494271230444044e-10,
      # (log_nu, log_eta)
      0, -0.017182572893641725, -0.030268383592785394, 0.025456147362160489,
      0.048476150945907289, 0.047667772098585913, 0.000062792022854467856,
      1.1257376444880078e-08,
      # (log_nu, log_phi)
      0, -0.018752315499794457, -0.13564534316371215, -0.13962225756403018,
      -0.10313995026740982, -0.068770058416855193, -0.000020960492167070209,
      -2.1247078395031084e-09,
      # (log_nu, log_nu)
      0, -0.067822547241232213, -0.027599647678290063, -0.024877759849780688,
      -0.020147478298309754, -0.014429008321484338, -5.2292360752176112e-06,
      -5.3116589399342653e-10
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
      0, 0.27216552697590868, 0.49236596391733093, 0.64699663922063049,
      0.74926864926535517, 0.81649658092772603, 0.99750933610763290,
      0.99997500093746094,
      # log_eta
      0, 0.19627346057233323, 0.070748519904957342, -0.058294894957957942,
      -0.13128734068840000, -0.16035315251861841, -0.012636966540530814,
      -0.00022516160074714828,
      # log_phi
      0, 0.21420434993474294, 0.31705384040131158, 0.31973671124275344,
      0.27933261047173329, 0.23134069792952238, 0.0042183230631417312,
      0.000042496812699207130,
      # log_nu
      0, -0.19395184919496493, -0.13800198742521154, -0.079586428232621532,
      -0.044173523862038895, -0.025030522678956797, -5.2641375606191973e-06,
      -5.3120130545877438e-10
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
      0, -0.23090995361450968, -0.083233552829361579, 0.068582229362303462,
      0.15445569492752941, 0.18865076766896284, 0.014867019459448016,
      0.00026489600087899798,
      # (delta, log_phi)
      0, -0.25200511757028581, -0.37300451811919010, -0.37616083675618052,
      -0.32862660055498034, -0.27216552697590868, -0.0049627330154608602,
      -0.000049996250234361329,
      # (delta, log_nu)
      0, 0.22817864611172344, 0.16235527932377829, 0.093631092038378273,
      0.051968851602398700, 0.029447673739949173, 6.1931030124931733e-06,
      6.2494271230444044e-10,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.23090995361450968, -0.083233552829361579, 0.068582229362303462,
      0.15445569492752941, 0.18865076766896284, 0.014867019459448016,
      0.00026489600087899798,
      # (log_eta, log_eta)
      0, 0.056395141700436081, 0.066442953728098367, -0.055575997841199235,
      -0.089067769858484146, -0.049204816956439167, 0.062511941089373389,
      0.0021606147005935510,
      # (log_eta, log_phi)
      0, 0.061547213934039318, 0.29775878951814140, 0.30482406369536885,
      0.18950443000072276, 0.070987545410903964, -0.020866998577016452,
      -0.00040779261624360493,
      # (log_eta, log_nu)
      0, -0.017182572893641725, -0.030268383592785394, 0.025456147362160489,
      0.048476150945907289, 0.047667772098585913, 0.000062792022854467856,
      1.1257376444880078e-08,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.25200511757028581, -0.37300451811919010, -0.37616083675618052,
      -0.32862660055498034, -0.27216552697590868, -0.0049627330154608602,
      -0.000049996250234361329,
      # (log_phi, log_eta)
      0, 0.061547213934039318, 0.29775878951814140, 0.30482406369536885,
      0.18950443000072276, 0.070987545410903964, -0.020866998577016452,
      -0.00040779261624360493,
      # (log_phi, log_phi)
      0, -0.16660338328257784, -0.086469229200357705, 0.081793112178378787,
      0.19112231242802804, 0.23134069792952238, 0.0083736860805649291,
      0.000084987251195219538,
      # (log_phi, log_nu)
      0, -0.018752315499794457, -0.13564534316371215, -0.13962225756403018,
      -0.10313995026740982, -0.068770058416855193, -0.000020960492167070209,
      -2.1247078395031084e-09,
      # (log_nu, alpha)
      rep(0, m),
      # (log_nu, delta)
      0, 0.22817864611172344, 0.16235527932377829, 0.093631092038378273,
      0.051968851602398700, 0.029447673739949173, 6.1931030124931733e-06,
      6.2494271230444044e-10,
      # (log_nu, log_eta)
      0, -0.017182572893641725, -0.030268383592785394, 0.025456147362160489,
      0.048476150945907289, 0.047667772098585913, 0.000062792022854467856,
      1.1257376444880078e-08,
      # (log_nu, log_phi)
      0, -0.018752315499794457, -0.13564534316371215, -0.13962225756403018,
      -0.10313995026740982, -0.068770058416855193, -0.000020960492167070209,
      -2.1247078395031084e-09,
      # (log_nu, log_nu)
      0, -0.067822547241232213, -0.027599647678290063, -0.024877759849780688,
      -0.020147478298309754, -0.014429008321484338, -5.2292360752176112e-06,
      -5.3116589399342653e-10
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

  object <- structure(list(stats = lltd$stats_1), class = "loglogistic5")

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

  true_value <- 0.13550858237432013

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
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

  value <- rss_fn(theta[c(2, 3)])

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  theta <- lltd$theta_5
  theta[3:5] <- log(theta[3:5])

  true_gradient <- c(
    0.37626542160678606, 0.13115216002236047, -0.068495061846838105,
    -0.077974141956824872, 0.058349828827166034
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 10.448511267521428, -0.082576427570544590, 3.6334852478615707,
      -1.3337668584372,
      # delta
      10.448511267521428, 7.7223662584479821, -0.29183144912487565,
      2.2315512750678686, -0.68861279711690616,
      # log_eta
      -0.08257642757054459, -0.29183144912487565, 0.20646906549981534,
      -0.13212057310461277, -0.081620488535810301,
      # log_phi
      3.6334852478615707, 2.2315512750678686, -0.13212057310461277,
      1.1102331695379725, -0.32982499549780055,
      # log_nu
      -1.3337668584372, -0.68861279711690616, -0.081620488535810301,
      -0.32982499549780055, 0.1992660473652395
    ),
    nrow = 5,
    ncol = 5
  )

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
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

  gh <- rss_gh(theta[c(2, 3)])

  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, 2)
  expect_length(gh$H, 2 * 2)

  expect_equal(gh$G, true_gradient[c(2, 3)])
  expect_equal(gh$H, true_hessian[c(2, 3), c(2, 3)])
})

test_that("mle_asy", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- rep(1, length(y))

  max_iter <- 10000

  theta <- c(
    0, 1, 3.3501562141129222, 2.2334446014256944, 3.0041670174896128
  )

  true_value <- c(
    0.86567490749559893, -0.77298635934282346, 3.3501562141129222,
    2.2334446014256944, 3.0041670174896128
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
    alpha = 0.86567490749559893, delta = -0.77298635934282346,
    eta = exp(3.3501562141129222), phi = exp(2.2334446014256944),
    nu = exp(3.0041670174896128)
  )

  rss_value <- 0.058276476180507351

  fitted_values <- rep(
    c(
      0.86567490749559893, 0.7901639477776979, 0.6645425348388443,
      0.508922320974665, 0.329951365602280, 0.142289339340291,
      0.092688548152775, 0.092688548152775
    ),
    k
  )

  residuals <- c(
    -0.01297490749559893, -0.10857490749559893, 0.07132509250440107,
    0.0373360522223021, -0.0123639477776979, 0.0944360522223021,
    -0.1084425348388443, 0.0402574651611557, 0.036377679025335,
    -0.022922320974665, 0.026177679025335, -0.069022320974665,
    0.025348634397720, -0.016751365602280, 0.019948634397720,
    -0.000189339340291, -0.075888548152775, 0.045911451847225, 0.030011451847225
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.86567490749559893, delta = -0.77298635934282346,
    eta = exp(3.3501562141129222), phi = exp(2.2334446014256944),
    nu = exp(3.0041670174896128)
  )

  rss_value <- 0.058276476180507351

  fitted_values <- rep(
    c(
      0.86567490749559893, 0.7901639477776979, 0.6645425348388443,
      0.508922320974665, 0.329951365602280, 0.142289339340291,
      0.092688548152775, 0.092688548152775
    ),
    k
  )

  residuals <- c(
    -0.01297490749559893, -0.10857490749559893, 0.07132509250440107,
    0.0373360522223021, -0.0123639477776979, 0.0944360522223021,
    -0.1084425348388443, 0.0402574651611557, 0.036377679025335,
    -0.022922320974665, 0.026177679025335, -0.069022320974665,
    0.025348634397720, -0.016751365602280, 0.019948634397720,
    -0.000189339340291, -0.075888548152775, 0.045911451847225, 0.030011451847225
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter, c(0.5, -1, 25, 8, 15), c(1, -0.5, 30, 12, 30)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.7, -0.6, 29, 11, 16), max_iter,
    c(0.5, -1, 25, 8, 15), c(1, -0.5, 30, 12, 30)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, 2, 1, 1, 1), max_iter,
    c(0.5, -1, 25, 8, 15), c(1, -0.5, 30, 12, 30)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta = c(
    alpha = 0.8, delta = -0.8, eta = exp(1.2607766612631776),
    phi = exp(1.9912508109520300), nu = exp(0.29620502796118245)
  )

  rss_value <- 0.098711571468416316

  fitted_values <- rep(
    c(
      0.8, 0.77882319799854372, 0.67672318574486147, 0.49873167092601701,
      0.31956127658455147, 0.19261518024345144, 0.000079037608724045848,
      2.3427366607360274e-08
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.1370, 0.048676802001456279, -0.0010231979985437210,
    0.10577680200145628, -0.12062318574486147, 0.028076814255138528,
    0.046568329073982991, -0.012731670926017009, 0.036368329073982991,
    -0.058831670926017009, 0.035738723415448528, -0.0063612765845514721,
    0.030338723415448528, -0.050515180243451441, 0.016720962391275954,
    0.13852096239127595, 0.12269997657263339
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.8, rep(-Inf, 3)), c(0.8, -0.8, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(0.8, -0.8, 1, 1, 1), max_iter,
    c(0.8, -0.8, rep(-Inf, 3)), c(0.8, -0.8, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), max_iter,
    c(0.8, -0.8, rep(-Inf, 3)), c(0.8, -0.8, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta = c(
    alpha = 0.8, delta = -0.8, eta = exp(1.2607766612631776),
    phi = exp(1.9912508109520300), nu = exp(0.29620502796118245)
  )

  rss_value <- 0.098711571468416316

  fitted_values <- rep(
    c(
      0.8, 0.77882319799854372, 0.67672318574486147, 0.49873167092601701,
      0.31956127658455147, 0.19261518024345144, 0.000079037608724045848,
      2.3427366607360274e-08
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.1370, 0.048676802001456279, -0.0010231979985437210,
    0.10577680200145628, -0.12062318574486147, 0.028076814255138528,
    0.046568329073982991, -0.012731670926017009, 0.036368329073982991,
    -0.058831670926017009, 0.035738723415448528, -0.0063612765845514721,
    0.030338723415448528, -0.050515180243451441, 0.016720962391275954,
    0.13852096239127595, 0.12269997657263339
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.8, 3, 5, 1), c(0.8, -0.8, 8, 10, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.8, -0.8, 7, 9, 2.1), max_iter,
    c(0.8, -0.8, 3, 5, 1), c(0.8, -0.8, 8, 10, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 10), max_iter,
    c(0.8, -0.8, 3, 5, 1), c(0.8, -0.8, 8, 10, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta = c(
    alpha = 0.90681779366690136, delta = -0.81378790717255029,
    eta = exp(2.7080648141767070), phi = exp(2.1486558471371228),
    nu = exp(2.3522233117187749)
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
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loglogistic5_new(
    x, y, w, c(1, -1, 1, 1, 1), max_iter, NULL, NULL
  )

  result <- fit(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta = c(
    alpha = 0.90681779366690136, delta = -0.81378790717255029,
    eta = exp(2.7080648141767070), phi = exp(2.1486558471371228),
    nu = exp(2.3522233117187749)
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
    c(0.5, -1, 10, 8, 5), c(1, -0.5, 20, 10, 15)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.7, -0.6, 18, 9.5, 7), max_iter,
    c(0.5, -1, 10, 8, 5), c(1, -0.5, 20, 10, 15)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(-2, -5, 0.5, 20, 0.1), 10000,
    c(0.5, -1, 10, 8, 5), c(1, -0.5, 20, 10, 15)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 5)
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta = c(
    alpha = 0.9, delta = -0.9, eta = exp(1.7810821447704601),
    phi = exp(2.0909375125759564), nu = exp(1.3229883342089993)
  )

  rss_value <- 0.041203263957338706

  fitted_values <- rep(
    c(
      0.9, 0.83059019593271452, 0.69255393105735286, 0.51034283205585752,
      0.31437042735601886, 0.15843123090622139, 2.9670052902024567e-07,
      3.4359102105208779e-13
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.0030901959327145152, -0.052790195932714515,
    0.054009804067285485, -0.13645393105735286, 0.012246068942647139,
    0.034957167944142476, -0.024342832055857524, 0.024757167944142476,
    -0.070442832055857524, 0.040929572643981141, -0.0011704273560188594,
    0.035529572643981141, -0.016331230906221388, 0.016799703299470980,
    0.13859970329947098, 0.12269999999965641
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with same equalities
  object <- loglogistic5_new(
    x, y, w, c(0.9, -0.9, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values with different equalities
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 3)), c(0.9, -0.9, rep(Inf, 3))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
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
    alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE
  )

  theta = c(
    alpha = 0.9, delta = -0.9, eta = exp(1.7810821447704601),
    phi = exp(2.0909375125759564), nu = exp(1.3229883342089993)
  )

  rss_value <- 0.041203263957338706

  fitted_values <- rep(
    c(
      0.9, 0.83059019593271452, 0.69255393105735286, 0.51034283205585752,
      0.31437042735601886, 0.15843123090622139, 2.9670052902024567e-07,
      3.4359102105208779e-13
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.0030901959327145152, -0.052790195932714515,
    0.054009804067285485, -0.13645393105735286, 0.012246068942647139,
    0.034957167944142476, -0.024342832055857524, 0.024757167944142476,
    -0.070442832055857524, 0.040929572643981141, -0.0011704273560188594,
    0.035529572643981141, -0.016331230906221388, 0.016799703299470980,
    0.13859970329947098, 0.12269999999965641
  )

  object <- loglogistic5_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 5, 5, 1), c(0.9, -0.9, 15, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0.9, -0.9, 8, 7, 2), max_iter,
    c(0.9, -0.9, 5, 5, 1), c(0.9, -0.9, 15, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loglogistic5_new(
    x, y, w, c(0, 1, 0.5, 0.5, 0.5), max_iter,
    c(0.9, -0.9, 5, 5, 1), c(0.9, -0.9, 15, 10, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loglogistic5_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)
})

test_that("fisher_info", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  max_iter <- 10000

  theta <- lltd$theta_5
  names(theta) <- c("alpha", "delta", "eta", "phi", "nu")

  sigma <- lltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.9600000000000, 3572.2954245320960, 9.9745955351015971,
      258.77353469936811, -255.55983146098816, -361.21156590873761,
      # delta
      3572.2954245320960, 2591.2144943994337, -40.923062458555474,
      156.43965827917267, -128.14086175964824, -2053.1823914615566,
      # eta
      9.9745955351015971, -40.923062458555474, 30.120442992389692,
      -4.0134045057114794, -9.3731945220393315, 781.70402665038452,
      # phi
      258.77353469936811, 156.43965827917267, -4.0134045057114794,
      17.311396456393088, -12.162137948260707, 244.66467743210642,
      # nu
      -255.55983146098816, -128.14086175964824, -9.3731945220393315,
      -12.162137948260707, 13.350580942202232, -574.23470712498419,
      # sigma
      -361.21156590873761, -2053.1823914615566, 781.70402665038452,
      244.66467743210642, -574.23470712498419, 52994.534857391555
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "nu", "sigma"
  )

  object <- loglogistic5_new(x, y, w, NULL, max_iter, NULL, NULL)

  fim <- fisher_info(object, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
  expect_equal(fim, true_value)
})

test_that("drda: 'lower_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c("a", "b", "c", "d", "e")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = matrix(-Inf, nrow = 5, ncol = 2),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      upper_bound = c("a", "b", "c", "d", "e")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 5),
      upper_bound = matrix(Inf, nrow = 5, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      lower_bound = rep(-Inf, 6),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be of length 5"
  )
})

test_that("drda: 'start' argument errors", {
  x <- lltd$D$x
  y <- lltd$D$y

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c("a", "b", "c", "d", "e")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(-Inf, 1, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(1, 1, 1, 1, 1, 1)
    ),
    "'start' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, -1, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 0, 1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, -1, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 0, 1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, -1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loglogistic5",
      start = c(0, 1, 1, 1, 0)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.097909681560146038)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.87325260271259105)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.0061191350953578714)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.41629510972061593
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(naac(result), 1 - 0.097909681560146038)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.87325260271259105)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.0061191350953578714)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.41629510972061593
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(nauc(result), 0.84987748063201568)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.16335338880480354)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.99497853656838019)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.79822066734103837
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loglogistic5")

  expect_equal(naac(result), 1 - 0.84987748063201568)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.16335338880480354)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.99497853656838019)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.79822066734103837
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
