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

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  i <- c(1, 2)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "gompertz"))
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

  object <- gompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "gompertz"))
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
    gompertz_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    gompertz_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    gompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

  m <- length(x)

  true_value <- c(
    0.9, 0.89999900489042372, 0.69362857548958928, 0.59130628183991215,
    0.49565490901875908, 0.41545956071125755, 0.20315447757070948,
    0.20002130241396156
  )

  value <- gompertz_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "gompertz")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = ltd$stats_1), class = "gompertz_fit")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 1.4215851089662344e-06, 0.29481632072915817, 0.44099102594298265,
      0.57763584425891561, 0.69220062755534635, 0.99549360347041503,
      0.99996956798005491,
      # eta
      0, 0.00034834526092876898, 0.50412525420491786, -0.50547408067930202,
      -1.3314557455680168, -1.7825246603050775, -0.16995739715365151,
      -0.0022154173413352887,
      # phi
      0, 1.3397894651106499e-06, 0.025206262710245893, 0.025273704033965101,
      0.022190929092800280, 0.017825246603050775, 0.00031473592065491021,
      2.1302089820531622e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- gompertz_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

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
      0, -0.00049763608704109854, -0.72017893457845408, 0.72210582954186003,
      1.9020796365257383, 2.5464638004358250, 0.24279628164807359,
      0.0031648819161932696,
      # (delta, phi)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.00049763608704109854, -0.72017893457845408, 0.72210582954186003,
      1.9020796365257383, 2.5464638004358250, 0.24279628164807359,
      0.0031648819161932696,
      # (eta, eta)
      0, -0.11288378602671264, -0.22322944347833114, 0.18325381188667304,
      3.6044240367362504, 11.267704843977305, 9.1362476238810089,
      0.23039639175120797,
      # (eta, phi)
      0, -0.00042077051314394212, 0.24090115492854237, 0.24357434974531736,
      0.16183555698239863, 0.065575417590734696, -0.013771617874712026,
      -0.00020023290224793758,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (phi, eta)
      0, -0.00042077051314394212, 0.24090115492854237, 0.24357434974531736,
      0.16183555698239863, 0.065575417590734696, -0.013771617874712026,
      -0.00020023290224793758,
      # (phi, phi)
      0, -1.6698784915194178e-06, -0.00055807360869582784,
      0.00045813452971668259, 0.0010012288990934029, 0.0011267704843977305,
      0.000031331439039372459, 2.1301441545045115e-07
    ),
    dim = c(m, 4, 4)
  )

  H <- gompertz_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 1.4215851089662344e-06, 0.29481632072915817, 0.44099102594298265,
      0.57763584425891561, 0.69220062755534635, 0.99549360347041503,
      0.99996956798005491,
      # eta
      0, 0.00034834526092876898, 0.50412525420491786, -0.50547408067930202,
      -1.3314557455680168, -1.7825246603050775, -0.16995739715365151,
      -0.0022154173413352887,
      # phi
      0, 1.3397894651106499e-06, 0.025206262710245893, 0.025273704033965101,
      0.022190929092800280, 0.017825246603050775, 0.00031473592065491021,
      2.1302089820531622e-06
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
      0, -0.00049763608704109854, -0.72017893457845408, 0.72210582954186003,
      1.9020796365257383, 2.5464638004358250, 0.24279628164807359,
      0.0031648819161932696,
      # (delta, phi)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.00049763608704109854, -0.72017893457845408, 0.72210582954186003,
      1.9020796365257383, 2.5464638004358250, 0.24279628164807359,
      0.0031648819161932696,
      # (eta, eta)
      0, -0.11288378602671264, -0.22322944347833114, 0.18325381188667304,
      3.6044240367362504, 11.267704843977305, 9.1362476238810089,
      0.23039639175120797,
      # (eta, phi)
      0, -0.00042077051314394212, 0.24090115492854237, 0.24357434974531736,
      0.16183555698239863, 0.065575417590734696, -0.013771617874712026,
      -0.00020023290224793758,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (phi, eta)
      0, -0.00042077051314394212, 0.24090115492854237, 0.24357434974531736,
      0.16183555698239863, 0.065575417590734696, -0.013771617874712026,
      -0.00020023290224793758,
      # (phi, phi)
      0, -1.6698784915194178e-06, -0.00055807360869582784,
      0.00045813452971668259, 0.0010012288990934029, 0.0011267704843977305,
      0.000031331439039372459, 2.1301441545045115e-07
    ),
    dim = c(m, 4, 4)
  )

  gh <- gompertz_gradient_hessian(x, theta)

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
  theta <- ltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0, 1.4215851089662344e-06, 0.29481632072915817, 0.44099102594298265,
      0.57763584425891561, 0.69220062755534635, 0.99549360347041503,
      0.99996956798005491,
      # log_eta
      0.0, 0.000034834526092876898, 0.050412525420491786, -0.050547408067930202,
      -0.13314557455680168, -0.17825246603050775, -0.016995739715365151,
      -0.00022154173413352887,
      # phi
      0.0, 1.3397894651106499e-06, 0.025206262710245893, 0.025273704033965101,
      0.022190929092800280, 0.017825246603050775, 0.00031473592065491021,
      2.1302089820531622e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- gompertz_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

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
      0, -0.000049763608704109854, -0.072017893457845408, 0.072210582954186003,
      0.19020796365257383, 0.25464638004358250, 0.024279628164807359,
      0.00031648819161932696,
      # (delta, phi)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.000049763608704109854, -0.072017893457845408, 0.072210582954186003,
      0.19020796365257383, 0.25464638004358250, 0.024279628164807359,
      0.00031648819161932696,
      # (log_eta, log_eta)
      0, -0.0010940033341742495, 0.048180230985708474, -0.048714869949063472,
      -0.097101334189439179, -0.065575417590734696, 0.074366736523444938,
      0.0020824221833785508,
      # (log_eta, phi)
      0, -0.000042077051314394212, 0.024090115492854237, 0.024357434974531736,
      0.016183555698239863, 0.0065575417590734696, -0.0013771617874712026,
      -0.000020023290224793758,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (phi, log_eta)
      0, -0.000042077051314394212, 0.024090115492854237, 0.024357434974531736,
      0.016183555698239863, 0.0065575417590734696, -0.0013771617874712026,
      -0.000020023290224793758,
      # (phi, phi)
      0, -1.6698784915194178e-06, -0.00055807360869582784,
      0.00045813452971668259, 0.0010012288990934029, 0.0011267704843977305,
      0.000031331439039372459, 2.1301441545045115e-07
    ),
    dim = c(m, 4, 4)
  )

  H <- gompertz_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- ltd$stats_1[, 1]
  theta <- ltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0.0, 1.4215851089662344e-06, 0.29481632072915817, 0.44099102594298265,
      0.57763584425891561, 0.69220062755534635, 0.99549360347041503,
      0.99996956798005491,
      # log_eta
      0.0, 0.000034834526092876898, 0.050412525420491786, -0.050547408067930202,
      -0.13314557455680168, -0.17825246603050775, -0.016995739715365151,
      -0.00022154173413352887,
      # phi
      0.0, 1.3397894651106499e-06, 0.025206262710245893, 0.025273704033965101,
      0.022190929092800280, 0.017825246603050775, 0.00031473592065491021,
      2.1302089820531622e-06
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
      0, -0.000049763608704109854, -0.072017893457845408, 0.072210582954186003,
      0.19020796365257383, 0.25464638004358250, 0.024279628164807359,
      0.00031648819161932696,
      # (delta, phi)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.000049763608704109854, -0.072017893457845408, 0.072210582954186003,
      0.19020796365257383, 0.25464638004358250, 0.024279628164807359,
      0.00031648819161932696,
      # (log_eta, log_eta)
      0, -0.0010940033341742495, 0.048180230985708474, -0.048714869949063472,
      -0.097101334189439179, -0.065575417590734696, 0.074366736523444938,
      0.0020824221833785508,
      # (log_eta, phi)
      0, -0.000042077051314394212, 0.024090115492854237, 0.024357434974531736,
      0.016183555698239863, 0.0065575417590734696, -0.0013771617874712026,
      -0.000020023290224793758,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -1.9139849501580713e-06, -0.036008946728922704, -0.036105291477093002,
      -0.031701327275428972, -0.025464638004358250, -0.00044962274379272887,
      -3.0431556886473746e-06,
      # (phi, log_eta)
      0, -0.000042077051314394212, 0.024090115492854237, 0.024357434974531736,
      0.016183555698239863, 0.0065575417590734696, -0.0013771617874712026,
      -0.000020023290224793758,
      # (phi, phi)
      0, -1.6698784915194178e-06, -0.00055807360869582784,
      0.00045813452971668259, 0.0010012288990934029, 0.0011267704843977305,
      0.000031331439039372459, 2.1301441545045115e-07
    ),
    dim = c(m, 4, 4)
  )

  gh <- gompertz_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(list(stats = ltd$stats_1), class = "gompertz")

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
  theta <- ltd$theta_g
  theta[3] <- log(theta[3])

  true_value <- 0.24809286101003823

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "gompertz"
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
  theta <- ltd$theta_g
  theta[3] <- log(theta[3])

  true_gradient <- c(
    1.9187338373330136, 0.98251185903665707, -0.12716128182964986,
    0.027606363431717869
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 7.7696659452385520, -0.71316228871823523, 0.23654099685649101,
      # delta
      7.7696659452385520, 5.4138110462849834, -0.26595524194628816,
      0.071428920313740434,
      # log_eta
      -0.71316228871823523, -0.26595524194628816, 0.044607391404419511,
      0.0061855101193171994,
      # phi
      0.23654099685649101, 0.071428920313740434, 0.0061855101193171994,
      0.0064996017731646237
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = ltd$stats_1, m = nrow(ltd$stats_1)),
    class = "gompertz"
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

  theta <- c(0, 1, -1.8198935068871131, -3.7323345311514175)

  true_value <- c(
    0.83750398486290947, -0.75377221534475194, -1.8198935068871131,
    -3.7323345311514175
  )

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

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

  rss_value <- 0.0607778298867072

  theta <- c(
    alpha = 0.83750398486290947, delta = -0.75377221534475194,
    eta = exp(-1.8198935068871131), phi = -3.7323345311514175
  )

  fitted_values <- rep(
    c(
      0.83750398486290947, 0.83750398486290947, 0.65963821695088722,
      0.48330993332321737, 0.32969963891713439, 0.22441628476824283,
      0.083856457410307348, 0.083731807287638583
    ),
    k
  )

  residuals <- c(
    0.015196015137090534, -0.080403984862909466, 0.099496015137090534,
    -0.010003984862909466, -0.059703984862909466, 0.047096015137090534,
    -0.10353821695088722, 0.045161783049112785, 0.061990066676782626,
    0.0026900666767826257, 0.051790066676782626, -0.043409933323217374,
    0.025600361082865611, -0.016499638917134389, 0.020200361082865611,
    -0.082316284768242828, -0.067056457410307348, 0.054743542589692652,
    0.038968192712361417
  )

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  object <- gompertz_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.0607778298867072

  theta <- c(
    alpha = 0.83750398486290947, delta = -0.75377221534475194,
    eta = exp(-1.8198935068871131), phi = -3.7323345311514175
  )

  fitted_values <- rep(
    c(
      0.83750398486290947, 0.83750398486290947, 0.65963821695088722,
      0.48330993332321737, 0.32969963891713439, 0.22441628476824283,
      0.083856457410307348, 0.083731807287638583
    ),
    k
  )

  residuals <- c(
    0.015196015137090534, -0.080403984862909466, 0.099496015137090534,
    -0.010003984862909466, -0.059703984862909466, 0.047096015137090534,
    -0.10353821695088722, 0.045161783049112785, 0.061990066676782626,
    0.0026900666767826257, 0.051790066676782626, -0.043409933323217374,
    0.025600361082865611, -0.016499638917134389, 0.020200361082865611,
    -0.082316284768242828, -0.067056457410307348, 0.054743542589692652,
    0.038968192712361417
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.6, -0.6, 2, 2), max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(-2, 2, 7, -8), max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.17387087229710121

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-2.0754477675542197),
    phi = -1.4340713249001234
  )

  fitted_values <- rep(
    c(
      0.8, 0.79999999999999980, 0.64725490118103268, 0.49240389185498950,
      0.33009961989054187, 0.19270533171641992, -0.098585719005589583,
      -0.099997335143268946
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.027500000000000198, -0.022199999999999802,
    0.084600000000000198, -0.091154901181032677, 0.057545098818967323,
    0.0528961081450105, -0.0064038918549894998, 0.042696108145010500,
    -0.0525038918549895, 0.025200380109458134, -0.016899619890541866,
    0.019800380109458134, -0.050605331716419919, 0.11538571900558958,
    0.23718571900558958, 0.22269733514326895
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.8, -0.9, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.17387087229710121

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(-2.0754477675542197),
    phi = -1.4340713249001234
  )

  fitted_values <- rep(
    c(
      0.8, 0.79999999999999980, 0.64725490118103268, 0.49240389185498950,
      0.33009961989054187, 0.19270533171641992, -0.098585719005589583,
      -0.099997335143268946
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.027500000000000198, -0.022199999999999802,
    0.084600000000000198, -0.091154901181032677, 0.057545098818967323,
    0.0528961081450105, -0.0064038918549894998, 0.042696108145010500,
    -0.0525038918549895, 0.025200380109458134, -0.016899619890541866,
    0.019800380109458134, -0.050605331716419919, 0.11538571900558958,
    0.23718571900558958, 0.22269733514326895
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.8, -0.9, 0.5, 2), max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, 8, -5), max_iter,
    c(0.8, -0.9, 0.05, -3), c(0.8, -0.9, 5, 3)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.027430943413179750

  theta <- c(
    alpha = 0.84880506428954864, delta = -0.76445273232771366,
    eta = exp(-1.6109893031741691), phi = -3.2943267184810845
  )

  fitted_values <- rep(
    c(
      0.84880506428954864, 0.84880506428954864, 0.71143950746621771,
      0.49564428958059577, 0.30871211456949977, 0.19496610079954203,
      0.084370589196723306, 0.084352332803671482
    ),
    k
  )

  residuals <- c(
    0.0038949357104513572, -0.091705064289548643, 0.088194935710451357,
    -0.021305064289548643, -0.071005064289548643, 0.035794935710451357,
    -0.15533950746621771, -0.0066395074662177070, 0.049655710419404225,
    -0.0096442895805957749, 0.039455710419404225, -0.055744289580595775,
    0.046587885430500226, 0.0044878854305002261, 0.041187885430500226,
    -0.052866100799542033, -0.067570589196723306, 0.054229410803276694,
    0.038347667196328518
  )

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  object <- gompertz_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.027430943413179750

  theta <- c(
    alpha = 0.84880506428954864, delta = -0.76445273232771366,
    eta = exp(-1.6109893031741691), phi = -3.2943267184810845
  )

  fitted_values <- rep(
    c(
      0.84880506428954864, 0.84880506428954864, 0.71143950746621771,
      0.49564428958059577, 0.30871211456949977, 0.19496610079954203,
      0.084370589196723306, 0.084352332803671482
    ),
    k
  )

  residuals <- c(
    0.0038949357104513572, -0.091705064289548643, 0.088194935710451357,
    -0.021305064289548643, -0.071005064289548643, 0.035794935710451357,
    -0.15533950746621771, -0.0066395074662177070, 0.049655710419404225,
    -0.0096442895805957749, 0.039455710419404225, -0.055744289580595775,
    0.046587885430500226, 0.0044878854305002261, 0.041187885430500226,
    -0.052866100799542033, -0.067570589196723306, 0.054229410803276694,
    0.038347667196328518
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.6, -0.6, 2, 2), max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(2, -2, 7, -8), max_iter,
    c(0.5, -1, 0.05, -5), c(1, -0.5, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.056975135420412259

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-1.8247274007834970),
    phi = -3.3265712508687191
  )

  fitted_values <- rep(
    c(
      0.9, 0.9, 0.70686056203863183, 0.49858900619492442, 0.31078143736730623,
      0.17934588458173208, 0.00016575815950820939, 5.2211214259512707e-08
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.0725, -0.1222, -0.0154, -0.15076056203863183,
    -0.0020605620386318339, 0.046710993805075579, -0.012589006194924421,
    0.036510993805075579, -0.058689006194924421, 0.044518562632693765,
    0.0024185626326937652, 0.039118562632693765, -0.03724588458173208,
    0.016634241840491791, 0.13843424184049179, 0.12269994778878574
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.9, -0.9, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.9, -0.9, rep(-Inf, 2)), c(0.9, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  rss_value <- 0.056975135420412259

  theta <- c(
    alpha = 0.9, delta = -0.9, eta = exp(-1.8247274007834970),
    phi = -3.3265712508687191
  )

  fitted_values <- rep(
    c(
      0.9, 0.9, 0.70686056203863183, 0.49858900619492442, 0.31078143736730623,
      0.17934588458173208, 0.00016575815950820939, 5.2211214259512707e-08
    ),
    k
  )

  residuals <- c(
    -0.0473, -0.1429, 0.037, -0.0725, -0.1222, -0.0154, -0.15076056203863183,
    -0.0020605620386318339, 0.046710993805075579, -0.012589006194924421,
    0.036510993805075579, -0.058689006194924421, 0.044518562632693765,
    0.0024185626326937652, 0.039118562632693765, -0.03724588458173208,
    0.016634241840491791, 0.13843424184049179, 0.12269994778878574
  )

  object <- gompertz_new(
    x, y, w, NULL, max_iter,
    c(0.9, -0.9, 0.05, -5), c(0.9, -0.9, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0.9, -0.9, 0.5, 2), max_iter,
    c(0.9, -0.9, 0.05, -5), c(0.9, -0.9, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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
  object <- gompertz_new(
    x, y, w, c(0, 1, 8, -8), max_iter,
    c(0.9, -0.9, 0.05, -5), c(0.9, -0.9, 5, 5)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "gompertz_fit"))
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

  theta <- ltd$theta_g
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- ltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.96, 2575.7736547102818, -2729.0666750822329, 78.333784600134089,
      -24869.75366811211,
      # delta
      2575.7736547102818, 1783.9750097718574, -919.06108253631968,
      24.812675215544203, -13840.186165313782,
      # eta
      -2729.0666750822329, -919.06108253631968, 6872.4198199637835,
      7.7338509890991883, 22399.260352877255,
      # phi
      78.333784600134089, 24.812675215544203, 7.7338509890991883,
      2.2266959647072065, -368.96905053777609,
      # sigma
      -24869.753668112110, -13840.186165313782, 22399.260352877255,
      -368.96905053777609, 117592.66277720903
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- gompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

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
      y ~ x, mean_function = "gompertz",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
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
      y ~ x, mean_function = "gompertz",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
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
      y ~ x, mean_function = "gompertz",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = rep(1, 5)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "gompertz",
      start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "gompertz")

  expect_equal(nauc(result), 0.43882466615261643)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.39636244305025756)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.40109615704980747)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 0.018304282204149767
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- ltd$D$x
  y <- ltd$D$y
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "gompertz")

  expect_equal(naac(result), 1 - 0.43882466615261643)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.39636244305025756)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.40109615704980747)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(1, 5), ylim = c(0.3, 0.7)), 1 - 0.018304282204149767
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "gompertz")

  expect_equal(nauc(result), 0.61829885821624581)
  expect_equal(nauc(result, xlim = c(-2, 2)), 0.67083050518516769)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.69668945965796646)
  expect_equal(nauc(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 0.62072202939460955
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- ltd$D$x
  y <- rev(ltd$D$y)
  w <- ltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "gompertz")

  expect_equal(naac(result), 1 - 0.61829885821624581)
  expect_equal(naac(result, xlim = c(-2, 2)), 1 - 0.67083050518516769)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.69668945965796646)
  expect_equal(naac(result, xlim = c(-15, -10), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(-5, -1), ylim = c(0.3, 0.7)), 1 - 0.62072202939460955
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})
