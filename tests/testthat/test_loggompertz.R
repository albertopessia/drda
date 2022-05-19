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

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  i <- c(1, 2)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  expect_true(inherits(object, "loggompertz"))
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

  object <- loggompertz_new(x, y, w, start, max_iter, lower_bound, upper_bound)

  expect_true(inherits(object, "loggompertz"))
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
    loggompertz_new(x, y, w, c(0, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 1, 1, 1), max_iter, NULL, NULL),
    "'start' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 0, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, -1, 1), max_iter, NULL, NULL),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 1, 0), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, c(0, 1, 1, -1), max_iter, NULL, NULL),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 3)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 3), rep(Inf, 4)),
    "'lower_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), rep(Inf, 3)),
    "'upper_bound' must be of length 4"
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, 0, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, -1, Inf)),
    "'upper_bound[3]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, 0)),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )

  expect_error(
    loggompertz_new(x, y, w, NULL, max_iter, rep(-Inf, 4), c(1, 1, Inf, -1)),
    "'upper_bound[4]' cannot be negative nor zero",
    fixed = TRUE
  )
})

test_that("Function value", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

  m <- length(x)

  true_value <- c(
    1, 0.99835911398420645, 0.82183032092156685, 0.57555097969061526,
    0.42486123076253040, 0.33801933438930586, 0.15212234596215889,
    0.15002124973437721
  )

  value <- loggompertz_fn(x, theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loggompertz")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)

  object <- structure(list(stats = lltd$stats_1), class = "loggompertz_fit")

  value <- fn(object, object$stats[, 1], theta)

  expect_type(value, "double")
  expect_length(value, m)
  expect_equal(value, true_value)
})

test_that("Gradient (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.0019304541362277092, 0.20961138715109782, 0.49935178859927617,
      0.67663384616172894, 0.77880078307140487, 0.99750312239746012,
      0.99997500031249740,
      # eta
      0, 0.0093970540520846307, 0.062120960821991612, -0.053740420946163813,
      -0.10559269877403762, -0.11471250798831215, -0.0063500361305660149,
      -0.00011258642934322865,
      # phi
      0, 0.0041022150394838821, 0.11135604942402072, 0.11790250564149576,
      0.089865432693354624, 0.066198066561069414, 0.00084787765403784111,
      8.4997875026562279e-06
    ),
    nrow = m,
    ncol = 4
  )

  G <- loggompertz_gradient(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

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
      0, -0.011055357708334860, -0.073083483319990132, 0.063224024642545663,
      0.12422670444004426, 0.13495589175095547, 0.0074706307418423704,
      0.00013245462275673958,
      # (delta, phi)
      0, -0.0048261353405692731, -0.13100711696943614, -0.13870883016646560,
      -0.10572403846277015, -0.077880078307140487, -0.00099750312239746012,
      -9.9997500031249740e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.011055357708334860, -0.073083483319990132, 0.063224024642545663,
      0.12422670444004426, 0.13495589175095547, 0.0074706307418423704,
      0.00013245462275673958,
      # (eta, eta)
      0, -0.045204776057939509, -0.0077973141424894308, 0.0029938447029538161,
      0.030242642409067612, 0.059634488615294076, 0.018975450654134089,
      0.00059650372086101613,
      # (eta, phi)
      0, -0.017682705989635783, 0.041700808527062246, 0.052382979149550082,
      0.019194496020505641, -0.0013147191159589372, -0.0021097255890769194,
      -0.000040783552121669912,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.0048261353405692731, -0.13100711696943614, -0.13870883016646560,
      -0.10572403846277015, -0.077880078307140487, -0.00099750312239746012,
      -9.9997500031249740e-06,
      # (phi, eta)
      0, -0.017682705989635783, 0.041700808527062246, 0.052382979149550082,
      0.019194496020505641, -0.0013147191159589372, -0.0021097255890769194,
      -0.000040783552121669912,
      # (phi, phi)
      0, -0.0094350945908129289, -0.047326321005208805, -0.0091701948832274482,
      0.0039316126803342648, 0.0066198066561069414, 0.00016872765315353038,
      1.6998725026562190e-06
    ),
    dim = c(m, 4, 4)
  )

  H <- loggompertz_hessian(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (1)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.0019304541362277092, 0.20961138715109782, 0.49935178859927617,
      0.67663384616172894, 0.77880078307140487, 0.99750312239746012,
      0.99997500031249740,
      # eta
      0, 0.0093970540520846307, 0.062120960821991612, -0.053740420946163813,
      -0.10559269877403762, -0.11471250798831215, -0.0063500361305660149,
      -0.00011258642934322865,
      # phi
      0, 0.0041022150394838821, 0.11135604942402072, 0.11790250564149576,
      0.089865432693354624, 0.066198066561069414, 0.00084787765403784111,
      8.4997875026562279e-06
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
      0, -0.011055357708334860, -0.073083483319990132, 0.063224024642545663,
      0.12422670444004426, 0.13495589175095547, 0.0074706307418423704,
      0.00013245462275673958,
      # (delta, phi)
      0, -0.0048261353405692731, -0.13100711696943614, -0.13870883016646560,
      -0.10572403846277015, -0.077880078307140487, -0.00099750312239746012,
      -9.9997500031249740e-06,
      # (eta, alpha)
      rep(0, m),
      # (eta, delta)
      0, -0.011055357708334860, -0.073083483319990132, 0.063224024642545663,
      0.12422670444004426, 0.13495589175095547, 0.0074706307418423704,
      0.00013245462275673958,
      # (eta, eta)
      0, -0.045204776057939509, -0.0077973141424894308, 0.0029938447029538161,
      0.030242642409067612, 0.059634488615294076, 0.018975450654134089,
      0.00059650372086101613,
      # (eta, phi)
      0, -0.017682705989635783, 0.041700808527062246, 0.052382979149550082,
      0.019194496020505641, -0.0013147191159589372, -0.0021097255890769194,
      -0.000040783552121669912,
      # (phi, alpha)
      rep(0, m),
      # (phi, delta)
      0, -0.0048261353405692731, -0.13100711696943614, -0.13870883016646560,
      -0.10572403846277015, -0.077880078307140487, -0.00099750312239746012,
      -9.9997500031249740e-06,
      # (phi, eta)
      0, -0.017682705989635783, 0.041700808527062246, 0.052382979149550082,
      0.019194496020505641, -0.0013147191159589372, -0.0021097255890769194,
      -0.000040783552121669912,
      # (phi, phi)
      0, -0.0094350945908129289, -0.047326321005208805, -0.0091701948832274482,
      0.0039316126803342648, 0.0066198066561069414, 0.00016872765315353038,
      1.6998725026562190e-06
    ),
    dim = c(m, 4, 4)
  )

  gh <- loggompertz_gradient_hessian(x, theta)

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
  theta <- lltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.0019304541362277092, 0.20961138715109782, 0.49935178859927617,
      0.67663384616172894, 0.77880078307140487, 0.99750312239746012,
      0.99997500031249740,
      # log_eta
      0, 0.018794108104169261, 0.12424192164398322, -0.10748084189232763,
      -0.21118539754807525, -0.22942501597662429, -0.012700072261132030,
      -0.00022517285868645729,
      # log_phi
      0, 0.020511075197419411, 0.55678024712010359, 0.58951252820747881,
      0.44932716346677312, 0.33099033280534707, 0.0042393882701892055,
      0.000042498937513281139
    ),
    nrow = m,
    ncol = 4
  )

  G <- loggompertz_gradient_2(x, theta)

  expect_type(G, "double")
  expect_length(G, m * 4)
  expect_equal(G, true_gradient)
})

test_that("Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

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
      0, -0.022110715416669719, -0.14616696663998026, 0.12644804928509133,
      0.24845340888008853, 0.26991178350191093, 0.014941261483684741,
      0.00026490924551347917,
      # (delta, log_phi)
      0, -0.024130676702846366, -0.65503558484718070, -0.69354415083232801,
      -0.52862019231385073, -0.38940039153570243, -0.0049875156119873006,
      -0.000049998750015624870,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.022110715416669719, -0.14616696663998026, 0.12644804928509133,
      0.24845340888008853, 0.26991178350191093, 0.014941261483684741,
      0.00026490924551347917,
      # (log_eta, log_eta)
      0, -0.16202499612758877, 0.093052665074025501, -0.095505463080512362,
      -0.090214827911804800, 0.0091129384845520113, 0.063201730355404326,
      0.0021608420247576072,
      # (log_eta, log_phi)
      0, -0.17682705989635783, 0.41700808527062246, 0.52382979149550082,
      0.19194496020505641, -0.013147191159589372, -0.021097255890769194,
      -0.00040783552121669912,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.024130676702846366, -0.65503558484718070, -0.69354415083232801,
      -0.52862019231385073, -0.38940039153570243, -0.0049875156119873006,
      -0.000049998750015624870,
      # (log_phi, log_eta)
      0, -0.17682705989635783, 0.41700808527062246, 0.52382979149550082,
      0.19194496020505641, -0.013147191159589372, -0.021097255890769194,
      -0.00040783552121669912,
      # (log_phi, log_phi)
      0, -0.21536628957290381, -0.62637777801011654, 0.36025765612679261,
      0.54761748047512974, 0.49648549920802060, 0.0084575795990274650,
      0.000084995750079686615
    ),
    dim = c(m, 4, 4)
  )

  H <- loggompertz_hessian_2(x, theta)

  expect_type(H, "double")
  expect_length(H, m * 4 * 4)
  expect_equal(H, true_hessian)
})

test_that("Gradient and Hessian (2)", {
  x <- lltd$stats_1[, 1]
  theta <- lltd$theta_g

  m <- length(x)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, m),
      # delta
      0, 0.0019304541362277092, 0.20961138715109782, 0.49935178859927617,
      0.67663384616172894, 0.77880078307140487, 0.99750312239746012,
      0.99997500031249740,
      # log_eta
      0, 0.018794108104169261, 0.12424192164398322, -0.10748084189232763,
      -0.21118539754807525, -0.22942501597662429, -0.012700072261132030,
      -0.00022517285868645729,
      # log_phi
      0, 0.020511075197419411, 0.55678024712010359, 0.58951252820747881,
      0.44932716346677312, 0.33099033280534707, 0.0042393882701892055,
      0.000042498937513281139
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
      0, -0.022110715416669719, -0.14616696663998026, 0.12644804928509133,
      0.24845340888008853, 0.26991178350191093, 0.014941261483684741,
      0.00026490924551347917,
      # (delta, log_phi)
      0, -0.024130676702846366, -0.65503558484718070, -0.69354415083232801,
      -0.52862019231385073, -0.38940039153570243, -0.0049875156119873006,
      -0.000049998750015624870,
      # (log_eta, alpha)
      rep(0, m),
      # (log_eta, delta)
      0, -0.022110715416669719, -0.14616696663998026, 0.12644804928509133,
      0.24845340888008853, 0.26991178350191093, 0.014941261483684741,
      0.00026490924551347917,
      # (log_eta, log_eta)
      0, -0.16202499612758877, 0.093052665074025501, -0.095505463080512362,
      -0.090214827911804800, 0.0091129384845520113, 0.063201730355404326,
      0.0021608420247576072,
      # (log_eta, log_phi)
      0, -0.17682705989635783, 0.41700808527062246, 0.52382979149550082,
      0.19194496020505641, -0.013147191159589372, -0.021097255890769194,
      -0.00040783552121669912,
      # (log_phi, alpha)
      rep(0, m),
      # (log_phi, delta)
      0, -0.024130676702846366, -0.65503558484718070, -0.69354415083232801,
      -0.52862019231385073, -0.38940039153570243, -0.0049875156119873006,
      -0.000049998750015624870,
      # (log_phi, log_eta)
      0, -0.17682705989635783, 0.41700808527062246, 0.52382979149550082,
      0.19194496020505641, -0.013147191159589372, -0.021097255890769194,
      -0.00040783552121669912,
      # (log_phi, log_phi)
      0, -0.21536628957290381, -0.62637777801011654, 0.36025765612679261,
      0.54761748047512974, 0.49648549920802060, 0.0084575795990274650,
      0.000084995750079686615
    ),
    dim = c(m, 4, 4)
  )

  gh <- loggompertz_gradient_hessian_2(x, theta)

  expect_type(gh, "list")
  expect_type(gh$G, "double")
  expect_type(gh$H, "double")

  expect_length(gh$G, m * 4)
  expect_length(gh$H, m * 4 * 4)

  expect_equal(gh$G, true_gradient)
  expect_equal(gh$H, true_hessian)

  object <- structure(list(stats = lltd$stats_1), class = "loggompertz")

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
  theta <- lltd$theta_g
  theta[3:4] <- log(theta[3:4])

  true_value <- 0.32075900098848013

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
    class = "loggompertz"
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
  theta <- lltd$theta_g
  theta[3:4] <- log(theta[3:4])

  true_gradient <- c(
    2.2653108698938061, 0.73068466091804248, -0.075702268537242359,
    0.57850423102806821
  )

  true_hessian <- matrix(
    c(
      # alpha
      19, 8.2261048577719928, -1.0136637259706368, 5.2206369313459388,
      # delta
      8.2261048577719928, 6.0552995869783359, -0.70635145977914647,
      1.9088026101488051,
      # log_eta
      -1.0136637259706368, -0.70635145977914647, 0.17854179243952868,
      -0.20592244468846133,
      # log_phi
      5.2206369313459388, 1.9088026101488051, -0.20592244468846133,
      2.7235196164555188
    ),
    nrow = 4,
    ncol = 4
  )

  object <- structure(
    list(stats = lltd$stats_1, m = nrow(lltd$stats_1)),
    class = "loggompertz"
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

  theta <- c(0, 1, 0.63524784841453301, 1.6101961281319852)

  true_value <- c(
    0.83777215113750696, -0.75729507682865963, 0.63524784841453301,
    1.6101961281319852
  )

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

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

  rss_value <- 0.073630385904407914

  theta <- c(
    alpha = 0.83777215113750696, delta = -0.75729507682865963,
    eta = exp(0.63524784841453301), phi = exp(1.6101961281319852)
  )

  fitted_values <- rep(
    c(
      0.83777215113750696, 0.83509727314890451, 0.67312406791826807,
      0.46539726317063078, 0.33641138140141234, 0.26005233684193214,
      0.083128285791155688, 0.080511485993225808
    ),
    k
  )

  residuals <- c(
    0.014927848862493037, -0.080672151137506963, 0.099227848862493037,
    -0.0075972731489045085, -0.057297273148904508, 0.049502726851095492,
    -0.11702406791826807, 0.031675932081731933, 0.079902736829369223,
    0.020602736829369223, 0.069702736829369223, -0.025497263170630777,
    0.018888618598587661, -0.023211381401412339, 0.013488618598587661,
    -0.11795233684193214, -0.066328285791155688, 0.055471714208844312,
    0.042188514006774192
  )

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loggompertz_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.073630385904407914

  theta <- c(
    alpha = 0.83777215113750696, delta = -0.75729507682865963,
    eta = exp(0.63524784841453301), phi = exp(1.6101961281319852)
  )

  fitted_values <- rep(
    c(
      0.83777215113750696, 0.83509727314890451, 0.67312406791826807,
      0.46539726317063078, 0.33641138140141234, 0.26005233684193214,
      0.083128285791155688, 0.080511485993225808
    ),
    k
  )

  residuals <- c(
    0.014927848862493037, -0.080672151137506963, 0.099227848862493037,
    -0.0075972731489045085, -0.057297273148904508, 0.049502726851095492,
    -0.11702406791826807, 0.031675932081731933, 0.079902736829369223,
    0.020602736829369223, 0.069702736829369223, -0.025497263170630777,
    0.018888618598587661, -0.023211381401412339, 0.013488618598587661,
    -0.11795233684193214, -0.066328285791155688, 0.055471714208844312,
    0.042188514006774192
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loggompertz_new(
    x, y, w, c(0.6, -0.6, 4, 8), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loggompertz_new(
    x, y, w, c(-2, 2, 7, 1), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.17654034036023363

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(0.34256135717772634),
    phi = exp(1.8046122081510387)
  )

  fitted_values <- rep(
    c(
      0.8, 0.79248376718906222, 0.65161512055560499, 0.47490216690684499,
      0.34359567627384339, 0.25186908386433202, -0.082746689754940315,
      -0.099320242477951822
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.035016232810937784, -0.014683767189062216,
    0.092116232810937784, -0.095515120555604989, 0.053184879444395011,
    0.070397833093155009, 0.011097833093155009, 0.060197833093155009,
    -0.035002166906844991, 0.011704323726156605, -0.030395676273843395,
    0.0063043237261566051, -0.10976908386433202, 0.099546689754940315,
    0.22134668975494032, 0.22202024247795182
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  # initial values with same equalities
  object <- loggompertz_new(
    x, y, w, c(0.8, -0.9, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  # initial values with different equalities
  object <- loggompertz_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

test_that("fit_constrained: equalities and inequalities", {
  x <- lltd$D$x
  y <- lltd$D$y

  n <- length(y)
  w <- rep(1, n)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.17654034036023363

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(0.34256135717772634),
    phi = exp(1.8046122081510387)
  )

  fitted_values <- rep(
    c(
      0.8, 0.79248376718906222, 0.65161512055560499, 0.47490216690684499,
      0.34359567627384339, 0.25186908386433202, -0.082746689754940315,
      -0.099320242477951822
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.035016232810937784, -0.014683767189062216,
    0.092116232810937784, -0.095515120555604989, 0.053184879444395011,
    0.070397833093155009, 0.011097833093155009, 0.060197833093155009,
    -0.035002166906844991, 0.011704323726156605, -0.030395676273843395,
    0.0063043237261566051, -0.10976908386433202, 0.099546689754940315,
    0.22134668975494032, 0.22202024247795182
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0.8, -0.9, 3, 7), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 8, 1), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

test_that("fit (weighted)", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.035138156284360913

  theta <- c(
    alpha = 0.84575428432937011, delta = -0.76656134360156587,
    eta = exp(0.89238310459026989), phi = exp(1.6622042309903930)
  )

  fitted_values <- rep(
    c(
      0.84575428432937011, 0.84573608713048233, 0.73789083646643889,
      0.47593010582855891, 0.31155933858113756, 0.22406801104864969,
      0.079774485527624186, 0.079195048422139125
    ),
    k
  )

  residuals <- c(
    0.0069457156706298950, -0.088654284329370105, 0.091245715670629895,
    -0.018236087130482331, -0.067936087130482331, 0.038863912869517669,
    -0.18179083646643889, -0.033090836466438890, 0.069369894171441088,
    0.010069894171441088, 0.059169894171441088, -0.036030105828558912,
    0.043740661418862440, 0.0016406614188624399, 0.038340661418862440,
    -0.081968011048649687, -0.062974485527624186, 0.058825514472375814,
    0.043504951577860875
  )

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  object <- loggompertz_new(x, y, w, c(0, 1, 1, 1), max_iter, NULL, NULL)

  result <- fit(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = TRUE, delta = TRUE, eta = TRUE, phi = TRUE)

  rss_value <- 0.035138156284360913

  theta <- c(
    alpha = 0.84575428432937011, delta = -0.76656134360156587,
    eta = exp(0.89238310459026989), phi = exp(1.6622042309903930)
  )

  fitted_values <- rep(
    c(
      0.84575428432937011, 0.84573608713048233, 0.73789083646643889,
      0.47593010582855891, 0.31155933858113756, 0.22406801104864969,
      0.079774485527624186, 0.079195048422139125
    ),
    k
  )

  residuals <- c(
    0.0069457156706298950, -0.088654284329370105, 0.091245715670629895,
    -0.018236087130482331, -0.067936087130482331, 0.038863912869517669,
    -0.18179083646643889, -0.033090836466438890, 0.069369894171441088,
    0.010069894171441088, 0.059169894171441088, -0.036030105828558912,
    0.043740661418862440, 0.0016406614188624399, 0.038340661418862440,
    -0.081968011048649687, -0.062974485527624186, 0.058825514472375814,
    0.043504951577860875
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values within the boundaries
  object <- loggompertz_new(
    x, y, w, c(0.6, -0.6, 4, 8), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w)

  # initial values outside the boundaries
  object <- loggompertz_new(
    x, y, w, c(-2, 2, 7, 1), max_iter,
    c(0.5, -1, 1, 3), c(1, -0.5, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
  expect_true(inherits(result, "loglogistic"))
  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss_value)
  expect_equal(result$df.residual, object$n - 4)
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

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.12152790634234742

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(0.67344700327430744),
    phi = exp(1.8267797417585254)
  )

  fitted_values <- rep(
    c(
      0.8, 0.79991222146934487, 0.71604576658544182, 0.49162778848529489,
      0.31063536742142226, 0.19268868846616760, -0.096135387777275743,
      -0.099957631352575174
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.027587778530655130, -0.022112221469344870,
    0.084687778530655130, -0.15994576658544182, -0.011245766585441821,
    0.053672211514705107, -0.0056277884852948934, 0.043472211514705107,
    -0.051727788485294893, 0.044664632578577739, 0.0025646325785777394,
    0.039264632578577739, -0.050588688466167597, 0.11293538777727574,
    0.23473538777727574, 0.22265763135257517
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  # initial values with same equalities
  object <- loggompertz_new(
    x, y, w, c(0.8, -0.9, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

  # initial values with different equalities
  object <- loggompertz_new(
    x, y, w, c(0, 1, 1, 1), max_iter,
    c(0.8, -0.9, rep(-Inf, 2)), c(0.8, -0.9, rep(Inf, 2))
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

test_that("fit_constrained (weighted): equalities and inequalities", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  n <- length(y)

  k <- as.numeric(table(x))

  max_iter <- 10000

  estimated <- c(alpha = FALSE, delta = FALSE, eta = TRUE, phi = TRUE)

  rss_value <- 0.12152790634234742

  theta <- c(
    alpha = 0.8, delta = -0.9, eta = exp(0.67344700327430744),
    phi = exp(1.8267797417585254)
  )

  fitted_values <- rep(
    c(
      0.8, 0.79991222146934487, 0.71604576658544182, 0.49162778848529489,
      0.31063536742142226, 0.19268868846616760, -0.096135387777275743,
      -0.099957631352575174
    ),
    k
  )

  residuals <- c(
    0.0527, -0.0429, 0.137, 0.027587778530655130, -0.022112221469344870,
    0.084687778530655130, -0.15994576658544182, -0.011245766585441821,
    0.053672211514705107, -0.0056277884852948934, 0.043472211514705107,
    -0.051727788485294893, 0.044664632578577739, 0.0025646325785777394,
    0.039264632578577739, -0.050588688466167597, 0.11293538777727574,
    0.23473538777727574, 0.22265763135257517
  )

  object <- loggompertz_new(
    x, y, w, NULL, max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0.8, -0.9, 3, 7), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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
  object <- loggompertz_new(
    x, y, w, c(0, 1, 8, 1), max_iter,
    c(0.8, -0.9, 1, 3), c(0.8, -0.9, 5, 12)
  )

  result <- fit_constrained(object)

  expect_true(inherits(result, "loggompertz_fit"))
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

test_that("fisher_info", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  max_iter <- 10000

  theta <- lltd$theta_g
  names(theta) <- c("alpha", "delta", "eta", "phi")

  sigma <- lltd$sigma

  true_value <- matrix(c(
      # alpha
      6206.96, 2748.9449629599373, -189.18362192800753, 344.04436069543435,
      -28355.127259362132,
      # delta
      2748.9449629599373, 2023.9485432497464, -124.69711327516109,
      134.06347723771206, -10502.299351663864,
      # eta
      -189.18362192800753, -124.69711327516109, 23.773250144905426,
      -11.783571918074014, 850.51618712697142,
      # phi
      344.04436069543435, 134.06347723771206, -11.783571918074014,
      29.686894742510176, -1441.3142462178932,
      # sigma
      -28355.127259362132, -10502.299351663864, 850.51618712697142,
      -1441.3142462178932, 113909.14953385276
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "delta", "eta", "phi", "sigma"
  )

  object <- loggompertz_new(x, y, w, NULL, max_iter, NULL, NULL)

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
      y ~ x, mean_function = "loggompertz",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
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
      y ~ x, mean_function = "loggompertz",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
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
      y ~ x, mean_function = "loggompertz",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(0, Inf, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(-Inf, 1, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = rep(1, 5)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(0, 1, -1, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(0, 1, 0, 1)
    ),
    "parameter 'eta' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(0, 1, 1, -1)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "loggompertz",
      start = c(0, 1, 1, 0)
    ),
    "parameter 'phi' cannot be negative nor zero"
  )
})

test_that("nauc: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loggompertz")

  expect_equal(nauc(result), 0.085299967613769906, tolerance = 1.0e-7)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.84575366076886644)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.0059283502711447546)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.33844687600885550
  )
  expect_equal(nauc(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 0.0)
})

test_that("naac: decreasing", {
  x <- lltd$D$x
  y <- lltd$D$y
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loggompertz")

  expect_equal(naac(result), 1 - 0.085299967613769906)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.84575366076886644)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.0059283502711447546)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.33844687600885550
  )
  expect_equal(naac(result, xlim = c(10, 15), ylim = c(0.3, 0.7)), 1.0)
})

test_that("nauc: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loggompertz")

  expect_equal(nauc(result), 0.88478280795251264)
  expect_equal(nauc(result, xlim = c(0, 2)), 0.129049490857619206)
  expect_equal(nauc(result, ylim = c(0.3, 0.7)), 0.99551714173529925)
  expect_equal(nauc(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 0.0)
  expect_equal(
    nauc(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 0.85655631688941921
  )
  expect_equal(nauc(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 1.0)
})

test_that("naac: increasing", {
  x <- lltd$D$x
  y <- rev(lltd$D$y)
  w <- lltd$D$w

  result <- drda(y ~ x, weights = w, mean_function = "loggompertz")

  expect_equal(naac(result), 1 - 0.88478280795251264)
  expect_equal(naac(result, xlim = c(0, 2)), 1 - 0.129049490857619206)
  expect_equal(naac(result, ylim = c(0.3, 0.7)), 1 - 0.99551714173529925)
  expect_equal(naac(result, xlim = c(0, 2), ylim = c(0.3, 0.7)), 1.0)
  expect_equal(
    naac(result, xlim = c(5, 8), ylim = c(0.3, 0.7)), 1 - 0.85655631688941921
  )
  expect_equal(naac(result, xlim = c(9, 12), ylim = c(0.3, 0.7)), 0.0)
})
