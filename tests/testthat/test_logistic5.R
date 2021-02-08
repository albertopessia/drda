context("5-parameter logistic - core functions")

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2)

  true_value <- c(
    0.89998272669845415, 0.89827524246036657, 0.75019144346736942,
    0.047052490645792413, 0.040000850995162837, 0.040000000085267377
  )

  value <- logistic5_function(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # log_omega
      0.85998272669845415, 0.85827524246036657, 0.71019144346736942,
      0.0070524906457924134, 8.5099516283666758e-07, 8.5267376783781469e-11,
      # eta
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # phi
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # log_nu
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09
    ),
    nrow = 6,
    ncol = 5
  )

  true_hessian <- array(
    c(
      # (alpha, ...)
      rep(0, 30),
      # (log_omega, alpha)
      rep(0, 6),
      # (log_omega, log_omega)
      0.85998272669845415, 0.85827524246036657, 0.71019144346736942,
      0.0070524906457924134, 8.5099516283666758e-07, 8.5267376783781469e-11,
      # (log_omega, eta)
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # (log_omega, phi)
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # (log_omega, log_nu)
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09,
      # (eta, alpha)
      rep(0, 6),
      # (eta, log_omega)
      -0.000093408380497224582, -0.0053476072761497761, -0.10403715378083063,
      0.019241514719677547, 6.4655250500714644e-06, 1.0411333261602936e-09,
      # (eta, eta)
      -0.00050511444418713689, -0.016555252126485643, -0.060637799043344759,
      0.049883501712141417, 0.000049098048382061985, 1.2712402410105175e-08,
      # (eta, phi)
      0.00016953809123729721, 0.0089408616320290739, 0.021478649997785072,
      0.053683659136403524, 0.000024123213411596926, 3.9939380477898473e-09,
      # (eta, log_nu)
      9.3805989611171948e-10, 5.3597039009804577e-06, 0.0085715606355158565,
      0.039930425361417561, 0.000070015304961175913, 2.0858519166902048e-08,
      # (phi, alpha)
      rep(0, 6),
      # (phi, log_omega)
      0.000034546082682501184, 0.0034443247589330582, 0.25925513615688965,
      0.025655352959570063, 3.4005945386908866e-06, 3.4106611099876862e-10,
      # (phi, eta)
      0.00016953809123729721, 0.0089408616320290739, 0.021478649997785072,
      0.053683659136403524, 0.000024123213411596926, 3.9939380477898473e-09,
      # (phi, phi)
      -0.000069090083756049672, -0.0068679160064153061, -0.37654877818008748,
      0.088681780821584741, 0.000013582081688859556, 1.3642440673798294e-09,
      # (phi, log_nu)
      -3.4693134127485286e-10, -3.4521160387059138e-06, -0.021359879993633144,
      0.053240567148556747, 0.000036825108839864453, 6.8330672303858801e-09,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, log_omega)
      8.6734287058557399e-11, 8.6447455933317324e-07, 0.0063015237621949961,
      0.021049325905316226, 0.000010065592915696157, 1.7935503452654302e-09,
      # (log_nu, eta)
      9.3805989611171948e-10, 5.3597039009804577e-06, 0.0085715606355158565,
      0.039930425361417561, 0.000070015304961175913, 2.0858519166902048e-08,
      # (log_nu, phi)
      -3.4693134127485286e-10, -3.4521160387059138e-06, -0.021359879993633144,
      0.053240567148556747, 0.000036825108839864453, 6.8330672303858801e-09,
      # (log_nu, log_nu)
      8.6733125674846536e-11, 8.6331893228632990e-07, 0.0055845141256710526,
      0.053441912635589904, 0.00011068910773858665, 3.6103283407515775e-08
    ),
    dim = c(6, 5, 5)
  )

  gradient_hessian <- logistic5_gradient_hessian(x, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 5)
  expect_length(gradient_hessian$H, 6 * 5 * 5)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

context("5-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2, -log(2), log(3 / 2))

  n <- c(3, 3, 2, 4, 3)

  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_value <- 0.13844046588658472

  rss <- logistic5_rss(stats)

  expect_type(rss, "closure")

  value <- rss(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2))
  rss <- logistic5_rss_fixed(stats, known_param)

  expect_type(rss, "closure")

  value <- rss(c(log(43 / 50), -2))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2, -log(2), log(3 / 2))

  n <- c(3, 3, 2, 4, 3)

  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_gradient <- c(
    -0.92903069002014085, -0.28833791237257394, 0.022267352061200450,
    -0.086374079772717690, -0.010097206230551137
  )

  true_hessian <- matrix(
    c(
      # alpha
      15.000000000000000, 6.6033693099798592, -0.14741189907774185,
      0.63157849846052230, 0.096833141608282843,
      # log_omega
      6.6033693099798592, 5.1492248569956934, -0.13897258392352741,
      0.29154887863609586, -0.00055060290012352080,
      # eta
      -0.14741189907774185, -0.13897258392352741, 0.018237239090142242,
      -0.077452292586773564, -0.017846532854813774,
      # phi
      0.63157849846052230, 0.29154887863609586, -0.077452292586773564,
      0.21294298805991181, -0.0090213900412146364,
      # log_nu
      0.096833141608282843, -0.00055060290012352080, -0.017846532854813774,
      -0.0090213900412146364, -0.020700058407993170
    ),
    nrow = 5,
    ncol = 5
  )

  gh <- logistic5_rss_gradient_hessian(stats)

  expect_type(gh, "closure")

  gradient_hessian <- gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 5)
  expect_length(gradient_hessian$H, 5 * 5)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2))
  gh <- logistic5_rss_gradient_hessian_fixed(stats, known_param)

  expect_type(gh, "closure")

  gradient_hessian <- gh(c(log(43 / 50), -2))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[2:3])
  expect_equal(gradient_hessian$H, true_hessian[2:3, 2:3])
})

context("5-parameter logistic - general functions")

test_that("logistic5_init", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(
    4 / 100, log(43 / 50), -2, -3 / 2, -log(2), log(3 / 2)
  )

  n <- c(3, 3, 2, 4, 3)

  m <- c(
    376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120
  )

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  true_value <- c(
    0.11297166830415549, -0.14005886043338885, -1.5025, -2.3025850929940460, 0.5
  )

  stats <- cbind(x, n, m, v)

  start <- logistic5_init(stats)

  expect_type(start, "double")
  expect_length(start, 5)
  expect_equal(start, true_value)
})

test_that("logistic5_fisher_info_normal", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  stats <- suff_stats(x, y)

  n <- length(y)

  theta <- c(
    alpha = 4 / 100,
    beta = 9 / 10,
    eta = -2,
    phi = -3 / 2,
    nu = 1 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5191.5736239671087, 132.90987565500003, 23.679974158167395,
      86.980812626015735, 85.283559876424630, 0,
      # beta
      132.90987565500003, 2542.6066247228912, -72.803115749724175,
      174.41094047950246, 9.0075218161025453, 0,
      # eta
      23.679974158167395, -72.803115749724175, 9.4208692315701512,
      -20.590332645351424, 0.57269249241964757, 0,
      # phi
      86.980812626015735, 174.41094047950246, -20.590332645351424,
      55.031696576965175, 4.7687936450855367, 0,
      # nu
      85.283559876424630, 9.0075218161025453, 0.57269249241964757,
      4.7687936450855367, 3.6713059622893367, 0,
      # sigma
      rep(0, 5), 6800
    ),
    nrow = 6,
    ncol = 6
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "nu", "sigma"
  )

  fim <- logistic5_fisher_info_normal(stats, n, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 6 * 6)
  expect_equal(fim, true_value)
})

context("5-parameter logistic - fit")

test_that("logistic5_fit: fit", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.093212121358460132,
    beta = 0.093212121358460132 + exp(-0.18985173265400376),
    eta = -1.7617462932355769,
    phi = -0.47836972294568200,
    nu = exp(1.3765334489390746)
  )

  rss <- 0.024883087882351184

  fitted_values <- c(
    rep(0.92028391866853350, 3), rep(0.91971916938223621, 2),
    rep(0.89002742082078101, 2), rep(0.55337929341255560, 5),
    rep(0.26271803120343567, 3), rep(0.15412982230606348, 4),
    0.11508521369102614
  )

  residuals <- c(
    0.00771608133146650, -0.03228391866853350, 0.05971608133146650,
    0.02828083061776379, -0.06371916938223621, 0.00697257917921899,
    -0.00702742082078101, -0.06537929341255560, -0.02137929341255560,
    0.03262070658744440, 0.01262070658744440, 0.04562070658744440,
    -0.00371803120343567, 0.00228196879656433, -0.01971803120343567,
    -0.03712982230606348, -0.01112982230606348, 0.02387017769393652,
    0.06487017769393652, -0.02308521369102614
  )

  result <- logistic5_fit(x, y, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  result <- logistic5_fit(x, y, c(0, 0, -1, 0, 0), max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic5_fit_constrained: inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.093212121358461020,
    beta = 0.093212121358461020 + exp(-0.18985173265400473),
    eta = -1.7617462932355648,
    phi = -0.47836972294567659,
    nu = exp(1.3765334489390619)
  )

  rss <- 0.024883087882351184

  fitted_values <- c(
    rep(0.92028391866853359, 3), rep(0.91971916938223627, 2),
    rep(0.89002742082078076, 2), rep(0.55337929341255556, 5),
    rep(0.26271803120343537, 3), rep(0.15412982230606358, 4),
    0.11508521369102661
  )

  residuals <- c(
    0.00771608133146641, -0.03228391866853359, 0.05971608133146641,
    0.02828083061776373, -0.06371916938223627, 0.00697257917921924,
    -0.00702742082078076, -0.06537929341255556, -0.02137929341255556,
    0.03262070658744444, 0.01262070658744444, 0.04562070658744444,
    -0.00371803120343537, 0.00228196879656463, -0.01971803120343537,
    -0.03712982230606358, -0.01112982230606358, 0.02387017769393642,
    0.06487017769393642, -0.02308521369102661
  )

  result <- logistic5_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, -0.5, -3, -2, 0),
    upper_bound = c(1,  0.5,  0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)

  # initial values within the boundaries
  result <- logistic5_fit_constrained(
    x, y, c(0.5, 0, -1.5, -1, 1), max_iter,
    lower_bound = c(0, -0.5, -3, -2, 0),
    upper_bound = c(1,  0.5,  0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)

  # initial values outside the boundaries
  result <- logistic5_fit_constrained(
    x, y, c(2, -1, 1, 1, -1), max_iter,
    lower_bound = c(0, -0.5, -3, -2, 0),
    upper_bound = c(1,  0.5,  0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
})

test_that("logistic5_fit_constrained: equalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605569,
    phi = -0.30897155961772666,
    nu = exp(0.85096464960136723)
  )

  rss <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706217, 3), rep(0.96390948874813836, 2),
    rep(0.83729075842125315, 2), rep(0.55558345438725122, 5),
    rep(0.29108572543250873, 3), rep(0.14085756611083844, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706217, -0.105387436096706217, -0.013387436096706217,
    -0.01590948874813836, -0.10790948874813836, 0.05970924157874685,
    0.04570924157874685, -0.06758345438725122, -0.02358345438725122,
    0.03041654561274878, 0.01041654561274878, 0.04341654561274878,
    -0.03208572543250873, -0.02608572543250873, -0.04808572543250873,
    -0.02385756611083844, 0.00214243388916156, 0.03714243388916156,
    0.07814243388916156, 0.02497284742870274
  )

  result <- logistic5_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  result <- logistic5_fit_constrained(
    x, y, c(0, 0, -1, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with different equalities
  result <- logistic5_fit_constrained(
    x, y, c(1, 3, -1, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic5_fit_constrained: equalities and inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.75885255907605569,
    phi = -0.30897155961772666,
    nu = exp(0.85096464960136723)
  )

  rss <- 0.053861132351488352

  fitted_values <- c(
    rep(0.993387436096706217, 3), rep(0.96390948874813836, 2),
    rep(0.83729075842125315, 2), rep(0.55558345438725122, 5),
    rep(0.29108572543250873, 3), rep(0.14085756611083844, 4),
    0.06702715257129726
  )

  residuals <- c(
    -0.065387436096706217, -0.105387436096706217, -0.013387436096706217,
    -0.01590948874813836, -0.10790948874813836, 0.05970924157874685,
    0.04570924157874685, -0.06758345438725122, -0.02358345438725122,
    0.03041654561274878, 0.01041654561274878, 0.04341654561274878,
    -0.03208572543250873, -0.02608572543250873, -0.04808572543250873,
    -0.02385756611083844, 0.00214243388916156, 0.03714243388916156,
    0.07814243388916156, 0.02497284742870274
  )

  result <- logistic5_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -3, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic5_fit_constrained(
    x, y, c(0, 0, -2, -1, 1), max_iter,
    lower_bound = c(0, 0, -3, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic5_fit_constrained(
    x, y, c(1, 2, -5, 2, -1), max_iter,
    lower_bound = c(0, 0, -3, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("5-parameter logistic - weighted fit")

test_that("logistic5_fit_weighted: fit", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  w <- c(
    1.46, 1.385, 1.704,
    0.96, 0,
    0.055, 1.071,
    0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269,
    0, 1.294, 0.038, 1.278,
    0.157
  )

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.14021510699415423,
    beta = 0.14021510699415423 + exp(-0.22629299711936212),
    eta = -1.3532016035342649,
    phi = -0.36746911119363781,
    nu = exp(0.88825453177528800)
  )

  rss <- 0.014141550871844299

  fitted_values <- c(
    rep(0.93758527225933412, 3), rep(0.93513517257012652, 1),
    rep(0.88595597201927951, 2), rep(0.55164531553540003, 4),
    rep(0.26479534106274697, 3), rep(0.17495286982311324, 3),
    0.14985594300774285
  )

  residuals <- c(
    -0.00958527225933412, -0.04958527225933412, 0.04241472774066588,
    0.01286482742987348, 0.01104402798072049, -0.00295597201927951,
    -0.06364531553540003, -0.01964531553540003, 0.01435468446459997,
    0.04735468446459997, -0.00579534106274697, 0.00020465893725303,
    -0.02179534106274697, -0.03195286982311324, 0.00304713017688676,
    0.04404713017688676, -0.05785594300774285
  )

  result <- logistic5_fit_weighted(x, y, w, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  result <- logistic5_fit_weighted(x, y, w, c(0, 0, -1, 0, 0), max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 5)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic5_fit_weighted_constrained: inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  w <- c(
    1.46, 1.385, 1.704,
    0.96, 0,
    0.055, 1.071,
    0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269,
    0, 1.294, 0.038, 1.278,
    0.157
  )

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0.046928317026565997,
    beta = 0.046928317026565997 + exp(-0.1),
    eta = -1.0006613313287804,
    phi = -0.32207641903196782,
    nu = exp(1)
  )

  rss <- 0.019185779018228248

  fitted_values <- c(
    rep(0.95052596409995778, 3), rep(0.93962179909607343, 1),
    rep(0.85184882797766230, 2), rep(0.55691822191298501, 4),
    rep(0.28293674769204406, 3), rep(0.14894441582152249, 3),
    0.09066655775367701
  )

  residuals <- c(
    -0.02252596409995778, -0.06252596409995778, 0.02947403590004222,
    0.00837820090392657, 0.04515117202233770, 0.03115117202233770,
    -0.06891822191298501, -0.02491822191298501, 0.00908177808701499,
    0.04208177808701499, -0.02393674769204406, -0.01793674769204406,
    -0.03993674769204406, -0.00594441582152249, 0.02905558417847751,
    0.07005558417847751, 0.00133344224632299
  )

  result <- logistic5_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0),
    upper_bound = c(1,  0.1,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values within the boundaries
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(0.5, 0, -1, -0.5, 0.5), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0),
    upper_bound = c(1,  0.1,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values outside the boundaries
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(2, -1, -3, 1, 2), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0),
    upper_bound = c(1,  0.1,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 5)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic5_fit_weighted_constrained: equalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  w <- c(
    1.46, 1.385, 1.704,
    0.96, 0,
    0.055, 1.071,
    0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269,
    0, 1.294, 0.038, 1.278,
    0.157
  )

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011306,
    phi = -0.33982624766042563,
    nu = exp(0.86063297591301602)
  )

  rss <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829322, 3), rep(0.95779593343449219, 1),
    rep(0.82640795793649349, 2), rep(0.55492427795777896, 4),
    rep(0.30125409492277097, 3), rep(0.15182795902080914, 3),
    0.07523599271751830
  )

  residuals <- c(
    -0.063572996129829322, -0.103572996129829322, -0.011572996129829322,
    -0.00979593343449219, 0.07059204206350651, 0.05659204206350651,
    -0.06692427795777896, -0.02292427795777896, 0.01107572204222104,
    0.04407572204222104, -0.04225409492277097, -0.03625409492277097,
    -0.05825409492277097, -0.00882795902080914, 0.02617204097919086,
    0.06717204097919086, 0.01676400728248170
  )

  result <- logistic5_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with same equalities
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(0, 0, -1, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with different equalities
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(1, 3, -1, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic5_fit_weighted_constrained: equalities and inequalities", {
  max_iter <- 10000

  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

  w <- c(
    1.46, 1.385, 1.704,
    0.96, 0,
    0.055, 1.071,
    0.134, 1.825, 0, 1.169, 0.628,
    0.327, 1.201, 0.269,
    0, 1.294, 0.038, 1.278,
    0.157
  )

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.72502089617011306,
    phi = -0.33982624766042563,
    nu = exp(0.86063297591301602)
  )

  rss <- 0.036717820758155562

  fitted_values <- c(
    rep(0.991572996129829322, 3), rep(0.95779593343449219, 1),
    rep(0.82640795793649349, 2), rep(0.55492427795777896, 4),
    rep(0.30125409492277097, 3), rep(0.15182795902080914, 3),
    0.07523599271751830
  )

  residuals <- c(
    -0.063572996129829322, -0.103572996129829322, -0.011572996129829322,
    -0.00979593343449219, 0.07059204206350651, 0.05659204206350651,
    -0.06692427795777896, -0.02292427795777896, 0.01107572204222104,
    0.04407572204222104, -0.04225409492277097, -0.03625409492277097,
    -0.05825409492277097, -0.00882795902080914, 0.02617204097919086,
    0.06717204097919086, 0.01676400728248170
  )

  result <- logistic5_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -2, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)

  # initial values within the boundaries
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(0, 0, -1, -1, 0.5), max_iter,
    lower_bound = c(0, 0, -2, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)

  # initial values outside the boundaries
  result <- logistic5_fit_weighted_constrained(
    x, y, w, c(1, 2, -5, 2, -1), max_iter,
    lower_bound = c(0, 0, -2, -2, 0),
    upper_bound = c(0, 0,  0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 3)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
})

context("5-parameter logistic - drda fit")

test_that("drda: 'lower_bound' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

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
    "'lower_bound' and 'upper_bound' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = c(0, -1, -Inf, -Inf, 0),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound[2]' cannot be smaller than 'lower_bound[1]'",
    fixed = TRUE
  )
})

test_that("drda: 'upper_bound' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

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
    "'lower_bound' and 'upper_bound' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      lower_bound = rep(-Inf, 5),
      upper_bound = c(1, 0, Inf, Inf, 0)
    ),
    "'upper_bound[2]' cannot be smaller than 'upper_bound[1]'",
    fixed = TRUE
  )
})

test_that("drda: 'start' argument errors", {
  x <- round(
    rep(
      -log(c(1000, 100, 10, 1, 0.1, 0.01, 0.001)),
      times = c(3, 2, 2, 5, 3, 4, 1)
    ),
    digits = 3
  )

  y <- c(
    0.928, 0.888, 0.98,
    0.948, 0.856,
    0.897, 0.883,
    0.488, 0.532, 0.586, 0.566, 0.599,
    0.259, 0.265, 0.243,
    0.117, 0.143, 0.178, 0.219,
    0.092
  )

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
      start = c(0, Inf, -1, 0, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(-Inf, 0, -1, 0, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 0, -1, 0, 1, 1)
    ),
    "'start' must be of length 5"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, -1, -1, 0, 1)
    ),
    "parameter 'beta' is smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, 0, 0, 1)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, -1, 0, 0)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic5",
      start = c(0, 1, -1, 0, -1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )
})
