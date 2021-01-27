context("four-parameter logistic fit")

test_that("logistic4_fisher_info", {
  stats <- matrix(
    c(
      -log(c(10000, 1000, 100, 10, 1, 0.1, 0.01)),
      rep(3, 7),
      c(0.85777, 0.90623, 0.82218, 0.51285, 0.11214, 0.07049, 0.0325),
      c(0.00101, 0.00447, 0.00412, 7e-05, 0.01114, 3e-04, 0.0028)
    ),
    nrow = 7,
    ncol = 4
  )
  colnames(stats) <- c("x", "n", "m", "v")
  param <- c(0.05, 0.9, -1.5, -2)
  sigma <- 0.07

  true_value <- matrix(
    c(
      1870.6712073826437, 186.31033844604499, 33.423679753536860,
      107.18330814627171, 0, 186.31033844604499, 2042.4224014395520,
      -47.989631045699760, 130.36237337243564, 0, 33.423679753536860,
      -47.989631045699760, 7.0397756134557709, -9.2598379346301116, 0,
      107.18330814627171, 130.36237337243564, -9.2598379346301116,
      58.570645838394806, 0, 0, 0, 0, 0, 8571.4285714285714
    ),
    nrow = 5,
    ncol = 5
  )

  fim <- logistic4_fisher_info(stats, param, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

test_that("logistic4_fit_unconstrained", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0.20363741036518857,
    "maximum" = 2.4851366482692220,
    "growth_rate" = -1.9444176076033746,
    "logx_midpoint" = -1.0185895274452179
  )
  sigma <- 0.050499678048221869
  loglik <- 41.350661959289913

  fitted_values <- c(
    rep(2.4851123759859562, 5), rep(2.4830029818125129, 5),
    rep(2.3115254653996470, 5), rep(0.48028924884678951, 5),
    rep(0.20720992799323027, 5)
  )

  residuals <- c(
    0.030127624014043775, -0.022402375985956225, -0.015272375985956225,
    0.015027624014043775, -0.0047823759859562254, -0.075142981812512868,
    0.078707018187487132, 0.015787018187487132, -0.035092981812512868,
    0.012977018187487132, -0.0095954653996469986, -0.078805465399646999,
    -0.052495465399646999, 0.063444534600353001, 0.077524534600353001,
    0.047424751153210485, 0.045497751153210485, -0.084873248846789515,
    -0.018497248846789515, 0.010423751153210485, -0.056137927993230272,
    0.037973072006769728, 0.014551072006769728, -0.032051927993230272,
    0.035684072006769728
  )

  result <- logistic4_fit_unconstrained(x, y, NULL, max_iter)

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  result <- logistic4_fit_unconstrained(x, y, c(0, 1, -1, 0), max_iter)

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: inequalities I", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0.20363741036518857,
    "maximum" = 2.4851366482692220,
    "growth_rate" = -1.9444176076033746,
    "logx_midpoint" = -1.0185895274452179
  )
  sigma <- 0.050499678048221869
  loglik <- 41.350661959289913

  fitted_values <- c(
    rep(2.4851123759859562, 5), rep(2.4830029818125129, 5),
    rep(2.3115254653996470, 5), rep(0.48028924884678951, 5),
    rep(0.20720992799323027, 5)
  )

  residuals <- c(
    0.030127624014043775, -0.022402375985956225, -0.015272375985956225,
    0.015027624014043775, -0.0047823759859562254, -0.075142981812512868,
    0.078707018187487132, 0.015787018187487132, -0.035092981812512868,
    0.012977018187487132, -0.0095954653996469986, -0.078805465399646999,
    -0.052495465399646999, 0.063444534600353001, 0.077524534600353001,
    0.047424751153210485, 0.045497751153210485, -0.084873248846789515,
    -0.018497248846789515, 0.010423751153210485, -0.056137927993230272,
    0.037973072006769728, 0.014551072006769728, -0.032051927993230272,
    0.035684072006769728
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 2, -2, -2),
    upper_bound = c(1, 3, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(0.5, 2.5, -1.5, -1), max_iter,
    lower_bound = c(0, 2, -2, -2),
    upper_bound = c(1, 3, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(2, 1, -0.5, 1), max_iter,
    lower_bound = c(0, 2, -2, -2),
    upper_bound = c(1, 3, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: inequalities II", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 2,
    "growth_rate" = -2,
    "logx_midpoint" = -0.50584907039612223
  )
  sigma <- 0.39519075737116401
  loglik <- -10.084378475791494

  fitted_values <- c(
    rep(1.9999944994803418, 5), rep(1.9994500977592080, 5),
    rep(1.9464669540461577, 5), rep(0.53329530680944199, 5),
    rep(0.0072456748394291598, 5)
  )

  residuals <- c(
    0.51524550051965822, 0.46271550051965822, 0.46984550051965822,
    0.50014550051965822, 0.48033550051965822, 0.40840990224079201,
    0.56225990224079201, 0.49933990224079201, 0.44845990224079201,
    0.49652990224079201, 0.35546304595384232, 0.28625304595384232,
    0.31256304595384232, 0.42850304595384232, 0.44258304595384232,
    -0.0055813068094419874, -0.0075083068094419874, -0.13787930680944199,
    -0.071503306809441987, -0.042582306809441987, 0.14382632516057084,
    0.23793732516057084, 0.21451532516057084, 0.16791232516057084,
    0.23564832516057084
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(-1, 1, -2, -1),
    upper_bound = c(0, 2, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(-0.5, 1.5, -1.5, -0.5), max_iter,
    lower_bound = c(-1, 1, -2, -1),
    upper_bound = c(0, 2, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(-2, 3, -3, 1), max_iter,
    lower_bound = c(-1, 1, -2, -1),
    upper_bound = c(0, 2, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: equalities", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 3,
    "growth_rate" = -0.89810058089211118,
    "logx_midpoint" = -1.3988129648706804
  )
  sigma <- 0.32357815255700665
  loglik <- -6.2233279732074851

  fitted_values <- c(
    rep(2.9788485634678064, 5), rep(2.8404909133496183, 5),
    rep(2.0774010564297024, 5), rep(0.66484752644712418, 5),
    rep(0.10424796908351360, 5)
  )

  residuals <- c(
    -0.46360856346780644, -0.51613856346780644, -0.50900856346780644,
    -0.47870856346780644, -0.49851856346780644, -0.43263091334961827,
    -0.27878091334961827, -0.34170091334961827, -0.39258091334961827,
    -0.34451091334961827, 0.22452894357029763, 0.15531894357029763,
    0.18162894357029763, 0.29756894357029763, 0.31164894357029763,
    -0.13713352644712418, -0.13906052644712418, -0.26943152644712418,
    -0.20305552644712418, -0.17413452644712418, 0.046824030916486399,
    0.14093503091648640, 0.11751303091648640, 0.070910030916486399,
    0.13864603091648640
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 3, -Inf, -Inf),
    upper_bound = c(0, 3, Inf, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  result <- logistic4_fit_constrained(
    x, y, c(0, 3, 1, 1), max_iter,
    lower_bound = c(0, 3, -Inf, -Inf),
    upper_bound = c(0, 3, Inf, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with different equalities
  result <- logistic4_fit_constrained(
    x, y, c(-1, 2, 1, 1), max_iter,
    lower_bound = c(0, 3, -Inf, -Inf),
    upper_bound = c(0, 3, Inf, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: equalities and inequalities", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 3,
    "growth_rate" = -1,
    "logx_midpoint" = -2
  )
  sigma <- 0.40018269077966560
  loglik <- -11.535340488553747

  fitted_values <- c(
    rep(2.9779954247442063, 5), rep(2.7935807511300708, 5),
    rep(1.7252230270189805, 5), rep(0.35760876606635267, 5),
    rep(0.040058452765874794, 5)
  )

  residuals <- c(
    -0.46275542474420625, -0.51528542474420625, -0.50815542474420625,
    -0.47785542474420625, -0.49766542474420625, -0.38572075113007083,
    -0.23187075113007083, -0.29479075113007083, -0.34567075113007083,
    -0.29760075113007083, 0.57670697298101955, 0.50749697298101955,
    0.53380697298101955, 0.64974697298101955, 0.66382697298101955,
    0.17010523393364733, 0.16817823393364733, 0.037807233933647332,
    0.10418323393364733, 0.13310423393364733, 0.11101354723412521,
    0.20512454723412521, 0.18170254723412521, 0.13509954723412521,
    0.20283554723412521
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 3, -2, -3),
    upper_bound = c(0, 3, -1, -2)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(0, 3, -1.5, -2.5), max_iter,
    lower_bound = c(0, 3, -2, -3),
    upper_bound = c(0, 3, -1, -2)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(1, 1, -3, -1), max_iter,
    lower_bound = c(0, 3, -2, -3),
    upper_bound = c(0, 3, -1, -2)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("four-parameter logistic fit - difficult")

test_that("logistic4_fit_constrained", {
  max_iter <- 10000

  x <- rep(-log(c(10000, 1000, 100, 10, 1, 0.1, 0.01)), each = 3)
  y <- c(
    0.877362, 0.812841, 0.883113,
    0.873494, 0.845769, 0.999422,
    0.888961, 0.735539, 0.842040,
    0.518041, 0.519261, 0.501252,
    0.253209, 0.083937, -0.000719,
    0.049249, 0.070804, 0.091425,
    0.041096, -0.036148, 0.092564
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 1,
    "growth_rate" = -0.64478026194068988,
    "logx_midpoint" = -2.4237429873134359
  )
  sigma <- 0.087699027767824044
  loglik <- 22.363900894124325

  fitted_values <- c(
    rep(0.98757895258110592, 3), rep(0.94740983368278412, 3),
    rep(0.80321982692343109, 3), rep(0.48047987151771678, 3),
    rep(0.17324786231650001, 3), rep(0.045327992438758544, 3),
    rep(0.010643509240940189, 3)
  )

  residuals <- c(
    -0.11021695258110593, -0.17473795258110593, -0.10446595258110593,
    -0.073915833682784112, -0.10164083368278411, 0.052012166317215888,
    0.085741173076568845, -0.067680826923431155, 0.038820173076568845,
    0.037561128482283166, 0.038781128482283166, 0.020772128482283166,
    0.079961137683499994, -0.089310862316500006, -0.17396686231650001,
    0.0039210075612414536, 0.025476007561241454, 0.046097007561241454,
    0.030452490759059810, -0.046791509240940190, 0.081920490759059810
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 1, -Inf, -Inf),
    upper_bound = c(0, 1, 0, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 19)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # try starting from another point
  result <- logistic4_fit_constrained(
    x, y, c(0, 1, -0.1, 0), max_iter,
    lower_bound = c(0, 1, -Inf, -Inf),
    upper_bound = c(0, 1, 0, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 19)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})
