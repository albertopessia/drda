context("four-parameter logistic weighted fit")

test_that("logistic4_weighted_fisher_info", {
  x <- rep(-log(c(10000, 1000, 100, 10, 1, 0.1, 0.01)), each = 3)
  w <- c(
    0.99087, 1.09524, 0.97454, 0.97332, 1.107, 1.01284, 1.05281, 1.01943,
    1.03254, 0.91983, 0.97139, 0.95902, 1.03779, 1.00684, 0.96938, 0.93563,
    1.0166, 1.01108, 0.98231, 1.06603, 0.95987
  )

  param <- c(0.05, 0.9, -1.5, -2)
  sigma <- 0.07

  true_value <- matrix(
    c(
      1862.8734926102293, 179.59247102183111, 34.297737907655957,
      103.73990413324456, 0, 179.59247102183111, 2082.9129939175371,
      -47.785087466830068, 125.24049641959011, 0, 34.297737907655957,
      -47.785087466830068, 6.9815011762936831, -8.7042204037096264, 0,
      103.73990413324456, 125.24049641959011, -8.7042204037096264,
      55.789275832253070, 0, 0, 0, 0, 0, 8629.2
    ),
    nrow = 5,
    ncol = 5
  )

  fim <- logistic4_weighted_fisher_info(x, w, param, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

test_that("logistic4_weighted_fit_unconstrained", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )
  w <- c(
    1.218965, 2.157773, 0.053758, 1.349539, 0.123137,
    1.369847, 2.098511, 0.061778, 1.268822, 0.507726,
    0.437259, 1.597353, 1.040318, 0.595828, 0.370306,
    0.013929, 1.060632, 0.002076, 0.010446, 0.131549,
    0.081876, 0.208199, 0.364845, 0.918672, 0.751648
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0.20395213202042778,
    "maximum" = 2.4890781444733631,
    "growth_rate" = -1.7951682327440798,
    "logx_midpoint" = -1.0168009501939123
  )
  sigma <- 0.044489144688820571
  loglik <- 29.427697512015472

  fitted_values <- c(
    rep(2.4890197832248375, 5), rep(2.4854422410831057, 5),
    rep(2.2823971606745102, 5), rep(0.52111777627602181, 5),
    rep(0.20983909674711762, 5)
  )

  residuals <- c(
    0.026220216775162453, -0.026309783224837547, -0.019179783224837547,
    0.011120216775162453, -0.0086897832248375469, -0.077582241083105730,
    0.076267758916894270, 0.013347758916894270, -0.037532241083105730,
    0.010537758916894270, 0.019532839325489807, -0.049677160674510193,
    -0.023367160674510193, 0.092572839325489807, 0.10665283932548981,
    0.0065962237239781905, 0.0046692237239781905, -0.12570177627602181,
    -0.059325776276021809, -0.030404776276021809, -0.058767096747117623,
    0.035343903252882377, 0.011921903252882377, -0.034681096747117623,
    0.033054903252882377
  )

  result <- logistic4_weighted_fit_unconstrained(x, y, w, NULL, max_iter)

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  result <- logistic4_weighted_fit_unconstrained(
    x, y, w, c(0, 1, -1, 0), max_iter
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_weighted_fit_constrained: inequalities I", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )
  w <- c(
    1.218965, 2.157773, 0.053758, 1.349539, 0.123137,
    1.369847, 2.098511, 0.061778, 1.268822, 0.507726,
    0.437259, 1.597353, 1.040318, 0.595828, 0.370306,
    0.013929, 1.060632, 0.002076, 0.010446, 0.131549,
    0.081876, 0.208199, 0.364845, 0.918672, 0.751648
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0.20395213202042778,
    "maximum" = 2.4890781444733631,
    "growth_rate" = -1.7951682327440798,
    "logx_midpoint" = -1.0168009501939123
  )
  sigma <- 0.044489144688820571
  loglik <- 29.427697512015472

  fitted_values <- c(
    rep(2.4890197832248375, 5), rep(2.4854422410831057, 5),
    rep(2.2823971606745102, 5), rep(0.52111777627602181, 5),
    rep(0.20983909674711762, 5)
  )

  residuals <- c(
    0.026220216775162453, -0.026309783224837547, -0.019179783224837547,
    0.011120216775162453, -0.0086897832248375469, -0.077582241083105730,
    0.076267758916894270, 0.013347758916894270, -0.037532241083105730,
    0.010537758916894270, 0.019532839325489807, -0.049677160674510193,
    -0.023367160674510193, 0.092572839325489807, 0.10665283932548981,
    0.0065962237239781905, 0.0046692237239781905, -0.12570177627602181,
    -0.059325776276021809, -0.030404776276021809, -0.058767096747117623,
    0.035343903252882377, 0.011921903252882377, -0.034681096747117623,
    0.033054903252882377
  )

  result <- logistic4_weighted_fit_constrained(
    x, y, w, NULL, 1000,
    lower_bound = c(0, 2, -2, -2),
    upper_bound = c(1, 3, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(0.5, 2.5, -1.5, -1), 1000,
    lower_bound = c(0, 2, -2, -2),
    upper_bound = c(1, 3, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(2, 1, -0.5, 1), max_iter,
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

test_that("logistic4_weighted_fit_constrained: inequalities II", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )
  w <- c(
    1.218965, 2.157773, 0.053758, 1.349539, 0.123137,
    1.369847, 2.098511, 0.061778, 1.268822, 0.507726,
    0.437259, 1.597353, 1.040318, 0.595828, 0.370306,
    0.013929, 1.060632, 0.002076, 0.010446, 0.131549,
    0.081876, 0.208199, 0.364845, 0.918672, 0.751648
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 2,
    "growth_rate" = -2,
    "logx_midpoint" = -0.38369690927574105
  )
  sigma <- 0.37806983476770585
  loglik <- -24.068145157520762

  fitted_values <- c(
    rep(1.9999956917196018, 5), rep(1.9995692638189248, 5),
    rep(1.9578256008540269, 5), rep(0.63408629020251394, 5),
    rep(0.0092415256296815808, 5)
  )

  residuals <- c(
    0.51524430828039821, 0.46271430828039821, 0.46984430828039821,
    0.50014430828039821, 0.48033430828039821, 0.40829073618107523,
    0.56214073618107523, 0.49922073618107523, 0.44834073618107523,
    0.49641073618107523, 0.34410439914597311, 0.27489439914597311,
    0.30120439914597311, 0.41714439914597311, 0.43122439914597311,
    -0.10637229020251394, -0.10829929020251394, -0.23867029020251394,
    -0.17229429020251394, -0.14337329020251394, 0.14183047437031842,
    0.23594147437031842, 0.21251947437031842, 0.16591647437031842,
    0.23365247437031842
  )

  result <- logistic4_weighted_fit_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(-1, 1, -2, -1),
    upper_bound = c(0, 2, -1, 0)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 21)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(-0.5, 1.5, -1.5, -0.5), max_iter,
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
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(-2, 3, -3, 1), max_iter,
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

test_that("logistic4_weighted_fit_constrained: equalities", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )
  w <- c(
    1.218965, 2.157773, 0.053758, 1.349539, 0.123137,
    1.369847, 2.098511, 0.061778, 1.268822, 0.507726,
    0.437259, 1.597353, 1.040318, 0.595828, 0.370306,
    0.013929, 1.060632, 0.002076, 0.010446, 0.131549,
    0.081876, 0.208199, 0.364845, 0.918672, 0.751648
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 3,
    "growth_rate" = -0.55090343806284802,
    "logx_midpoint" = -1.2184442720398583
  )
  sigma <- 0.27137679970069061
  loglik <- -16.916025718670329

  fitted_values <- c(
    rep(2.8748521867525742, 5), rep(2.5979003685164140, 5),
    rep(1.9350854370143040, 5), rep(1.0146542686053593, 5),
    rep(0.37702690131062104, 5)
  )

  residuals <- c(
    -0.35961218675257422, -0.41214218675257422, -0.40501218675257422,
    -0.37471218675257422, -0.39452218675257422, -0.19004036851641402,
    -0.036190368516414022, -0.099110368516414022, -0.14999036851641402,
    -0.10192036851641402, 0.36684456298569605, 0.29763456298569605,
    0.32394456298569605, 0.43988456298569605, 0.45396456298569605,
    -0.48694026860535935, -0.48886726860535935, -0.61923826860535935,
    -0.55286226860535935, -0.52394126860535935, -0.22595490131062104,
    -0.13184390131062104, -0.15526590131062104, -0.20186890131062104,
    -0.13413290131062104
  )

  result <- logistic4_weighted_fit_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 3, -Inf, -Inf),
    upper_bound = c(0, 3, Inf, Inf)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(0, 3, 1, 1), max_iter,
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
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(-1, 2, 1, 1), max_iter,
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

test_that("logistic4_weighted_fit_constrained: equalities and inequalities", {
  max_iter <- 10000

  x <- rep(-log(c(1000, 100, 10, 1, 0.1)), each = 5)
  y <- c(
    2.515240, 2.462710, 2.469840, 2.500140, 2.480330,
    2.407860, 2.561710, 2.498790, 2.447910, 2.495980,
    2.301930, 2.232720, 2.259030, 2.374970, 2.389050,
    0.527714, 0.525787, 0.395416, 0.461792, 0.490713,
    0.151072, 0.245183, 0.221761, 0.175158, 0.242894
  )
  w <- c(
    1.218965, 2.157773, 0.053758, 1.349539, 0.123137,
    1.369847, 2.098511, 0.061778, 1.268822, 0.507726,
    0.437259, 1.597353, 1.040318, 0.595828, 0.370306,
    0.013929, 1.060632, 0.002076, 0.010446, 0.131549,
    0.081876, 0.208199, 0.364845, 0.918672, 0.751648
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 3,
    "growth_rate" = -1,
    "logx_midpoint" = -2
  )
  sigma <- 0.36566958617509652
  loglik <- -24.371573110772934

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

  result <- logistic4_weighted_fit_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 3, -2, -3),
    upper_bound = c(0, 3, -1, -2)
  )

  expect_true(result$converged)
  expect_equal(result$coefficients, param_est)
  expect_equal(result$sigma, sigma)
  expect_equal(result$loglik, loglik)
  expect_equal(result$df.residual, 23)
  expect_identical(result$weights, w)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(0, 3, -1.5, -2.5), max_iter,
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
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(1, 1, -3, -1), max_iter,
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

context("four-parameter logistic weighted fit - difficult")

test_that("logistic4_weighted_fit_constrained", {
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
  w <- c(
    0.990868, 1.095238, 0.974544,
    0.973318, 1.107001, 1.012844,
    1.052806, 1.019427, 1.032544,
    0.919827, 0.971385, 0.959019,
    1.037789, 1.006835, 0.969383,
    0.935633, 1.016597, 1.011085,
    0.982307, 1.066032, 0.959870
  )

  # values computed with Mathematica
  param_est <- c(
    "minimum" = 0,
    "maximum" = 1,
    "growth_rate" = -0.63998970490154362,
    "logx_midpoint" = -2.4252577287646781
  )
  sigma <- 0.088627511775198154
  loglik <- 22.178018857136147

  fitted_values <- c(
    rep(0.98716147280684045, 3), rep(0.94627995667364474, 3),
    rep(0.80140860770236510, 3), rep(0.48038276919850788, 3),
    rep(0.17477740371273640, 3), rep(0.046275122654902929, 3),
    rep(0.010993448119827933, 3)
  )

  residuals <- c(
    -0.10979947280684047, -0.17432047280684047, -0.10404847280684047,
    -0.072785956673644729, -0.10051095667364473, 0.053142043326355271,
    0.087552392297634831, -0.065869607702365169, 0.040631392297634831,
    0.037658230801492069, 0.038878230801492069, 0.020869230801492069,
    0.078431596287263600, -0.090840403712736400, -0.17549640371273640,
    0.0029738773450970691, 0.024528877345097069, 0.045149877345097069,
    0.030102551880172066, -0.047141448119827934, 0.081570551880172066
  )

  # automatic initial values are unfortunately bad and algorithm converges to
  # a local optimum or not at all
  result <- logistic4_weighted_fit_constrained(
    x, y, w, NULL, max_iter,
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
  result <- logistic4_weighted_fit_constrained(
    x, y, w, c(0, 1, -0.1, 0), max_iter,
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
