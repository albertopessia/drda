context("6-parameter logistic - core functions")

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2, 3 / 2)

  true_value <- c(
    0.42221710418125004, 0.42171092652477733, 0.37575797785755273,
    0.046454735980600713, 0.040000850149265818, 0.040000000085266528
  )

  value <- logistic6_function(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2, 1 / 2, 3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # log_omega
      0.38221710418125004, 0.38171092652477733, 0.33575797785755273,
      0.0064547359806007129, 8.5014926581760944e-07, 8.5266527756032689e-11,
      # eta
      -0.000027676835132687997, -0.0015860669501348013, -0.033819316025770067,
      0.016847800200509503, 6.4558872590255383e-06, 1.0411177759771462e-09,
      # phi
      0.000010235979146562789, 0.0010215652316203735, 0.084275963560716097,
      0.022463733600679338, 3.3955254655155444e-06, 3.4106101689568101e-10,
      # log_nu
      0.30995139895241304, 0.30954109512977705, 0.27365642866850772,
      0.020345465001078442, 0.000010057277525685746, 1.7935341844903725e-09,
      # log_xi
      -0.76442909037292680, -0.76291107043374447, -0.62937797393474741,
      -0.0016776051608617570, -2.5357988774466817e-09, -2.5470642248697705e-15
    ),
    nrow = 6,
    ncol = 6
  )

  true_hessian <- array(
    c(
      # (alpha, ...)
      rep(0, 36),
      # (log_omega, alpha)
      rep(0, 6),
      # (log_omega, log_omega)
      0.38221710418125004, 0.38171092652477733, 0.33575797785755273,
      0.0064547359806007129, 8.5014926581760944e-07, 8.5266527756032689e-11,
      # (log_omega, eta)
      -0.000027676835132687997, -0.0015860669501348013, -0.033819316025770067,
      0.016847800200509503, 6.4558872590255383e-06, 1.0411177759771462e-09,
      # (log_omega, phi)
      0.000010235979146562789, 0.0010215652316203735, 0.084275963560716097,
      0.022463733600679338, 3.3955254655155444e-06, 3.4106101689568101e-10,
      # (log_omega, log_nu)
      0.30995139895241304, 0.30954109512977705, 0.27365642866850772,
      0.020345465001078442, 0.000010057277525685746, 1.7935341844903725e-09,
      # (log_omega, log_xi)
      -0.76442909037292680, -0.76291107043374447, -0.62937797393474741,
      -0.0016776051608617570, -2.5357988774466817e-09, -2.5470642248697705e-15,
      # (eta, alpha)
      rep(0, 6),
      # (eta, log_omega)
      -0.000027676835132687997, -0.0015860669501348013, -0.033819316025770067,
      0.016847800200509503, 6.4558872590255383e-06, 1.0411177759771462e-09,
      # (eta, eta)
      -0.00014966654512113986, -0.0049151222824604186, -0.022033188829639366,
      0.040691115014078006, 0.000048988285040711751, 1.2712117605288108e-08,
      # (eta, phi)
      0.000050234568891350354, 0.0026549841461145773, 0.012767570346070067,
      0.043022953218431006, 0.000024068017042147825, 3.9938472952908821e-09,
      # (eta, log_nu)
      -0.000022443796599450935, -0.0012851291354092528, -0.025441940164379255,
      0.038446220915158790, 0.000069926976157246812, 2.0858233544335177e-08,
      # (eta, log_xi)
      0.000083029949497691868, 0.0047550172812319847, 0.095091408114754192,
      -0.0065681903916336688, -2.8884630598278052e-08, -4.6650084931125022e-14,
      # (phi, alpha)
      rep(0, 6),
      # (phi, log_omega)
      0.000010235979146562789, 0.0010215652316203735, 0.084275963560716097,
      0.022463733600679338, 3.3955254655155444e-06, 3.4106101689568101e-10,
      # (phi, eta)
      0.000050234568891350354, 0.0026549841461145773, 0.012767570346070067,
      0.043022953218431006, 0.000024068017042147825, 3.9938472952908821e-09,
      # (phi, phi)
      -0.000020471547105604366, -0.0020390294716921548, -0.13682175910245932,
      0.072339760025027567, 0.000013551717657746544, 1.3642135032685379e-09,
      # (phi, log_nu)
      8.3005962517132876e-06, 0.00082773507307800837, 0.063399981849819896,
      0.051261627886878386, 0.000036778651599976862, 6.8329736629477352e-09,
      # (phi, log_xi)
      -0.000030707731845927761, -0.0030626451990868243, -0.23696280667266186,
      -0.0087575871888448917, -1.5192102157816607e-08, -1.5282157093078623e-14,
      # (log_nu, alpha)
      rep(0, 6),
      # (log_nu, log_omega)
      0.30995139895241304, 0.30954109512977705, 0.27365642866850772,
      0.020345465001078442, 0.000010057277525685746, 1.7935341844903725e-09,
      # (log_nu, eta)
      -0.000022443796599450935, -0.0012851291354092528, -0.025441940164379255,
      0.038446220915158790, 0.000069926976157246812, 2.0858233544335177e-08,
      # (log_nu, phi)
      8.3005962517132876e-06, 0.00082773507307800837, 0.063399981849819896,
      0.051261627886878386, 0.000036778651599976862, 6.8329736629477352e-09,
      # (log_nu, log_nu)
      -0.058602443935219952, -0.058524387549721785, -0.047971111555920137,
      0.053556142044358048, 0.00011061568486094975, 3.6102991956890671e-08,
      # (log_nu, log_xi)
      0.14453556075579985, 0.14475353032653957, 0.15590308930299996,
      -0.0021506440975162657, -2.4930719747841768e-08, -4.8481994812079120e-14,
      # (log_xi, alpha)
      rep(0, 6),
      # (log_xi, log_omega)
      -0.76442909037292680, -0.76291107043374447, -0.62937797393474741,
      -0.0016776051608617570, -2.5357988774466817e-09, -2.5470642248697705e-15,
      # (log_xi, eta)
      0.000083029949497691868, 0.0047550172812319847, 0.095091408114754192,
      -0.0065681903916336688, -2.8884630598278052e-08, -4.6650084931125022e-14,
      # (log_xi, phi)
      -0.000030707731845927761, -0.0030626451990868243, -0.23696280667266186,
      -0.0087575871888448917, -1.5192102157816607e-08, -1.5282157093078623e-14,
      # (log_xi, log_nu)
      0.14453556075579985, 0.14475353032653957, 0.15590308930299996,
      -0.0021506440975162657, -2.4930719747841768e-08, -4.8481994812079120e-14,
      # (log_xi, log_xi)
      1.5288428268799306, 1.5242908182679455, 1.1402745445331639,
      -0.0010235832726989319, -2.5244533240149399e-09, -2.5469500967997703e-15
    ),
    dim = c(6, 6, 6)
  )

  gradient_hessian <- logistic6_gradient_hessian(x, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 6)
  expect_length(gradient_hessian$H, 6 * 6 * 6)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

context("6-parameter logistic - RSS functions")

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

  true_value <- 2.0908902035928604

  rss <- logistic6_rss(stats)

  expect_type(rss, "closure")

  value <- rss(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2), NA)
  rss <- logistic6_rss_fixed(stats, known_param)

  expect_type(rss, "closure")

  value <- rss(c(log(43 / 50), -2, log(3 / 2)))

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
    -4.5432784577966121, -1.4803879620386775, 0.030726379870560891,
    -0.098256479294546431, -1.2078753918106346, 2.9116476844300818
  )

  true_hessian <- matrix(
    c(
      # alpha
      15.000000000000000, 2.9891215422033879, -0.0050692949435275129,
      0.26151245173284690, 2.4872023714204765, -5.8474868585403523,
      # log_omega
      2.9891215422033879, -0.37937502970889297, 0.0066031690384121262,
      -0.039902270261351402, -0.31371360203813016, 0.73879892428711024,
      # eta
      -0.0050692949435275129, 0.0066031690384121262, 0.016681432220990243,
      -0.037754568939522141, -0.0055461232427009715, -0.056647742370221930,
      # phi
      0.26151245173284690, -0.039902270261351402, -0.037754568939522141,
      0.13344876935609788, -0.038587395248642860, 0.14652892520075174,
      # log_nu
      2.4872023714204765, -0.31371360203813016, -0.0055461232427009715,
      -0.038587395248642860, 0.92894109931730305, -2.3526424176254161,
      # log_xi
      -5.8474868585403523, 0.73879892428711024, -0.056647742370221930,
      0.14652892520075174, -2.3526424176254161, -1.4043333859738446
    ),
    nrow = 6,
    ncol = 6
  )

  gh <- logistic6_rss_gradient_hessian(stats)

  expect_type(gh, "closure")

  gradient_hessian <- gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6)
  expect_length(gradient_hessian$H, 6 * 6)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, -3 / 2, -log(2), NA)
  gh <- logistic6_rss_gradient_hessian_fixed(stats, known_param)

  expect_type(gh, "closure")

  gradient_hessian <- gh(c(log(43 / 50), -2, log(3 / 2)))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 3)
  expect_length(gradient_hessian$H, 3 * 3)

  expect_equal(gradient_hessian$G, true_gradient[c(2, 3, 6)])
  expect_equal(gradient_hessian$H, true_hessian[c(2, 3, 6), c(2, 3, 6)])
})

context("6-parameter logistic - general functions")

test_that("logistic6_init", {
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
    -0.017959172658911955, -2.7745983376934573, -1.00499999999999989,
    0.57564627324851103, -1, -1
  )

  stats <- cbind(x, n, m, v)

  start <- logistic6_init(stats)

  expect_type(start, "double")
  expect_length(start, 6)
  expect_equal(start, true_value)
})

test_that("logistic6_fisher_info_normal", {
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
    nu = 1 / 2,
    xi = 3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      6085.1433179765774, 699.07299953734511, 16.237238579181247,
      86.124012799896435, 1036.4244469696333, -772.89570186952358, 0,
      # beta
      699.07299953734511, 516.71068294873238, -10.887523968805600,
      27.012697248651210, 722.02276982000374, -583.49068403166273, 0,
      # eta
      16.237238579181247, -10.887523968805600, 1.4844212990353402,
      -1.5225998656546614, -14.242119593664837, 11.976827529012201, 0,
      # phi
      86.124012799896435, 27.012697248651210, -1.5225998656546614,
      6.6843744037638038, 39.220817474330747, -28.746434710532119, 0,
      # nu
      1036.4244469696333, 722.02276982000374, -14.242119593664837,
      39.220817474330747, 1010.7439620495048, -814.84583690216071, 0,
      # xi
      -772.89570186952358, -583.49068403166273, 11.976827529012201,
      -28.746434710532119, -814.84583690216071, 659.48717278612922, 0,
      # sigma
      rep(0, 6), 6800
    ),
    nrow = 7,
    ncol = 7
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "nu", "xi", "sigma"
  )

  fim <- logistic6_fisher_info_normal(stats, n, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 7 * 7)
  expect_equal(fim, true_value)
})

context("6-parameter logistic - fit")

test_that("logistic6_fit: fit", {
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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss <- 0.024883087882351184

  fitted_values <- c(
    rep(0.9202839186685336, 3), rep(0.9197191693822363, 2),
    rep(0.8900274208207811, 2), rep(0.5533792934125557, 5),
    rep(0.2627180312034356, 3), rep(0.1541298223060635, 4),
    0.11508521369102623
  )

  residuals <- c(
    0.0077160813314664, -0.0322839186685336, 0.0597160813314664,
    0.0282808306177637, -0.0637191693822363, 0.0069725791792189,
    -0.0070274208207811, -0.0653792934125557, -0.0213792934125557,
    0.0326207065874443, 0.0126207065874443, 0.0456207065874443,
    -0.0037180312034356, 0.0022819687965644, -0.0197180312034356,
    -0.03712982230606349, -0.01112982230606349, 0.02387017769393651,
    0.06487017769393651, -0.02308521369102623
  )

  result <- logistic6_fit(x, y, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  result <- logistic6_fit(x, y, c(0, 0, -1, 0, 0, 0), max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic6_fit_constrained: inequalities", {
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

  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  theta <- c(
    alpha = 0.11015274566316652,
    beta = 0.110152745663166526 + exp(-0.015957337613761089),
    eta = -1.4048937624493838,
    phi = -0.66709123288456063,
    nu = exp(1),
    xi = exp(0.51669269729838785)
  )

  rss <- 0.025326988257878513

  fitted_values <- c(
    rep(0.92387938085890750, 3), rep(0.92204263235590002, 2),
    rep(0.88001856248146377, 2), rep(0.55580425395287030, 5),
    rep(0.25641505220964583, 3), rep(0.15480915223246012, 4),
    0.12373639849382178
  )

  residuals <- c(
    0.00412061914109250, -0.03587938085890750, 0.05612061914109250,
    0.02595736764409998, -0.06604263235590002, 0.01698143751853623,
    0.00298143751853623, -0.06780425395287030, -0.02380425395287030,
    0.03019574604712970, 0.01019574604712970, 0.04319574604712970,
    0.00258494779035417, 0.00858494779035417, -0.01341505220964583,
    -0.03780915223246012, -0.01180915223246012, 0.02319084776753988,
    0.06419084776753988, -0.03173639849382178
  )

  result <- logistic6_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-5)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic6_fit_constrained(
    x, y, c(0.5, 0, -1, -0.5, 0.5, 0.5), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-5)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic6_fit_constrained(
    x, y, c(2, -1, -3, 1, 2, -1), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-5)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic6_fit_constrained: equalities", {
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -2.3448538310034034,
    phi = -1.1869812279392416,
    nu = exp(2.1013467432775474),
    xi = exp(0.68268209514681377)
  )

  rss <- 0.027858205619881242

  fitted_values <- c(
    rep(0.91990299490116622, 3), rep(0.91975017005499427, 2),
    rep(0.89071256487960993, 2), rep(0.5492685788947333, 5),
    rep(0.2842905938017210, 3), rep(0.14692068718566469, 4),
    0.07590586225363996
  )

  residuals <- c(
    0.00809700509883378, -0.03190299490116622, 0.06009700509883378,
    0.02824982994500573, -0.06375017005499427, 0.00628743512039007,
    -0.00771256487960993, -0.0612685788947333, -0.0172685788947333,
    0.0367314211052667, 0.0167314211052667, 0.0497314211052667,
    -0.0252905938017210, -0.0192905938017210, -0.0412905938017210,
    -0.02992068718566469, -0.00392068718566469, 0.03107931281433531,
    0.07207931281433531, 0.01609413774636004
  )

  result <- logistic6_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  result <- logistic6_fit_constrained(
    x, y, c(0, 0, -1, 0, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with different equalities
  result <- logistic6_fit_constrained(
    x, y, c(1, 3, -1, 0, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic6_fit_constrained: equalities and inequalities", {
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -2.3448538310034034,
    phi = -1.1869812279392416,
    nu = exp(2.1013467432775474),
    xi = exp(0.68268209514681377)
  )

  rss <- 0.027858205619881242

  fitted_values <- c(
    rep(0.91990299490116622, 3), rep(0.91975017005499427, 2),
    rep(0.89071256487960993, 2), rep(0.5492685788947333, 5),
    rep(0.2842905938017210, 3), rep(0.14692068718566469, 4),
    0.07590586225363996
  )

  residuals <- c(
    0.00809700509883378, -0.03190299490116622, 0.06009700509883378,
    0.02824982994500573, -0.06375017005499427, 0.00628743512039007,
    -0.00771256487960993, -0.0612685788947333, -0.0172685788947333,
    0.0367314211052667, 0.0167314211052667, 0.0497314211052667,
    -0.0252905938017210, -0.0192905938017210, -0.0412905938017210,
    -0.02992068718566469, -0.00392068718566469, 0.03107931281433531,
    0.07207931281433531, 0.01609413774636004
  )

  result <- logistic6_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic6_fit_constrained(
    x, y, c(0, 0, -2, -1, 1, 0.5), max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic6_fit_constrained(
    x, y, c(1, 2, -5, 2, -1, 2), max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("6-parameter logistic - weighted fit")

test_that("logistic6_fit_weighted: fit", {
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

  # logistic6 model is basically unidentifiable: many parameters are associated
  # with the same residual sum of squares
  # there is no point in testing the values of `result$coefficients`
  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss <- 0.014141550871844299

  fitted_values <- c(
    rep(0.9375852722593341, 3), rep(0.9351351725701265, 1),
    rep(0.8859559720192795, 2), rep(0.5516453155354000, 4),
    rep(0.2647953410627469, 3), rep(0.17495286982311317, 3),
    0.14985594300774278
  )

  residuals <- c(
    -0.0095852722593341, -0.0495852722593341, 0.0424147277406659,
    0.0128648274298735, 0.0110440279807205, -0.0029559720192795,
    -0.0636453155354000, -0.0196453155354000, 0.0143546844646000,
    0.0473546844646000, -0.0057953410627469, 0.0002046589372531,
    -0.0217953410627469, -0.03195286982311317, 0.00304713017688683,
    0.04404713017688683, -0.05785594300774278
  )

  result <- logistic6_fit_weighted(x, y, w, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  result <- logistic6_fit_weighted(x, y, w, c(0, 0, -1, 0, 0, 0), max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 6)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic6_fit_weighted_constrained: inequalities", {
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

  estimated <- c(
    alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  rss <- 0.0141415508718443

  fitted_values <- c(
    rep(0.93758527235398976, 3), rep(0.93513517255536545, 1),
    rep(0.88595597148284482, 2), rep(0.55164531568757359, 4),
    rep(0.26479534053012690, 3), rep(0.17495287000603935, 3),
    0.14985594381005773
  )

  residuals <- c(
    -0.00958527235398976, -0.04958527235398976, 0.04241472764601024,
    0.01286482744463455, 0.01104402851715518, -0.00295597148284482,
    -0.06364531568757359, -0.01964531568757359, 0.01435468431242641,
    0.04735468431242641, -0.00579534053012690, 0.00020465946987310,
    -0.02179534053012690, -0.03195287000603935, 0.00304712999396065,
    0.04404712999396065, -0.05785594381005773
  )

  result <- logistic6_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 6)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values within the boundaries
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(0.5, 0, -1, -0.5, 0.5, 0.5), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 6)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values outside the boundaries
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(2, -1, -3, 1, 2, -1), max_iter,
    lower_bound = c(0, -0.1, -2, -1, 0, 0),
    upper_bound = c(1,  0.1,  0,  0, 1, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 6)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic6_fit_weighted_constrained: equalities", {
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -2.0281792165635454,
    phi = -1.2409276005266773,
    nu = exp(2.0224177185354462),
    xi = exp(0.48833621411130149)
  )

  rss <- 0.016942492408652008

  fitted_values <- c(
    rep(0.93741397025419554, 3), rep(0.93679545474110372, 1),
    rep(0.88551343724525149, 2), rep(0.5471760024605463, 4),
    rep(0.2955709497045542, 3), rep(0.15934552408476684, 3),
    0.08588003931261228
  )

  residuals <- c(
    -0.00941397025419554, -0.04941397025419554, 0.04258602974580446,
    0.01120454525889628, 0.01148656275474851, -0.00251343724525149,
    -0.0591760024605463, -0.0151760024605463, 0.0188239975394537,
    0.0518239975394537, -0.0365709497045542, -0.0305709497045542,
    -0.0525709497045542, -0.01634552408476684, 0.01865447591523316,
    0.05965447591523316, 0.00611996068738772
  )

  result <- logistic6_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with same equalities
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(0, 0, -1, 0, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with different equalities
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(1, 3, -1, 0, 0, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic6_fit_weighted_constrained: equalities and inequalities", {
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

  estimated <- c(
    alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE, nu = TRUE, xi = TRUE
  )

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -2.0281792165635454,
    phi = -1.2409276005266773,
    nu = exp(2.0224177185354462),
    xi = exp(0.48833621411130149)
  )

  rss <- 0.016942492408652008

  fitted_values <- c(
    rep(0.93741397025419554, 3), rep(0.93679545474110372, 1),
    rep(0.88551343724525149, 2), rep(0.5471760024605463, 4),
    rep(0.2955709497045542, 3), rep(0.15934552408476684, 3),
    0.08588003931261228
  )

  residuals <- c(
    -0.00941397025419554, -0.04941397025419554, 0.04258602974580446,
    0.01120454525889628, 0.01148656275474851, -0.00251343724525149,
    -0.0591760024605463, -0.0151760024605463, 0.0188239975394537,
    0.0518239975394537, -0.0365709497045542, -0.0305709497045542,
    -0.0525709497045542, -0.01634552408476684, 0.01865447591523316,
    0.05965447591523316, 0.00611996068738772
  )

  result <- logistic6_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(0, 0, -2, -1, 1, 0.5), max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic6_fit_weighted_constrained(
    x, y, w, c(1, 2, -5, 2, -1, 2), max_iter,
    lower_bound = c(0, 0, -3, -2, 0, 0),
    upper_bound = c(0, 0, -1,  0, 3, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-6)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("6-parameter logistic - drda fit")

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
      y ~ x, mean_function = "logistic6",
      lower_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = matrix(-Inf, nrow = 6, ncol = 2),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 6)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c( 0, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c(Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' and 'upper_bound' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c(0, -1, -Inf, -Inf, 0, 0),
      upper_bound = rep(Inf, 6)
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
      y ~ x, mean_function = "logistic6",
      upper_bound = c("a", "b", "c", "d", "e", "f")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 6),
      upper_bound = matrix(Inf, nrow = 6, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 7),
      upper_bound = rep(Inf, 7)
    ),
    "'lower_bound' and 'upper_bound' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      lower_bound = rep(-Inf, 6),
      upper_bound = c(1, 0, Inf, Inf, 0, 0)
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
      y ~ x, mean_function = "logistic6",
      start = c("a", "b", "c", "d", "e", "f")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, Inf, -1, 0, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(-Inf, 0, -1, 0, 1, 1)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 0, -1, 0, 1, 1, 1)
    ),
    "'start' must be of length 6"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, -1, -1, 0, 1, 1)
    ),
    "parameter 'beta' is smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, 0, 0, 1, 1)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 0, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, -1, 1)
    ),
    "parameter 'nu' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 1, 0)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic6",
      start = c(0, 1, -1, 0, 1, -1)
    ),
    "parameter 'xi' cannot be negative nor zero"
  )
})
