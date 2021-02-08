context("4-parameter logistic - core functions")

test_that("Function value", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_value <- c(
    0.89998272678518785, 0.89827610635754658, 0.75615618502842375,
    0.080786250932707432, 0.040427955721356054, 0.040004281666562474
  )

  value <- logistic4_function(x, theta)

  expect_type(value, "double")
  expect_length(value, 6)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian", {
  x <- -log(c(1000, 100, 10, 1, 0.1, 0.01))
  theta <- c(4 / 100, 9 / 10, -2, -3 / 2)

  true_gradient <- matrix(
    c(
      # alpha
      rep(1, 6),
      # log_omega
      0.85998272678518785, 0.85827610635754658, 0.71615618502842375,
      0.040786250932707432, 0.00042795572135605419, 4.2816665624737117e-06,
      # eta
      -0.000093407442446748936, -0.0053422529404918635, -0.096137223993448991,
      0.058277891052876651, 0.0016265282457603059, 0.000026140172899299094,
      # phi
      0.000034545735754643967, 0.0034408761005103580, 0.23956892504646165,
      0.077703854737168868, 0.00085548552155061679, 8.5632904908345112e-06
    ),
    nrow = 6,
    ncol = 4
  )

  true_hessian <- array(
    c(
      # (alpha, ...)
      rep(0, 24),
      # (log_omega, alpha)
      rep(0, 6),
      # (log_omega, log_omega)
      0.85998272678518785, 0.85827610635754658, 0.71615618502842375,
      0.040786250932707432, 0.00042795572135605419, 4.2816665624737117e-06,
      # (log_omega, eta)
      -0.000093407442446748936, -0.0053422529404918635, -0.096137223993448991,
      0.058277891052876651, 0.0016265282457603059, 0.000026140172899299094,
      # (log_omega, phi)
      0.000034545735754643967, 0.0034408761005103580, 0.23956892504646165,
      0.077703854737168868, 0.00085548552155061679, 8.5632904908345112e-06,
      # (eta, alpha)
      rep(0, 6),
      # (eta, log_omega)
      -0.000093407442446748936, -0.0053422529404918635, -0.096137223993448991,
      0.058277891052876651, 0.0016265282457603059, 0.000026140172899299094,
      # (eta, eta)
      -0.00050510429899797611, -0.016522099929072128, -0.051347268812537329,
      0.079125196968925631, 0.0061788564529920891, 0.00015958861514360448,
      # (eta, phi)
      0.00016953451261236921, 0.0089212330591913757, 0.0081702409985757300,
      0.066648335256649741, 0.0028220761397279995, 0.000047998179978742605,
      # (phi, alpha)
      rep(0, 6),
      # (phi, log_omega)
      0.000034545735754643967, 0.0034408761005103580, 0.23956892504646165,
      0.077703854737168868, 0.00085548552155061679, 8.5632904908345112e-06,
      # (phi, eta)
      0.00016953451261236921, 0.0089212330591913757, 0.0081702409985757300,
      0.066648335256649741, 0.0028220761397279995, 0.000047998179978742605,
      # (phi, phi)
      -0.000069088696086429638, -0.0068541628780712289, -0.31885641694258540,
      0.14066701683364557, 0.0017092682062478156, 0.000017126410446066421
    ),
    dim = c(6, 4, 4)
  )

  gradient_hessian <- logistic4_gradient_hessian(x, theta)

  expect_type(gradient_hessian, "list")
  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 6 * 4)
  expect_length(gradient_hessian$H, 6 * 4 * 4)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)
})

context("4-parameter logistic - RSS functions")

test_that("Value of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2)

  n <- c(3, 3, 2, 4, 3)

  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_value <- 0.11303184522146127

  rss <- logistic4_rss(stats)

  expect_type(rss, "closure")

  value <- rss(theta)

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)

  known_param <- c(4 / 100, NA, NA, -3 / 2)
  rss <- logistic4_rss_fixed(stats, known_param)

  expect_type(rss, "closure")

  value <- rss(c(log(43 / 50), -2))

  expect_type(value, "double")
  expect_length(value, 1)
  expect_equal(value, true_value)
})

test_that("Gradient and Hessian of the RSS", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2)

  n <- c(3, 3, 2, 4, 3)

  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000,
    177363 / 25000000, 560629 / 112500000
  )

  stats <- cbind(x, n, m, v)

  true_gradient <- c(
    -0.78088225962005132, -0.28930436711562624, 0.011139573467794527,
    -0.087637515178570666
  )

  true_hessian <- matrix(
    c(
      # alpha
      15.000000000000000, 6.7515177403799487, 0.029409719813073703,
      0.80294599111504562,
      # log_omega
      6.7515177403799487, 5.1717341134316600, -0.13104549720149550,
      0.27712690508157734,
      # eta
      0.029409719813073703, -0.13104549720149550, 0.022131256419080663,
      -0.045877037358683589,
      # phi
      0.80294599111504562, 0.27712690508157734, -0.045877037358683589,
      0.19227987353030561
    ),
    nrow = 4,
    ncol = 4
  )

  gh <- logistic4_rss_gradient_hessian(stats)

  expect_type(gh, "closure")

  gradient_hessian <- gh(theta)

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 4)
  expect_length(gradient_hessian$H, 4 * 4)

  expect_equal(gradient_hessian$G, true_gradient)
  expect_equal(gradient_hessian$H, true_hessian)

  known_param <- c(4 / 100, NA, NA, -3 / 2)
  gh <- logistic4_rss_gradient_hessian_fixed(stats, known_param)

  expect_type(gh, "closure")

  gradient_hessian <- gh(c(log(43 / 50), -2))

  expect_type(gradient_hessian$G, "double")
  expect_type(gradient_hessian$H, "double")

  expect_length(gradient_hessian$G, 2)
  expect_length(gradient_hessian$H, 2 * 2)

  expect_equal(gradient_hessian$G, true_gradient[2:3])
  expect_equal(gradient_hessian$H, true_hessian[2:3, 2:3])
})

context("4-parameter logistic - general functions")

test_that("logistic4_init", {
  x <- -log(c(1000, 100, 10, 1, 0.1))
  theta <- c(4 / 100, log(43 / 50), -2, -3 / 2)

  n <- c(3, 3, 2, 4, 3)

  m <- c(376 / 375, 3091 / 3750, 8989 / 10000, 1447 / 10000, 11 / 120)

  v <- c(
    643663 / 450000000, 31087 / 112500000, 961 / 160000, 177363 / 25000000,
    560629 / 112500000
  )

  true_value <- c(
    0.10241913285726964, -0.17645537182325855, -2, -1.4802332740676008
  )

  stats <- cbind(x, n, m, v)

  start <- logistic4_init(stats)

  expect_type(start, "double")
  expect_length(start, 4)
  expect_equal(start, true_value)
})

test_that("logistic4_fisher_info_normal", {
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
    phi = -3 / 2
  )

  sigma <- 0.05

  true_value <- matrix(c(
      # alpha
      5035.9390472283091, 203.94858412428799, 100.15628253319072,
      181.09830066706291, 0,
      # beta
      203.94858412428799, 2556.1637845231149, -62.902195809034419,
      169.69326402671243, 0,
      # eta
      100.15628253319072, -62.902195809034419, 14.212019449141211,
      -9.3705833463297613, 0,
      # phi
      181.09830066706291, 169.69326402671243, -9.3705833463297613,
      57.950059278822026, 0,
      # sigma
      rep(0, 4), 6800
    ),
    nrow = 5,
    ncol = 5
  )

  rownames(true_value) <- colnames(true_value) <- c(
    "alpha", "beta", "eta", "phi", "sigma"
  )

  fim <- logistic4_fisher_info_normal(stats, n, theta, sigma)

  expect_type(fim, "double")
  expect_length(fim, 5 * 5)
  expect_equal(fim, true_value)
})

context("4-parameter logistic - fit")

test_that("logistic4_fit: fit", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.14236056369990992,
    beta = 0.14236056369990992 + exp(-0.23272371082134968),
    eta = -0.89523277708641233,
    phi = 0.13714073896752300
  )

  rss <- 0.030080292750857811

  fitted_values <- c(
    rep(0.93329083211408324, 3), rep(0.92353786225971087, 2),
    rep(0.85458328542760035, 2), rep(0.56283675776095414, 5),
    rep(0.24201206136017816, 3), rep(0.15661550382373062, 4),
    0.14420322010290019
  )

  residuals <- c(
    -0.00529083211408324, -0.04529083211408324, 0.04670916788591676,
    0.02446213774028913, -0.06753786225971087, 0.04241671457239965,
    0.02841671457239965, -0.07483675776095414, -0.03083675776095414,
    0.02316324223904586, 0.00316324223904586, 0.03616324223904586,
    0.01698793863982184, 0.02298793863982184, 0.00098793863982184,
    -0.03961550382373062, -0.01361550382373062, 0.02138449617626938,
    0.06238449617626938, -0.05220322010290019
  )

  result <- logistic4_fit(x, y, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  result <- logistic4_fit(x, y, c(0, 0, -1, 0), max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: inequalities", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.12654836760564175,
    beta = 0.12654836760564175 + exp(-0.2),
    eta = -0.81946721811752066,
    phi = 0.16558638357709105
  )

  rss <- 0.031021869952591774

  fitted_values <- c(
    rep(0.94279939044917876, 3), rep(0.92918328145981286, 2),
    rep(0.84963760209486180, 2), rep(0.56364513570595485, 5),
    rep(0.24760017802035200, 3), rep(0.14753394255411784, 4),
    0.12979816009612679
  )

  residuals <- c(
    -0.01479939044917876, -0.05479939044917876, 0.03720060955082124,
    0.01881671854018714, -0.07318328145981286, 0.04736239790513820,
    0.03336239790513820, -0.07564513570595485, -0.03164513570595485,
    0.02235486429404515, 0.00235486429404515, 0.03535486429404515,
    0.01139982197964800, 0.01739982197964800, -0.00460017802035200,
    -0.03053394255411784, -0.00453394255411784, 0.03046605744588216,
    0.07146605744588216, -0.03779816009612679
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(0.5, 0, -0.5, 0.5), max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(2, -1, -3, -1), max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 1.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: equalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.48361565993858148,
    phi = 0.55079105348629506
  )

  rss <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552128, 3), rep(0.92367933936267506, 2),
    rep(0.79901316975086137, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936762, 3), rep(0.12339358998664541, 4),
    0.04417373430922942
  )

  residuals <- c(
    -0.045588472981552128, -0.085588472981552128, 0.006411527018447872,
    0.02432066063732494, -0.06767933936267506, 0.09798683024913863,
    0.08398683024913863, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936762, -0.03497945701936762, -0.05697945701936762,
    -0.00639358998664541, 0.01960641001335459, 0.05460641001335459,
    0.09560641001335459, 0.04782626569077058
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with same equalities
  result <- logistic4_fit_constrained(
    x, y, c(0, 0, -1, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values with different equalities
  result <- logistic4_fit_constrained(
    x, y, c(1, 3, -1, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

test_that("logistic4_fit_constrained: equalities and inequalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.48361565993858148,
    phi = 0.55079105348629506
  )

  rss <- 0.061009639061112450

  fitted_values <- c(
    rep(0.973588472981552128, 3), rep(0.92367933936267506, 2),
    rep(0.79901316975086137, 2), rep(0.56620181868293205, 5),
    rep(0.29997945701936762, 3), rep(0.12339358998664541, 4),
    0.04417373430922942
  )

  residuals <- c(
    -0.045588472981552128, -0.085588472981552128, 0.006411527018447872,
    0.02432066063732494, -0.06767933936267506, 0.09798683024913863,
    0.08398683024913863, -0.07820181868293205, -0.03420181868293205,
    0.01979818131706795, -0.00020181868293205, 0.03279818131706795,
    -0.04097945701936762, -0.03497945701936762, -0.05697945701936762,
    -0.00639358998664541, 0.01960641001335459, 0.05460641001335459,
    0.09560641001335459, 0.04782626569077058
  )

  result <- logistic4_fit_constrained(
    x, y, NULL, max_iter,
    lower_bound = c(0, 0, -1, 0),
    upper_bound = c(0, 0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(0, 0, -0.8, 0.5), max_iter,
    lower_bound = c(0, 0, -1, 0),
    upper_bound = c(0, 0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_constrained(
    x, y, c(1, 2, -5, -1), max_iter,
    lower_bound = c(0, 0, -1, 0),
    upper_bound = c(0, 0,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, length(y) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("4-parameter logistic - weighted fit")

test_that("logistic4_fit_weighted: fit", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.17286261579329588,
    beta = 0.17286261579329588 + exp(-0.26243436674316789),
    eta = -0.96536781199483179,
    phi = -0.0061161452508493088
  )

  rss <- 0.015426894657174924

  fitted_values <- c(
    rep(0.94105809852012520, 3), rep(0.93306866024846045, 1),
    rep(0.86650297102635021, 2), rep(0.55631567380252901, 4),
    rep(0.24759852669347848, 3), rep(0.18172932146179900, 3),
    0.17383247345211985
  )

  residuals <- c(
    -0.01305809852012520, -0.05305809852012520, 0.03894190147987480,
    0.01493133975153955, 0.03049702897364979, 0.01649702897364979,
    -0.06831567380252901, -0.02431567380252901, 0.00968432619747099,
    0.04268432619747099, 0.01140147330652152, 0.01740147330652152,
    -0.00459852669347848, -0.03872932146179900, -0.00372932146179900,
    0.03727067853820100, -0.08183247345211985
  )

  result <- logistic4_fit_weighted(x, y, w, NULL, max_iter)

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  result <- logistic4_fit_weighted(x, y, w, c(0, 0, -1, 0), max_iter)

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

test_that("logistic4_fit_weighted_constrained: inequalities", {
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

  estimated <- c(alpha = TRUE, beta = TRUE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0.13727313366883312,
    beta = 0.13727313366883312 + exp(-0.2),
    eta = -0.79858420471931035,
    phi = 0.078372236750389013
  )

  rss <- 0.018015759502831172

  fitted_values <- c(
    rep(0.95292434651090841, 3), rep(0.93700817267641609, 1),
    rep(0.84963993729419210, 2), rep(0.55944477087279272, 4),
    rep(0.25576749618362561, 3), rep(0.15873543600035093, 3),
    0.14076156632218356
  )

  residuals <- c(
    -0.02492434651090841, -0.06492434651090841, 0.02707565348909159,
    0.01099182732358391, 0.04736006270580790, 0.03336006270580790,
    -0.07144477087279272, -0.02744477087279272, 0.00655522912720728,
    0.03955522912720728, 0.00323250381637439, 0.00923250381637439,
    -0.01276749618362561, -0.01573543600035093, 0.01926456399964907,
    0.06026456399964907, -0.04876156632218356
  )

  result <- logistic4_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 2.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values within the boundaries
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(0.5, 0, -0.5, 0.5), max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 2.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])

  # initial values outside the boundaries
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(2, -1, -3, -1), max_iter,
    lower_bound = c(0, -0.2, -1, 0),
    upper_bound = c(1,  0.2,  0, 1)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta, tolerance = 2.0e-7)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 4)
  expect_equal(result$fitted.values, fitted_values, tolerance = 1.0e-7)
  expect_equal(result$residuals, residuals, tolerance = 1.0e-7)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic4_fit_weighted_constrained: equalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.45801428680160602,
    phi = 0.56172839498959692
  )

  rss <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275455, 3), rep(0.91423176703152800, 1),
    rep(0.78786209884225038, 2), rep(0.56396744745889192, 4),
    rep(0.31055470762448156, 3), rep(0.13565276472766772, 3),
    0.05182418247944707
  )

  residuals <- c(
    -0.04036185505275455, -0.08036185505275455, 0.01163814494724545,
    0.03376823296847200, 0.10913790115774962, 0.09513790115774962,
    -0.07596744745889192, -0.03196744745889192, 0.00203255254110808,
    0.03503255254110808, -0.05155470762448156, -0.04555470762448156,
    -0.06755470762448156, 0.00734723527233228, 0.04234723527233228,
    0.08334723527233228, 0.04017581752055293
  )

  result <- logistic4_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with same equalities
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(0, 0, -1, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])

  # initial values with different equalities
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(1, 3, -1, 0), max_iter,
    lower_bound = c(0, 0, -Inf, -Inf),
    upper_bound = c(0, 0,  Inf,  Inf)
  )

  expect_true(result$converged)
  expect_false(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
  expect_equal(result$weights, w[w > 0])
})

test_that("logistic4_fit_weighted_constrained: equalities and inequalities", {
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

  estimated <- c(alpha = FALSE, beta = FALSE, eta = TRUE, phi = TRUE)

  theta <- c(
    alpha = 0,
    beta = 1,
    eta = -0.45801428680160602,
    phi = 0.56172839498959692
  )

  rss <- 0.040269628142735194

  fitted_values <- c(
    rep(0.96836185505275455, 3), rep(0.91423176703152800, 1),
    rep(0.78786209884225038, 2), rep(0.56396744745889192, 4),
    rep(0.31055470762448156, 3), rep(0.13565276472766772, 3),
    0.05182418247944707
  )

  residuals <- c(
    -0.04036185505275455, -0.08036185505275455, 0.01163814494724545,
    0.03376823296847200, 0.10913790115774962, 0.09513790115774962,
    -0.07596744745889192, -0.03196744745889192, 0.00203255254110808,
    0.03503255254110808, -0.05155470762448156, -0.04555470762448156,
    -0.06755470762448156, 0.00734723527233228, 0.04234723527233228,
    0.08334723527233228, 0.04017581752055293
  )

  result <- logistic4_fit_weighted_constrained(
    x, y, w, NULL, max_iter,
    lower_bound = c(0, 0, -2, 0),
    upper_bound = c(0, 0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values within the boundaries
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(0, 0, -0.8, 0.5), max_iter,
    lower_bound = c(0, 0, -2, 0),
    upper_bound = c(0, 0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)

  # initial values outside the boundaries
  result <- logistic4_fit_weighted_constrained(
    x, y, w, c(1, 2, -5, -1), max_iter,
    lower_bound = c(0, 0, -2, 0),
    upper_bound = c(0, 0,  0, 2)
  )

  expect_true(result$converged)
  expect_true(result$constrained)
  expect_equal(result$estimated, estimated)
  expect_equal(result$coefficients, theta)
  expect_equal(result$rss, rss)
  expect_equal(result$df.residual, sum(w > 0) - 2)
  expect_equal(result$fitted.values, fitted_values)
  expect_equal(result$residuals, residuals)
})

context("4-parameter logistic - drda fit")

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
      y ~ x, mean_function = "logistic4",
      lower_bound = c("a", "b", "c", "d")
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = matrix(-Inf, nrow = 4, ncol = 2),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 4)
    ),
    "'lower_bound' and 'upper_bound' must have the same length"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c( 0, -Inf, -Inf, -Inf),
      upper_bound = c(-1, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be larger than 'upper_bound'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c(Inf, -Inf, -Inf, -Inf),
      upper_bound = c(Inf, Inf, Inf, Inf)
    ),
    "'lower_bound' cannot be equal to infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' and 'upper_bound' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c(0, -1, -Inf, -Inf),
      upper_bound = rep(Inf, 4)
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
      y ~ x, mean_function = "logistic4",
      upper_bound = c("a", "b", "c", "d")
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 4),
      upper_bound = matrix(Inf, nrow = 4, ncol = 2)
    ),
    "'upper_bound' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = c(-Inf, -Inf, -Inf, -Inf),
      upper_bound = c(-Inf, Inf, Inf, Inf)
    ),
    "'upper_bound' cannot be equal to -infinity"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 5),
      upper_bound = rep(Inf, 5)
    ),
    "'lower_bound' and 'upper_bound' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      lower_bound = rep(-Inf, 4),
      upper_bound = c(1, 0, Inf, Inf)
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
      y ~ x, mean_function = "logistic4",
      start = c("a", "b", "c", "d")
    ),
    "'start' must be a numeric vector"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, Inf, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(-Inf, 0, -1, 0)
    ),
    "'start' must be finite"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 0, -1, 0, 1)
    ),
    "'start' must be of length 4"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, -1, -1, 0)
    ),
    "parameter 'beta' is smaller than 'alpha'"
  )

  expect_error(
    drda(
      y ~ x, mean_function = "logistic4",
      start = c(0, 1, 0, 0)
    ),
    "parameter 'eta' cannot be initialized to zero"
  )
})
