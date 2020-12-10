context("Basic auxiliary functions")

test_that("ntrm_max_abs", {
  x <- c(
    -0.186298, -0.370595, -2.902886, -1.392832, -0.024334, 0.174643, -0.817747,
    -0.555250, 0.161866, 1.167719
  )

  expect_equal(ntrm_max_abs(x), 2.902886)
})

test_that("ntrm_norm", {
  x <- c(
    -0.186298, -0.370595, -2.902886, -1.392832, -0.024334, 0.174643, -0.817747,
    -0.555250, 0.161866, 1.167719
  )

  expect_equal(ntrm_norm(x), 3.5967553450964664)

  x <- c(
    -0.0517961279, -0.1030359211, -0.8070846420, -0.3872468006, -0.0067655422,
    0.0485557074, -0.2273568596, -0.1543752484, 0.0450033390, 0.3246590018
  )

  expect_equal(ntrm_norm(x), 1)
})

test_that("ntrm_solve_quadratic_equation", {
  # simple equations
  result <- ntrm_solve_quadratic_equation(1, -1, -1)
  expect_equal(result, c(-0.618033988749895, 1.61803398874989))

  result <- ntrm_solve_quadratic_equation(pi, -1, -1)
  expect_equal(result, c(-0.427053366380407, 0.745363252564198))

  # numerically unstable equations
  result <- ntrm_solve_quadratic_equation(1, 200, -0.000015)
  expect_equal(result, c(-200.000000075, 0.000000075))

  result <- ntrm_solve_quadratic_equation(
    1, -1.786737601482363, 2.054360090947453e-8
  )
  expect_equal(result, c(1.149782767465722e-8, 1.786737589984535))

  # zero discriminant
  result <- ntrm_solve_quadratic_equation(1, 2, 1)
  expect_equal(result, rep(-1, 2))

  result <- ntrm_solve_quadratic_equation(pi^2, 2 * pi, 1)
  expect_equal(result, rep(-1 / pi, 2))

  # zero coefficients
  result <- ntrm_solve_quadratic_equation(1, 0, -pi^2)
  expect_equal(result, c(-pi, pi))

  result <- ntrm_solve_quadratic_equation(pi, 0, -pi^2)
  expect_equal(result, c(-sqrt(pi), sqrt(pi)))

  result <- ntrm_solve_quadratic_equation(1, 1, 0)
  expect_equal(result, c(-1, 0))

  result <- ntrm_solve_quadratic_equation(pi, 1, 0)
  expect_equal(result, c(-1 / pi, 0))

  result <- ntrm_solve_quadratic_equation(1, 0, 0)
  expect_equal(result, c(0, 0))

  # almost zero coefficients
  result <- ntrm_solve_quadratic_equation(1, 1.0e-13, -pi^2)
  expect_equal(result, c(-pi, pi))

  result <- ntrm_solve_quadratic_equation(pi, 1.0e-13, -pi^2)
  expect_equal(result, c(-sqrt(pi), sqrt(pi)))

  result <- ntrm_solve_quadratic_equation(1, 1, 1.0e-13)
  expect_equal(result, c(-1, 0))

  result <- ntrm_solve_quadratic_equation(pi, 1, 1.0e-13)
  expect_equal(result, c(-1 / pi, 0))

  result <- ntrm_solve_quadratic_equation(1, 1.0e-13, 1.0e-13)
  expect_equal(result, c(0, 0))

  # negative discriminant
  expect_error(ntrm_solve_quadratic_equation(1, 1, pi))
  expect_error(ntrm_solve_quadratic_equation(1, 0, 1))
  expect_error(ntrm_solve_quadratic_equation(1, 1.0e-13, 1))
})

test_that("ntrm_line_search", {
  expect_equal(ntrm_line_search(function(x) 2 * x^2 - x + 3, 0, 1), 0.25)
  expect_equal(ntrm_line_search(function(x) -dnorm(x), -10, 10), 0)
})
