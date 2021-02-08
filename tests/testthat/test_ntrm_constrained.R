context("Interior point Newton trust region method")

test_that("ntrm_constrained: Rosenbrock function", {
  set.seed(3470830)

  max_iter <- 10000

  fn <- function(x) {
    (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  }

  gradient <- function(x) {
    c(2 * (200 * x[1]^3 - 200 * x[1] * x[2] + x[1] - 1), 200 * (x[2] - x[1]^2))
  }

  hessian <- function(x) {
    matrix(
      c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1], -400 * x[1], 200),
      nrow = 2, ncol = 2
    )
  }

  gh <- function(x) {
    list(G = gradient(x), H = hessian(x))
  }

  # feasible region contains the global minimum

  # initial point is within the feasible region
  result <- ntrm_constrained(fn, gh, c(2.5, -0.8), max_iter, c(0, -1), c(3,  5))
  expect_equal(result$optimum, c(1, 1))

  result <- ntrm_constrained(
    fn, gh, c(2.5, -0.8), max_iter, c(-Inf, -3), c(3,  Inf)
  )
  expect_equal(result$optimum, c(1, 1))

  # initial point is outside the feasible region
  result <- ntrm_constrained(fn, gh, c(-5, 10), max_iter, c(0, -1), c(3,  5))
  expect_equal(result$optimum, c(1, 1))

  result <- ntrm_constrained(
    fn, gh, c(10, -10), max_iter, c(-Inf, -3), c(3,  Inf)
  )
  expect_equal(result$optimum, c(1, 1))

  # try different random starting points
  for (i in 1:5) {
    result <- ntrm_constrained(
      fn, gh, rnorm(2, sd = 3), max_iter, c(0, -1), c(3,  5)
    )
    expect_equal(result$optimum, c(1, 1))

    result <- ntrm_constrained(
      fn, gh, rnorm(2, sd = 3), max_iter, c(-Inf, -3), c(3,  Inf)
    )
    expect_equal(result$optimum, c(1, 1))
  }

  # feasible region does not contain the global minimum

  # initial point is within the feasible region
  result <- ntrm_constrained(
    fn, gh, c(0, 0), max_iter, c(-0.5, -0.5), c(0.5, 0.5)
  )
  expect_equal(result$optimum, c(0.5, 0.25))

  result <- ntrm_constrained(
    fn, gh, c(-10, 10), max_iter, c(-Inf, 2), c(0,  Inf)
  )
  expect_equal(result$optimum, c(-1.4111899, 2))

  # initial point is outside the feasible region
  result <- ntrm_constrained(
    fn, gh, c(1, -1), max_iter, c(-0.5, -0.5), c(0.5, 0.5)
  )
  expect_equal(result$optimum, c(0.5, 0.25))

  result <- ntrm_constrained(
    fn, gh, c(10, -10), max_iter, c(-Inf, 2), c(0,  Inf)
  )
  expect_equal(result$optimum, c(-1.4111899, 2))

  # try different random starting points
  for (i in 1:5) {
    result <- ntrm_constrained(
      fn, gh, rnorm(2, sd = 3), max_iter, c(-0.5, -0.5), c(0.5, 0.5)
    )
    expect_equal(result$optimum, c(0.5, 0.25))

    result <- ntrm_constrained(
      fn, gh, rnorm(2, sd = 3), max_iter, c(-Inf, 2), c(0,  Inf)
    )
    expect_equal(result$optimum, c(-1.4111899, 2))
  }
})
