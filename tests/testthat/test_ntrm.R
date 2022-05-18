test_that("ntrm_solve_tr_subproblem: random vector", {
  set.seed(4705968)

  G <- c(-6.26752849, -1.40272211, -0.02658175, -1.85259931)
  H <- matrix(
    c(
      15, 21.188528, -0.5116276, 3.120284,
      21.1885279, 48.405099, -1.4542101, 3.925457,
      -0.5116276, -1.454210, -0.3024076, -1.191577,
      3.1202843, 3.925457, -1.1915771, 3.543116
    ),
    nrow = 4,
    ncol = 4
  )

  m <- function(x) {
    sum(G * x) + sum(x * (H %*% x)) / 2
  }

  result <- ntrm_solve_tr_subproblem(G, H, 1)

  expect_equal(result$m, m(result$p))

  for (j in seq_len(10)) {
    # generate a random vector and set it on the boundary
    bad_p <- runif(4)
    bad_p <- bad_p / sqrt(sum(bad_p^2))

    expect_lte(result$m, m(bad_p) + 1.0e-8)
  }
})

test_that("ntrm_solve_tr_subproblem: random Hessian", {
  set.seed(4306006)

  for (i in seq_len(1000)) {
    G <- rnorm(4)
    H <- matrix(rnorm(16), nrow = 4, ncol = 4)
    H <- H + t(H)

    result <- ntrm_solve_tr_subproblem(G, H, 1)

    m <- function(x) {
      sum(G * x) + sum(x * (H %*% x)) / 2
    }

    expect_equal(result$m, m(result$p))
    expect_lte(result$m, m(rep(0, 4)) + 1.0e-8)

    for (j in seq_len(10)) {
      # generate a random vector and set it on the boundary
      bad_p <- runif(4)
      bad_p <- bad_p / sqrt(sum(bad_p^2))

      expect_lte(result$m, m(bad_p) + 1.0e-8)

      bad_p <- bad_p * runif(1)

      expect_lte(result$m, m(bad_p) + 1.0e-8)
    }
  }
})

test_that("ntrm_solve_tr_subproblem: various situations", {
  set.seed(7160982)

  H <- matrix(runif(16), nrow = 4, ncol = 4)
  H <- t(H) * H + 4 * diag(4)

  H_eig <- eigen(H, symmetric = TRUE)
  U <- H_eig$vectors

  G <- rep(0, 4)
  G[1] <- 1

  true_p <- -solve(H, G)
  true_p_sn <- sum(true_p^2)

  true_m <- sum(G * true_p) + sum(true_p * (H %*% true_p)) / 2

  # An interior solution
  delta <- sqrt(true_p_sn) + 1
  result <- ntrm_solve_tr_subproblem(G, H, delta)

  expect_lte(sum(result$p^2), delta^2)
  expect_equal(result$m, true_m)
  expect_equal(result$p, true_p)

  # A boundary solution
  delta <- sqrt(true_p_sn) / 2
  result <- ntrm_solve_tr_subproblem(G, H, delta)

  expect_gt(result$m, true_m)
  expect_equal(sqrt(sum(result$p^2)), delta)

  # hard case problem
  L <- rep(0.1, 4)
  L[1] <- -1

  H <- U %*% diag(L) %*% t(U)
  H <- (t(H) + H) / 2
  G <- U[, 2]

  expect_true(isSymmetric(H))
  expect_lt(abs(sum(G * U[, 1])), 1.0e-12)

  true_p <- -solve(H, G)
  true_p_sn <- sum(true_p^2)
  true_m <- sum(G * true_p) + sum(true_p * (H %*% true_p)) / 2

  delta <- sqrt(true_p_sn) / 2
  result <- ntrm_solve_tr_subproblem(G, H, delta)

  expect_equal(sqrt(sum(result$p^2)), delta)
})

test_that("ntrm: Rosenbrock function", {
  set.seed(2823441)

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

  # try different random starting points
  for (i in 1:10) {
    result <- ntrm(fn, gh, rnorm(2, sd = 3), 1000)
    expect_equal(result$optimum, c(1, 1))
  }
})
