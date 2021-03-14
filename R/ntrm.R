# Solution to the (unconstrained) trust-region sub-problem
#
# Compute the current best step of the trust-region sub-problem.
#
# @param Q numeric matrix of dot products between eigenvectors and gradient.
# @param X numeric matrix of eigenvectors.
# @param L numeric vector of eigenvalues.
# @param lambda lambda parameter of the trust-region method.
#
# @return Current best solution to the trust-region sub-problem.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_calc_p <- function(Q, X, L, lambda) {
  # Equation (4.38) at page 84 from Nocedal and Wright (2006)
  - colSums((t(X) * Q) / (L + lambda))
}

# Safeguard parameter computation
#
# Update the current value of `lambda` to hopefully obtain a non-singular
# matrix.
#
# @param lambda lambda parameter of the trust-region method.
# @param G gradient vector.
# @param H hessian matrix.
# @param delta radius of the trust region.
#
# @return Updated scalar value `lambda`.
#
# @references
# Jorge J Moré and D C Sorensen. Computing a Trust Region Step.
# **SIAM Journal on Scientific and Statistical Computing**, 4(3):553-572, 1983.
# doi: 10.1137/0904038.
ntrm_safeguard <- function(lambda, G, H, delta) {
  # Equations are on page 560 of Moré and Sorensen (1983)
  lambda_S <- max(-diag(H))

  # squared Euclidean norm
  G_norm <- sum(G^2)

  # operator norm: maximum absolute column sum
  H_norm <- max(colSums(abs(H)))

  lambda_L <- max(0, lambda_S, G_norm / delta - H_norm)
  lambda_U <- G_norm / delta + H_norm

  # page 558
  lambda <- min(max(lambda, lambda_L), lambda_U)
  if (lambda <= lambda_S) {
    lambda <- max(lambda_U / 1000, sqrt(lambda_L * lambda_U))
  }

  lambda
}

# Iterative solution to the (unconstrained) trust-region sub-problem
#
# Compute the best step of the trust-region sub-problem.
#
# @param G gradient vector.
# @param H hessian matrix.
# @param delta radius of the trust region.
#
# @return List with the best step solution `p` to the trust-region sub-problem
#   and the value `m` of the objective function `m(p)`.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_solve_tr_subproblem <- function(G, H, delta) {
  eps <- sqrt(.Machine$double.eps)
  delta_2 <- delta^2

  k <- length(G)

  H_eig <- eigen(H, symmetric = TRUE)
  H_ev_min <- H_eig$values[k]

  # dot products between eigenvectors and gradient
  Q <- colSums(H_eig$vectors * G)

  p <- rep(0, k)

  # Unconstrained solution
  if (H_ev_min > eps) {
    p <- ntrm_calc_p(Q, H_eig$vectors, H_eig$values, 0)
  }

  if ((H_ev_min < eps) || (sum(p^2) > delta_2)) {
    # either hard case or we went outside the trust region
    hard_case <- FALSE

    if (abs(Q[k]) < eps) {
      # this might be a hard case because the gradient is orthogonal to all
      # eigenvectors associated with the lowest eigenvalue.
      #
      # lambda is taken to be -H_ev_min and we only need to find a multiple of
      # an orthogonal eigenvector that lands on the boundary

      # Equation (4.45) at page 88 from Nocedal and Wright (2006)
      j <- 1:max(which(H_eig$values > H_ev_min))
      p <- ntrm_calc_p(
        Q[j], H_eig$vectors[, j], H_eig$values[j], -H_ev_min
      )

      p_norm_2 <- sum(p^2)

      if (p_norm_2 <= delta_2) {
        hard_case <- TRUE
        tau <- sqrt(delta_2 - p_norm_2)
        p <- tau * H_eig$vectors[, k] - p
      }
    }

    if (!hard_case) {
      # not hard case, we can apply the root finding technique
      # Algorithm 4.3 at page 87 in Nocedal and Wright (2006)
      B <- H

      # lambda cannot be lower than this lower bound
      lambda_lb <- -H_ev_min + 1.0e-12
      lambda <- ntrm_safeguard(lambda_lb, G, H, delta)

      for (i in seq_len(10)) {
        diag(B) <- diag(H) + lambda

        R <- tryCatch(chol(B), error = function(e) NULL)

        if (is.null(R)) {
          lambda <- 10 * lambda
          next
        }

        p <- -solve(R, solve(t(R), G))
        q <- solve(t(R), p)

        p_norm_2 <- sum(p^2)
        q_norm_2 <- sum(q^2)

        lambda_old <- lambda
        lambda <- lambda +
          (p_norm_2 / q_norm_2) *
            (sqrt(p_norm_2) - delta) /
            delta

        if (lambda < lambda_lb) {
          lambda <- (lambda_lb + lambda_old) / 2
        }

        if (abs(lambda - lambda_old) < eps) {
          break
        }
      }
    }
  }

  # Equation (4.2) at page 68 in Nocedal and Wright (2006)
  # function f is dropped because it does not depend on p and also because it
  # disappears when computing m(0) - m(p)
  m <- sum(G * p) + 0.5 * sum(p * (H %*% p))

  list(
    p = p,
    m = m
  )
}

# Newton with Trust Region Method
#
# Find the minimum of a function using the Newton's method, combined with the
# Trust Region method.
#
# @details
# The code is adapted from the Julia implementation of
# [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl). See
# 'src/multivariate/solvers/second_order/newton_trust_region.jl'
# file for the original source code.
#
# @param fn function handle to evaluate the function to minimize.
# @param gh function handle to evaluate both the gradient and Hessian of the
#   function to minimize.
# @param init numeric vector of starting values.
# @param max_iter integer value for the maximum number of iterations.
# @param update_fn function handle to update a subset of parameters with the
#   actual optimum.
#
# @return A list containing the value `optimum` for which `f(optimum)` is
#   minimum, the value `minimum` containing `f(optimum)`, and a boolean
#   value `converged` stating if the algorithm converged or not.
#
# @references
# Jorge J Moré and D C Sorensen. Computing a Trust Region Step.
# **SIAM Journal on Scientific and Statistical Computing**, 4(3):553-572, 1983.
# doi: 10.1137/0904038.
#
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
#
# Patrick Kofod Mogensen and Asbjørn Nilsen Riseth.
# Optim: A mathematical optimization package for Julia.
# **Journal of Open Source Software**, 3(24):615, 2018.
# doi: 10.21105/joss.00615.
ntrm <- function(fn, gh, init, max_iter, update_fn = NULL) {
  eps <- 1.0e-10
  converged <- FALSE

  delta <- 1

  cur_optimum <- init
  cur_minimum <- fn(cur_optimum)

  f_converged <- 0
  x_converged <- FALSE
  g_converged <- FALSE
  for (i in seq_len(max_iter)) {
    if ((i %% 10 == 0) && !is.null(update_fn)) {
      cur_optimum <- update_fn(cur_optimum)
      cur_minimum <- fn(cur_optimum)
    }

    gradient_hessian <- gh(cur_optimum)

    g_min <- ntrm_max_abs(gradient_hessian$G)
    g_converged <- g_min <= eps

    if (!g_converged && any(is.infinite(gradient_hessian$H))) {
      # when the model is close to non-identifiability the gradient might be
      # close to zero, but not enough for our previous check to succeed
      #
      # if the Hessian contains infinite values, we relax the condition for
      # convergence
      g_converged <- g_min <= sqrt(.Machine$double.eps)

      if (!g_converged) {
        # this was our last try because we cannot continue our search
        break
      }
    }

    # we use a counter for f_converged because the objective function might
    # be very flat around the optimum. Once we hit a flat region, we explore the
    # function 10 more times until we give up.
    if (x_converged || g_converged || (f_converged > 10)) {
      # this solution is good enough
      converged <- TRUE
      break
    }

    # perform a correction in case the Hessian is not positive definite
    B <- ntrm_correct_hessian(gradient_hessian$H)
    result <- ntrm_solve_tr_subproblem(gradient_hessian$G, B, delta)

    candidate <- cur_optimum + result$p

    cur_value <- fn(candidate)
    cur_diff <- cur_minimum - cur_value

    rho <- if (result$m > 0) {
      # This can happen if the trust region radius is too large and the
      # Hessian is not positive definite. We should shrink the trust region.
      -1
    } else if (abs(result$m) <= eps) {
      # This should only happen when the step is very small, in which case
      # we should accept the step and assess convergence
      1
    } else if (is.nan(cur_value)) {
      # We went somewhere wrong
      -1
    } else {
      # m(0) = f
      # m(p) = f + G' p + 0.5 p' H p
      # m(0) - m(p) = -(G'p + 0.5 p' H p) = -result$m
      -cur_diff / result$m
    }

    if (rho < 0.25) {
      delta <- delta / 4
    } else if ((rho > 0.75) && (ntrm_norm(result$p) >= delta)) {
      delta <- min(2 * delta, 100)
    }

    if (rho > 0.1) {
      if ((abs(cur_value - cur_minimum) / abs(cur_value)) <= eps) {
        f_converged <- f_converged + 1
      }

      x_converged <- {
        (max(abs(candidate - cur_optimum)) / max(abs(candidate))) <= eps
      }

      cur_optimum <- candidate
      cur_minimum <- cur_value
    }
  }

  # did it converge at the very last iteration?
  if (i == max_iter) {
    g_converged <- ntrm_max_abs(gh(cur_optimum)$G) <= eps
    if (x_converged || g_converged || (f_converged > 10)) {
      converged <- TRUE
    }
  }

  list(
    optimum = cur_optimum,
    minimum = cur_minimum,
    converged = converged,
    iterations = i
  )
}
