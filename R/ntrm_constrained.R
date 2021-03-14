# Actual reduction
#
# Evaluate the actual reduction `ared(p)` for the proposed step `p`.
#
# @details
# The actual reduction is defined as the difference of merit functions
# `phi(x, s) - phi(x + p_x, s + p_s)` where
# `phi(x, s) = f(x) - mu sum(log(s)) + nu |C(x) - s|`.
#
# @param nu penalty parameter.
# @param f_current current value of the original function, i.e. `f(x)`.
# @param f_candidate candidate value of original function, i.e. `f(x + p_x)`.
# @param s_current current value of the slack variables.
# @param s_candidate candidate value of the slack variables, i.e. `s + p_s`.
# @param m0_current norm of `C(x) - s`.
# @param m0_candidate norm of `C(x + p_x) - (s + p_s)`.
# @param mu barrier penalty parameter.
#
# @return Value of the actual reduction `ared(p)` for the proposed step `p`.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_ared <- function(
  nu, f_current, f_candidate, s_current, s_candidate, m0_current, m0_candidate,
  mu
) {
  # Equation (19.40) at page 582 of Nocedal and Wright (2006)
  ntrm_merit(nu, f_current, s_current, m0_current, mu) -
  ntrm_merit(nu, f_candidate, s_candidate, m0_candidate, mu)
}

# Process constraint requirements
#
# Given a set of lower bounds and upper bounds on the parameters, construct
# the Jacobian of the constraint inequalities.
#
# @param lower_bound numeric vector of lower bounds.
# @param upper_bound numeric vector of upper bounds.
#
# @return A list containing the Jacobian matrix `A` and the indices of the
#   finite lower bounds `idx_lb` and finite upper bounds `idx_ub` constraints.
ntrm_create_constraints <- function(lower_bound, upper_bound) {
  idx_lb <- which(!is.infinite(lower_bound))
  idx_ub <- which(!is.infinite(upper_bound))

  n <- length(lower_bound)

  m_lb <- length(idx_lb)
  m_ub <- length(idx_ub)
  m <- m_lb + m_ub

  # A is the Jacobian of the constraint equations
  # the first m_lb rows are associated with functions f(x_i) = (x_i - l_i)
  # the last m_ub rows are associated with functions f(x_i) = (u_i - x_i)
  A <- matrix(0, nrow = m, ncol = n)
  A[seq_len(m_lb) + m * (idx_lb - 1)] <- 1
  A[((m_lb + 1):m) + m * (idx_ub - 1)] <- -1

  list(
    idx_lb = idx_lb,
    idx_ub = idx_ub,
    A = A
  )
}

# Maximum error of the current solution
#
# Compute a measure of how far the Lagrangian gradient is from zero.
#
# @details
# The solution to the minimization problem is such that the gradient of our
# Lagrangian function is equal to zero. After splitting the full Lagrangian
# gradient into three vectors, the system of equations is the following:
#
# G - A' z = 0
# S Z - mu e = 0
# C(x) - s = 0
#
# where `G` is the gradient of the original function that must be minimized,
# `C(.)` are functions evaluating the inequality constraints (in our case
# they are the distances from the boundary), `A` is the Jacobian matrix of
# `C(.)`, `z` are the Lagrangian multipliers, `s` are the slack variables,
# `Z` is the diagonal matrix of Lagrange multipliers, `S` is the diagonal
# matrix of slack variables,`e` is a vector of ones, and `mu` is a penalty
# barrier parameter.
#
# We compute the maximum norm of the previous three vectors.
#
# @param obj list containing all relevant variables.
# @param mu barrier penalty parameter.
#
# @return Maximum of the three norms.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_error <- function(obj, mu) {
  # See equation (19.10) at page 567 of Nocedal and Wright (2006)
  max(
    ntrm_max_abs(obj$G - crossprod(obj$A, obj$z)),
    ntrm_max_abs(obj$s * obj$z - mu),
    obj$m0
  )
}

# Inertia correction and regularization
#
# Update the Hessian matrix so that the primal-dual matrix is positive
# definite.
#
# @param H Hessian matrix.
# @param A Jacobian matrix of the inequality constraints.
# @param slack_variable current value of the slack variables.
# @param lagrange_multiplier current value of the Lagrange multipliers.
# @param kappa perturbation used in the previous interior-point iteration.
#
# @return Updated value of the `kappa` parameter.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_inertia_correction <- function(
  H, A, slack_variable, lagrange_multiplier, kappa
) {
  eps <- sqrt(.Machine$double.eps)

  n_variables <- ncol(A)
  n_constraints <- nrow(A)
  n_slacks <- length(slack_variable)
  n_total <- n_variables + n_constraints

  idx <- seq_len(n_total)

  # X = | H  0  A' |
  #     | 0  V -I  |
  #     | A -I  0  |
  #
  # where V = Z / S is the diagonal matrix of lagrange multipliers divided by
  # the slack variables
  X <- rbind(
    cbind(
      H,
      matrix(0, nrow = n_variables, ncol = n_slacks),
      t(A)
    ),
    cbind(
      matrix(0, nrow = n_slacks, ncol = n_variables),
      diag(lagrange_multiplier / slack_variable, nrow = n_slacks),
      -diag(n_slacks)
    ),
    cbind(
      A,
      -diag(n_slacks),
      matrix(0, nrow = n_constraints, ncol = n_constraints)
    )
  )

  eigval <- eigen(X, symmetric = TRUE, only.values = TRUE)$values

  some_zero <- any(abs(eigval) < eps)
  some_negative <- any(eigval[idx] < -eps)
  some_positive <- any(eigval[-idx] > eps)

  if (some_zero || some_negative || some_positive) {
    if (kappa == 0) {
      kappa <- 1.0e-4
    } else {
      kappa <- kappa / 2
    }

    repeat {
      X <- rbind(
        cbind(
          H + kappa * diag(n_variables),
          matrix(0, nrow = n_variables, ncol = n_slacks),
          t(A)
        ),
        cbind(
          matrix(0, nrow = n_slacks, ncol = n_variables),
          diag(lagrange_multiplier / slack_variable, nrow = n_slacks),
          -diag(n_slacks)
        ),
        cbind(
          A,
          -diag(n_slacks),
          matrix(0, nrow = n_constraints, ncol = n_constraints)
        )
      )

      eigval <- eigen(X, symmetric = TRUE, only.values = TRUE)$values

      some_zero <- any(abs(eigval) < eps)
      some_negative <- any(eigval[idx] < -eps)
      some_positive <- any(eigval[-idx] > eps)

      if (some_zero || some_negative || some_positive) {
        kappa <- 10 * kappa
      } else {
        break
      }
    }
  }

  kappa
}

# Initialize solution object
#
# @param fn function handle to evaluate the function to minimize.
# @param gh function handle to evaluate both the gradient and Hessian of the
#   function to minimize.
# @param C function handle to evaluate the distances of the current solution
#   from the boundary.
# @param x initial values of the solution.
# @param s initial values of slack variables.
# @param A Jacobian of the constraints
# @param mu barrier penalty parameter.
#
# @return List with all important variables initialized.
ntrm_init_obj <- function(fn, gh, C, x, s, A, mu) {
  n_var <- ncol(A)
  n_slk <- nrow(A)
  n_tot <- n_var + n_slk

  obj <- list(
    # Jabobian of constraints
    A = A,
    # number of variables
    k = n_var,
    # number of constraints (slacks)
    m = n_slk,
    # number of unknowns
    n = n_var + n_slk,
    # index of original variables in the joint vector
    idx_x = seq_len(n_var),
    # index of slack variables in the joint vector
    idx_s = (n_var + 1):n_tot,
    # optimum
    x = x,
    # function value
    f = fn(x),
    # slack variables
    s = s,
    # slack residuals (distances from boundary)
    r = C(x) - s
  )

  # Equation (19.41+) at page 580
  obj$m0 <- ntrm_norm(obj$r)

  # Lagrangian Hessian (with respect to cur_optimum) is the same as the
  # Hessian of the function to optimize. This is because our constraints
  # are linear (either `x - lb - s` or `ub - x - s`) and the second
  # derivatives are all equal to zero
  gradient_hessian <- gh(obj$x)
  obj$G <- gradient_hessian$G
  obj$H <- gradient_hessian$H

  obj$B <- cbind(A, -diag(s, nrow = n_slk))
  obj$Y <- rbind(
    cbind(diag(n_tot), t(obj$B)),
    cbind(obj$B, matrix(0, nrow = n_slk, ncol = n_slk))
  )

  obj$u <- c(obj$G, -rep(mu, obj$m))

  # compute multipliers with equations (19.36)-(19.38) at page 581
  obj$z <- ntrm_lagrange_multiplier(obj, mu)

  obj
}

# Compute Lagrange multipliers
#
# Update the Lagrange multipliers given the current optimum and slack
# variables.
#
# @details
# Given the current value of the optimum and slack variables, find the
# least-squares Lagrange multipliers. Define the matrix `B = cbind(A, -S)`
# where `A` is the Jacobian of the constraint functions and `S` is the diagonal
# matrix of slack variables. Define the vector `u = c(G, -mu e)` where `G` is
# the gradient of the original function and `e` is a vector of ones.
# The least-squares multipliers are `z = (B B')^(-1) B u`. If `z[i]` is
# negative, it is substituted by the formula
# `min(0.001, mu / slack_variable[i])`.
#
# @param obj list containing all the relevant variables.
# @param mu barrier penalty parameter.
#
# @return Numeric vector of least-squares Lagrange multipliers.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_lagrange_multiplier <- function(obj, mu) {
  # Compute Lagrange multipliers with equations (19.36)-(19.38) at page 581
  # of Nocedal and Wright (2006)
  #
  # z = (B B')^(-1) B u => B B' z = B u
  #
  # Note that the system can be solved by our Y matrix:
  #
  # Y = | I   B' | | x | = |    0 |
  #     | B   0  | | z |   | -B u |
  idx_z <- (obj$n + 1):(obj$m + obj$n)
  b <- rep(0, obj$n)
  b[idx_z] <- -obj$B %*% obj$u
  z <- solve(obj$Y, b)[idx_z]

  # Lagrange multipliers cannot be negative
  idx_neg <- z < 0
  if (any(idx_neg)) {
    z[idx_neg] <- pmin(1.0e-3, mu / obj$s[idx_neg])
  }

  z
}

# Merit function
#
# Evaluate the merit function
# `phi(x, s) = f(x) - mu sum(log(s)) + nu |C(x) - s|`.
#
# @param nu penalty parameter.
# @param f value of the original function `f(x)` we are minimizing.
# @param s vector of slack variables.
# @param m0 norm of `C(x) - s`.
# @param mu barrier penalty parameter.
#
# @return Value of the merit function evaluated at the provided parameter.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_merit <- function(nu, f, s, m0, mu) {
  # Equation (19.26) at page 575
  f - mu * sum(log(s)) + nu * m0
}

# Objective function value
#
# Compute the value of the sub-problem objective function.
#
# @details
# The function evaluated here is
#
# `h(p) = 0.5 p' Q p + u' p`
#
# @param p numeric vector at which the function is to be evaluated.
# @param Q numeric matrix representing the quadratic part.
# @param u numeric vector representing the linear part.
#
# @return Value of the function.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_objective_function <- function(p, Q, u) {
  0.5 * sum(p * (Q %*% p)) + sum(p * u)
}

# Predicted reduction
#
# Evaluate the predicted reduction `pred(p)` for the proposed step `p`.
#
# @details
# The predicted reduction is defined as the difference `q(0, nu) - q(p, nu)`
# where `q(p, nu) = h(p_x, p_s) + nu m(p)`, `h(p_x, p_s)` is equal to
# `0.5 (p_x' H p_x + p_s' S^(-1) Z p_s) + (G' p_x - mu e' S^(-1) p_s)`, and
# `m(p) = |A p_x - p_s + C(x) - s|`. `G` is the gradient of the original
# function that must be minimized, `C(.)` are functions evaluating the
# inequality constraints (in our case they are the distances from the
# boundary), `A` is the Jacobian matrix of `C(.)`, `z` are the Lagrangian
# multipliers, `s` are the slack variables, `Z` is the diagonal matrix of
# Lagrange multipliers, `S` is the diagonal matrix of slack variables,`e` is a
# vector of ones, and `nu` and `mu` are penalty parameters.
#
# Note that `pred(p)` is simply equal to `nu (m(0) - m(p)) - h(p_x, p_s)`.
#
# @param nu penalty parameter.
# @param m0 norm of `C(x) - s`.
# @param mp norm of `A p_x - p_s + C(x) - s`.
# @param h value of the function `h(p_x, p_s)` (see details).
#
# @return Value of the predicted reduction `pred(p)` for the proposed step `p`.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_pred <- function(nu, m0, mp, h) {
  nu * (m0 - mp) - h
}

# Evaluate the proposed step
#
# Compute the actual reduction, predicted reduction, and their ratio.
#
# @details
# The actual reduction is defined as the difference of merit functions
# `phi(x, s) - phi(x + p_x, s + p_s)` where
# `phi(x, s) = f(x) - mu sum(log(s)) + nu |C(x) - s|`.
# The predicted reduction is instead simply equal to
# `nu (m(0) - m(p)) - h(p_x, p_s)`.
#
# @param obj list containing all relevant variables of the current best
#   solution.
# @param candidate list containing all relevant variables of the proposed new
#   solution.
# @param mu barrier penalty parameter.
# @param nu constraints penalty parameter.
#
# @return A list containing the updated parameter `nu`, the actual reduction,
#   the predicted reduction, and their ratio.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_reductions <- function(obj, candidate, mu, nu) {
  nu <- ntrm_update_nu(nu, candidate$h, obj$m0, candidate$mp)

  ared <- ntrm_ared(
    nu,
    obj$f, candidate$f,
    obj$s, candidate$s,
    obj$m0, candidate$m0,
    mu
  )

  pred <- ntrm_pred(
    nu, obj$m0, candidate$mp, candidate$h
  )

  rho <- ared / pred

  if (is.nan(rho)) {
    # both reductions are zero and we reject this step
    rho <- 0
  }

  list(
    nu = nu,
    ared = ared,
    pred = pred,
    rho = rho
  )
}

# Step computation
#
# Compute the step of the interior point trust-region method
#
# @details
# This is the core function of the interior point trust-region method. We seek
# the vector `p = (p_x, p_s)` minimizing
# `h(p_x, p_s) = G' p_x + 0.5 p_x' H p_x - mu e' p_s + 0.5 p_s' Z S p_s`
# subject to the constraints `A p_x - S p_s + (C(x) - s) = k`, `|p| <= delta`,
# and `min(p_s) >= -tau`.
#
# `G` is the gradient of the original function, `H` is the Hessian matrix of
# the Lagrangian function, `mu` is the barrier penalty parameter, `Z` is the
# diagonal matrix of Lagrange multipliers, `S` is the diagonal matrix of slack
# variables `s`, `C(.)` are functions evaluating the inequality constraints
# (in our case they are the distances from the boundary), `A` is the
# Jacobian matrix of `C(.)`, and `e` is a vector of ones.
#
# Our objective function is made of two separate quadratic programs. The first
# part is associated with the original variables:
#
# QP1(p_x) = 0.5 p_x' H p_x + G' p_x
#
# The second part is associated with the slack variables:
#
# QP2(p_s) = 0.5 p_s' V p_s - mu' p_s
#
# where `V = Z S`.
#
# Set `w(x) = C(x) - s`, `p = c(p_x, p_s)`,
# `Q = rbind(cbind(H, 0), cbind(0, V))`, `B = cbind(A, -S)`, `u = c(G, -mu e)`,
# `y(x) = k - w(x)`. The problem becomes:
#
# min_{p} 0.5 p' Q p + p' u
# subject to B p = y(x)
#
# which we can solve with the projected conjugate gradient method.
#
# @param obj list containing all relevant variables.
# @param delta radius of the trust region.
# @param tau fraction to the boundary parameter.
#
# @return List with the step of the interior point trust-region method,
#   satisfying the boundary constraints.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_step <- function(obj, delta, tau) {
  # Section 19.5, 'Step computation' at page 580
  v <- ntrm_step_normal(obj$Y, obj$B, obj$r, delta, tau, obj$idx_s)

  # Equation (19.41+b) at page 582
  mp <- ntrm_norm(obj$B %*% v + obj$r)

  # apply correction if Hessian matrix is (almost) singular
  kappa <- ntrm_inertia_correction(obj$H, obj$A, obj$s, obj$z, obj$kappa)

  # See Algorithm 16.2 at page 461 of Nocedal and Wright (2006)
  V <- diag(obj$z * obj$s, nrow = obj$m)

  Q <- matrix(0, nrow = obj$n, ncol = obj$n)
  Q[seq_len(obj$k), seq_len(obj$k)] <- obj$H + kappa * diag(obj$k)
  Q[(obj$k + 1):obj$n, (obj$k + 1):obj$n] <- V

  # Algorithm 16.2 at page 461
  step <- ntrm_step_tilde(Q, obj$u, obj$B, v, delta, tau, obj$idx_s)

  # Equation (19.32) at page 580
  p <- step$p
  p[obj$idx_s] <- obj$s * p[obj$idx_s]

  list(
    p_x = p[obj$idx_x],
    p_s = p[obj$idx_s],
    p_norm = ntrm_norm(p),
    h = step$h,
    mp = mp,
    kappa = kappa
  )
}

# Normal step computation
#
# Compute the normal step of the interior point trust-region method.
#
# @details
# We want to find the minimum of `(B v + w)'(B v + w)` where `B = cbind(A, -S)`
# and `w = C(x) - s`. `C(.)` are functions evaluating the inequality
# constraints (in our case they are the distances from the boundary), `A` is
# the Jacobian matrix of `C(.)`, `s` are the slack variables, and `S` is the
# diagonal matrix of slack variables.
#
# The solution must also satisfy the constraints `|v| <= 0.8 delta` and
# `min(v_s) >= -tau / 2`, where `v_s` are the values of `v` corresponding to
# the slack variables.
#
# The solution of the problem is computed with the dogleg method.
#
# @param Y numeric matrix used to solve the linear system.
# @param B numeric matrix defined as `B = cbind(A, -S)`.
# @param slack_residual numeric vector of distances from the boundary.
# @param delta radius of the trust region.
# @param tau fraction to the boundary parameter.
# @param idx_s integer vector indicating the position of the slack variables in
#   the `initial_value` starting vector.
#
# @return Numeric vector approximately minimizing `|B v + C(x) - s|^2` subject
#   to the constraints.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_step_normal <- function(Y, B, slack_residual, delta, tau, idx_s) {
  # See Section 4.1 (pages 73-76), Section 18.5 (pages 547-548), and
  # Section 19.5 (page 580) of Nocedal and Wright (2006)
  rho <- 0.8 * delta
  omega <- - tau / 2

  # f(v) = 0.5 v' B' B v + v' B' w + 0.5 w' w
  #
  # The gradient is g = B' w + B' B v and the steepest descent direction is -g.
  # Start at a step of length 0, i.e. v = 0, to obtain g = B' w.
  #
  # The unbounded Gauss-Newton step (with the smallest norm) is equal to
  #
  # v_gn = -B' (B B')^(-1) w
  #
  # while the steepest descent step is
  #
  # v_sd = -(G' G / (G' B' B G)) G
  #
  # Gauss-Newton step v_gn can be obtained by our Y matrix:
  #
  # Y = | I   B' | | v_gn | = |  0 |
  #     | B   0  | |    x |   | -w |
  n <- ncol(B)
  b <- c(rep(0, n), -slack_residual)
  v_gn <- solve(Y, b)[seq_len(n)]
  v_gn_norm <- ntrm_norm(v_gn)

  v <- if (v_gn_norm <= rho) {
    # Gauss-Newton step is within the boundary -> accept it
    v_gn
  } else {
    # Steepest Descent step
    g <- as.numeric(crossprod(B, slack_residual))
    v_sd <- -(sum(g^2) / sum((B %*% g)^2)) * g
    v_sd_norm <- ntrm_norm(v_sd)

    if (v_sd_norm >= rho) {
      # both Gauss-Newton step and Steepest Descent step are outside the
      # boundary. We return the Steepest Descent shrunk to the boundary
      (rho / v_sd_norm) * v_sd
    } else {
      # dogleg method
      #
      # solve | theta * v_gn + (1 - theta) * v_sd |^2 = delta^2 for theta, that
      # is
      #
      # theta^2 x'x + 2 * theta * x'v_sd + v_sd'v_sd - delta^2 = 0
      #
      # where x = v_gn - v_sd
      x <- v_gn - v_sd
      a <- ntrm_solve_quadratic_equation(
        sum(x^2),
        2 * sum(x * v_sd),
        v_sd_norm^2 - delta^2
      )

      a[2] * v_gn + (1 - a[2]) * v_sd
    }
  }

  small_v <- min(v[idx_s])
  if (small_v < omega) {
    v <- (omega / small_v) * v
  }

  v
}

# Shrink a candidate step
#
# Shrink a candidate step so that `min(p_s) >= -tau`.
#
# @param p previous step in the iteration.
# @param a magnitude of the new proposed direction move.
# @param d direction towards the new step solution.
# @param tau fraction to the boundary parameter.
# @param idx_s integer vector indicating the position of the slack variables in
#   the vector `p`.
#
# @return New step obtained by starting at `p` and moving along direction `d`
#   such that `min(p_new[idx_s]) >= -tau`.
ntrm_step_shrink <- function(p, a, d, tau, idx_s) {
  p_new <- p + a * d

  i <- which.min(p_new[idx_s])
  small_p <- p_new[idx_s][i]

  if (small_p >= -tau) {
    # no need to shrink the solution
    a
  } else {
    small_d <- d[idx_s][i]

    a_shrink <- - (tau + p[idx_s][i]) / small_d

    if (abs(a_shrink) > abs(a)) {
      # the direction is not shrunk but enlarged
      0
    } else {
      a_shrink
    }
  }
}

# Step computation
#
# Apply the projected gradient method to find the candidate step of the
# interior point trust-region method.
#
# @details
# The problem we want to solve is
#
# min_{p} 0.5 p' Q p + p' u
# subject to B p = y(x)
#
# If the initial value satisfies the linear constraints, the projected gradient
# method is guaranteed to automatically iterate among solutions that also
# satisfy the constraints.
#
# @param Q numeric matrix representing the quadratic part.
# @param u numeric vector representing the linear part.
# @param B numeric matrix representing the linear constraints.
# @param initial_value initial vector satisfying the equality constraints.
# @param delta radius of the trust region.
# @param tau fraction to the boundary parameter.
# @param idx_s integer vector indicating the position of the slack variables in
#   the `initial_value` starting vector.
#
# @return List containing the step of the interior point trust-region method
#   satisfying the boundary constraints, and the minimum value of the objective
#   function.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
ntrm_step_tilde <- function(Q, u, B, initial_value, delta, tau, idx_s) {
  eps <- sqrt(.Machine$double.eps)

  m <- nrow(B)
  n <- ncol(B)

  # The projection, equation (16.28d), is done by the augmented system (16.34)
  #
  # K = | Q B' | | g | = | r |
  #     | B 0  | | y |   | 0 |
  K <- rbind(
    cbind(Q, t(B)),
    cbind(B, matrix(0, nrow = m, ncol = m))
  )

  p <- p_new <- p_old <- initial_value

  h <- ntrm_objective_function(p, Q, u)
  r <- Q %*% p_old + u

  # Equation (16.34) at page 463
  b <- c(r, rep(0, m))
  g <- solve(K, b)[seq_len(n)]
  dot <- sum(r * g)
  d <- -g

  repeat {
    W <- Q %*% d
    dqf <- sum(d * W)
    a <- dot / dqf

    p_new <- p_old + a * d

    # new proposed step is `p + a d` but `a` is computed according to the
    # unrestricted problem. We want instead that | p + a d | <= delta and
    # min(p_s) >= -tau
    #
    # the best coefficient satisfying (p + a d)'(p + a d) <= delta^2 is within
    # the interval delimited by the solutions of
    # d'd a^2 + 2 d'p a + p'p - delta^2  = 0.
    #
    # the function to minimize is h(a) = 0.5 d' Q d a^2 + d' (Q p + u) a, which
    # is a quadratic equation. Its minimum is obviously the one we just
    # computed.
    #
    # If the global minimum is outside the feasible set, being a quadratic
    # equation, the best coefficient for a feasible solution must be on the
    # boundary of the interval.
    if (dqf > 0) {
      p_new_norm <- ntrm_norm(p_new)

      if (p_new_norm > delta) {
        k_sol <- ntrm_solve_quadratic_equation(
          sum(d^2),
          2 * sum(d * p_old),
          sum(p_old^2) - delta^2
        )

        k <- if (a > 0) {
          k_sol[2]
        } else {
          k_sol[1]
        }

        # check we can shrink this solution for a better step?
        a_shrink <- ntrm_step_shrink(p_old, k, d, tau, idx_s)

        if (a_shrink != 0) {
          p_tmp <- p_old + a_shrink * d
          hs <- ntrm_objective_function(p_tmp, Q, u)

          if (hs < h) {
            p <- p_tmp
            h <- hs
          }
        }

        # stop algorithm because the proposed step went outside the boundary
        break
      }
    } else {
      # current search direction is a direction of non-positive curvature
      # we scan the whole interval, update the best solution and return
      f <- function(x) {
        x * (0.5 * dqf * x + sum(d * (Q %*% p_old + u)))
      }

      k <- ntrm_solve_quadratic_equation(
        sum(d^2),
        2 * sum(d * p_old),
        sum(p_old^2) - delta^2
      )

      a <- ntrm_line_search(f, k[1], k[2])
      a <- ntrm_step_shrink(p_old, a, d, tau, idx_s)

      if (a != 0) {
        p_tmp <- p_old + a * d
        hs <- ntrm_objective_function(p_tmp, Q, u)

        if (hs < h) {
          p <- p_tmp
          h <- hs
        }
      }

      break
    }

    # norm of p_new is less then delta, so coefficient `a` is surely in
    # [k[1], k[2]], but is min(p_s) >= -tau?
    a_shrink <- ntrm_step_shrink(p_old, a, d, tau, idx_s)

    if (a_shrink == a) {
      p <- p_new
      h <- ntrm_objective_function(p, Q, u)
    } else if (a_shrink != 0) {
      p_tmp <- p_old + a_shrink * d
      hs <- ntrm_objective_function(p_tmp, Q, u)

      if (hs < h) {
        p <- p_tmp
        h <- hs
      }
    }

    converged <- (ntrm_max_abs(p_new - p_old) / ntrm_max_abs(p_new)) < eps

    if (converged) {
      break
    }

    r <- r + a * W
    b <- c(r, rep(0, m))
    g <- solve(K, b)[seq_len(n)]

    dot_old <- dot
    dot <- sum(r * g)

    if (abs(dot) < eps) {
      break
    }

    p_old <- p_new

    d <- -g + (dot / dot_old) * d
  }

  list(
    p = p,
    h = h
  )
}

# Update penalty parameter
#
# Update the penalty parameter of the merit function.
#
# @param nu current value of the penalty parameter.
# @param h value of the interior point trust-region objective function.
# @param m0 norm of `C(x) - s`.
# @param mp norm of `A p_x - p_s + C(x) - s`.
#
# @return Updated penalty parameter such that `pred(d) >= 0.3 nu (m0 - mp)`.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
#
# Richard H. Byrd, Mary E. Hribar, and Jorge Nocedal. An interior point
# algorithm for large scale nonlinear programming.
# **SIAM Journal on Optimization**, 9(4):877–900, 1999.
# doi: 10.1137/S1052623497325107.
ntrm_update_nu <- function(nu, h, m0, mp) {
  # We use Equation (19.41+) at page 582 to get an inequality similar to
  # Equation (18.36) at page 542 of Nocedal and Wright (2006)
  #
  # We choose rho = 0.3 as in NITRO (Byrd et al. 1998)
  # 1 / (1 - rho) = 1 / (1 - 0.3) = 1 / 0.7
  denominator <- (0.7 * (m0 - mp))

  if (abs(denominator) <= .Machine$double.eps) {
    # basically m0 and mp are equal
    nu
  } else {
    nu_lb <- h / denominator

    if (is.nan(nu_lb) || is.infinite(nu_lb) || (nu >= nu_lb)) {
      # old value satisfies the inequality or the new lower bound is invalid
      nu
    } else {
      nu_lb + 0.5
    }
  }
}

# Update the solution
#
# Update the old solution with a better one.
#
# @param obj list containing all relevant variables of the current best
#   solution.
# @param candidate list containing all relevant variables of the proposed new
#   solution.
# @param gh function handle to evaluate both the gradient and Hessian of the
#   function to minimize.
# @param mu barrier penalty parameter.
ntrm_update_solution <- function(obj, candidate, gh, mu) {
  obj$x <- candidate$x
  obj$f <- candidate$f
  obj$s <- candidate$s
  obj$r <- candidate$r
  obj$m0 <- candidate$m0

  gradient_hessian <- gh(obj$x)
  obj$G <- gradient_hessian$G
  obj$H <- gradient_hessian$H

  obj$B <- cbind(obj$A, -diag(obj$s, nrow = obj$m))
  obj$Y <- rbind(
    cbind(diag(obj$n), t(obj$B)),
    cbind(obj$B, matrix(0, nrow = obj$m, ncol = obj$m))
  )
  obj$u <- c(obj$G, -rep(mu, obj$m))

  obj$z <- ntrm_lagrange_multiplier(obj, mu)

  obj
}

# Newton Trust-Region Interior-Point Method
#
# Find the minimum of a function using the Newton method, combined with the
# Trust-Region Interior-Point method.
#
# @param fn function handle to evaluate the function to minimize.
# @param gh function handle to evaluate both the gradient and Hessian of the
#   function to minimize.
# @param init numeric vector of starting values.
# @param max_iter integer value for the maximum number of iterations.
# @param lower_bound numeric vector of lower bounds for each variable.
# @param upper_bound numeric vector of upper bounds for each variable.
#
# @return Numeric vector `x` for which `fn(x)` is the minimum within the
#   boundaries specified by `lower_bound` and `upper_bound`.
#
# @references
# Jorge Nocedal and Stephen J Wright. **Numerical optimization**. Springer,
# New York, NY, USA, second edition, 2006. ISBN 978-0-387-30303-1.
#
# Richard H. Byrd, Mary E. Hribar, and Jorge Nocedal. An interior point
# algorithm for large scale nonlinear programming.
# **SIAM Journal on Optimization**, 9(4):877–900, 1999.
# doi: 10.1137/S1052623497325107.
#
# Todd Dennis Plantenga.
# **Large-scale nonlinear constrained optimization using trust regions**.
# Doctoral dissertation, Northwestern University, 1994.
ntrm_constrained <- function(fn, gh, init, max_iter, lower_bound, upper_bound) {
  # ---------
  # For details see Chapter 19 of Nocedal and Wright (2006).
  #
  # Constraints on the parameters are of the form
  # l_1 <= x_1 <= u_1
  # ...
  # l_k <= x_k <= u_k
  # where the lower bounds can be -Inf and the upper bounds can be Inf.
  #
  # In this particular application we don't have equality constraints.
  #
  # Inequalities are transformed into equalities by the introduction of slack
  # variables s = (s_l, s_u)^{T} such that
  #
  # a_1(x_1) - s_l_1 = x_1 - l_1 - s_l_1 = 0
  # ...
  # a_k(x_k) - s_l_k = x_k - l_k - s_l_k = 0
  #
  # b_1(x_1) - s_u_1 = u_1 - x_1 - s_u_1 = 0
  # ...
  # b_k(x_k) - s_u_k = u_k - x_k - s_u_k = 0
  #
  # Define the function C(x) = (a(x), b(x))^{T}.
  #
  # Define the Langrangian function as
  #
  # L(x, s, z) = fn(x) - z' (C(x) - s)
  #
  # where z are the Lagrangian multipliers (at most of length 2k).
  #
  # Define now the Jacobian matrix A of C(x)
  # (`|` represent the matrix delimiter, not the determinant):
  #
  # |  1,  0, ...,  0,  0|
  # |  0,  1, ...,  0,  0|
  # ...
  # |  0,  0, ...,  1,  0|
  # |  0,  0, ...,  0,  1|
  # | -1,  0, ...,  0,  0|
  # |  0, -1, ...,  0,  0|
  # ...
  # |  0,  0, ..., -1,  0|
  # |  0,  0, ...,  0, -1|
  #
  # A is a (n_constraints)-by-(n_variables) rectangular matrix.
  #
  # The primal-dual system is defined as in equation (19.12):
  #
  # | Hessian(L),  0,  A' | |  p_x | = - | Gradient(x) - A' z |
  # |          0,  V, -I  | |  p_s |     |    z - mu S^{-1} e |
  # |          A, -I,  0  | | -p_z |     |           C(x) - s |
  #
  # where e = (1, 1, ..., 1)', Z and S are diagonal matrices made from z and s
  # respectively, and V = S^{-1} Z.
  #
  # Parameter mu is updated during the iterations to prevent z and s to approach
  # zero too fast.
  #
  # Starting from an initial feasible solution l <= x <= u and slack variables
  # s >= 0, the algorithm will do basically the following:
  # - Find a new value of the Lagrange multipliers z.
  # - Find step p = (p_x, p_s)' of at most length delta such that fn(x + p) is
  #   sufficiently decreased under the constraints.
  # - Update parameter mu if necessary.
  # - Repeat.
  # -------
  eps_1 <- 1.0e-10
  eps_2 <- 1.0e-8
  eps_counter <- 0

  delta <- 1
  mu <- 0.1
  eta <- 1.0e-8
  tau <- 0.995
  nu <- 0.5

  converged <- FALSE

  # setup constraints
  constraints <- ntrm_create_constraints(lower_bound, upper_bound)

  il <- constraints$idx_lb
  iu <- constraints$idx_ub

  lb <- lower_bound[il]
  ub <- upper_bound[iu]

  C <- function(x) {
    c(x[il] - lb, ub - x[iu])
  }

  cur_optimum <- init

  # initial value must be within the constraint region
  is_out <- (init < lower_bound) | (init > upper_bound)
  for (i in seq_len(ncol(constraints$A))) {
    if (is_out[i]) {
      if (is.infinite(lower_bound[i])) {
        # value is greater than the upper bound
        cur_optimum[i] <- upper_bound[i] - 1
      } else if (is.infinite(upper_bound[i])) {
        # value is smaller than the lower bound
        cur_optimum[i] <- lower_bound[i] + 1
      } else {
        # value is outside the finite interval
        cur_optimum[i] <- (lower_bound[i] + upper_bound[i]) / 2
      }
    }
  }

  # the algorithm is very sensitive to the initial slack variables
  #
  # we try to find a value for which the error is already small
  s <- tryCatch(
    {
      ntrm_line_search(
        function(x) {
          slacks <- rep(x, nrow(constraints$A))
          ntrm_error(
            ntrm_init_obj(fn, gh, C, cur_optimum, slacks, constraints$A, mu),
            mu
          )
        },
        0.01,
        100
      )
    },
    error = function(e) {
      1
    }
  )

  # if the starting point is way too far from the true solution, the previous
  # strategy might not work
  obj <- ntrm_init_obj(
    fn, gh, C, cur_optimum, rep(s, nrow(constraints$A)), constraints$A, mu
  )

  error <- ntrm_error(obj, mu)
  if (error > 100) {
    obj <- ntrm_init_obj(
      fn, gh, C, cur_optimum, rep(1, nrow(constraints$A)), constraints$A, mu
    )
  }

  i <- 0
  repeat {
    error <- ntrm_error(obj, 0)

    if (error < eps_1) {
      converged <- TRUE
      break
    } else if (error < eps_2) {
      # the algorithm might be stuck, but the error is small enough for
      # declaring it converged
      if (eps_counter >= 100) {
        converged <- TRUE
        break
      } else {
        eps_counter <- eps_counter + 1
      }
    }

    if (i == max_iter) {
      break
    }

    obj$kappa <- 0

    repeat {
      error <- ntrm_error(obj, mu)

      if (error < mu || (delta < eps_1 && error < eps_2) || i == max_iter) {
        break
      }

      candidate <- ntrm_step(obj, delta, tau)
      candidate$x <- obj$x + candidate$p_x
      candidate$f <- fn(candidate$x)
      candidate$s <- obj$s + candidate$p_s
      candidate$r <- C(candidate$x) - candidate$s
      candidate$m0 <- ntrm_norm(candidate$r)

      # penalty parameter might have been updated
      obj$kappa <- candidate$kappa

      acceptance <- ntrm_reductions(obj, candidate, mu, nu)
      nu <- acceptance$nu

      if (acceptance$rho >= eta) {
        obj <- ntrm_update_solution(obj, candidate, gh, mu)

        # Equation (3.55) of Byrd et al. (1999)
        if (acceptance$rho >= 0.9) {
          delta <- max(delta, 7 * candidate$p_norm)
        } else if (acceptance$rho >= 0.3) {
          delta <- max(delta, 2 * candidate$p_norm)
        }
      } else {
        delta <- 0.25 * delta
      }

      i <- i + 1
    }

    delta <- 1
    nu <- 0.5
    mu <- 0.2 * mu

    # since we updated the penalty parameter, recompute the associated variables
    obj$u <- c(obj$G, -rep(mu, obj$m))
    obj$z <- ntrm_lagrange_multiplier(obj, mu)
  }

  list(
    optimum = obj$x,
    minimum = obj$f,
    converged = converged,
    iterations = i
  )
}
