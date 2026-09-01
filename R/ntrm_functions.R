# Hessian matrix correction
#
# Correct the Hessian matrix to be positive definite
#
# @param H Hessian matrix to correct.
#
# @return Matrix `H` if it is already positive definite, otherwise a new matrix
#  `H + w I` that is positive definite.
ntrm_correct_hessian <- function(H) {
  z <- diag(H)
  h <- min(z)
  b <- 1.0e-3
  w <- 0

  if (h <= 0) {
    w <- b - h
  }

  diag(H) <- z + w
  H_chol <- tryCatch(chol(H), error = function(e) NULL)

  i <- 0
  while (is.null(H_chol) && (i < 100)) {
    w <- max(2 * w, b)

    diag(H) <- z + w
    H_chol <- tryCatch(chol(H), error = function(e) NULL)

    i <- i + 1
  }

  if (is.null(H_chol)) {
    stop("Hessian matrix is ill-conditioned", call. = FALSE)
  }

  H
}

# Maximum absolute value
#
# Find the maximum absolute value in a vector, i.e. `max(abs(x))`, calling
# `abs` only once.
#
# @param x numeric vector.
#
# @return Maximum absolute value in `x`.
ntrm_max_abs <- function(x) {
  min_max <- range(x)
  max(abs(min_max[1]), min_max[2])
}

# Euclidean norm of a vector
#
# Compute the Euclidean norm of a vector `x`, i.e. `|x| = sqrt(sum(x^2))`.
#
# @param x numeric vector.
#
# @return Euclidean norm of `x`.
ntrm_norm <- function(x) {
  sqrt(sum(x^2))
}

# Solutions of a quadratic equation
#
# Solve the equation `a x^2 + b x + c = 0` for x, in a numerically stable
# manner. However, the discriminant is still evaluated using the actual
# formula. This solver is not a general quadratic equation solver, but only
# for our specific case of NTRM line-sphere intersection.
#
# @param a quadratic coefficient.
# @param b linear coefficient.
# @param c constant.
#
# @return Numeric vector of length 2 with the solutions of the quadratic
#   equation.
ntrm_solve_quadratic_equation <- function(a, b, c) {
  # ~100 ulps
  u <- 100 * .Machine$double.eps

  # Check if it is a degenerate quadratic equation, i.e. at most linear.
  #
  # When a is very close to zero, the equation is basically `b x + c = 0` with
  # solution `x = -c / b`, unless also b is very close to zero.
  #
  # There is also a problem of scale. `a` and `b` can be very small or very
  # large, therefore we need to use a relative tolerance. If `b` is very large,
  # we use a relative tolerance based on its value. If `b` is very small, we
  # use an absolute tolerance based on the machine epsilon (this is the reason
  # why we put a 1 in the `max` function).
  if (abs(a) <= u * max(abs(b), 1)) {
    # `a` is close to zero, so we check if `b` is also close to zero.
    if (abs(b) <= u) {
      # `b` is close to zero, so we check if `c` is also close to zero.
      if (abs(c) <= u) {
        # all coefficients are basically zero, so the equation is effectively
        # 0 = 0. Since there is no x, we return two zeros by convention.
        return(c(0, 0))
      }

      # `a` is close to zero (compared to `b`), `b` is close to zero, but `c`
      # is not. The equation is effectively `c = 0`, so no solution is possible.
      stop(
        "no finite solutions while solving quadratic equation",
        call. = FALSE
      )
    }

    # `a` is close to zero (compared to `b`) but `b` is not, so the equation is
    # effectively `b x + c = 0` with solution `x = -c / b`.
    return(rep(-c / b, 2))
  }

  # In general, we should check for a negative discriminant. However, in our
  # application, we know that the discriminant is always non-negative. Negative
  # discriminants are likely numerical errors, so we set them to zero.
  discriminant <- max(b^2 - 4 * a * c, 0)

  if (discriminant == 0) {
    return(rep(-b / (2 * a), 2))
  }

  s <- 1
  if (b < 0) {
    s <- -1
  }

  t <- b + s * sqrt(discriminant)
  x <- -c(t / (2 * a), 2 * c / t)

  if (x[1] <= x[2]) {
    x
  } else {
    c(x[2], x[1])
  }
}

# Find the minimum of a function by line search
#
# Using the Brent's method, find the value of `x` that minimize a function
# `f(x)` in the interval `[a, b]`.
#
# @details
# This is a direct port of the algorithm from section 10.3 at page 498 of Press
# et al. (2007).
#
# @param f one-dimensional function to minimize.
# @param a lower bound.
# @param b upper bound.
#
# @return Scalar value `a <= x <= b` such that `f(x)` is minimum.
#
# @references
# William H. Press et al. **Numerical recipes**. Cambridge University Press,
# Cambridge, UK, third edition, 2007. ISBN 978-0-511-33555-6.
ntrm_line_search <- function(f, a, b) {
  # 1 - golden_ratio
  k <- 0.381966

  tol_0 <- 1.0e-10
  eps_zero <- 1.0e-10

  d <- 0
  e <- 0

  fa <- f(a)
  fb <- f(b)

  y <- a
  fy <- fa

  if (fb < fa) {
    y <- b
    fy <- fb
  }

  x <- w <- v <- b
  fx <- fw <- fv <- fb

  for (i in seq_len(101)) {
    xm <- (a + b) / 2

    tol_1 <- tol_0 * abs(x) + eps_zero
    tol_2 <- 2 * tol_1

    if (abs(x - xm) <= tol_2) {
      if (fx < fy) {
        y <- x
        fy <- fx
      }
      break
    }

    if (abs(e) > tol_1) {
      r <- (x - w) * (fx - fv)
      q <- (x - v) * (fx - fw)
      p <- (x - v) * q - (x - w) * r
      q <- 2 * (q - r)

      if (q > 0) {
        p <- -p
      }

      q <- abs(q)
      e_tmp <- e
      e <- d

      if (
        (abs(p) >= abs(0.5 * q * e_tmp)) ||
          (p <= q * (a - x)) ||
          (p >= q * (b - x))
      ) {
        e <- if (x >= xm) {
          a - x
        } else {
          b - x
        }

        d <- k * e
      } else {
        d <- p / q
        u <- x + d

        if (((u - a) < tol_2) || ((b - u) < tol_2)) {
          d <- sign(xm - x) * tol_1
        }
      }
    } else {
      e <- if (x >= xm) {
        a - x
      } else {
        b - x
      }

      d <- k * e
    }

    u <- if (abs(d) >= tol_1) {
      x + d
    } else {
      x + sign(d) * tol_1
    }

    fu <- f(u)

    if (fu <= fx) {
      if (u >= x) {
        a <- x
      } else {
        b <- x
      }

      v <- w
      w <- x
      x <- u

      fv <- fw
      fw <- fx
      fx <- fu
    } else {
      if (u < x) {
        a <- u
      } else {
        b <- u
      }

      if ((fu <= fw) || (w == x)) {
        v <- w
        w <- u

        fv <- fw
        fw <- fu
      } else if ((fu <= fv) || (v == x) || (v == w)) {
        v <- u
        fv <- fu
      }
    }
  }

  if (i == 101) {
    stop("line search not converged", call. = FALSE)
  }

  y
}
