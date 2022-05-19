#' @importFrom stats integrate
#'
#' @export
nauc.drda <- function(object, xlim = NULL, ylim = c(0, 1)) {
  if (is.null(xlim)) {
    xlim <- if (inherits(object, "loglogistic")) {
      c(0, 1000)
    } else {
      c(-10, 10)
    }
  } else {
    if (length(xlim) != 2) {
      stop("'xlim' must be of length 2", call. = FALSE)
    }

    if (!is.numeric(xlim)) {
      stop("'xlim' must be a numeric vector of length 2", call. = FALSE)
    }

    if (xlim[1] >= xlim[2]) {
      stop("'xlim[1]' cannot be larger or equal to 'xlim[2]'", call. = FALSE)
    }

    if (inherits(object, "loglogistic")) {
      # predictor cannot be negative, so we cannot integrate for x < 0
      if (xlim[1] < 0) {
        xlim[1] <- 0
      }

      if (xlim[2] < 0) {
        return(0.0)
      }
    }
  }

  if (length(ylim) != 2) {
    stop("'ylim' must be of length 2", call. = FALSE)
  }

  if (!is.numeric(ylim)) {
    stop("'ylim' must be a numeric vector of length 2", call. = FALSE)
  }

  if (ylim[1] >= ylim[2]) {
    stop("'ylim[1]' cannot be larger or equal to 'ylim[2]'", call. = FALSE)
  }

  k <- 1
  if (
    inherits(object, "loglogistic6_fit") || inherits(object, "logistic6_fit")
  ) {
    k <- 1 / object$coefficients[6]^(1 / object$coefficients[5])
  }

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]

  # in case the curve intersect `ylim`, these are the values at which it happens
  tmp <- inverse_fn(object, ylim)

  # value of the integral
  I <- 0

  # check if we really need to perform an integration
  flag <- TRUE

  # we might change the range of integration
  xlim_new <- xlim

  if (delta >= 0) {
    # curve is monotonically increasing
    lb <- alpha
    ub <- alpha + delta * k

    if (lb < ylim[1]) {
      # the curve in `c(0, tmp[1])` is to be considered zero
      if (tmp[1] > xlim[2]) {
        # the integral is simply zero
        flag <- FALSE
      } else if (tmp[1] > xlim[1]) {
        xlim_new[1] <- tmp[1]
      }
    }

    if (ub > ylim[2]) {
      # the curve in `c(tmp[2], Inf)` is equal to `ylim[2]`
      if (tmp[2] < xlim[1]) {
        # not much to do in this case, the curve in the requested range is flat
        I <- I + (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

        # stop here
        flag <- FALSE
      } else if (tmp[2] < xlim[2]) {
        # standard before `tmp[2]` and flat in c(tmp[2], xlim[2])
        I <- I + (xlim[2] - tmp[2]) * (ylim[2] - ylim[1])

        # modify the range of integration
        xlim_new[2] <- tmp[2]
      }
    }
  } else {
    # curve is monotonically decreasing
    lb <- alpha + delta * k
    ub <- alpha

    if (ub > ylim[2]) {
      # the first part of the curve in `c(0, tmp[2])` is equal to `ylim[2]`
      if (tmp[2] > xlim[2]) {
        # not much to do in this case, the curve in the requested range is flat
        I <- I + (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])

        # stop here
        flag <- FALSE
      } else if (tmp[2] > xlim[1]) {
        # flat in c(xlim[1], tmp[2]) and standard after `tmp[2]`
        I <- I + (tmp[2] - xlim[1]) * (ylim[2] - ylim[1])

        # modify the range of integration
        xlim_new[1] <- tmp[2]
      }
    }

    if (lb < ylim[1]) {
      # the curve after `tmp[1]` is to be considered zero
      if (tmp[1] < xlim[1]) {
        # the integral is simply zero
        flag <- FALSE
      } else if (tmp[1] < xlim[2]) {
        xlim_new[2] <- tmp[1]
      }
    }
  }

  if (flag) {
    # we remove `ylim[1]` to shift the curve to zero
    f <- function(x) {
      fn(object, x, object$coefficients) - ylim[1]
    }

    I <- I + integrate(
      f, lower = xlim_new[1], upper = xlim_new[2],
      rel.tol = sqrt(.Machine$double.eps)
    )$value
  }

  nauc <- I / ((xlim[2] - xlim[1]) * (ylim[2] - ylim[1]))
  names(nauc) <- NULL

  nauc
}
