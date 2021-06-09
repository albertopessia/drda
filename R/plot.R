#' Model fit plotting
#'
#' Plot maximum likelihood curves fitted with `drda`.
#'
#' @param x `drda` object as returned by the \code{link[drda]{drda}} function.
#' @param ... other `drda` objects or parameters to be passed to the plotting
#'   functions. See 'Details'.
#'
#' @details
#'
#' This function provides a scatter plot of the observed data, overlaid with the
#' maximum likelihood curve fit.
#' If multiple fit objects are given, they will all be placed in the same plot.
#'
#' Accepted plotting arguments are:
#'
#' \describe{
#'   \item{base}{character string with the base used to print the values on the
#'     x axis. Accepted values are `e` for the natural logarithm (the default),
#'     `10` for base 10, and `2` for base 2.}
#'   \item{col}{curve color(s). By default, up to 9 color-blind friendly
#'     colors are provided.}
#'   \item{xlab, ylab}{axis labels.}
#'   \item{xlim, ylim}{the range of x and y values with sensible defaults.}
#'   \item{level}{level of confidence intervals. Set to zero or a negative value
#'     to disable confidence intervals.}
#'   \item{legend}{custom labels for the legend model names.}
#' }
#'
#' @return No return value.
#'
#' @importFrom graphics axis box curve legend lines par polygon points
#' @importFrom graphics plot.default
#' @importFrom grDevices col2rgb dev.flush dev.hold extendrange rgb
#' @importFrom stats qchisq
#'
#' @export
plot.drda <- function(x, ...) {
  dotargs <- list(...)

  named <- if (is.null(names(dotargs))) {
    rep_len(FALSE, length(dotargs))
  } else {
    names(dotargs) != ""
  }

  fit_objects <- dotargs[!named]

  is_drda <- vapply(fit_objects, function(x) inherits(x, "drda"), FALSE)
  fit_objects <- fit_objects[is_drda]

  if (length(fit_objects)) {
    return(plot.drdalist(c(list(x), fit_objects), ...))
  }

  col <- dotargs$col
  if (is.null(col)) {
    col <- "#EE6677FF"
  }

  base <- dotargs$base
  k <- 1
  if (!is.null(base)) {
    if (base == "10") {
      k <- log(10)
    } else if (base == "2") {
      k <- log(2)
    } else {
      # if it's not 2 nor 10, we set by default "e"
      base <- "e"
    }
  } else {
    if (x$is_log) {
      base <- "e"
    } else {
      # by default we use a log10 scale
      base <- "10"
      k <- log(10)
    }
  }

  xlab <- dotargs$xlab
  if (is.null(xlab)) {
    xlab <- if (base == "e") {
      "log(Predictor)"
    } else {
      "Predictor"
    }
  }

  ylab <- dotargs$ylab
  if (is.null(ylab)) {
    ylab <- "Response"
  }

  theta <- x$coefficients

  alpha <- 0
  beta <- 1
  eta <- -1
  phi <- 0

  if (x$mean_function != "logistic2") {
    alpha <- theta[1]
    beta <- theta[2]
    eta <- theta[3]
    phi <- theta[4]
  } else {
    eta <- theta[1]
    phi <- theta[2]
  }

  xv <- x$model[, 2]
  yv <- x$model[, 1]
  wv <- x$weights

  idx <- !is.na(yv) & !is.na(xv) & !is.na(wv) & !(wv == 0)

  if (sum(idx) != length(yv)) {
    yv <- yv[idx]
    xv <- xv[idx]
    wv <- wv[idx]
  }

  if (!x$is_log) {
    # our functions are defined with the natural logarithm in mind
    # we will take care of the base later
    xv <- log(xv)
  }

  xlim <- dotargs$xlim
  if (is.null(xlim)) {
    xlim <- extendrange(xv, f = 0.08)

    if (xlim[1] > phi) {
      xlim[1] <- phi - 50
    } else if (xlim[2] < phi) {
      xlim[2] <- phi + 50
    }
  }

  ylim <- dotargs$ylim
  if (is.null(ylim)) {
    ylim <- extendrange(yv, f = 0.08)

    if (ylim[1] > alpha) {
      ylim[1] <- alpha - 0.5
    }

    if (ylim[2] < beta) {
      ylim[2] < beta + 0.5
    }

    if ((ylim[1] > 0) && (ylim[1] < 1)) {
      ylim[1] <- 0
    }

    if ((ylim[2] > 0) && (ylim[2] < 1)) {
      ylim[2] <- 1
    }
  }

  xx <- seq(xlim[1], xlim[2], length.out = 500)
  mu <- fn(x, xx, theta)

  level <- dotargs$level
  if (is.null(level)) {
    level <- 0.95
  }

  if ((level > 0) && (level < 1)) {
    q <- qchisq(level, sum(x$estimated))
    cv <- curve_variance(x, xx)
    cs <- sqrt(q * cv)

    upper_bound <- mu + cs
    lower_bound <- mu - cs
  }

  if (base == "e") {
    x_axis_ticks <- pretty(xlim)
    x_axis_labels <- TRUE
  } else {
    xlim <- xlim / k

    xv <- xv / k
    xx <- xx / k

    x1 <- floor(xv[1])
    x2 <- ceiling(xv[length(xv)])

    x_axis_ticks <- seq(x1, x2, by = ceiling((x2 - x1) / 6))

    x_axis_labels <- FALSE
    if (base == "10" || base == "2") {
      x_axis_labels <- str2expression(paste(base, "^", x_axis_ticks, sep = ""))
    }

    # the following is based on https://stackoverflow.com/a/6956596/7073122
    x_axis_minor <- if (base == "10") {
      log10(pretty(10^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else {
      log2(pretty(2^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    }

    x_axis_minor <- x_axis_minor[-c(1, length(x_axis_minor))]
    x_axis_minor <- c(outer(x_axis_minor, x_axis_ticks, `+`))
    x_axis_minor <- x_axis_minor[
      x_axis_minor > xlim[1] & x_axis_minor < xlim[2]
    ]
  }

  dev.hold()

  plot.default(
    xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    axes = FALSE
  )

  axis(1, at = x_axis_ticks, labels = x_axis_labels)
  if (base != "e") {
    axis(1, at = x_axis_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  }
  axis(2, at = pretty(ylim))
  box()

  points(xv, yv, col = col)
  lines(xx, mu, lty = 2, lwd = 2, col = col)

  if ((level > 0) && (level < 1)) {
    xci <- c(xx, rev(xx))
    yci <- c(upper_bound, rev(lower_bound))

    cc <- col2rgb(col) / 255
    cc <- rgb(cc[1], cc[2], cc[3], 0.08)

    polygon(xci, yci, col = cc, border = FALSE)
  }

  if (x$mean_function == "logistic2" || x$mean_function == "logistic4") {
    f <- fn(x, phi, theta)

    if (base != "e") {
      phi <- phi / k
    }

    if (phi > xlim[1] && phi < xlim[2]) {
      lines(
        x = rep(phi, 2), y = c(par("usr")[3], f), lty = 3, col = col
      )

      lines(
        x = c(par("usr")[1], phi), y = rep(f, 2), lty = 3, col = col
      )
    }
  }

  location <- if (eta <= 0) {
    "bottomleft"
  } else {
    "topleft"
  }

  legend_labels <- dotargs$legend
  if (is.null(legend_labels)) {
    legend_labels <- x$mean_function
  }

  legend(
    location, col = col, lty = 2, lwd = 2, bg = "white", legend = legend_labels
  )

  dev.flush()

  invisible(NULL)
}

#' @export
plot.drdalist <- function(x, ...) {
  n_curves <- length(x)

  if (n_curves == 1) {
    return(plot.drda(x[[1L]], ...))
  }

  M <- length(x)

  # non-log fits get the precedence
  is_log <- all(vapply(x, function(y) y$is_log, TRUE))

  dotargs <- list(...)

  col <- dotargs$col
  if (is.null(col)) {
    # choose a colorblind friendly palette
    # https://personal.sron.nl/~pault/
    col <- if (M <= 6) {
      # Paul Tol's bright palette
      c(
        "#EE6677FF", "#4477AAFF", "#228833FF",
        "#AA3377FF", "#66CCEEFF", "#CCBB44FF"
      )
    } else {
      # Paul Tol's muted palette
      c(
        "#CC6677FF", "#332288FF", "#117733FF",
        "#882255FF", "#88CCEEFF", "#999933FF",
        "#AA4499FF", "#44AA99FF", "#DDCC77FF"
      )
    }
  }

  if (length(col) < M) {
    col <- c(col, rep("#000000FF", M - length(col)))
  }

  base <- dotargs$base
  k <- 1
  if (!is.null(base)) {
    if (base == "10") {
      k <- log(10)
    } else if (base == "2") {
      k <- log(2)
    } else {
      # if it's not 2 nor 10, we set by default "e"
      base <- "e"
    }
  } else {
    if (is_log) {
      base <- "e"
    } else {
      # by default we use a log10 scale
      base <- "10"
      k <- log(10)
    }
  }

  xlab <- dotargs$xlab
  if (is.null(xlab)) {
    xlab <- if (base == "e") {
      "log(Predictor)"
    } else {
      "Predictor"
    }
  }

  ylab <- dotargs$ylab
  if (is.null(ylab)) {
    ylab <- "Response"
  }

  # these values we need to decide the best xlim and ylim
  alpha <- 0
  beta <- 1
  eta <- 0
  phi_min <- 0
  phi_max <- 0

  xv <- vector("list", M)
  yv <- vector("list", M)
  wv <- vector("list", M)

  for (i in seq_len(M)) {
    theta <- x[[i]]$coefficients

    if (x[[i]]$mean_function != "logistic2") {
      if (theta[1] < alpha) {
        alpha <- theta[1]
      }

      if (theta[2] > beta) {
        beta <- theta[2]
      }

      if (theta[4] < phi_min) {
        phi_min <- theta[4]
      }

      if (theta[3] <= 0) {
        eta <- eta - 1
      } else {
        eta <- eta + 1
      }

      if (theta[4] > phi_max) {
        phi_max <- theta[4]
      }
    } else {
      if (theta[1] <= 0) {
        eta <- eta - 1
      } else {
        eta <- eta + 1
      }

      if (theta[2] < phi_min) {
        phi_min <- theta[2]
      }

      if (theta[2] > phi_max) {
        phi_max <- theta[2]
      }
    }

    xv[[i]] <- x[[i]]$model[, 2]
    yv[[i]] <- x[[i]]$model[, 1]
    wv[[i]] <- x[[i]]$weights

    idx <- !is.na(yv[[i]]) & !is.na(xv[[i]]) & !is.na(wv[[i]]) & !(wv[[i]] == 0)

    if (sum(idx) != length(yv[[i]])) {
      yv[[i]] <- yv[[i]][idx]
      xv[[i]] <- xv[[i]][idx]
      wv[[i]] <- wv[[i]][idx]
    }

    if (!x[[i]]$is_log) {
      # our functions are defined with the natural logarithm in mind
      # we will take care of the base later
      xv[[i]] <- log(xv[[i]])
    }
  }

  # we use the average of eta to get an idea of the
  xlim <- dotargs$xlim
  if (is.null(xlim)) {
    # use all values to decide the best xlim
    xlim <- extendrange(unlist(xv), f = 0.08)

    # check if some of the fits are bad
    if (xlim[1] > phi_min) {
      xlim[1] <- phi_min - 50
    } else if (xlim[2] < phi_max) {
      xlim[2] <- phi_max + 50
    }
  }

  ylim <- dotargs$ylim
  if (is.null(ylim)) {
    ylim <- extendrange(unlist(yv), f = 0.08)

    if (ylim[1] > alpha) {
      ylim[1] <- alpha - 0.5
    }

    if (ylim[2] < beta) {
      ylim[2] < beta + 0.5
    }

    if ((ylim[1] > 0) && (ylim[1] < 1)) {
      ylim[1] <- 0
    }

    if ((ylim[2] > 0) && (ylim[2] < 1)) {
      ylim[2] <- 1
    }
  }

  xx <- seq(xlim[1], xlim[2], length.out = 500)
  mu <- lapply(x, function(y) fn(y, xx, y$coefficients))

  level <- dotargs$level
  if (is.null(level)) {
    level <- 0.95
  }

  if ((level > 0) && (level < 1)) {
    upper_bound <- vector("list", M)
    lower_bound <- vector("list", M)

    for (i in seq_len(M)) {
      q <- qchisq(level, sum(x[[i]]$estimated))
      cv <- curve_variance(x[[i]], xx)
      cs <- sqrt(q * cv)

      upper_bound[[i]] <- mu[[i]] + cs
      lower_bound[[i]] <- mu[[i]] - cs
    }
  }

  if (base == "e") {
    x_axis_ticks <- pretty(xlim)
    x_axis_labels <- TRUE
  } else {
    xlim <- xlim / k

    xv <- lapply(xv, function(y) y / k)
    xx <- xx / k

    x1 <- min(vapply(xv, function(y) floor(y[1]), 0))
    x2 <- max(vapply(xv, function(y) ceiling(y[length(y)]), 0))

    x_axis_ticks <- seq(x1, x2, by = ceiling((x2 - x1) / 6))

    x_axis_labels <- FALSE
    if (base == "10" || base == "2") {
      x_axis_labels <- str2expression(paste(base, "^", x_axis_ticks, sep = ""))
    }

    # the following is based on https://stackoverflow.com/a/6956596/7073122
    x_axis_minor <- if (base == "10") {
      log10(pretty(10^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else {
      log2(pretty(2^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    }

    x_axis_minor <- x_axis_minor[-c(1, length(x_axis_minor))]
    x_axis_minor <- c(outer(x_axis_minor, x_axis_ticks, `+`))
    x_axis_minor <- x_axis_minor[
      x_axis_minor > xlim[1] & x_axis_minor < xlim[2]
    ]
  }

  dev.hold()

  plot.default(
    xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    axes = FALSE
  )

  axis(1, at = x_axis_ticks, labels = x_axis_labels)
  if (base != "e") {
    axis(1, at = x_axis_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  }
  axis(2, at = pretty(ylim))
  box()

  for (i in seq_len(M)) {
    points(xv[[i]], yv[[i]], col = col[i])
    lines(xx, mu[[i]], lty = 2, lwd = 2, col = col[i])

    if ((level > 0) && (level < 1)) {
      xci <- c(xx, rev(xx))
      yci <- c(upper_bound[[i]], rev(lower_bound[[i]]))

      cc <- col2rgb(col[i]) / 255
      cc <- rgb(cc[1], cc[2], cc[3], 0.08)

      polygon(xci, yci, col = cc, border = FALSE)
    }

    mean_fn <- x[[i]]$mean_function
    phi <- if (mean_fn == "logistic4") {
      x[[i]]$coefficients[4]
    } else  if (mean_fn == "logistic2") {
      x[[i]]$coefficients[2]
    } else {
      NA_real_
    }

    if (!is.na(phi)) {
      f <- fn(x[[i]], phi, x[[i]]$coefficients)

      if (base != "e") {
        phi <- phi / k
      }

      if (phi > xlim[1] && phi < xlim[2]) {
        lines(
          x = rep(phi, 2), y = c(par("usr")[3], f), lty = 3, col = col[i]
        )

        lines(
          x = c(par("usr")[1], phi), y = rep(f, 2), lty = 3, col = col[i]
        )
      }
    }
  }

  location <- if (eta <= 0) {
    "bottomleft"
  } else {
    "topleft"
  }

  legend_labels <- dotargs$legend
  if (is.null(legend_labels)) {
    legend_labels <- vapply(x, function(y) y$mean_function, "a")
  }

  legend(
    location, col = col[1:M], lty = 2, lwd = 2, bg = "white",
    legend = legend_labels
  )

  dev.flush()

  invisible(NULL)
}
