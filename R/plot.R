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
#' If multiple fit objects from the same family of models are given, they will
#' all be placed in the same plot.
#'
#' Accepted plotting arguments are:
#'
#' \describe{
#'   \item{base}{character string with the base used for printing the values on
#'     the `x` axis. Accepted values are `10` for base 10, `2` for base 2, `e`
#'     for base e, or `n` (default) for no log-scale printing.}
#'   \item{col}{curve color(s). By default, up to 9 color-blind friendly
#'     colors are provided.}
#'   \item{xlab, ylab}{axis labels.}
#'   \item{xlim, ylim}{the range of x and y values with sensible defaults.}
#'   \item{level}{level of confidence intervals. Set to zero or a negative value
#'     to disable confidence intervals.}
#'   \item{midpoint}{if `FALSE` do not show guidelines associated with the
#'     curve mid-point.}
#'   \item{plot_data}{if `FALSE` do not show data points used for fitting in the
#'     plot.}
#'   \item{legend_show}{if `FALSE` do not show the legend.}
#'   \item{legend_location}{character string with custom legend position. See
#'     \code{link[graphics]{legend}} for possible keywords.}
#'   \item{legend}{custom labels for the legend model names.}
#' }
#'
#' @return No return value.
#'
#' @importFrom graphics axis box curve legend lines par polygon points
#' @importFrom graphics plot.default
#' @importFrom grDevices adjustcolor dev.flush dev.hold extendrange
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

  if (length(fit_objects) > 0) {
    return(plot.drdalist(c(list(x), fit_objects), ...))
  }

  col <- dotargs[["col"]]
  if (is.null(col)) {
    col <- "#EE6677FF"
  }

  xlab <- dotargs[["xlab"]]
  if (is.null(xlab)) {
    xlab <- "Predictor"
  }

  ylab <- dotargs[["ylab"]]
  if (is.null(ylab)) {
    ylab <- "Response"
  }

  main <- dotargs[["main"]]
  if (is.null(main)) {
    main <- ""
  }

  level <- dotargs[["level"]]
  if (is.null(level)) {
    level <- 0.95
  }

  midpoint <- dotargs[["midpoint"]]
  if (is.null(midpoint)) {
    midpoint <- TRUE
  }

  plot_data <- dotargs[["plot_data"]]
  if (is.null(plot_data)) {
    plot_data <- TRUE
  }

  params <- plot_params(
    x, dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]]
  )

  dev.hold()

  plot.default(
    params$xlim, params$ylim, type = "n", xlim = params$xlim,
    ylim = params$ylim, xlab = xlab, ylab = ylab, axes = FALSE, main = main
  )

  if (is.null(params$x_axis_ticks_1)) {
    axis(1, at = params$x_axis_ticks, labels = params$x_axis_labels)
  } else {
    axis(1, at = params$x_axis_ticks_1, labels = params$x_axis_labels_1)
    axis(1, at = params$x_axis_ticks_2, labels = params$x_axis_labels_2)

    axis(1, at = params$x_axis_ticks_1[2], labels = FALSE, tcl = -par("tcl"))
    axis(1, at = params$x_axis_ticks_2[1], labels = FALSE, tcl = -par("tcl"))

    axis(1, at = params$x_axis_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  }

  axis(2, at = pretty(params$ylim))

  box_x <- par("usr")[params$box$x]
  box_y <- par("usr")[params$box$y]

  if (!is.null(params$box$z)) {
    box_x[1] <- params$box$z[1]
    box_x[10] <- params$box$z[2]
  }

  lines(x = box_x, y = box_y)

  if ((level > 0) && (level < 1)) {
    q <- qchisq(level, sum(x$estimated))
    cs <- sqrt(q * params$cv)
    upper_bound <- params$mu + cs
    lower_bound <- params$mu - cs

    xci <- c(params$xx, rev(params$xx))
    yci <- c(upper_bound, rev(lower_bound))

    yci[yci < par("usr")[3]] <- par("usr")[3]
    yci[yci > par("usr")[4]] <- par("usr")[4]

    polygon(xci, yci, col = adjustcolor(col, 0.08), border = FALSE)
  }

  if (plot_data) {
    points(params$xv, params$yv, col = col)
  }
  lines(params$xx, params$mu, lty = 2, lwd = 2, col = col)

  if (midpoint && !is.null(params$midpoint_x)) {
    lines(x = params$midpoint_x, y = params$midpoint_y, lty = 3, col = col)
  }

  legend_show <- dotargs[["legend_show"]]
  if (is.null(legend_show) || legend_show) {
    legend_location <- dotargs[["legend_location"]]

    if (is.null(legend_location)) {
      legend_location <- if (params$mu[1] < params$mu[length(params$mu)]) {
        "bottomright"
      } else {
        "topright"
      }
    }

    legend_labels <- dotargs[["legend"]]
    if (is.null(legend_labels)) {
      legend_labels <- x$mean_function
    }

    legend(
      legend_location, col = col, lty = 2, lwd = 2, bg = "white",
      legend = legend_labels
    )
  }

  dev.flush()

  invisible(NULL)
}

#' @export
plot.drdalist <- function(x, ...) {
  n_curves <- length(x)

  if (n_curves == 1) {
    return(plot.drda(x[[1L]], ...))
  }

  dotargs <- list(...)

  col <- dotargs[["col"]]
  if (is.null(col)) {
    # choose a colorblind friendly palette
    # https://personal.sron.nl/~pault/
    col <- if (n_curves <= 6) {
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

  if (length(col) < n_curves) {
    col <- c(col, rep("#000000FF", n_curves - length(col)))
  }

  xlab <- dotargs[["xlab"]]
  if (is.null(xlab)) {
    xlab <- "Predictor"
  }

  ylab <- dotargs[["ylab"]]
  if (is.null(ylab)) {
    ylab <- "Response"
  }

  main <- dotargs[["main"]]
  if (is.null(main)) {
    main <- ""
  }

  level <- dotargs[["level"]]
  if (is.null(level)) {
    level <- 0.95
  }

  midpoint <- dotargs[["midpoint"]]
  if (is.null(midpoint)) {
    midpoint <- TRUE
  }

  plot_data <- dotargs[["plot_data"]]
  if (is.null(plot_data)) {
    plot_data <- TRUE
  }

  params <- vector("list", n_curves)

  plot_type <- 1
  if (inherits(x[[n_curves]], "loglogistic")) {
    plot_type <- 2
  }

  params[[n_curves]] <- plot_params(
    x[[n_curves]], dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]]
  )

  for (i in seq_len(n_curves - 1)) {
    if (
      (inherits(x[[i]], "logistic") && plot_type != 1) ||
      (inherits(x[[i]], "loglogistic") && plot_type != 2)
    ) {
      stop("curves defined on different domains", call. = FALSE)
    }

    params[[i]] <- plot_params(
      x[[i]], dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]]
    )
  }

  tmp <- vapply(params, function(w) w$xlim, numeric(2))

  j1 <- which.min(tmp[1, ])
  j2 <- which.max(tmp[2, ])
  xlim <- c(tmp[1, j1], tmp[2, j2])

  tmp <- vapply(params, function(w) w$ylim, numeric(2))
  ylim <- c(min(tmp[1, ]), max(tmp[2, ]))

  dev.hold()

  plot.default(
    xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    axes = FALSE
  )

  if (plot_type == 1) {
    axis(1, at = params[[j2]]$x_axis_ticks, labels = params[[j2]]$x_axis_labels)
  } else if (plot_type == 2) {
    # on the left side of the plot we must use the smallest values
    axis(
      1, at = params[[j1]]$x_axis_ticks_1, labels = params[[j1]]$x_axis_labels_1
    )

    # default ticks are those of "j1"
    tks_1 <- params[[j1]]$x_axis_ticks_2
    lbl_1 <- params[[j1]]$x_axis_labels_2

    # does "j2" have extra ticks to add?
    # remove the "gap" tick because we are using that of "j1"
    tks_2 <- params[[j2]]$x_axis_ticks_2[-1]
    lbl_2 <- params[[j2]]$x_axis_labels_2[-1]

    idx <- !(tks_2 %in% tks_1)
    if (any(idx)) {
      tks_1 <- c(tks_1, tks_2[idx])
      lbl_1 <- c(lbl_1, lbl_2[idx])

      ord <- order(tks_1)
      tks_1 <- tks_1[ord]
      lbl_1 <- lbl_1[ord]
    }

    axis(1, at = tks_1, labels = lbl_1)

    axis(
      1, at = params[[j1]]$x_axis_ticks_1[2], labels = FALSE, tcl = -par("tcl")
    )
    axis(
      1, at = params[[j1]]$x_axis_ticks_2[1], labels = FALSE, tcl = -par("tcl")
    )
  }

  axis(
    1,
    at = sort(unique(c(params[[j1]]$x_axis_minor, params[[j2]]$x_axis_minor))),
    labels = FALSE, tcl = par("tcl") * 0.5
  )

  axis(2, at = pretty(ylim))

  box_x <- par("usr")[params[[j1]]$box$x]
  box_y <- par("usr")[params[[j1]]$box$y]

  if (!is.null(params[[j1]]$box$z)) {
    box_x[1] <- params[[j1]]$box$z[1]
    box_x[10] <- params[[j1]]$box$z[2]
  }

  lines(x = box_x, y = box_y)

  for (i in seq_len(n_curves)) {
    if ((level > 0) && (level < 1)) {
      q <- qchisq(level, sum(x[[i]]$estimated))
      cs <- sqrt(q * params[[i]]$cv)
      upper_bound <- params[[i]]$mu + cs
      lower_bound <- params[[i]]$mu - cs

      xci <- c(params[[i]]$xx, rev(params[[i]]$xx))
      yci <- c(upper_bound, rev(lower_bound))

      yci[yci < par("usr")[3]] <- par("usr")[3]
      yci[yci > par("usr")[4]] <- par("usr")[4]

      polygon(xci, yci, col = adjustcolor(col[i], 0.08), border = FALSE)
    }

    if (plot_data) {
      points(params[[i]]$xv, params[[i]]$yv, col = col[i])
    }
    lines(params[[i]]$xx, params[[i]]$mu, lty = 2, lwd = 2, col = col[i])

    if (midpoint && !is.null(params[[i]]$midpoint_x)) {
      midpoint_x <- c(xlim[1], params[[i]]$midpoint_x[-1])
      midpoint_y <- c(params[[i]]$midpoint_y[-3], ylim[1])
      lines(x = midpoint_x, y = midpoint_y, lty = 3, col = col[i])
    }
  }

  legend_show <- dotargs[["legend_show"]]
  if (is.null(legend_show) || legend_show) {
    legend_location <- dotargs[["legend_location"]]

    if (is.null(legend_location)) {
      v <- vapply(params, function(w) w$mu[1] < w$mu[length(w$mu)], FALSE)

      legend_location <- if (mean(v) > 0.5) {
        "bottomright"
      } else {
        "topright"
      }
    }

    legend_labels <- dotargs[["legend"]]
    if (is.null(legend_labels)) {
      legend_labels <- vapply(x, function(w) w$mean_function, "a")
    }

    legend(
      legend_location, col = col, lty = 2, lwd = 2, bg = "white",
      legend = legend_labels
    )
  }

  dev.flush()

  invisible(NULL)
}

# Initialize graphical parameters for a curve defined over the whole real line.
#
# @param x `drda` object.
# @param base character string with the base used for printing the values on the
#   `x` axis.
# @param xlim the range of `x` values.
# @param ylim the range of `y` values.
#
# @return List with processed graphical parameters.
plot_params.logistic <- function(x, base, xlim, ylim) {
  # these constants are used for proper scaling on the requested base
  k <- 1

  if (!is.null(base)) {
    if (base == "10") {
      k <- log(10)
    } else if (base == "2") {
      k <- log(2)
    } else if (base != "e" && base != "n") {
      # only "n", e", "2", and "10" are supported.
      stop("base value not supported", call. = FALSE)
    }
  } else {
    # by default we do not change the scale
    base <- "n"
  }

  theta <- x$coefficients

  lb <- theta[1]
  ub <- theta[1]
  mp <- 0

  if (theta[2] >= 0) {
    ub <- theta[1] + theta[2]
  } else {
    lb <- theta[1] + theta[2]
  }

  if (x$mean_function == "logistic4") {
    mp <- theta[4]
  } else if (x$mean_function == "logistic2") {
    mp <- theta[4]
  } else if (x$mean_function == "logistic5") {
    mp <- theta[4] + (log(theta[5]) - log(2^theta[5] - 1)) / theta[3]
  } else if (x$mean_function == "gompertz") {
    mp <- theta[4] - log(log(2)) / theta[3]
  } else if (x$mean_function == "logistic6") {
    q <- theta[6]^(-1 / theta[5])

    if (theta[2] >= 0) {
      ub <- theta[1] + theta[2] * q
    } else {
      lb <- theta[1] + theta[2] * q
    }

    mp <- theta[4] + (
      log(theta[5]) - log(theta[6]) - log(2^theta[5] - 1)
    ) / theta[3]
  } else {
    stop("model not supported", call. = FALSE)
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

  # make sure data is sorted according to `xv`
  ord <- order(xv, yv, wv)
  xv <- xv[ord]
  yv <- yv[ord]
  wv <- wv[ord]

  if (is.null(xlim)) {
    xlim <- extendrange(xv, f = 0.08)

    if (xlim[1] > mp) {
      xlim[1] <- mp - 50
    } else if (xlim[2] < mp) {
      xlim[2] <- mp + 50
    }
  }

  if (is.null(ylim)) {
    ylim <- if (x$mean_function != "logistic2") {
      extendrange(yv, f = 0.08)
    } else {
      c(0, 1)
    }

    if (ylim[1] > lb) {
      ylim[1] <- lb - 0.5
    }

    if (ylim[2] < ub) {
      ylim[2] <- ub + 0.5
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
  cv <- curve_variance(x, xx)

  x_axis_ticks <- pretty(xlim)
  x_axis_labels <- TRUE

  x_axis_minor <- NULL

  x_axis_ticks_1 <- NULL
  x_axis_labels_1 <- NULL
  x_axis_ticks_2 <- NULL
  x_axis_labels_2 <- NULL

  box <- list(
    x = c(
      # bottom-left border
      2, 1,
      # left border
      1, 1,
      # top border
      1, 2,
      # right border
      2, 2,
      # bottom-right border
      2, 2
    ),
    y = c(
      # bottom-left border
      3, 3,
      # left border
      3, 4,
      # top border
      4, 4,
      # right border
      3, 4,
      # bottom-right border
      3, 3
    ),
    z = NULL
  )

  if (base != "n") {
    xlim <- xlim / k

    xv <- xv / k
    xx <- xx / k

    x1 <- floor(xv[1])
    x2 <- ceiling(xv[length(xv)])

    x_axis_ticks <- seq(x1, x2, by = ceiling((x2 - x1) / 6))
    x_axis_labels <- str2expression(paste(base, "^", x_axis_ticks, sep = ""))

    # the following is based on https://stackoverflow.com/a/6956596/7073122
    x_axis_minor <- if (base == "10") {
      log10(pretty(10^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else if (base == "2") {
      log2(pretty(2^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else {
      log(pretty(exp(x_axis_ticks[1:2]), 3)) - x_axis_ticks[1]
    }

    x_axis_minor <- x_axis_minor[-c(1, length(x_axis_minor))]
    x_axis_minor <- c(outer(x_axis_minor, x_axis_ticks, `+`))
    x_axis_minor <- x_axis_minor[
      x_axis_minor > xlim[1] & x_axis_minor < xlim[2]
    ]
  }

  midpoint_x <- NULL
  midpoint_y <- NULL

  f <- fn(x, mp, theta)

  if (base != "n") {
    mp <- mp / k
  }

  if (mp > xlim[1] && mp < xlim[2]) {
    midpoint_x <- c(xlim[1], mp, mp)
    midpoint_y <- c(f, f, ylim[1])
  }

  list(
    xlim = xlim,
    ylim = ylim,
    xv = xv,
    yv = yv,
    xx = xx,
    zz = xx,
    mu = mu,
    cv = cv,
    x_axis_ticks = x_axis_ticks,
    x_axis_labels = x_axis_labels,
    x_axis_ticks_1 = x_axis_ticks_1,
    x_axis_labels_1 = x_axis_labels_1,
    x_axis_ticks_2 = x_axis_ticks_2,
    x_axis_labels_2 = x_axis_labels_2,
    x_axis_minor = x_axis_minor,
    box = box,
    midpoint_x = midpoint_x,
    midpoint_y = midpoint_y
  )
}

# Initialize graphical parameters for a curve defined only for positive values.
#
# @param x `drda` object.
# @param base character string with the base used for printing the values on the
#   `x` axis.
# @param xlim the range of `x` values.
# @param ylim the range of `y` values.
#
# @return List with processed graphical parameters.
plot_params.loglogistic <- function(x, base, xlim, ylim) {
  # these constants are used for proper scaling on the requested base
  k <- 1
  h <- exp(3)

  if (!is.null(base)) {
    if (base == "10") {
      k <- log(10)
      h <- exp(5)
    } else if (base == "2") {
      k <- log(2)
      h <- exp(2)
    } else if (base != "e" && base != "n") {
      # only "n", e", "2", and "10" are supported.
      stop("base value not supported", call. = FALSE)
    }
  } else {
    # by default we do not change the scale
    base <- "n"
  }

  theta <- x$coefficients

  lb <- theta[1]
  ub <- theta[1]
  mp <- 0

  if (theta[2] >= 0) {
    ub <- theta[1] + theta[2]
  } else {
    lb <- theta[1] + theta[2]
  }

  if (x$mean_function == "loglogistic4") {
    mp <- log(theta[4])
  } else if (x$mean_function == "loglogistic2") {
    mp <- log(theta[4])
  } else if (x$mean_function == "loglogistic5") {
    mp <- log(theta[4]) + (log(theta[5]) - log(2^theta[5] - 1)) / theta[3]
  } else if (x$mean_function == "loggompertz") {
    mp <- log(theta[4]) - log(log(2)) / theta[3]
  } else if (x$mean_function == "loglogistic6") {
    q <- theta[6]^(-1 / theta[5])

    if (theta[2] >= 0) {
      ub <- theta[1] + theta[2] * q
    } else {
      lb <- theta[1] + theta[2] * q
    }

    mp <- log(theta[4]) + (
      log(theta[5]) - log(theta[6]) - log(2^theta[5] - 1)
    ) / theta[3]
  } else {
    stop("model not supported", call. = FALSE)
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

  # make sure data is sorted according to `xv`
  ord <- order(xv, yv, wv)
  xv <- xv[ord]
  yv <- yv[ord]
  wv <- wv[ord]

  len <- length(xv)
  zero_x <- xv == 0

  if (base != "n") {
    xv[zero_x] <- min(xv[!zero_x]) / h
    xv <- log(xv)
  }

  if (is.null(xlim)) {
    xlim <- if (base == "n") {
      c(0, xv[len] * 1.08)
    } else {
      c(xv[1], xv[len] * 1.08)
    }

    if (base == "n") {
      if (xlim[1] > exp(mp)) {
        xlim[1] <- exp(mp - 3)
      } else if (xlim[2] < exp(mp)) {
        xlim[2] <- exp(mp + 3)
      }
    } else {
      if (xlim[1] > mp) {
        xlim[1] <- mp - 3
      } else if (xlim[2] < mp) {
        xlim[2] <- mp + 3
      }
    }
  } else {
    if (xlim[1] < 0 || xlim[2] <= 0 || xlim[1] >= xlim[2]) {
      stop("wrong 'xlim' values", call. = FALSE)
    }

    if (base != "n") {
      if (xlim[1] == 0) {
        xlim <- c(xv[1], log(xlim[2]))
      } else {
        xlim <- log(xlim)
      }
    }
  }

  if (is.null(ylim)) {
    ylim <- if (x$mean_function != "loglogistic2") {
      extendrange(yv, f = 0.08)
    } else {
      c(0, 1)
    }

    if (ylim[1] > lb) {
      ylim[1] <- lb - 0.5
    }

    if (ylim[2] < ub) {
      ylim[2] <- ub + 0.5
    }

    if ((ylim[1] > 0) && (ylim[1] < 1)) {
      ylim[1] <- 0
    }

    if ((ylim[2] > 0) && (ylim[2] < 1)) {
      ylim[2] <- 1
    }
  }

  xx <- seq(xlim[1], xlim[2], length.out = 500)
  zz <- if (base == "n") {
    xx
  } else {
    exp(xx)
  }

  mu <- fn(x, zz, theta)
  cv <- curve_variance(x, zz)

  x_axis_ticks <- pretty(xlim)
  x_axis_labels <- TRUE

  x_axis_minor <- NULL

  x_axis_ticks_1 <- NULL
  x_axis_labels_1 <- NULL
  x_axis_ticks_2 <- NULL
  x_axis_labels_2 <- NULL

  box <- list(
    x = c(
      # bottom-left border
      2, 1,
      # left border
      1, 1,
      # top border
      1, 2,
      # right border
      2, 2,
      # bottom-right border
      2, 2
    ),
    y = c(
      # bottom-left border
      3, 3,
      # left border
      3, 4,
      # top border
      4, 4,
      # right border
      3, 4,
      # bottom-right border
      3, 3
    ),
    z = NULL
  )

  f <- 0.0

  if (base != "n") {
    f <- fn(x, exp(mp), theta)

    mp <- mp / k

    xlim <- xlim / k

    xv <- xv / k
    xx <- xx / k

    x1 <- floor(xv[1])
    x2 <- ceiling(xv[len])

    x_axis_ticks <- seq(x1, x2, by = ceiling((x2 - x1) / 6))
    x_axis_labels <- str2expression(
      if (any(zero_x)) {
        c("0", paste(base, "^", x_axis_ticks[-1], sep = ""))
      } else {
        paste(base, "^", x_axis_ticks, sep = "")
      }
    )

    # the following is based on https://stackoverflow.com/a/6956596/7073122
    x_axis_minor <- if (base == "10") {
      log10(pretty(10^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else if (base == "2") {
      log2(pretty(2^x_axis_ticks[1:2], 9)) - x_axis_ticks[1]
    } else {
      log(pretty(exp(x_axis_ticks[1:2]), 3)) - x_axis_ticks[1]
    }

    x_axis_minor <- x_axis_minor[-c(1, length(x_axis_minor))]
    x_axis_minor <- c(outer(x_axis_minor, x_axis_ticks, `+`))
    x_axis_minor <- x_axis_minor[
      x_axis_minor > xlim[1] & x_axis_minor < xlim[2]
    ]

    if (any(zero_x)) {
      # when there is a zero dose we must show a gap in the x-axis
      x_axis_ticks[1] <- xv[1]

      p1 <- 0.25 * (x_axis_ticks[2] - x_axis_ticks[1])
      p2 <- 0.35 * (x_axis_ticks[2] - x_axis_ticks[1])

      x1 <- x_axis_ticks[1] + p1
      x2 <- x_axis_ticks[1] + p2

      x_axis_minor <- x_axis_minor[x_axis_minor >= x2]

      x_axis_ticks_1 <- c(x_axis_ticks[1], x1)
      x_axis_labels_1 <- c(x_axis_labels[1], "")
      x_axis_ticks_2 <- c(x2, x_axis_ticks[-1])
      x_axis_labels_2 <- c("", x_axis_labels[-1])

      box$z <- c(x1, x2)
    }
  } else {
    mp <- exp(mp)
    f <- fn(x, mp, theta)
  }

  midpoint_x <- NULL
  midpoint_y <- NULL

  if (mp > xlim[1] && mp < xlim[2]) {
    midpoint_x <- c(xlim[1], mp, mp)
    midpoint_y <- c(f, f, ylim[1])
  }

  list(
    xlim = xlim,
    ylim = ylim,
    xv = xv,
    yv = yv,
    xx = xx,
    zz = zz,
    mu = mu,
    cv = cv,
    x_axis_ticks = x_axis_ticks,
    x_axis_labels = x_axis_labels,
    x_axis_ticks_1 = x_axis_ticks_1,
    x_axis_labels_1 = x_axis_labels_1,
    x_axis_ticks_2 = x_axis_ticks_2,
    x_axis_labels_2 = x_axis_labels_2,
    x_axis_minor = x_axis_minor,
    box = box,
    midpoint_x = midpoint_x,
    midpoint_y = midpoint_y
  )
}
