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
#'   \item{level}{level of confidence intervals (default is 0.95). Set to zero
#'     or a negative value to disable confidence intervals.}
#'   \item{midpoint}{if `TRUE` (default) shows guidelines associated with the
#'     curve mid-point.}
#'   \item{plot_data}{if `TRUE` (default) shows in the plot the data points used
#'     for fitting.}
#'   \item{legend_show}{if `TRUE` (default) shows the legend.}
#'   \item{legend_location}{character string with custom legend position. See
#'     \code{link[graphics]{legend}} for possible keywords.}
#'   \item{legend}{custom labels for the legend model names.}
#'   \item{show}{If `TRUE` (default) a figure is plotted, otherwise the function
#'     returns a list with values to create the figure manually.}
#' }
#'
#' @return If argument `show` is `TRUE`, no return value. If argument `show` is
#'   `FALSE`, a list with all plotting data.
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

  display_plot <- dotargs[["show"]]
  if (is.null(display_plot)) {
    display_plot <- TRUE
  }

  params <- plot_params(
    x, dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]], level
  )

  if (!display_plot) {
    return(output_params(params))
  }

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

  axis(2, at = params$y_axis_ticks)

  box_x <- par("usr")[params$box$x]
  box_y <- par("usr")[params$box$y]

  if (!is.null(params$box$z)) {
    box_x[1] <- params$box$z[1]
    box_x[10] <- params$box$z[2]
  }

  lines(x = box_x, y = box_y)

  if ((level > 0) && (level < 1)) {
    xci <- params$xci
    yci <- params$yci

    yci[yci < par("usr")[3]] <- par("usr")[3]
    yci[yci > par("usr")[4]] <- par("usr")[4]

    polygon(xci, yci, col = adjustcolor(col, 0.08), border = FALSE)
  }

  if (plot_data) {
    points(params$xv, params$yv, col = col)
  }
  lines(params$xf, params$yf, lty = 2, lwd = 2, col = col)

  if (midpoint && !is.null(params$midpoint_x)) {
    lines(x = params$midpoint_x, y = params$midpoint_y, lty = 3, col = col)
  }

  legend_show <- dotargs[["legend_show"]]
  if (is.null(legend_show) || legend_show) {
    legend_location <- dotargs[["legend_location"]]

    if (is.null(legend_location)) {
      legend_location <- if (params$yf[1] < params$yf[length(params$yf)]) {
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

  display_plot <- dotargs[["show"]]
  if (is.null(display_plot)) {
    display_plot <- TRUE
  }

  params <- vector("list", n_curves)

  plot_type <- 1
  if (inherits(x[[n_curves]], "loglogistic")) {
    plot_type <- 2
  }

  params[[n_curves]] <- plot_params(
    x[[n_curves]], dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]],
    level
  )

  for (i in seq_len(n_curves - 1)) {
    if (
      (inherits(x[[i]], "logistic") && plot_type != 1) ||
      (inherits(x[[i]], "loglogistic") && plot_type != 2)
    ) {
      stop("curves defined on different domains", call. = FALSE)
    }

    params[[i]] <- plot_params(
      x[[i]], dotargs[["base"]], dotargs[["xlim"]], dotargs[["ylim"]], level
    )
  }

  # find common graphical parameters to fit all plots into one figure
  tmp <- vapply(params, function(w) w$xlim, numeric(2))

  j1 <- which.min(tmp[1, ])
  j2 <- which.max(tmp[2, ])
  xlim <- c(tmp[1, j1], tmp[2, j2])

  tmp <- vapply(params, function(w) w$ylim, numeric(2))
  ylim <- c(min(tmp[1, ]), max(tmp[2, ]))

  common <- list()

  common$xlim <- xlim
  common$ylim <- ylim

  if (plot_type == 1) {
    common$x_axis_ticks_1 <- params[[j2]]$x_axis_ticks
    common$x_axis_labels_1 <- params[[j2]]$x_axis_labels
  } else if (plot_type == 2) {
    if (is.null(params[[j1]]$x_axis_ticks_1)) {
      # there is no zero to plot, that is no gap to plot
      # default ticks are those of "j1"
      tks_1 <- params[[j1]]$x_axis_ticks
      lbl_1 <- params[[j1]]$x_axis_labels

      # does "j2" have extra ticks to add?
      tks_2 <- params[[j2]]$x_axis_ticks
      lbl_2 <- params[[j2]]$x_axis_labels

      idx <- !(tks_2 %in% tks_1)
      if (any(idx)) {
        tks_1 <- c(tks_1, tks_2[idx])
        lbl_1 <- c(lbl_1, lbl_2[idx])

        ord <- order(tks_1)
        tks_1 <- tks_1[ord]
        lbl_1 <- lbl_1[ord]
      }

      common$x_axis_ticks_1 <- tks_1
      common$x_axis_labels_1 <- lbl_1
    } else {
      # we got the "gap" to plot
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

      common$x_axis_ticks_1 <- params[[j1]]$x_axis_ticks_1
      common$x_axis_labels_1 <- params[[j1]]$x_axis_labels_1
      common$x_axis_ticks_2 <- tks_1
      common$x_axis_labels_2 <- lbl_1
    }
  }

  common$x_axis_minor <- sort(
    unique(c(params[[j1]]$x_axis_minor, params[[j2]]$x_axis_minor))
  )

  common$y_axis_ticks <- pretty(ylim)

  common$box <- params[[j1]]$box

  if (!display_plot) {
    return(output_params(params, common))
  }

  dev.hold()

  plot.default(
    xlim, ylim, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
    axes = FALSE
  )

  axis(1, at = common$x_axis_ticks_1, labels = common$x_axis_labels_1)

  if (!is.null(common$x_axis_ticks_2)) {
    axis(1, at = common$x_axis_ticks_2, labels = common$x_axis_labels_2)
    axis(
      1, at = tail(common$x_axis_ticks_1, 1), labels = FALSE, tcl = -par("tcl")
    )
    axis(
      1, at = head(common$x_axis_ticks_2, 1), labels = FALSE, tcl = -par("tcl")
    )
  }

  axis(1, at = common$x_axis_minor, labels = FALSE, tcl = par("tcl") * 0.5)

  axis(2, at = common$y_axis_ticks)

  box_x <- par("usr")[common$box$x]
  box_y <- par("usr")[common$box$y]

  if (!is.null(common$box$z)) {
    box_x[1] <- common$box$z[1]
    box_x[10] <- common$box$z[2]
  }

  lines(x = box_x, y = box_y)

  for (i in seq_len(n_curves)) {
    if ((level > 0) && (level < 1)) {
      xci <- params[[i]]$xci
      yci <- params[[i]]$yci

      yci[yci < par("usr")[3]] <- par("usr")[3]
      yci[yci > par("usr")[4]] <- par("usr")[4]

      polygon(xci, yci, col = adjustcolor(col[i], 0.08), border = FALSE)
    }

    if (plot_data) {
      points(params[[i]]$xv, params[[i]]$yv, col = col[i])
    }
    lines(params[[i]]$xf, params[[i]]$yf, lty = 2, lwd = 2, col = col[i])

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
      v <- vapply(params, function(w) w$yf[1] < w$yf[length(w$yf)], FALSE)

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
# @param object `drda` object.
# @param base character string with the base used for printing the values on the
#   `x` axis.
# @param xlim the range of `x` values.
# @param ylim the range of `y` values.
# @param level level of confidence intervals.
#
# @return List with processed graphical parameters.
#
#' @export
plot_params.logistic <- function(object, base, xlim, ylim, level) {
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

  theta <- object$coefficients

  lb <- theta[1]
  ub <- theta[1]
  mp <- 0

  if (theta[2] >= 0) {
    ub <- theta[1] + theta[2]
  } else {
    lb <- theta[1] + theta[2]
  }

  if (object$mean_function == "logistic4") {
    mp <- theta[4]
  } else if (object$mean_function == "logistic2") {
    mp <- theta[4]
  } else if (object$mean_function == "logistic5") {
    mp <- theta[4] + (log(theta[5]) - log(2^theta[5] - 1)) / theta[3]
  } else if (object$mean_function == "gompertz") {
    mp <- theta[4] - log(log(2)) / theta[3]
  } else if (object$mean_function == "logistic6") {
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

  xv <- object$model[, 2]
  yv <- object$model[, 1]
  wv <- object$weights

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
    ylim <- if (object$mean_function != "logistic2") {
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

  xf <- seq(xlim[1], xlim[2], length.out = 500)
  yf <- fn(object, xf, theta)
  cv <- curve_variance(object, xf)

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
    xf <- xf / k

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

  f <- fn(object, mp, theta)

  if (base != "n") {
    mp <- mp / k
  }

  if (mp > xlim[1] && mp < xlim[2]) {
    midpoint_x <- c(xlim[1], mp, mp, use.names = FALSE)
    midpoint_y <- c(f, f, ylim[1], use.names = FALSE)
  }

  xci <- NULL
  yci <- NULL

  if ((level > 0) && (level < 1)) {
    q <- qchisq(level, sum(object$estimated))
    cs <- sqrt(q * cv)
    upper_bound <- yf + cs
    lower_bound <- yf - cs

    xci <- c(xf, rev(xf))
    yci <- c(upper_bound, rev(lower_bound))
  }

  list(
    xlim = xlim,
    ylim = ylim,
    xv = xv,
    yv = yv,
    xf = xf,
    yf = yf,
    xci = xci,
    yci = yci,
    x_axis_ticks = x_axis_ticks,
    x_axis_labels = x_axis_labels,
    x_axis_ticks_1 = x_axis_ticks_1,
    x_axis_labels_1 = x_axis_labels_1,
    x_axis_ticks_2 = x_axis_ticks_2,
    x_axis_labels_2 = x_axis_labels_2,
    x_axis_minor = x_axis_minor,
    y_axis_ticks = pretty(ylim),
    box = box,
    midpoint_x = midpoint_x,
    midpoint_y = midpoint_y
  )
}

# Initialize graphical parameters for a curve defined only for positive values.
#
# @param object `drda` object.
# @param base character string with the base used for printing the values on the
#   `x` axis.
# @param xlim the range of `x` values.
# @param ylim the range of `y` values.
# @param level level of confidence intervals.
#
# @return List with processed graphical parameters.
#' @export
plot_params.loglogistic <- function(object, base, xlim, ylim, level) {
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

  theta <- object$coefficients

  lb <- theta[1]
  ub <- theta[1]
  mp <- 0

  if (theta[2] >= 0) {
    ub <- theta[1] + theta[2]
  } else {
    lb <- theta[1] + theta[2]
  }

  if (object$mean_function == "loglogistic4") {
    mp <- log(theta[4])
  } else if (object$mean_function == "loglogistic2") {
    mp <- log(theta[4])
  } else if (object$mean_function == "loglogistic5") {
    mp <- log(theta[4]) + (log(theta[5]) - log(2^theta[5] - 1)) / theta[3]
  } else if (object$mean_function == "loggompertz") {
    mp <- log(theta[4]) - log(log(2)) / theta[3]
  } else if (object$mean_function == "loglogistic6") {
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

  xv <- object$model[, 2]
  yv <- object$model[, 1]
  wv <- object$weights

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
    ylim <- if (object$mean_function != "loglogistic2") {
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

  xf <- seq(xlim[1], xlim[2], length.out = 500)
  zz <- if (base == "n") {
    xf
  } else {
    exp(xf)
  }

  yf <- fn(object, zz, theta)
  cv <- curve_variance(object, zz)

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
    f <- fn(object, exp(mp), theta)

    mp <- mp / k

    xlim <- xlim / k

    xv <- xv / k
    xf <- xf / k

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
    f <- fn(object, mp, theta)
  }

  midpoint_x <- NULL
  midpoint_y <- NULL

  if (mp > xlim[1] && mp < xlim[2]) {
    midpoint_x <- c(xlim[1], mp, mp, use.names = FALSE)
    midpoint_y <- c(f, f, ylim[1], use.names = FALSE)
  }

  xci <- NULL
  yci <- NULL

  if ((level > 0) && (level < 1)) {
    q <- qchisq(level, sum(object$estimated))
    cs <- sqrt(q * cv)
    upper_bound <- yf + cs
    lower_bound <- yf - cs

    xci <- c(xf, rev(xf))
    yci <- c(upper_bound, rev(lower_bound))
  }

  list(
    xlim = xlim,
    ylim = ylim,
    xv = xv,
    yv = yv,
    xf = xf,
    yf = yf,
    xci = xci,
    yci = yci,
    x_axis_ticks = x_axis_ticks,
    x_axis_labels = x_axis_labels,
    x_axis_ticks_1 = x_axis_ticks_1,
    x_axis_labels_1 = x_axis_labels_1,
    x_axis_ticks_2 = x_axis_ticks_2,
    x_axis_labels_2 = x_axis_labels_2,
    x_axis_minor = x_axis_minor,
    y_axis_ticks = pretty(ylim),
    box = box,
    midpoint_x = midpoint_x,
    midpoint_y = midpoint_y
  )
}

# Plot parameters
#
# Convert a `params` object into a list of user-friendly graphical parameters.
#
# @param params List of graphical parameters as returned by `plot_params`.
# @param common List with common graphical parameters.
#
# @return List with user-friendly graphical parameters.
#
#' @importFrom utils head tail
output_params <- function(params, common = NULL) {
  f <- function(w, i) {
    list(
      data = data.frame(id = rep(i, length(w$xv)), x = w$xv, y = w$yv),
      fitted_curve = data.frame(id = rep(i, length(w$xf)), x = w$xf, y = w$yf),
      confidence_interval = data.frame(
        id = rep(i, length(w$xf)),
        x = w$xf,
        y_lower = rev(tail(w$yci, length(w$xf))),
        y_upper = head(w$yci, length(w$xf))
      ),
      midpoint = c(id = i, x = w$midpoint_x[2], y = w$midpoint_y[2]),
      limits = data.frame(x = w$xlim, y = w$ylim),
      x_axis_ticks = w$x_axis_ticks,
      x_axis_labels = w$x_axis_labels,
      x_axis_ticks_1 = w$x_axis_ticks_1,
      x_axis_labels_1 = w$x_axis_labels_1,
      x_axis_ticks_2 = w$x_axis_ticks_2,
      x_axis_labels_2 = w$x_axis_labels_2,
      x_axis_minor = w$x_axis_minor,
      y_axis_ticks = w$y_axis_ticks,
      box = w$box
    )
  }

  g <- function(params, common) {
    data <- vector("list", length(params))
    fitted_curve <- vector("list", length(params))
    confidence_interval <- vector("list", length(params))
    midpoint <- vector("list", length(params))

    for (i in seq_along(params)) {
      out <- f(params[[i]], i)
      data[[i]] <- out$data
      fitted_curve[[i]] <- out$fitted_curve
      confidence_interval[[i]] <- out$confidence_interval
      midpoint[[i]] <- out$midpoint
    }

    list(
      data = do.call("rbind", data),
      fitted_curve = do.call("rbind", fitted_curve),
      confidence_interval = do.call("rbind", confidence_interval),
      midpoint = do.call("rbind", midpoint),
      limits = data.frame(x = common$xlim, y = common$ylim),
      x_axis_ticks = common$x_axis_ticks,
      x_axis_labels = common$x_axis_labels,
      x_axis_ticks_1 = common$x_axis_ticks_1,
      x_axis_labels_1 = common$x_axis_labels_1,
      x_axis_ticks_2 = common$x_axis_ticks_2,
      x_axis_labels_2 = common$x_axis_labels_2,
      x_axis_minor = common$x_axis_minor,
      y_axis_ticks = common$y_axis_ticks,
      box = common$box
    )
  }

  if (is.null(common)) {
    f(params, 1)
  } else {
    g(params, common)
  }
}
