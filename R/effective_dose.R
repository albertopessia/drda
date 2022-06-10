#' @export
effective_dose.drda <- function(
    object, y = 0.5, type = "relative", level = 0.95
) {
  if (level <= 0 || level >= 1) {
    stop("Confidence level must be in the interval (0, 1)", call. = FALSE)
  }

  lbl <- as.character(round(y, digits = 2))

  alpha <- object$coefficients[1]
  delta <- object$coefficients[2]
  nu <- 1
  xi <- 1
  k <- 1

  if (
    inherits(object, "loglogistic5_fit") || inherits(object, "logistic5_fit")
  ) {
    nu <- object$coefficients[5]
  } else if (
    inherits(object, "loglogistic6_fit") || inherits(object, "logistic6_fit")
  ) {
    xi <- object$coefficients[6]
    k <- xi^(-1 / nu)
  }

  if (type == "relative") {
    y[y <= 0 | y >= 1] <- NA_real_
    y <- alpha + y * delta * k
  } else if (type == "absolute") {
    y1 <- alpha
    y2 <- alpha + delta * k

    if (delta > 0) {
      y[y < y1 | y > y2] <- NA_real_
    } else {
      y[y < y2 | y > y1] <- NA_real_
    }
  } else {
    stop("invalid value for `type`", call. = FALSE)
  }

  x <- inverse_fn(object, y)
  names(x) <- NULL

  idx <- which(object$estimated)
  if (
    inherits(object, "loglogistic2_fit") || inherits(object, "logistic2_fit")
  ) {
    idx <- idx - 2
  }

  V <- object$vcov[idx, idx, drop = FALSE]
  G <- inverse_fn_gradient(object, y)[, idx, drop = FALSE]

  std_err <- if (any(is.na(V))) {
    rep(NA_real_, length(y))
  } else{
    sqrt(diag(tcrossprod(crossprod(t(G), V), G)))
  }
  names(std_err) <- NULL

  q <- qnorm((1 - level) / 2)
  l <- round(level * 100)

  result <- matrix(
    c(
      x,
      x + q * std_err,
      x - q * std_err
    ),
    nrow = length(y),
    ncol = 3,
    dimnames = list(
      lbl, c("Estimate", paste0(c("Lower .", "Upper ."), c(l, l)))
    )
  )

  if (inherits(object, "loglogistic")) {
    result[result < 0] <- 0
  }

  result
}
