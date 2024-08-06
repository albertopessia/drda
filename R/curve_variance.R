# Variance evaluation
#
# Evaluate the variance of the maximum likelihood curve at different predictor
# values.
#
# @param object object fit of class `drda`.
# @param x numeric vector at which to evaluate the variance.
#
# @return Numeric vector with the variances of the maximum likelihood curve.
#
#' @export
curve_variance.drda <- function(object, x) {
  len <- length(x)

  idx <- which(object$estimated)
  if (
    inherits(object, "loglogistic2_fit") || inherits(object, "logistic2_fit")
  ) {
    idx <- idx - 2
  }

  V <- object$vcov[idx, idx, drop = FALSE]

  if (any(is.na(V))) {
    return(rep(NA_real_, len))
  }

  G <- gradient(object, x)

  variance <- rep(NA_real_, len)

  for (i in seq_len(len)) {
    variance[i] <- as.numeric(tcrossprod(crossprod(G[i, idx], V), G[i, idx]))
  }

  variance
}
