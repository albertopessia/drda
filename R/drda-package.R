#' Dose-response data analysis
#'
#' @description
#' `drda` is a package for fitting (log-)logistic curves and performing
#' dose-response data analysis.
#'
#' @section Available functions:
#'
#' Functions specific to `drda`:
#'
#' \itemize{
#'   \item{`drda`: main function for fitting observed data.}
#'   \item{`logistic2_fn`: 2-parameter logistic function.}
#'   \item{`logistic4_fn`: 4-parameter logistic function.}
#'   \item{`logistic5_fn`: 5-parameter logistic function.}
#'   \item{`logistic6_fn`: 6-parameter logistic function.}
#'   \item{`gompertz_fn`: Gompertz function.}
#'   \item{`loglogistic2_fn`: 2-parameter log-logistic function.}
#'   \item{`loglogistic4_fn`: 4-parameter log-logistic function.}
#'   \item{`loglogistic5_fn`: 5-parameter log-logistic function.}
#'   \item{`loglogistic6_fn`: 6-parameter log-logistic function.}
#'   \item{`loggompertz_fn`: log-Gompertz function.}
#'   \item{`nauc`: normalized area under the curve.}
#'   \item{`naac`: normalized area above the curve.}
#' }
#'
#' Functions expected for an object fit:
#'
#' \itemize{
#'   \item{`anova`: compare model fits.}
#'   \item{`deviance`: residual sum of squares of the model fit.}
#'   \item{`logLik`: value of the log-likelihood function associated to the
#'     model fit.}
#'   \item{`plot`: plotting function.}
#'   \item{`predict`: model predictions.}
#'   \item{`print`: basic model summaries.}
#'   \item{`residuals`: model residuals.}
#'   \item{`sigma`: residual standard deviation.}
#'   \item{`summary`: fit summaries.}
#'   \item{`vcov`: approximate variance-covariance matrix of model parameters.}
#'   \item{`weights`: model weights.}
#' }
#'
#' @references Malyutina A, Tang J, Pessia A (2023). drda: An R package for
#'   dose-response data analysis using logistic functions. Journal of
#'   Statistical Software, 106(4), 1-26. doi:10.18637/jss.v106.i04
#'
#' @name drda-package
#' @keywords internal
"_PACKAGE"
