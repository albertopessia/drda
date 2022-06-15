#' Vorinostat in OPM-2 cell-line dataset
#'
#' A dataset containing dose-response data of drug Vorinostat tested ex-vivo on
#' the OPM-2 cell-line.
#'
#' @format A data frame with 45 rows and 4 variables:
#' \describe{
#'   \item{response}{viability measures normalized using positive and
#'     negative controls}
#'   \item{dose}{drug concentrations (nM) used for testing}
#'   \item{log_dose}{natural logarithm of variable `dose`}
#'   \item{weight}{random weights included only for package demonstration}
#' }
"voropm2"
