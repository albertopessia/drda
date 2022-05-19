#' @export
naac.drda <- function(object, xlim = NULL, ylim = c(0, 1)) {
  1 - nauc(object, xlim, ylim)
}
