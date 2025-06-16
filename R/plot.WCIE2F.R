#' @method plot WCIE2F
#' @param object An object of class \code{WCIE2F} returned by the \code{WCIE2F} function.
#' @param ... Further arguments (currently not used).
#' @title summary.WCIE2F
#' @return Prints a plot of the exposition historic effects estimates.
#' @name plot.WCIE2F
#' @export
#'
#' @examples
#' \dontrun{
#' }
plot.WCIE2F <- function(object, ...) {
  print(object$effectplot)
  invisible(object)
}
