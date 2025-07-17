#' @method plot WCEland
#' @param x An object of class \code{WCEland} returned by the \code{WCEland} function.
#' @param ... Further arguments (currently not used).
#' @title summary.WCIE2F
#' @return Prints a plot of the exposition historic effects estimates.
#' @name plot.WCEland
#' @export
#'
plot.WCEland <- function(x, ...) {
  print(x$effectplot)
  invisible(x)
}
