#' @method summary WCIE2F
#' @export
summary.WCIE2F <- function(object, ...) {
  cat( "Call for",object$outcome_type,"model :\n")
  print(deparse(object$outcome_call))
  cat("\n")
  cat("Coefficient :\n")
  object$estimate
  cat("\n")
  cat("\u00e0 voir ce que l'on veut rajouter dans le summary ? \n")

}
