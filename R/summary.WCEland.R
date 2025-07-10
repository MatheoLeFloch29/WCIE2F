#' @method summary WCEland
#' @param object An object of class \code{WCEland} returned by the \code{WCEland} function.
#' @param ... Further arguments (currently not used).
#' @title summary.WCEland
#' @return Prints a summary of the WCEland analysis including model information, estimation results, and bootstrap details.
#' @name summary.WCEland
#' @export
#'
#' @examples
#' \dontrun{
#' res <- WCEland(mexpo = mexpo,
#'     var.time = "TIME",
#'     timerange = c(0, 10),
#'     data = data,
#'     model = outcome_model)
#' summary(res)
#' }
summary.WCEland <- function(object, ...) {

  cat("Weighted Cumulative Index Estimation (WCEland) \n")
  cat("     fitted by weighted bootstrap method \n")
  cat(" \n")
  cat(paste(deparse(object$call), collapse = " "), "\n")
  cat(" \n")
  if (object$reg.type=="logistic"|object$reg.type=="cox") {
    cat("Outcome model type:", object$reg.type, "\n")
    cat(" \n")
    cat("Statistical Model:", "\n")
    #cat(paste("     Dataset:", as.character(as.expression(object$call[]))),"\n")
    cat(paste("     Number of subjects:", object$n),"\n")
    #cat(paste("     Number of observations:", "wait"),"\n")
    #if(length(x$na.action))
    cat(paste("     Number of subjects deleted:", object$nb.subj.del),"\n")
    cat(" \n")
    cat("Weighting basis:", object$weightbasis, "\n")
    cat("Number of bootstrap replications:", object$nboot, "\n")
    cat("Time variable:", object$var.time, "\n")
    cat("Knots for weight basis:", length(object$knots.quantile), "\n\n")
    cat(" \n")
    cat("Parameter Estimates (bootstrap-based):\n")
    cat(" \n")

    parameters <- object$estimate[,-c(3:4)]
    colnames(parameters) <- c("coef", "Se", "Wald", "p-value")
    parameters$coef <- format(round(parameters$coef,5),scientific=FALSE)
    parameters$Se <- format(round(parameters$Se,5),scientific=FALSE)
    parameters$`p-value` <- format(round(parameters$`p-value`,5),scientific=FALSE)
    parameters$Wald <- format(round(parameters$Wald,3),scientific=FALSE)
    print(parameters, row.names = TRUE)

    cat(" \n")
    cat("Summary of Goodness-of-fit statistics (over all bootstrap samples):  \n")
    cat(paste("     maximum log-likelihood:", object$loglike)," \n")
    cat(paste("     AIC:", object$AIC)," \n")
    cat(" \n")
    cat("Variance-covariance matrix:\n")
    cat(" \n")
    print(object$V)
    cat(" \n")
    cat(" \n")
    cat("Summary of weighted cumulative exposure effect:\n")
    cat(" \n")
    cat("Weighted cumulative exposure effect per time observation:\n")
    cat(" \n")

    effects <- object$expositioneffect[,-c(4:5)]
    rownames(effects) <- NULL
    colnames(effects) <- c("Time","Effect", "Se")
    effects$Wald <- effects$Effect / (effects$Se)
    effects$`p-value` <- 2 * (1 - pnorm(abs(effects$Effect / (effects$Se))))
    effects$Effect <- format(round(effects$Effect,5),scientific=FALSE)
    effects$Se <- format(round(effects$Se,5),scientific=FALSE)
    effects$Wald <- format(round(effects$Wald,3),scientific=FALSE)
    effects$`p-value` <- format(round(effects$`p-value`,5),scientific=FALSE)
    print(effects, row.names = TRUE)

    cat(" \n")
    cat("     Mean effect:", sprintf("%.5f", object$mean.effect), "\n")
    cat("     Se of effect:", sprintf("%.5f", object$sd.mean.effect), "\n")
    cat("     Wald:", sprintf("%.3f", object$mean.effect / (object$sd.mean.effect)), "\n")
    cat("     p-value:", sprintf("%.5f", 2 * (1 - pnorm(abs(object$mean.effect / (object$sd.mean.effect))))), "\n")

  }
  invisible(object)
}


