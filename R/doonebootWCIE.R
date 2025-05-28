#' @title Bootstrap Estimation of WCIE for a Single Replicate
#'
#' @description
#' This internal function is used to perform one iteration of the bootstrap procedure
#' for estimating WCIE and fitting the corresponding outcome model using perturbed parameters.
#'
#' @details
#' The function takes a set of bootstrap parameters from \code{boot_params}, reconstructs the exposure model
#' using those parameters, and then estimates the WCIE values. These WCIE components are then used to
#' fit the outcome model. It returns the estimated regression coefficients, their variance-covariance matrix,
#' the AIC, and the log-likelihood of the outcome model.
#'
#' @param mexpo An object of class \code{hlme} from the \code{lcmm} package, used to model the exposure process.
#' The object must be created using the \code{hlme} function with specific arguments provided.
#' The fixed effects formula and the random effects formula must be specified. The \code{subject} argument
#' must indicate the subject ID, and the dataset must be provided via the \code{data} argument.
#' It is essential to include \code{returnData = TRUE} in the function call to ensure that the internal data can be accessed.
#' @param i Integer. The row index corresponding to a bootstrap sample in \code{boot_params}.
#' @param boot_params A matrix or data frame containing all bootstrap parameter vectors.
#' @param var.time character indicating the name of the time variable
#' in the model \code{mexpo}.
#' @param times Numeric vector of length 4
#' indicating the desired time window for exposure (min, max, step, alea).
#' @param weightbasis Type of temporal weighting function used to estimate the Weighted Cumulative Indirect Effects (WCIE).
#' This specifies the functional form used to model the influence of past exposures over time.
#' Currently, the following options are available: \code{"NS"} for natural splines (implemented),
#' \code{"BS"} for B-splines (to be developed), and \code{"PS"} for P-splines (to be developed).
#' @param knots number of internal knots
#' @param knots.vector list of 2 vector : one for the internal knots for the splines (used only
#' for splines temporal weighting function) and a second for the boundary knots.
#' (ex:knots.vector=list(knots=c(-15,-5),boundary.knots=c(-20,0)))
#' @param data A data frame containing the variables specified in the outcome model
#' \code{model} including the outcome variable. This dataset will be used to estimate the
#' outcome model, and the WCIE variables calculated previously will be added to this
#' dataset prior to model fitting.
#' @param reg.type Type of outcome model: \code{"logistic"}, \code{"cox"}, ...
#' @param model two-sided linear formula object for the outcome model. The response outcome is on
#' the left of ~ and the covariates are separated by + on the right of ~. To include the effect of past exposure,
#' you must explicitly add \code{WCIE} (or interaction terms such as \code{WCIE:sex}) to the formula.
#' For example, \code{Y ~ WCIE + age + sex} or \code{Y ~ WCIE:sex + age} are valid formulas.
#'
#'
#' @return A list containing:
#'   \item{coef}{Estimated coefficients of the outcome model.}
#'   \item{V}{Variance-covariance matrix (intra + inter) of the estimators.}
#'   \item{AIC}{AIC criteria of the outcome model.}
#'   \item{loglike}{log-likelihood of the outcome model.}
#'
#'
#' @name doOneBootWCIE
#'
#' @export
doOneBootWCIE <- function(boot_params,i,times,mexpo,
                          var.time,weightbasis,knots,knots.vector,
                          data, reg.type, model){
  call_hlme<-""
  # Transformer l'appel du modèle en chaîne de caractères et enlever data=nom du jeu de donnée de l'utilisateur
  mexpo$call$data <- NULL
  call_hlme <- deparse(mexpo$call)
  # remettre l'argument data = mexpo$data pour récupérer le jeu de donnée de l'utilisateur
  call_hlme <- gsub("\\)$", ", data=mexpo$data)", call_hlme)
  # Ajouter l'argument maxiter = 0 à la chaîne
  call_hlme <- gsub("\\)$", ", maxiter = 0)", call_hlme)
  #ajouter le B=boot_params[i,] pour utiliser les paramètres des n bootsrap
  call_hlme <- gsub("\\)$", ", B=boot_params[i,])", call_hlme)

  # Évaluer la nouvelle chaîne comme un appel de fonction
  m_expo_boot <- eval(parse(text = call_hlme))

  # obliger de faire sans passer par une fonction
  bpt_WCIE<-WCIEestimation(mexpo = m_expo_boot, times = times,
                           var.time = var.time,weightbasis = weightbasis,knots = knots,
                           knots.vector = knots.vector,
                           data = data, reg.type = reg.type, model = model)

  # changer la façon de récupérer la matrice de variance covariance suivant le reg.type


  return(list(coef=bpt_WCIE[[1]]$coefficients,V=vcov(bpt_WCIE[[1]]),
              AIC=bpt_WCIE[[8]],loglike=bpt_WCIE[[9]]))
}


