#'Estimation of the Weighted Cumulative Effect with Two-Level
#'
#' @description
#' This function estimates the effect of an exposure history on a health outcome
#' (such as a binary event, survival time, or repeated event),
#' using the WCIE (Weighted Cumulative Index of Exposure) approach.
#' It is designed for longitudinal data where exposure varies over time,
#' and where the goal is to model the effect of this time-varying exposure on an outcome.
#'
#'
#' A. Estimation follows the following steps:
#'
#' (1) Individual prediction of exposure:
#' Based on the mixed model \code{mexpo}, which must be an object of class \code{hlme},
#' individual exposure trajectories are predicted over a user-defined time window
#' (\code{timerange}) and frequency (\code{step}).
#' The presence or absence of a random intercept is taken into account, but random effects
#' must be explicitly specified.
#' The temporal structure (e.g., splines, polynomials) used in the original model is preserved,
#' but must be directly included in the \code{mexpo} model.
#'
#' (2) Reconstruction of the exposure history:
#' A temporal weighting basis is used to construct the WCIE, i.e., time-weighted exposure scores.
#' Currently, only natural splines (NS) are implemented. The spline knots are automatically
#'  computed based on time quantiles, according to the number of internal knots (\code{knots})
#'  specified by the user.
#'
#' (3) Computation of the WCIE:
#' The WCIE corresponds to the weighted sum of the products between time spline basis functions
#'  and predicted exposure values.
#' Each spline generates a separate WCIE component (e.g., WCIE1, WCIE2, ...).
#'
#' (4) Outcome modeling:
#' A model is then fitted to the outcome, incorporating the WCIE components as explanatory
#' variables.
#' Two types of models are currently available via the \code{reg.type} argument:
#' - "logistic": logistic regression (GLM model)
#' - "cox": Cox proportional hazards model (under development)
#' The final formula is automatically reconstructed based on the original model, replacing
#' the exposure term with the sum of WCIE components.
#'
#' The function output includes:
#' - the fitted outcome model (glm),
#' - the temporal predictions (Ypred) and random effects,
#' - the time quantiles used for splines,
#' - the final formula used in the outcome model.
#'
#'
#' \strong{B. USAGE PRECAUTIONS:}
#'
#' (1) The estimation of WCIE relies on the quality of the exposure model provided in \code{mexpo}.
#' Since this model is used to simulate exposure histories, incorrect specification of fixed
#' or random effects may bias the WCIE estimates. It is strongly recommended to use natural
#' splines (\code{NS}) for the time structure and ensure that returndata = TRUE is specified in the \code{hlme}
#' object.
#'
#' (2) The weighting function used to estimate the WCIE over time can
#' lead to instability or overfitting if the number of \code{knots} is too high or if the time window
#' (\code{timerange}) is poorly chosen. Users should inspect the shape of the weight function and
#' consider simplifying it if convergence issues or unreasonable effects are detected.
#'
#' (3) Convergence problems can occur in the outcome model if the WCIE variables are highly
#' collinear, especially when using a fine time grid (\code{step}) or too many \code{knots}. If the model
#' fails to converge, consider reducing the number of time points or the complexity of the
#' spline basis.
#'
#'
#' @param mexpo An object of class \code{hlme} from the \code{lcmm} package, used to model the exposure process.
#' The object must be created using the \code{hlme} function with specific arguments provided.
#' The fixed effects formula and the random effects formula must be specified. The \code{subject} argument
#' must indicate the subject ID, and the dataset must be provided via the \code{data} argument.
#' It is essential to include \code{returnData = TRUE} in the function call to ensure that the internal data can be accessed.
#' @param var.time character indicating the name of the time variable
#' in the model \code{mexpo}.
#' @param timerange Numeric vector of length 2
#' indicating the desired time window for exposure (min, max).
#' @param step Step between two consecutive time points in the time window indicating in \code{timerange}.
#' @param weightbasis Type of temporal weighting function used to estimate the Weighted Cumulative Indirect Effects (WCIE).
#' This specifies the functional form used to model the influence of past exposures over time.
#' Currently, the following options are available: \code{"NS"} for natural splines (implemented),
#' \code{"BS"} for B-splines (to be developed), and \code{"PS"} for P-splines (to be developed).
#' @param knots Number of internal knots for the splines (used only for splines temporal weighting function).
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
#' @return A list containing:
#'   \item{model}{outcome model object}
#'   \item{data_expo}{Intermediate data set with individual predictions.}
#'   \item{mexpo}{Exposure model \code{hlme} object provided at the beginning.}
#'   \item{data_outcome}{Data set used to fit the outcome model.}
#'   \item{call}{The matched call for the outcome model.}
#'   \item{splines.quantiles}{The internal quantiles used for the natural splines,
#'   which define the time-weighting function.}
#'   \item{boundary.quantiles}{The lower and upper bounds of the natural splines
#'   (5th and 95th percentiles), used to define the limits of the time-weighting function.}
#'   \item{AIC}{AIC criteria of the outcome model.}
#'   \item{loglike}{log-likelihood of the outcome model.}
#'
#'
#' @import dplyr
#' @importFrom splines ns
#' @importFrom stats glm quantile aggregate
#' @importFrom lcmm estimates VarCov predictY
#' @author Encore un giga beau gosse
#'
#' @seealso \code{\link{WCIE2F}}
#' @references
#' Maud Wagner et al. “Time-varying associations between an exposure history and a subsequent health
#' outcome : a landmark approach to identify critical windows”. In : BMC Med Res Methodol (2021).
#' doi : 10.1186/s12874-021-01403-w
#' @name WCIEestimation
#' @export
WCIEestimation <- function(mexpo,var.time, timerange, step,
                            weightbasis, knots, data, reg.type, model){

  #####################################################
  ##### 1) prediction individuelle de l'exposition ####
  #####################################################

  # fenêtre d'exposition souhaitée par l'utilisateur
  timerange_min <- timerange[1]
  timerange_max <- timerange[2]

  # Créer un jeu de données avec le même nb d'individu et le bon nombre de ligne par individu

  time_seq <- seq(from=timerange_min,to=timerange_max,by=step) # sequence de mesure dans la fenêtre choisis

  nb_ind <- mexpo$ns  # nb d'individu
  id_seq<-mexpo$pprob[[mexpo$call[[4]]]] # séquence d'identifiant utilisée dans le modèle (en vecteur)

  new_data <- data.frame(id= rep(id_seq, each = length(time_seq)) # reprendre la séquence d'ID données par les data de l'individu
  )

  # renommer l'identifiant comme celui de l'utilisateur
  colnames(new_data)[1]<- mexpo$call[[4]]

  # tire une séquence pour chaque individu de la fenêtre qu'il veut analyser
  new_data[var.time] <- rep(time_seq, nb_ind)

  #########################################################################################################################
  ###################### prendre en compte les différentes fonctions du temps possible que l'utilisateur peut rentrer #####
  #########################################################################################################################


  ######################## si la personne utilise des bs/ns/PF directement dans le modèle #################################
  # obligation de demander de remplir directementla fonction dans la formule du modèle hlmr sinon impossible de connaitre ce que la personne à
  # utiliser comme fonction

  # recompile les effets aléatoires pour la nouvelle fenêtre de données (pour n'importe quelle faction du temps ou pas)

  variable_RE <- model.matrix(as.formula(paste("~", mexpo$call[[3]][2])), data = new_data)

  # recup covariable utilisé dans le modèle
  covar<-mexpo$Xnames2[!mexpo$Xnames2 %in% var.time][-1]

  value_fixe_covar <- distinct(mexpo$data[c(covar,mexpo$call[[4]])])
  new_data <- merge(new_data, value_fixe_covar, by = mexpo$call[[4]])

  ######## predexpo ###############

  # prediction (du modèle)predictY au temps de la fenêtre donné
  predexpo <- predictY(mexpo,newdata = new_data,var.time = var.time)

  new_data <- cbind(new_data, predexpo$pred)

  ######### +  les predRE si il y  a des effets aléatoire dans le modèle ############

  # récupérer les paramêtres des effets aléatoires du model et les renommer (prendre en compte la présence d'un intercept aléatoire)
  predRE <- mexpo$predRE
  if((grepl("^-1", as.character(mexpo$call$random)[2])==F)){
    variable_RE <- variable_RE[,-1] # enlever l'intercept du varRE généré par le model.matrix
    for (h in 1:(ncol(predRE)-2)) {
      colnames(predRE)[h+2] <- paste0("predRE",h)
    }
  }else{ # si pas d'intercept aléatoire
    for (h in 1:(ncol(predRE)-1)) {
      colnames(predRE)[h+1] <- paste0("predRE",h)
    }
  }
  new_data <- merge(new_data, predRE,by =  mexpo$call[[4]]) #merge les predRE (effets aléatoires) par individu

  # si présence d'un intercept aléatoire alors le rajouter dans le calcul de la prédiction
  if(grepl("^-1", as.character(mexpo$call$random)[2])==F){
    new_data$Ypred <- new_data$Ypred + new_data$intercept
  }

  # predictY + effet aléatoire*variableRE
  # si au moins 2 colonnes alors faire la boucle sinon pas de boucle (else)
  if(is.vector(variable_RE)==F) {
    for (m in 1:(ncol(variable_RE))) {
      new_data$Ypred <- new_data$Ypred +
        (variable_RE[,m]*new_data[paste0("predRE",m)])
    }
  }else{
    new_data$Ypred <- new_data$Ypred +
      (variable_RE*new_data["predRE1"])
  }

  ################################################################################################################




  ##########################################################
  ###### 2) Récupérer l'historique d'exposition ############
  ##########################################################
  data_expo_pred <- new_data

  if (weightbasis=="NS") {

    # Créer un vecteur de probabilités (par exemple 5 quantiles => 0.2, 0.4, 0.6, 0.8, 1)
    probs <- seq(0, 1, length.out = knots+1)
    probs <- probs[-c(1, length(probs))]

    # Calculer les quantiles automatiquement
    quantiles <- quantile(data_expo_pred[var.time], probs = probs, na.rm = TRUE)

    b5  <- quantile(data_expo_pred[var.time],probs = c(0.05),na.rm=T)
    b95 <- quantile(data_expo_pred[var.time],probs = c(0.95),na.rm=T)

    ## splines recompile
    B2K <- as.matrix(ns(unlist(data_expo_pred[var.time]),knots = quantiles, Boundary.knots = c(b5, b95), intercept = T))
    data_expo_pred <-cbind(data_expo_pred,B2K)


    ## faire la prédiction * les nouveaux splines (Xi(Tu)*Bk(Tu))
    for (r in 1:length(colnames(B2K))) {
      data_expo_pred[paste0("COCO",r)] <- data_expo_pred$Ypred * data_expo_pred[colnames(B2K)[r]]
    } #Xi fois les Bk

    # somme cummulé pondéré des expositions (sum(Xi*Bk)=Fki)
    data_cum <- unique(data_expo_pred[mexpo$call[[4]]])
    for (r in 1:length(colnames(B2K))) {
      WCIE<-aggregate(data_expo_pred[[paste0("COCO", r)]] ~ data_expo_pred[[mexpo$call[[4]]]],
                      data = data_expo_pred, FUN = sum)
      colnames(WCIE)<-c(mexpo$call[[4]],paste0("WCIE", r)) #renommer les colonnes car aggregate ne renomme pas
      data_cum <-  cbind(data_cum,WCIE[paste0("WCIE", r)])
    }
  }
  if(weightbasis=="PS"){

    #to be develop
  }
  if(weightbasis=="BS"){
    #same

  }

  ##############################################
  ##### 3) Outcome model #######################
  ##############################################

  if (reg.type=="logistic"){
    data.outcome <- merge(data, data_cum, by=mexpo$call[[4]]) #récupère uniquement les individus utilisés dans le modèle d'exposition et également présent dans le data pour le modèle d'expo

    # remplacer les expo par les variables d'exposition dans la formule
    new_expo<-c(NULL)
    new_expo <- paste(paste0("WCIE", seq_len(ncol(B2K))), collapse = "+") ## donne "ns1+ns2+ns3+ns4"

    formdroite <- as.character(model[3]) ## la partie à droite du tilde
    formdroitebis <- gsub("\\bWCIE\\b", paste("(", new_expo, ")"), formdroite) # remplace "expo", par "(ns1+ns2+ns3+ns4)"

    new_formula <- as.formula(paste(as.character(model[2]),"~",formdroitebis))

    model_outcome <- glm(new_formula,family = binomial,data = data.outcome)


    # récupérer quelques statistiques pour le summary :

    #log likelihood and AIC
    AIC <- AIC(model_outcome)
    loglike <- as.numeric(logLik(model_outcome))

    # call
    glm_call <- deparse(model_outcome$call)
    real_formula <- deparse(new_formula)
    glm_call_updated <- gsub("new_formula", real_formula, glm_call)

  }
  if (reg.type=="cox"){

  }
  return(list(model=model_outcome,data_expo=new_data[,1:6],
              mexpo=mexpo, data_outcome=data_outcome,
              call=glm_call_updated,splines.quantiles=quantiles,
              boundary.quantiles=c(b5,b95),
              AIC = AIC, loglike=loglike
              )
         )
}







