#' Data simulation according to different outcome models from WCIEpackage
#'
#' @description This function simulates two sample, one exposition sample according to a model estimated
#' with hlme, functions and an outcome sample simulate in relation with the effect of the exposition.
#'
#'
#' @details to do
#'
#'
#'
#' @param object an object of class \code{hlme}
#' @param nsim simulation number
#' @param seed the random seed
#' @param times either a data frame with 2 columns containing IDs and measurement times, or a
#' vector of length 4 specifying the minimal and maximum measurement times, the spacing between
#' 2 consecutive visits and the margin around this spacing
#' @param tname the name of the variable representing the measurement times in \code{object}. Default
#' to the second column's name of times if it is a data frame, and to object$var.time otherwise.
#' @param n number of subjects to simulate*
#' @param internal.step Time increment used for high-resolution simulation before resampling to the user's desired time step.
#' @param Xbin an optional named list giving the probabilities of the binary covariates to
#' simulate. The list's names should match the binary covariate's names used in \code{object}.
#' @param Xcont an optional named list giving the mean and standard deviation of the Gaussian
#' covariates to simulate. The list's names should match the continuous covariate's names used
#' in object.
#' @param weightbasis Type of temporal weighting function used to estimate the Weighted Cumulative Indirect Effects (WCIE).
#' This specifies the functional form used to model the influence of past exposures over time.
#' Currently, the following options are available: \code{"NS"} for natural splines (implemented),
#' \code{"BS"} for B-splines (to be developed), and \code{"PS"} for P-splines (to be developed).
#' @param reg.type Type of outcome model: \code{"logistic"}, \code{"cox"}, ...
#' @param knots number of internal knots
#' @param knots.vector Vector of internal knots for the splines (used only for splines temporal weighting function).
#' @param model two-sided linear formula object for the outcome model. The response outcome is on
#' the left of ~ and the covariates are separated by + on the right of ~. To include the effect of past exposure,
#' you must explicitly add \code{WCIE} (or interaction terms such as \code{WCIE:sex}) to the formula.
#' For example, \code{Y ~ WCIE + age + sex} or \code{Y ~ WCIE:sex + age} are valid formulas.
#' @param coef.wcie named list giving the parameters of the WCIE variables to
#' simulate.
#' @param Xoutcome named list giving the parameters of the intercept and the covariables in the
#' logistic model to simulate de outcome variable. The list's names should match the intercept and
#' covariable names used in object
#'
#'
#'
#' @return
#' \item{exposition.data}{a data frame with one line per observation and one column per variable.
#'  Variables appears in the following order : subject id, measurement time, entry time, binary
#'  covariates, continuous covariates, longitudinal outcomes.}
#' \item{outcome.data}{a data frame with one line per individual and one column per variable.
#'  Variables appears in the following order : subject id, binary
#'  covariates, continuous covariates, binary outcomes.}
#'
#'
#' @import lcmm
#' @import ggplot2
#' @importFrom stats rnorm rbinom runif plogis median simulate
#'
#'
#'
#' @name simulateWCE
#'
#'
#' @export
simulateWCE <- function(object ,nsim=1, seed=NULL, times,internal.step, tname, n, Xbin=NULL, Xcont=NULL
                          ,weightbasis,reg.type,knots=NULL, knots.vector=NULL,model, coef.wcie, Xoutcome){

  if(!is.null(knots)==T & !is.null(knots.vector)==T) stop("You must have to specify knots or knots.vector")


  # parameters
  model_sigma_error <- object$best["stderr"] # sigma du modèle object
  fixe_sigma_error <- 0 # sigma fixé à 0 pour le simulate de hlme
  step_fixe <- internal.step # step fixé pour simulate
  n_after <- nchar(strsplit(as.character(step_fixe), "\\.")[[1]][2]) # chiffres après la virgule du step

  # variable name
  id <- object$call[[4]]
  y <- as.character(object$call$fixed[2])


  ########################################################
  ########## data exposition generation ##################
  ########################################################

  # fixer l'erreur de mesure manuellement à 0 dans le modèle
  object$best["stderr"] <- fixe_sigma_error

  args <- list(
    object = object,
    nsim = nsim,
    times = c(times[1], times[2], step_fixe, 0),
    tname = tname,
    n = n
  )



  # Ajouter seed si non nul
  if (!is.null(seed)) {
    args$seed <- seed
    set.seed(seed)
  }
  # Ajouter Xbin si non nul
  if (!is.null(Xbin)) {
    args$Xbin <- Xbin
  }

  # Ajouter Xcont si non nul
  if (!is.null(Xcont)) {
    args$Xcont <- Xcont
  }

  Sim_D <- do.call(simulate, args)

  # tirer la séquence souhaité pour chaque individu avec un aléa autour de chacune des visite
  t <- seq(times[1],times[2],times[3]) # sequence de temps pour tout les individus
  id_seq <- unique(Sim_D[[id]])

  # Génération de la table avec temps irréguliers pour chaque individu
  visit_data <- do.call(rbind, lapply(id_seq, function(id) {
    # Ajouter un bruit aléatoire à chaque point de temps "arrondir à n_after pour correspondre step_fixe
    t_ind <- t + round(runif(length(t), -abs(times[4]), abs(times[4])),n_after)

    # (optionnel) Forcer le dernier point à être exactement times[2] si il est = 0
     if (times[2] == 0) t_ind[length(t_ind)] <- times[2]
    data.frame(id = id, time = t_ind)
  }))
  colnames(visit_data) <- c(id,tname) #renommer les colonnes comme dans l'objet (model)

  # merge les simulations

  visit_data<-merge(visit_data, Sim_D[c(id,tname,y)],
                    by=c(id,tname))
  visit_data <- visit_data[order(visit_data[,id], visit_data[,tname]), ] #reordonner par id,time

  #appliquer à l'erreur de mesure et
  visit_data[y] <- visit_data[y] + model_sigma_error*rnorm(nrow(visit_data),0,1)

  # recup covariable utilisé dans le modèle
  covar<-object$Xnames2[!object$Xnames2 %in% tname][-1]

  #rajouter les covariables
  exposition_data <- merge(visit_data,
                           unique(Sim_D[c(id,covar)]),# dataset with only id and covariable
                           by=id
  )

  #######################################################
  ######### data outcome generation #####################
  #######################################################

  # pour chaque individu
  data_outcome <- data.frame(id=unique(exposition_data[id]))

      # générer les splines (prendre en compte les y sur le pas de 0.01 dans Sim_D)
      if (weightbasis=="NS") {
        # Créer un vecteur de probabilités (par exemple 5 quantiles => 0.2, 0.4, 0.6, 0.8, 1)


        if (is.null(knots)==F) {
          probs <- seq(0, 1, length.out = knots+2)
          probs <- probs[-c(1, length(probs))]

          # Calculer les quantiles automatiquement sur tout les temps observés
          quantiles <- quantile(exposition_data[,tname], probs = probs, na.rm = TRUE)
          b5  <- quantile(exposition_data[,tname],probs = c(0.05),na.rm=T)
          b95 <- quantile(exposition_data[,tname],probs = c(0.95),na.rm=T)
        }



        if (is.null(knots.vector)==F) {
          quantiles<- knots.vector$knots
          b5 <- knots.vector$boundary.knots[1]
          b95 <- knots.vector$boundary.knots[2]
        }

        if(length(quantiles)==0){
          ## splines recompile (Bk(u)) pour un nombre de noeuds internes à 0
          Zns <- as.matrix(ns(Sim_D[,tname], #knots = quantile,
                              Boundary.knots = c(b5, b95),
                              intercept = T))
        }else{
          ## splines recompile (Bk(u))
          Zns <- as.matrix(ns(Sim_D[,tname], knots = c(quantiles),
                              Boundary.knots = c(b5, b95),
                              intercept = T))
        }



        # récupérer les paramètres Sk donnée par l'utilisateur
        if(ncol(Zns) != length(coef.wcie)) stop("The number of WCIE parameters giving not correspond to the knots number")
        Sk <- coef.wcie

        weights<- Zns%*%Sk #somme(sk*Bk(u)) des k splines

        # proposer un graphique pour la forme de la fonction du poids donné à chaque exposition
        df_weights <- data.frame(
          Time = unique(Sim_D[, tname]),
          Weight = unique(as.numeric(weights))
        )

        # Graphique
        weight_graph<-ggplot(df_weights, aes(x = Time, y = Weight)) +
          geom_line() +
          labs(
            title = "Weight Function Over Time",
            x = "Time",
            y = "Weight"
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10)
          ) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray")

        #############################################################################

        # ressortir les données uniquement au temps d'intérêt times[1] à times[2]
        df_to_out <- unique(df_weights[df_weights$Time %in% t,])

        # somme(w(u)*Xi(true)*pas)
        tmpwcie <- weights*Sim_D[, y]*step_fixe #w(u)*Xi(true)*pas
        sommeWCIE <- tapply(tmpwcie, Sim_D[,id], sum)  #sum() par individu

      }
      #data_outcome$wcie <- sommeWCIE # pour verif glm
    if(reg.type=="logistic"){
      # Récupérer les paramètres (paramètres à donner) de
      # l'intercept et des covariables (présente dans le model d'expo)

      # si pas de covariable
      etai <- Xoutcome$intercept + sommeWCIE #ajouter le paramètre de l'intercept et la somme WCIE

      # récupérer les covariables à intégrer au modèle d'outcome
      outcome_covar <- Xoutcome[names(Xoutcome) !="intercept"]
      # si covariable faire alors ça
      if(length(outcome_covar)>0){
        # Ajouter les covariables de outcome_covar à data_outcome
        data_outcome <- merge(
          data_outcome,
          Sim_D[!duplicated(Sim_D[[id]]), c(id, names(outcome_covar))],
          by = id
        )

        # Initialiser une matrice pour stocker les produits pondérés
        covar_names <- names(outcome_covar)
        tmp_X <- matrix(NA, nrow = nrow(data_outcome), ncol = length(covar_names))
        colnames(tmp_X) <- covar_names

        # Boucle sur les covariables pour remplir tmp_X
        for (i in covar_names) {
          tmp_X[, i] <- data_outcome[[i]] * outcome_covar[[i]]
        }

        etai <- etai + rowSums(tmp_X[, names(outcome_covar),drop=F])
      }

      pi <- 1/(1+exp(-etai)) # probabilité prédite

      # Générer une valeur binaire selon cette proba (proba pi)
      Yobs <- rbinom(n = n, size = 1, prob = pi)

      data_outcome[,"y"] <-  Yobs # l'ajouter au dataoutcome
    }


  return(list(exposition.data = exposition_data,
              outcome.data=data_outcome,
              weight.plot=weight_graph,
              weight.data=df_to_out))
}

## glm(y ~ wcie, data=data_outcome) # -> le coef pour wcie devrait etre egal a 1, et intercept = Xlogisyic$interceot
