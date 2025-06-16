#' Estimation of the Weighted Cumulative Effect with Two-Level Bootstrap \code{WCIE2F}
#'
#'
#' @description This function estimates the effects of a past exposure on a health outcome
#' using the WCIE method, accounting for parameter uncertainty through a two-level bootstrap.
#' It uses a longitudinal exposure model \code{hlme} and an outcome model
#' to compute time-weighted effects of exposure.
#'
#'
#' @details \strong{A. The estimation follows these steps:}
#'
#' (1) Extraction of the parameters from the mixed model \code{mexpo} :
#' The function extracts the mean (\code{mu}) and the variance-covariance matrix (\code{Sigma})
#' estimated parameters using \code{estimates()} and \code{VarCov()}.
#'
#' (2) Generation of bootstrap parameter sets :
#' Using multivariate normal random sampling (\code{mvrnorm()}), \code{n_boot} samples of
#' parameters are simulated around the estimated mean and variance of the exposure model.
#'
#' (3) Cholesky transformation (if necessary) :
#' If the random effects are correlated (i.e., \code{idiag = 0}), the Cholesky components
#' are transformed into full variance-covariance matrices for each bootstrap.
#' If \code{idiag = 1}, only the variances are adjusted by squaring.
#'
#' (4) Repeated estimation of WCIE and the outcome model :
#' For each bootstrap parameter set, the WCIE is recalculated and the outcome model is fitted.
#'  This step uses the \code{doOneBoot()} function, which
#' returns for each iteration :
#'
#' \itemize{
#'  \item{The estimated coefficients of the outcome model,}
#'  \item{The intra-individual variance matrix.}
#' }
#'
#' (5) Calculation of bootstrap variances (Gelman-Rubin method) :
#' The total variance of the parameters is decomposed as follows :
#'
#' \itemize{
#'  \item\emph{Intra-bootstrap variance} : average of the variance matrices from
#'  each iteration.
#'  \item\emph{Inter-bootstrap variance} : dispersion of the estimators from the different
#'  simulated parameter sets, calculated using the Gelman-Rubin formula: \cr
#'  \deqn{\frac{n+1}{n(n-1)} \sum_{i=1}^n (\theta_i - \bar{\theta})(\theta_i - \bar{\theta})^T}
#'  \item\emph{Total variance} : sum of the two components above.
#' }
#'
#' (6) Construction of confidence intervals and p-values :
#' From the total variance, a corrected standard deviation is calculated. 95% confidence
#' intervals, z-values, and p-values are then produced for each parameter.
#'
#' (7) Temporal Exposure Effect Estimation :
#' Based on the fitted model, the WCIE effects are computed as a linear combination of splines weighted by their estimated coefficients.
#' This approach estimates the time-varying effect of past exposure over a specified time window (defined by \code{timerange}).
#'
#' Two scenarios are handled:
#' \itemize{
#'   \item \emph{Without interactions:} The effect at each time point is given by:
#'   \deqn{w(t) = B(t)^T \theta}
#'   where \eqn{B(t)} is the spline basis at time \eqn{t}, and \eqn{\theta} are the estimated WCIE coefficients.
#'
#'   \item \emph{With interactions:} The time-varying effect is modified according to covariates (e.g., sex, age) using the model’s design matrix.
#'   The spline × covariate interaction terms are incorporated as:
#'   \deqn{w(t,x) = B(t)^T \theta + \sum_j B(t)^T \theta_j x_j}
#' }
#'
#' (8) Uncertainty Quantification :
#' The variance of the WCIE effect at each time point is computed by using the Delta method
#' with an extended form for interaction cases. These variances allow the construction of 95% confidence intervals around the estimated effect.
#'
#' An overall average effect across the time window is also estimated
#' and its variance is derived based on the average of the spline basis across time.
#'
#' (9) Graphical Representation :
#' A graphical output is generated to visualize the estimated time-varying effect of exposure along with its confidence interval.
#' This helps identify critical periods where exposure has a significant impact on the outcome.
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
#' (3) The bootstrap procedure used to estimate standard errors can be computationally intensive.
#' With large datasets or many bootstrap replicates (\code{n_boot}), execution time may increase
#' significantly. It is advised to start with a smaller number of replicates to ensure model
#' convergence before increasing \code{n_boot} for final inference.
#'
#' (4) Convergence problems can occur in the outcome model if the WCIE variables are highly
#' collinear, especially when using a fine time grid (\code{step}) or too many \code{knots}. If the model
#' fails to converge, consider reducing the number of time points or the complexity of the
#' spline basis.
#'
#' (5) When the outcome model includes additional covariates or interaction terms, care should
#' be taken to avoid overfitting or misinterpretation of the WCIE effect. The WCIE component
#' should be interpreted as a marginal cumulative effect, conditional on other covariates in
#' the model
#'
#'
#' @param mexpo An object of class \code{hlme} from the \code{lcmm} package, used to model the exposure process.
#' The object must be created using the \code{hlme} function with specific arguments provided.
#' The fixed effects formula and the random effects formula must be specified. The \code{subject} argument
#' must indicate the subject ID, and the dataset must be provided via the \code{data} argument.
#' It is essential to include \code{returnData = TRUE} in the function call to ensure that the internal data can be accessed.
#' @param var.time character indicating the name of the time variable
#' in the model \code{mexpo}.
#' @param time.frame Numeric vector of length 3
#' indicating the desired time window for exposure (min, max, step).
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
#' @param n_boot Number of bootstrap replicates to perform (default = 500).
#' Bootstrapping is used to approximate the estimation uncertainty and compute the variance of the estimations.
#'
#'
#' @return A list containing:
#'
#' \item{estimate}{Table of final outcome model estimators: mean (Estimate),
#' standard error (Se), confidence intervals (con.low, conf.high), z-statistic,
#' and p-values.}
#' \item{data.expo}{Intermediate data set with individual predictions.}
#' \item{data.outcome}{Data set used to fit the outcome model.}
#' \item{effect.plot}{Graph representing the estimate effects of exposure history.}
#' \item{exposition.effect}{Table with the estimate effects of exposure history, is standard error, IC}
#' \item{mexpo}{Exposure model \code{hlme} object provided at the beginning.}
#' \item{reg.type}{Type of regression model used.}
#' \item{mean.effect}{Mean effect of exposure history.}
#' \item{sd.mean.effect}{Variance of the mean effect of exposure history over time.}
#' \item{nboot}{Number of bootstrap replicates.}
#' \item{call}{The matched call for the outcome model.}
#' \item{knots.quantile}{internal knots uses for splines (used only if \code{weightbasis = "NS"}).}
#' \item{V}{Variance-covariance matrix (intra + inter) of the estimators.}
#' \item{var.time}{Name of the time variable used in the model.}
#' \item{AIC}{Mean AIC of the outcome model across the \code{nboot} bootstrap replicates.}
#' \item{loglike}{Mean log-likelihood of the outcome model across the \code{nboot} bootstrap replicates.}
#' \item{n}{Number of subjects.}
#' \item{nb.subj.del}{Number of subjects removed.}
#' \item{time.processing}{Total computation time (proc.time() object).}
#'
#'
#' @import dplyr
#' @importFrom splines ns
#' @importFrom stats glm quantile aggregate pnorm qnorm as.formula binomial knots model.matrix step vcov formula logLik
#' @importFrom lcmm estimates VarCov predictY
#' @import ggplot2
#'
#' @author un super beau gosse
#'
#' @seealso
#' \code{\link{summary.WCIE2F}}
#' \code{\link{WCIEestimation}}
#' \code{\link{doOneBootWCIE}}
#'
#'
#' @references
#' Maud Wagner et al. “Time-varying associations between an exposure history and a subsequent health
#' outcome : a landmark approach to identify critical windows”. In : BMC Med Res Methodol (2021).
#' doi : 10.1186/s12874-021-01403-w
#'
#' @name WCIE2F
#'
#'
#' @export
WCIE2F <- function(mexpo,var.time, time.frame, weightbasis="NS", knots=NULL,knots.vector=NULL,
                   data, reg.type="RL", model,n_boot=500){

  ptm <- proc.time()

  if(is.null(mexpo$data)==T) stop("The argument mexpo need to specify returndata = T")
  if(!inherits(mexpo,"hlme")) stop("The argument mexpo must be a hlme object")
  if (is.null(data)==T) stop("the argument outcome_data is missing")
  if (is.null(model)==T) stop("the argument outcomeformula is missing")
  #if (timerange[1]<min(mexpo$data[var.time])) stop("the argument timerange must be equal or higher then the minimum time value")
  #if (timerange[2]>max(mexpo$data[var.time])) stop("the argument timerange must be equal or less then the maximum time value")
  if(is.null(knots)==T&is.null(knots.vector)==T) stop("You must have to specify knots or knots.vector")


  # Extraire la moyenne et la matrice de variance-covariance des paramètres estimés du modèle d'exposition
  mu <- as.matrix(estimates(mexpo))
  Sigma <- as.matrix(VarCov(mexpo))

  # Pour passer de la cholesky à la variance :

  # Générer les nouvelles valeurs des paramètres bootstrap
  boot_params <- MASS::mvrnorm(n = n_boot, mu = mu, Sigma = Sigma)

  NPROB = mexpo$N[1]
  NEF = mexpo$N[2]
  NVC = mexpo$N[3]
  idiag0 = mexpo$idiag
  nea0 = sum(mexpo$idea0)

  ######## passer de cholesky à varcov ##################

  if(idiag0==0 & NVC>0){

    for (i in 1:nrow(boot_params)) {
      # Extraire les paramètres Cholesky de l'échantillon i
      Cholesky <- boot_params[i, (NPROB+NEF+1):(NPROB+NEF+NVC)]

      # Construire la matrice triangulaire supérieure U
      U <- matrix(0, nrow = nea0, ncol = nea0)
      U[upper.tri(U, diag = TRUE)] <- Cholesky

      # Transformer en matrice de variance-covariance
      varcov <- t(U) %*% U

      # Remplacer dans boot_params
      boot_params[i, (NPROB+NEF+1):(NPROB+NEF+NVC)] <- varcov[upper.tri(varcov, diag = TRUE)]
    }

  }
  if(idiag0==1 & NVC>0){
    for (i in 1:dim(boot_params)[1]) {
      boot_params[i,(NPROB+NEF+1):(NPROB+NEF+NVC)] <- boot_params[i,(NPROB+NEF+1):(NPROB+NEF+NVC)]**2
    }
  }

  #  if(outcome_type=="RL"){

  ###########################################
  ########### start bootsrap ################
  ###########################################

  # save all the bootstrap dans boot_result

  boot_results <- lapply(1:n_boot, function(i) {
    doOneBootWCIE(i = i,boot_params=boot_params,
                  times = time.frame,mexpo=mexpo,knots.vector=knots.vector,
                  var.time = var.time,weightbasis = weightbasis,knots = knots,
                  data = data, reg.type = reg.type, model = model)}
  )

  # remplacer la boucle for par un replicate qui utilise la fonction doOneBoot et qui sort :
  # -une liste des paramètres estimés pour chaque bootstrap
  # -une liste des matrices de variancecovariance

  # parameters mean for n_boot bootstrap
  boot_est <- sapply(boot_results, function(x) x[[1]])
  # reprendre ici et faire mean
  ## Intra-individual variability
  var_intra <- Reduce("+",lapply(boot_results, function(x) x[[2]]))/n_boot

  ## Inter-individual variability (formule pour bootsrap de gelman rubin)
  #(M+1)/(M(M-1))sum(teta-mean(teta)²)
  mean_boot_est <- apply(boot_est, 1, function(x) mean(x))

  somme_var_var <- Reduce("+", lapply(1:ncol(boot_est), function(i) {
    (boot_est[, i] - mean_boot_est) %*% t(boot_est[, i] - mean_boot_est)# fait la somme des teta-mean(teta)²
  }))

  var_inter <- ((n_boot + 1)/(n_boot*(n_boot-1)))*somme_var_var #applique le (M+1)/(M(M-1))
  rownames(var_inter)<-colnames(var_intra)

  #essayer avec le M-1/M pour voir

  #######################################################
  ########## end bootstrap ##############################
  #######################################################

  ####################################################
  ############## parameters estimations ##############
  ####################################################

  # Var tot parameters
  var_tot <- var_inter + var_intra

  # récupérer la diagonale pour récupérer la variance corrigée avec le bootstrap
  var_estim<-as.matrix(diag(var_tot))

  # calcul écart type
  parameters_var<-cbind(mean_boot_est, sqrt(var_estim))
  colnames(parameters_var) <- c("Estimate","Se")
  parameters_var <- as.data.frame(parameters_var)

  #calculer les IC
  parameters_var$con.low <- parameters_var$Estimate - qnorm(0.975)  * (parameters_var$Se)
  parameters_var$conf.high <- parameters_var$Estimate + qnorm(0.975)  * (parameters_var$Se)

  #calculer la stat de test
  parameters_var$z_value <- parameters_var$Estimate / (parameters_var$Se)

  #calculer la p-valeur (pour une loi normal)
  parameters_var$p_values <- 2 * (1 - pnorm(abs(parameters_var$z_value)))
  #}



  #######################################################################################
  #################### calcul des effets de l'expositon #################################
  #######################################################################################


  # faire tourner un modèle pour quand même récupérer les nouvelles données (dans la nouvelle fenêtre)
  # les nouvelles data du modèle d'expo et avec l'outcome
  # le call du modèle d'exposition
  # utiliser pour le calcul des effets de l'exposition passée dans le temps

  WCIE <- WCEland(mexpo = mexpo,var.time = var.time, times = time.frame,
                         weightbasis = weightbasis, knots = knots,knots.vector=knots.vector,
                         data = data, reg.type = reg.type, model = model)

  new_data <- WCIE$data_expo #data exposition
  data_outcome <- WCIE$data_outcome #data outcome

  n <- n_distinct(data_outcome) # subject nb
  n_obs <- nrow(data_outcome) # nb d'observation total
  nb_subject_delete <- n_distinct(data)-n # nombre de sujet enlevee
  nb_obs_delete <- nrow(data)-n_obs # nb d'observation enlevee

  mean_AIC <- mean(sapply(boot_results, function(x) x$AIC)) #mean AIC bootstrap
  mean_loglike <- mean(sapply(boot_results, function(x) x$loglike)) #mean loglikelihood bootstrap

  # mettre une condition sur le type d'outcome
  if(reg.type=="logistic"){

    ##################### repérer si il y a des intéractions dans le modèle #####################

    # Identifier toutes les interactions contenant "WCIEX" (avant ou après le :)
    var <- grep("WCIE[0-9]+:", rownames(parameters_var), value = TRUE)
    var <- c(var, grep(":WCIE[0-9]+$", rownames(parameters_var), value = TRUE))

    ##################################################################
    ########### si pas d'interaction alors faire ça ##################
    ##################################################################

    # new matrice spline to compile the effect
    data_splines1 <- data.frame(unique(new_data[var.time]))

    data_splines <- seq(from=time.frame[1],to=time.frame[2],by=time.frame[3]) # sequence de mesure dans la fenêtre choisis

    ## splines recompile with the same parameters than put in the wcieestimation function
    new_splines <- as.matrix(ns(unlist(data_splines),knots = WCIE$splines.quantiles,
                                Boundary.knots = WCIE$boundary.quantiles,
                                intercept = T))
    # renommer WCIE1:WCIEk
    colnames(new_splines) <- paste0("WCIE", 1:ncol(new_splines))

    # faire une première boucle sur les WCIE sans intéraction

    # Calculer l'effet total avec une boucle
    effect <- data.frame(time=data_splines,eff_expo=0)
    effect <- effect %>% arrange(effect[1])
    # si il y a pas d'intéraction - faire juste une boucle sur les WCIE
    for (s in 1:dim(new_splines)[2]){
      effect$eff_expo <- effect$eff_expo + parameters_var[paste0("WCIE",s),"Estimate"]*new_splines[,s]
    }

    #################################################################
    ###### si présence d'intéraction alors faire ça : ###############
    #################################################################

    if(length(var)>0){

      # variable présente dans le modèle
      vars_model <- all.vars(formula(model))
      # Variables communes entre les deux
      vars_communes <- intersect(names(data_cov_int),vars_model)
      # recréer le jeu de donnée avec les mêmes valeurs

      # pour que model.matrix fonctionne il doit y avoir au moins deux lignes dans data_cov_int
      # donc créer une deuxième ligne factice et la supprimer après avoir récupérer le bon nom de variable

      # Extraire les variables catégorielles du modèle initial
      vars_categorielle <- names(Filter(is.factor, model$model))

      # Appliquer les mêmes niveaux aux variables catégorielles du nouveau dataset
      # obliger a ne pas mettre le as.factor dans le modèle
      for (var_cat in vars_categorielle) {

        if (var_cat %in% names(data_cov_int)) {
          niveaux <- levels(model$model[[var_cat]])  # Récupérer les niveaux du modèle
          data_cov_int[[var_cat]] <- factor(data_cov_int[[var_cat]], levels = niveaux)  # Appliquer
        }
      }

      # Générer la nouvelle matrice avec les mêmes variables que dans le modèle
      new_data_cov_int <- model.matrix(as.formula(paste("~", paste(colnames(data_cov_int), collapse = " + "))), data = data_cov_int)[,-1]
      new_data_cov_int <- as.data.frame(t(new_data_cov_int))

      # si une variable laisser new_data_cov_int comme celui de base
      if(ncol(data_cov_int)==1 & ncol(new_data_cov_int)==1){
        new_data_cov_int<-data_cov_int
      }

      doOneBoot_effectInt<-function(i){

        # Faire une boucle sur l'ensemble des variables WCIE sans interaction
        effect <- data.frame(time=data_splines,eff_expo=0)
        effect <- effect %>% arrange(effect[1])
        # si il y a pas d'intéraction - faire juste une boucle sur les WCIE
        for (s in 1:dim(new_splines)[2]){
          effect$eff_expo <- effect$eff_expo + boot_est[i,paste0("WCIE",s)]*new_splines[,s]
        }

        ######### + rajouter les intéractions #############

        for(q in 1:ncol(new_data_cov_int)){

          # repérer la place des intéraction
          repaire <- grep(paste0(colnames(new_data_cov_int)[q],"."), colnames(boot_est),value = T)
          repaire <- c(repaire, grep(paste0(".",colnames(new_data_cov_int)[q]), colnames(boot_est),value = T))

          #si 1 variable en interaction alors repasser sous forme de chiffre l'unique valeur pour pouvoir utiliser new_data_cov_int
          p<-1 #si plusieurs covariable intéraction alors réinitialise le compteur des splines à 1
          if (dim(new_data_cov_int)[2]==1) {
            new_data_cov_int2 <- data_cov_int[1,1]
            for (g in 1:length(repaire)) {
              #rajouter valeur
              effect$eff_expo <- effect$eff_expo + boot_est[i,repaire[g]] * new_splines[,p] * new_data_cov_int2[q]
              p <- ifelse(p == (knots + 1), 1, (p + 1))
            }
          }else{
            for (g in 1:length(repaire)) {
              #rajouter valeur
              effect$eff_expo <- effect$eff_expo + boot_est[i,repaire[g]] * new_splines[,p] * new_data_cov_int[,q]
              p <- ifelse(p == (knots + 1), 1, (p + 1))
            }
          }
        }
        return(effect)
      }

      boot_effect <- lapply(1:nrow(boot_est), doOneBoot_effectInt)

      mean_effect<- cbind(boot_effect[[1]][1],rowMeans(sapply(boot_effect, function(x) x[,2])))


    }


    ### si pas d'interaction alors le calcul de la variance des effets se fait ainsi :

    if(length(var)==0){

      # 1) récupérer la matrice de variance co-variances qui nous intéressent
      # Stocker les résultats dans la matrice
      eff_varCov_tot <- var_tot[grepl("WCIE\\d+", rownames(var_tot)),
                                grepl("WCIE\\d+", colnames(var_tot))]

      # 2) calculer la variance pour chaque temps v(w(t))=B(t)'v(teta)B(t)
      effect <- as.data.frame(effect)
      effect$var_eff <- 0
      for(z in 1:nrow(effect)){
        effect$var_eff[z] <- t(new_splines[z,]) %*% eff_varCov_tot %*% new_splines[z,]
      }



    }

    ### si présence d'intéraction alors faire ça : (à faire vérifier car pas sûr)
    if(length(var)>0){

      # 1) récupérer la matrice de variance co-variances qui nous intéressent (avec les intéractions)

      # Stocker les résultats dans la matrice
      eff_varCov_tot <- var_tot[grepl("WCIE\\d+", rownames(var_tot)),
                                grepl("WCIE\\d+", colnames(var_tot))]

      # 2) calculer la variance pour chaque temps v(w(t))=B(t)'v(teta)B(t) en prenant en compte l'intéraction

      #Créer la matrice des splines*les intéraction (autant de fois qu'il y a d'intéraction) B(t)
      # code à revoir (pas sûr mais semble ok)

      # sinon si présence de variable catégorielle alors créer la matrice splines*intéraction de cette manière
      new_matheo<-matrix(0,nrow = nrow(new_splines),ncol = nrow(eff_varCov_tot))
      colnames(new_matheo) <- colnames(eff_varCov_tot)

      # ajouter les valeurs des splines aux endroits correspondants
      for (k in 1:ncol(new_splines)) {
        for (f in 1:ncol(new_matheo)) {
          if(grepl(colnames(new_splines)[k], colnames(new_matheo)[f])){
            new_matheo[,f] <- new_splines[,k]
          }
        }
      }
      # multiplier les valeurs avec les covariables lorsqu'il y a une intéraction
      for (k in 1:ncol(new_data_cov_int)) {
        for (f in 1:ncol(new_matheo)) {
          if(grepl(colnames(new_data_cov_int)[k], colnames(new_matheo)[f])){
            new_matheo[,f] <- new_matheo[,f]*new_data_cov_int[,k]
          }
        }
      }


      # faire le calcul de V(w(t))
      effect <- as.data.frame(effect)
      effect$var_eff <- 0
      for(z in 1:nrow(effect)){
        effect$var_eff[z] <- t(new_matheo[z,]) %*% eff_varCov_tot %*% new_matheo[z,]
      }


      ######################## à verifier mais semble ok ##################################
      # calculer l'effet moyen de 0 à -T
      # calculer l'effet moyen wbarre = 1/T+1 somme(w(u))
      real_mean_effect <- 1/(nrow(effect)+1)*sum(effect[2])

      # sans interaction
      # calculer sa variance v(wbarre=(1/T+1 somme(B(t)'))*v(teta)*(1/T+1 somme(B(t)))
      col_means_splines <- colSums(new_matheo) / (nrow(effect)+1) #1/T+1(B(t))
      real_mean_var_effect <- t(col_means_splines) %*% eff_varCov_tot %*% col_means_splines #v(wbarre)
      #####################################################################################

    }



    # calculer l'interval de confiance des effets

    effect$var_eff <- sqrt(effect$var_eff)
    colnames(effect) <- c("Time","Effect","sd")
    effect <- as.data.frame(effect)

    #calculer les IC
    effect$conf.low <- effect$Effect - qnorm(0.975)  * (effect$sd)
    effect$conf.high <- effect$Effect + qnorm(0.975)  * (effect$sd)


    # calculer l'effet moyen de 0 à -T
    # calculer l'effet moyen wbarre = 1/T+1 somme(w(u))
    mean_effect <- 1/(nrow(effect)+1)*sum(effect[2])

    # sans interaction
    # calculer sa variance v(wbarre=(1/T somme(B(t)'))*v(teta)*(1/T somme(B(t))) confirmer qu'il n'y a pas de +1 ?
    means_splines <- colSums(new_splines) / (nrow(effect)) #1/T+1(B(t))
    real_mean_var_effect <- t(means_splines) %*% eff_varCov_tot %*% means_splines #v(wbarre)
    mean_variable_effect <- sqrt(real_mean_var_effect)

    # avec interaction


    ###########################################################
    ############ faire le graph si ok #########################
    ###########################################################


    ### BOX_PLOT

    # Effect - Intercept
    graph_effect <- ggplot(effect, aes(x=as.factor(Time), y=Effect, group = 1)) +
      geom_line(color="darkgray", linewidth=1)+  # Remplace le boxplot par une ligne
      geom_line(aes(y=conf.low), linetype="dashed", color="black") +  # Bord inférieur en pointillé
      geom_line(aes(y=conf.high), linetype="dashed", color="black") +  # Bord supérieur en pointillé
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="gray", alpha=0.3, linetype="dashed") +  # Ajout IC en pointillé
      theme(legend.position = "none") +
      xlab("Years preceding the outcome") +
      ylab("Estimate") +
      theme_classic() +
      theme(axis.text.x = element_text(size=6, face="bold"),
            plot.title = element_text(size=14, face="bold"),
            axis.title.x = element_text(size=12, face="bold"),
            axis.title.y = element_text(size=12, face="bold")) +
      geom_hline(yintercept=0, linetype="solid")

  }

  # time
  cost<-proc.time()-ptm

  result<- list(estimate=parameters_var,
                data.expo=WCIE[[2]],data.outcome=data_outcome,effectplot=graph_effect,
                expositioneffect=effect,
                mexpo=WCIE[[3]],reg.type=reg.type,mean.effect=mean_effect,
                sd.mean.effect=mean_variable_effect,nboot = n_boot,
                call=WCIE$call, #séparé les différentes parties du call
                knots.quantile=WCIE$splines.quantiles,V=var_tot,var.time=var.time,AIC=mean_AIC
                ,loglike=mean_loglike,n=n,nb.subj.del=nb_subject_delete,
                time.processing=cost)

  class(result) <- "WCIE2F"
  return(result)
}




