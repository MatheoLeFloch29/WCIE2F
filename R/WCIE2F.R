# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @param mexpo l'objet hlme contenant le modèle de toute l'exposition
#' @param var.time la variable temporelle utilisée
#' @param timerange fenêtre de temps d'exposition que l'utilisateur souhaite analyser
#' @param step pas de temps entre chaque observation dans la fenêtre de temps
#' @param weightbasis la fonction à attribuer au poids de chaque temps d"exposition
#' @param knots le nombre de noeuds internes si utilisation de splines
#' @param data les données avec la variable à étudier
#' @param reg.type le type de regression à utiliser pour étudier la variable d'intérêt
#' @param model la formule du modèle à utiliser pour étudier la variable
#' @param n_boot nombre d'échantillon bootsrap à tirer pour le calcul de la variance
#' @return estimation du modèle avec l'outcome + variance calculer avec bootstrap paramétrique
#' @import dplyr
#' @importFrom splines ns
#' @importFrom stats glm quantile aggregate pnorm qnorm as.formula binomial knots model.matrix step vcov
#' @importFrom lcmm estimates VarCov predictY
#' @importFrom utils data
#' @importFrom MASS mvrnorm
#' @name WCIE2F
#' @title Fonction WCIE2F
#' @export
WCIE2F <- function(mexpo,var.time, timerange, step=1, weightbasis="NS", knots=3,
                   data, reg.type="RL", model,n_boot=500){

  ptm <- proc.time()

  if(is.null(mexpo$data)==T) stop("The argument mexpo need to specify returndata = T")
  if(!inherits(mexpo,"hlme")) stop("The argument mexpo must be a hlme object")
  if (is.null(data)==T) stop("the argument outcome_data is missing")
  if (is.null(model)==T) stop("the argument outcomeformula is missing")
  if (timerange[1]<min(mexpo$data[var.time])) stop("the argument timerange must be equal or higher then the minimum time value")
  if (timerange[2]>max(mexpo$data[var.time])) stop("the argument timerange must be equal or less then the maximum time value")

  #####################################################
  ##### 1) prediction individuelle de l'exposition ####
  #####################################################


  #WCIE <- WCIE_estimation(mexpo = mexpo,windows_min = windows_min,windows_max = windows_max,
  #                        pas = pas,spline = spline,knots = knots,outcome_data = outcome_data,
  #                        outcome_type = outcome_type,outcomeformula = outcomeformula)

  ########################################################
  ####### 4) effect of the exposure ######################
  ########################################################

  # Nombre de tirages bootstrap dans n_boot

  # Extraire la moyenne et la matrice de variance-covariance des paramètres estimés
  mu <- as.matrix(estimates(mexpo))
  Sigma <- as.matrix(VarCov(mexpo))

  # Pour passer de la cholesky à la variance :


  # Générer les nouvelles valeurs des paramètres bootstrap
  boot_params <- mvrnorm(n = n_boot, mu = mu, Sigma = Sigma)


  NPROB = mexpo$N[1]
  NEF = mexpo$N[2]
  NVC = mexpo$N[3]
  idiag0 = mexpo$idiag
  nea0 = sum(mexpo$idea0)


  ######## passer de cholesky à varcov ##################

  if(idiag0==0 & NVC>0){

    for (i in 1:dim(boot_params)[1]) {
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
                          doOneBoot(i = i,boot_params=boot_params)}
                         )

  # remplacer la boucle for par un replicate qui utilise la fonction doOneBoot et qui sort :
  # -une liste des paramètres estimé pour chaque bootstrap
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

  #######################################################
  ########## end bootstrap ##############################
  #######################################################

  ####################################################
  ############## parameters estimations ##############
  ####################################################

  # faire mieux UE TOUT ça
  ## Inter-individual variability
  mean_estimate   <- as.matrix(rowMeans(boot_est))

  # Var tot parameters
  # regarder si la formule est exact (var_inter)
  var_tot <- var_inter + var_intra

  # récupérer la diagonale pour récupérer la variance corrigée avec le bootstrap
  var_estim<-as.matrix(diag(var_tot))

  parameters_var<-cbind(mean_estimate, sqrt(var_estim))
  colnames(parameters_var) <- c("Estimate","Sd")
  parameters_var <- as.data.frame(parameters_var)

  #calculer les IC
  parameters_var$con.low <- parameters_var$Estimate - qnorm(0.975)  * (parameters_var$Sd)
  parameters_var$conf.high <- parameters_var$Estimate + qnorm(0.975)  * (parameters_var$Sd)

  #calculer la stat de test
  parameters_var$z_value <- parameters_var$Estimate / (parameters_var$Sd)

  #calculer la p-valeur (pour une loi normal)
  parameters_var$p_values <- 2 * (1 - pnorm(abs(parameters_var$z_value)))
  #}



  # time
  cost<-proc.time()-ptm

  result<- list(estimate=parameters_var,boots_estimation=t(boot_est),
                new_windows_data=WCIE[[2]],outcome_model=WCIE[[1]],
                exposition_model=WCIE[[3]],outcome_type=outcome_type,
                outcome_call=WCIE[[4]],knots=knots,outcome_var_cov=var_tot,
                time_processing=cost)

  class(result) <- "WCIE2F"
  return(result)

}




