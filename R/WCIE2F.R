# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

WCIE2F <- function(mexpo,var.time, timerange, pas=1, weightbasis="NS", knots=3,
                   data, reg.type="RL", model,n_boot=500){

  ptm <- proc.time()

  if(is.null(mexpo$data)==T) stop("The argument mexpo need to specify returndata = T")
  if(!inherits(mexpo,"hlme")) stop("The argument mexpo must be a hlme object")
  if (any((windows_max-windows_min) %% pas != 0)) stop("the argument time must be a multiple of argument pas")
  if (is.null(outcome_data)==T) stop("the argument outcome_data is missing")
  if (is.null(outcomeformula)==T) stop("the argument outcomeformula is missing")
  if (windows_min<min(mexpo$data[sub(".*?\\((\\w+).*", "\\1", mexpo$call[[3]][2])])
  ) stop("the argument windows_min must be equal or higher then the minimum time value")
  if (windows_max> max(mexpo$data[sub(".*?\\((\\w+).*", "\\1", mexpo$call[[3]][2])])
  ) stop("the argument windows_max must be equal or less then the maximum time value")

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

  # Make a function to save all the bootstrap

  doOneBoot <- function(i){
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
    bpt_WCIE<-WCIE_estimation(mexpo = m_expo_boot, windows_min = windows_min,windows_max = windows_max,
                              pas = pas,spline = spline,knots = knots,outcome_data = outcome_data,
                              outcome_type = outcome_type,outcomeformula = outcomeformula)

    return(list(bpt_WCIE[[1]]$coefficients,vcov(bpt_WCIE[[1]])))
  }

  boot_results <- lapply(1:n_boot, doOneBoot)

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


summary.WCIE2F <- function(object) {
  cat( "Call for",object$outcome_type,"model :\n")
  print(deparse(object$outcome_call))
  cat("\n")
  cat("Coefficient :\n")
  object$estimate
  cat("\n")
  cat("à voir ce que l'on veut rajouter dans le summary ? \n")

}
