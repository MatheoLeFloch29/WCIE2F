# fonction utiliser pour faire i échantillon bootsrtap

#' @param i la ligne i correspondant aux paramètres dans boot_params
#' @return coefficient et la matrice de variance covariance
#' @export
#' @name dOneBoot
#' @title Fonction dOneBoot
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
  bpt_WCIE<-WCIE_estimation(mexpo = m_expo_boot, timerange = timerange,
                            step = step,var.time = var.time,weightbasis = weightbasis,knots = knots,
                            data = data, reg.type = re.type, model = model)

  return(list(bpt_WCIE[[1]]$coefficients,vcov(bpt_WCIE[[1]])))
}
