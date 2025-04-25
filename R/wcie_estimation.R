# description de la fonction WCIE_estimation
# Cette fonction permet d'obtenir les estimations beta_k de la mesure d'exposition : WCIE

#' @param mexpo l'objet hlme contenant le modèle de toute l'exposition
#' @param var.time la variable temporelle utilisée
#' @param timerange fenêtre de temps d'exposition que l'utilisateur souhaite analyser
#' @param step pas de temps entre chaque observation dans la fenêtre de temps
#' @param weightbasis la fonction à attribuer au poids de chaque temps d"exposition
#' @param knots le nombre de noeuds internes si utilisation de splines
#' @param data les données avec la variable à étudier
#' @param reg.type le type de regression à utiliser pour étudier la variable d'intérêt
#' @param model la formule du modèle à utiliser pour étudier la variable
#' @return estimation du modèle avec l'outcome
#' @import dplyr
#' @importFrom splines ns
#' @importFrom stats glm quantile aggregate
#' @importFrom lcmm estimates VarCov predictY
#' @export
#' @name WCIE_estimation
#' @title Fonction WCIE_estimation
WCIE_estimation <- function(mexpo,var.time, timerange, step,
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
    new_outcome_date <- merge(data, data_cum, by=mexpo$call[[4]]) #récupère uniquement les individus utilisés dans le modèle d'exposition

    # remplacer les expo par les variables d'exposition dans la formule
    new_expo<-c(NULL)
    new_expo <- paste(paste0("WCIE", seq_len(ncol(B2K))), collapse = "+") ## donne "ns1+ns2+ns3+ns4"

    formdroite <- as.character(model[3]) ## la partie à droite du tilde
    formdroitebis <- gsub("\\bWCIE\\b", paste("(", new_expo, ")"), formdroite) # remplace "expo", par "(ns1+ns2+ns3+ns4)"

    new_formula <- as.formula(paste(as.character(model[2]),"~",formdroitebis))

    model_outcome <- glm(new_formula,family = binomial,data = new_outcome_date)
  }
  if (reg.type=="cox"){

  }
  return(list(model=model_outcome,data=new_data[,1:6],mexpo=mexpo,
              call=new_formula,splines.quantiles=quantiles,
              boundary.quantiles=c(b5,b95)
              )
         )
}







