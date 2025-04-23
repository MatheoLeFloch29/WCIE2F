# description de la fonction WCIE_estimation
# Cette fonction permet d'obtenir les estimations beta_k de la mesure d'exposition : WCIE

# les paramètres à renseigner sont :
# mexpo : l'objet hlme contenant le modèle de toute l'exposition
# var.time : la variable temporelle utilisée
# timerange : la fenêtre de temps d'exposition que l'utilisateur souhaite analyser
# step : le pas de temps entre chaque observation dans la fenêtre de temps
# weightbasis : la fonction à attribuer au poids de chaque temps d"exposition
# knots : le nombre de noeuds internes si utilisation de splines
# data : les données avec la variable à étudier
# reg.type : le type de regression à utiliser pour étudier la variable d'intérêt
# model : la formule du modèle à utiliser pour étudier la variable


WCIE_estimation <- function(mexpo,var.time, timerange, step,
                            weightbasis, knots, data, reg.type, model){

  #####################################################
  ##### 1) prediction individuelle de l'exposition ####
  #####################################################

  # fenêtre d'exposition souhaitée par l'utilisateur
  timerange_min <- timerange[1]
  timerange_max <- timerange[2]

  # Créer un jeu de données avec le même nb d'individu et le bon nombre de ligne par individu

  time_seq <- seq(from=timerange_max,to=timerange_min,by=-step) # sequence de mesure dans la fenêtre choisis

  nb_ind <- mexpo$ns  # nb d'individu
  new_data <- data.frame(id= rep(seq(1,nb_ind),
                                 each = length(time_seq))
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
  var_RE2 <- model.matrix(as.formula(paste("~", mexpo$call[[3]][2])), data = new_data)

  for (h in 1:(ncol(var_RE2)-1)) {
    colnames(var_RE2)[h+1] <- paste0("var_RE",h)
  }
  var_RE1 <- var_RE2[,-1]

  new_data <- cbind(new_data, var_RE1)

  # recup covariable utilisé dans le modèle
  covar<-mexpo$Xnames2[!mexpo$Xnames2 %in% var.time][-1]

  value_fixe_covar <- distinct(mexpo$data[c(covar,mexpo$call[[4]])])
  new_data <- merge(new_data, value_fixe_covar, by = mexpo$call[[4]])

  ######## predexpo ###############

  # prediction (du modèle)predictY au temps de la fenêtre donné
  predexpo <- predictY(mexpo,newdata = new_data,var.time = var.time)

  new_data <- cbind(new_data, predexpo$pred)

  ######### +  les predRE si il y  a des effets aléatoire dans le modèle ############

  # si la deuxième ligne à le même nom alors renommer la deuxième (uniquement lorsque il y seulement la variable de temps prise en compte)
  if(colnames(mexpo$predRE)[ncol(mexpo$predRE)]==var.time){
    colnames(mexpo$predRE)[ncol(mexpo$predRE)]<- "Alea_effect"
  }

  # récupérer les paramêtres du model
  new_data <- merge(new_data,mexpo$predRE,by = mexpo$call[[4]])

  new_data$pred_final <- 0
  # si présence d'un intercept aléatoire alors le rajouter dans le calcul de la prédiction
  if(grepl("^-1", as.character(mexpo$call$random)[2])==F){
    new_data$pred_final <- new_data$Ypred + new_data$intercept
    # renommer les predRE
    new_mexpopredRE <- mexpo$predRE[-2]
    for (g in 1:(ncol(new_mexpopredRE)-1)) {
      colnames(new_mexpopredRE)[g+1] <- paste0("varRE_model",g)
    }
    new_data <- merge(new_data,new_mexpopredRE,by = mexpo$call[[4]])
    name_se<- "varRE_model"
  }else{
    new_mexpopredRE <- mexpo$predRE
    for (g in 1:(ncol(new_mexpopredRE)-1)) {
      colnames(new_mexpopredRE)[g+1] <- paste0("varRE_model",g)
    }
    new_data <- merge(new_data,new_mexpopredRE,by = mexpo$call[[4]])
    name_se<- "varRE_model"
  }

  # predictY + effet aléatoire
  for (m in 1:(ncol(var_RE2)-1)) {
    new_data$pred_final <- new_data$pred_final +
      (new_data[paste0("var_RE",m)]*new_data[paste0(name_se,m)])
  }
  ################################################################################################################




  ##########################################################
  ###### 2) Récupérer l'historique d'exposition ############
  ##########################################################
  data_expo_pred <- new_data
  data_expo_pred <- data_expo_pred %>%     #remettre dans l'ordre le  temps par individu
    arrange(.data[[mexpo$call[[4]]]],.data[[var.time]])


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
      data_expo_pred[paste0("COCO",r)] <- data_expo_pred$pred_final * data_expo_pred[colnames(B2K)[r]]
    } #Xi fois les Bk

    # somme cummulé pondéré des expositions (sum(Xi*Bk)=Fki)
    data_cum <- data_expo_pred[mexpo$call[[4]]]
    for (r in 1:length(colnames(B2K))) {

      WCIE<-aggregate(data_expo_pred[[paste0("COCO", r)]] ~ data_expo_pred[[mexpo$call[[4]]]],
                      data = data_expo_pred, FUN = sum) #
      colnames(WCIE)<-c(mexpo$call[[4]],paste0("WCIE", r))

      data_cum <-  merge(data_cum,WCIE,by = mexpo$call[[4]])
    }
    data_cum <- unique(data_cum)

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
    new_outcome_date <- merge(data, data_cum, by=mexpo$call[[4]])

    # remplacer les expo par les variables d'exposition dans la formule
    new_expo<-c(NULL)
    for(p in 1:length(colnames(B2K))) {
      new_expo <- paste(c(new_expo, paste0("WCIE",p)), collapse = "+") ## donne "ns1+ns2+ns3+ns4"
    }
    formdroite <- as.character(model[3]) ## la partie à droite du tilde
    formdroitebis <- gsub("\\bWCIE\\b", paste("(", new_expo, ")"), formdroite) # remplace "expo", par "(ns1+ns2+ns3+ns4)"

    new_formula <- as.formula(paste(as.character(model[2]),"~",formdroitebis))

    model_outcome <- glm(new_formula,family = binomial,data = new_outcome_date)
  }
  if (reg.type=="cox"){

  }
  return(list(model=model_outcome,data=data_expo_pred,mexpo=mexpo,call=new_formula))
}







