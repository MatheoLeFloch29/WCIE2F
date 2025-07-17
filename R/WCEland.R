#' Estimation of the Weighted Cumulative Effect with Two-Level Bootstrap \code{WCEland}
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
#' @importFrom survival coxph Surv
#' @import ggplot2
#'
#' @author un super beau gosse
#'
#' @seealso
#' \code{\link{summary.WCEland}}
#' \code{\link{WCIE2F}}
#' \code{\link{doOneBootWCIE}}
#'
#' @examples
#' \dontrun{
#' library(lcmm)
#' # Fit linear mixed-effects model for exposure
#' m_expo <- hlme(E ~ X + time,
#'                random = ~ time,
#'                subject = "ID",
#'                data = data_expo,
#'                returndata = TRUE)
#'
#' # Estimate weighted cumulative exposure effect in Cox model
#' m_outcome <- WCEland(mexpo = m_expo,
#'                      var.time = "time",
#'                      time.frame = c(-10, 0, 1),
#'                      weightbasis = "NS",
#'                      knots = 1,
#'                      data = data_outcome,
#'                      reg.type = "cox",
#'                      model = Surv(time, Y) ~ X + WCIE,
#'                      n_boot = 500)
#'
#'  summary(m_outcome)
#' }
#'
#' @references
#' Maud Wagner et al. “Time-varying associations between an exposure history and a subsequent health
#' outcome : a landmark approach to identify critical windows”. In : BMC Med Res Methodol (2021).
#' doi : 10.1186/s12874-021-01403-w
#'
#' @name WCEland
#'
#'
#' @export
WCEland <- function(mexpo,var.time, time.frame, weightbasis="NS", knots=NULL,knots.vector=NULL,
                   data, reg.type="RL", model,n_boot=500){

  ptm <- proc.time()

  if(is.null(mexpo$data)==T) stop("The argument mexpo need to specify returndata = T")
  if(!inherits(mexpo,"hlme")) stop("The argument mexpo must be a hlme object")
  if (is.null(data)==T) stop("the argument outcome_data is missing")
  if (is.null(model)==T) stop("the argument outcomeformula is missing")
  if (abs(time.frame[1])>abs(floor(min(mexpo$data[var.time])))) stop("The first time.frame argument must be equal to or greater than the minimum time value, rounded down to the nearest integer.")
  if (abs(time.frame[2])<abs(ceiling(max(mexpo$data[var.time])))) stop("the second argument time.frame must be equal or less then the maximum time value, rounded up to the nearest integer")
  if(is.null(knots)==T&is.null(knots.vector)==T) stop("You must have to specify knots or knots.vector")


  # Extract the mean and variance-covariance matrix of the estimated parameters from the exposure model
  mu <- as.matrix(estimates(mexpo))
  Sigma <- as.matrix(VarCov(mexpo))

  # To convert from the Cholesky decomposition to the variance:

  # Generate new bootstrap parameter values
  boot_params <- MASS::mvrnorm(n = n_boot, mu = mu, Sigma = Sigma)

  NPROB = mexpo$N[1]
  NEF = mexpo$N[2]
  NVC = mexpo$N[3]
  idiag0 = mexpo$idiag
  nea0 = sum(mexpo$idea0)

  ######## Convert from Cholesky to variance-covariance matrix ##################

  if(idiag0==0 & NVC>0){

    for (i in 1:nrow(boot_params)) {
      # Extract the Cholesky parameters of sample i
      Cholesky <- boot_params[i, (NPROB+NEF+1):(NPROB+NEF+NVC)]

      # Build the upper triangular matrix U
      U <- matrix(0, nrow = nea0, ncol = nea0)
      U[upper.tri(U, diag = TRUE)] <- Cholesky

      # Transform into the variance-covariance matrix
      varcov <- t(U) %*% U

      # Replace in boot_params
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

  # save all the bootstrap in boot_result

  boot_results <- lapply(1:n_boot, function(i) {
    doOneBootWCIE(i = i,boot_params=boot_params,
                  times = time.frame,mexpo=mexpo,knots.vector=knots.vector,
                  var.time = var.time,weightbasis = weightbasis,knots = knots,
                  data = data, reg.type = reg.type, model = model)}
  )

  # Replace the for loop with a replicate that uses the doOneBoot function and returns:

  # a list of estimated parameters for each bootstrap

  # a list of variance-covariance matrices

  # parameters mean for n_boot bootstrap
  boot_est <- sapply(boot_results, function(x) x[[1]])

  ## Intra-individual variability
  var_intra <- Reduce("+",lapply(boot_results, function(x) x[[2]]))/n_boot


  ## Inter-individual variability (formula for bootsrap from rubin)
  #(M+1)/(M(M-1))sum(teta-mean(teta)²)
  mean_boot_est <- apply(boot_est, 1, function(x) mean(x))

  somme_var_var <- Reduce("+", lapply(1:ncol(boot_est), function(i) {
    (boot_est[, i] - mean_boot_est) %*% t(boot_est[, i] - mean_boot_est)# do the sum of teta-mean(teta)²
  }))

  var_inter <- ((n_boot + 1)/(n_boot*(n_boot-1)))*somme_var_var #do the (M+1)/(M(M-1))
  rownames(var_inter)<-colnames(var_intra)


  #######################################################
  ########## end bootstrap ##############################
  #######################################################

  ####################################################
  ############## parameters estimations ##############
  ####################################################

  # Var tot parameters
  var_tot <- var_inter + var_intra

  # Extract the diagonal to obtain the bootstrap-corrected variance.
  var_estim<-as.matrix(diag(var_tot))

  # SE calcul
  parameters_var<-cbind(mean_boot_est, sqrt(var_estim))
  colnames(parameters_var) <- c("Estimate","Se")
  parameters_var <- as.data.frame(parameters_var)

  #IC calcul
  parameters_var$con.low <- parameters_var$Estimate - qnorm(0.975)  * (parameters_var$Se)
  parameters_var$conf.high <- parameters_var$Estimate + qnorm(0.975)  * (parameters_var$Se)

  # statistical test value
  parameters_var$z_value <- parameters_var$Estimate / (parameters_var$Se)

  #p-valeur (for normal law)
  parameters_var$p_values <- 2 * (1 - pnorm(abs(parameters_var$z_value)))
  #}



  #######################################################################################
  #################### Calculation of the exposure effects ##############################
  #######################################################################################


  # Run a model to retrieve the updated data (within the new time window)

  #Get the updated exposure model data along with the outcome data

  #Save the call of the exposure model

  #Use these for calculating the effects of past exposure over time

  WCIE <- WCIE2F(mexpo = mexpo,var.time = var.time, times = time.frame,
                         weightbasis = weightbasis, knots = knots,knots.vector=knots.vector,
                         data = data, reg.type = reg.type, model = model)

  new_data <- WCIE$data_expo #data exposition
  data_outcome <- WCIE$data_outcome #data outcome

  n <- n_distinct(data_outcome) # subject nb
  n_obs <- nrow(data_outcome) # nb total observation
  nb_subject_delete <- n_distinct(data)-n # Number of subjects removed
  nb_obs_delete <- nrow(data)-n_obs # Number of observations removed

  mean_AIC <- mean(sapply(boot_results, function(x) x$AIC)) #mean AIC bootstrap
  mean_loglike <- mean(sapply(boot_results, function(x) x$loglike)) #mean loglikelihood bootstrap

  # Put a condition based on the outcome type
  if(reg.type=="logistic"|reg.type=="cox"){

    ##################### Detect if there are any interactions in the model #####################

    # Identify all interactions containing "WCIEX" (either before or after the colon :)
    var <- grep("WCIE[0-9]+:", rownames(parameters_var), value = TRUE)
    var <- c(var, grep(":WCIE[0-9]+$", rownames(parameters_var), value = TRUE))

    ##################################################################
    ########### If there is no interaction, then do this. ############
    ##################################################################

    # new matrice spline to compile the effect
    data_splines1 <- data.frame(unique(new_data[var.time]))

    data_splines <- seq(from=time.frame[1],to=time.frame[2],by=time.frame[3]) # sequence de mesure dans la fenêtre choisis

    ## splines recompile with the same parameters than put in the WCEland function
    new_splines <- as.matrix(ns(unlist(data_splines),knots = WCIE$splines.quantiles,
                                Boundary.knots = WCIE$boundary.quantiles,
                                intercept = TRUE))
    # rename WCIE1:WCIEk
    colnames(new_splines) <- paste0("WCIE", 1:ncol(new_splines))

    # Loop first over the WCIE terms without interaction.

    # Calculate the total effect using a loop.
    effect <- data.frame(time=data_splines,eff_expo=0)

    effect <- effect %>% arrange(effect[1])
    # If there are no interactions, perform a single loop over the WCIE terms only.
    for (s in 1:dim(new_splines)[2]){
      effect$eff_expo <- effect$eff_expo + parameters_var[paste0("WCIE",s),"Estimate"]*new_splines[,s]
    }

    #################################################################
    ###### If interactions are present, proceed as follows : ########
    #################################################################

    if(length(var)>0){

      # variable present in the model
      vars_model <- all.vars(formula(model))
      # Common variables between the two
      vars_communes <- intersect(names(data_cov_int),vars_model)
      # recreate the dataset with the same values

      # for model.matrix to work, there must be at least two rows in data_cov_int
      # so create a second dummy row and delete it after retrieving the correct variable name

      # Extract the categorical variables from the initial model
      vars_categorielle <- names(Filter(is.factor, model$model))

      # Apply the same levels to the categorical variables in the new dataset
      #(must avoid using as.factor directly in the model)
      for (var_cat in vars_categorielle) {

        if (var_cat %in% names(data_cov_int)) {
          niveaux <- levels(model$model[[var_cat]])  # Retrieve the factor levels from the model
          data_cov_int[[var_cat]] <- factor(data_cov_int[[var_cat]], levels = niveaux)  # applicate
        }
      }

      # Generate the new matrix with the same variables as in the model
      new_data_cov_int <- model.matrix(as.formula(paste("~", paste(colnames(data_cov_int), collapse = " + "))), data = data_cov_int)[,-1]
      new_data_cov_int <- as.data.frame(t(new_data_cov_int))

      # If a variable, keep new_data_cov_int the same as the original
      if(ncol(data_cov_int)==1 & ncol(new_data_cov_int)==1){
        new_data_cov_int<-data_cov_int
      }

      doOneBoot_effectInt<-function(i){

        # Loop over all WCIE variables without interaction
        effect <- data.frame(time=data_splines,eff_expo=0)
        effect <- effect %>% arrange(effect[1])
        # If there is no interaction, just loop over the WCIE variables
        for (s in 1:dim(new_splines)[2]){
          effect$eff_expo <- effect$eff_expo + boot_est[i,paste0("WCIE",s)]*new_splines[,s]
        }

        ######### + ad the interactions #############

        for(q in 1:ncol(new_data_cov_int)){

          # Identify the position/index of the interactions
          repaire <- grep(paste0(colnames(new_data_cov_int)[q],"."), colnames(boot_est),value = T)
          repaire <- c(repaire, grep(paste0(".",colnames(new_data_cov_int)[q]), colnames(boot_est),value = T))

          # If there is 1 variable in interaction, convert its single value back to numeric
          # in order to use new_data_cov_int correctly
          p<-1 # If multiple covariates are in interaction, reset the spline counter to 1
          if (dim(new_data_cov_int)[2]==1) {
            new_data_cov_int2 <- data_cov_int[1,1]
            for (g in 1:length(repaire)) {
              #ad value
              effect$eff_expo <- effect$eff_expo + boot_est[i,repaire[g]] * new_splines[,p] * new_data_cov_int2[q]
              p <- ifelse(p == (knots + 1), 1, (p + 1))
            }
          }else{
            for (g in 1:length(repaire)) {
              #ad value
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


    ### If no interaction, then the variance of the effects is calculated as follows:

    if(length(var)==0){

      # 1) Retrieve the variance-covariance matrix of interest
      # Store the results in the matrix
      eff_varCov_tot <- var_tot[grepl("WCIE\\d+", rownames(var_tot)),
                                grepl("WCIE\\d+", colnames(var_tot))]

      # 2) Calculate the variance for each time point: v(w(t)) = B(t)' * V(theta) * B(t)
      effect <- as.data.frame(effect)
      effect$var_eff <- 0
      for(z in 1:nrow(effect)){
        effect$var_eff[z] <- t(new_splines[z,]) %*% eff_varCov_tot %*% new_splines[z,]
      }



    }

    ### If interactions are present, then do this: (to be verified as not certain)
    if(length(var)>0){

      # 1) Extract the relevant variance-covariance matrix (including interactions)

      # Store the results in the matrix
      eff_varCov_tot <- var_tot[grepl("WCIE\\d+", rownames(var_tot)),
                                grepl("WCIE\\d+", colnames(var_tot))]

      # 2) Calculate the variance for each time point v(w(t)) = B(t)' * v(theta) * B(t) taking interactions into account

      # Create the matrix of splines multiplied by interactions (as many times as there are interactions) B(t)
      # code to review (not certain but seems okay)

      # Otherwise, if there is a categorical variable, create the splines*interaction matrix as follows
      new_matheo<-matrix(0,nrow = nrow(new_splines),ncol = nrow(eff_varCov_tot))
      colnames(new_matheo) <- colnames(eff_varCov_tot)

      # Add the spline values at the corresponding positions
      for (k in 1:ncol(new_splines)) {
        for (f in 1:ncol(new_matheo)) {
          if(grepl(colnames(new_splines)[k], colnames(new_matheo)[f])){
            new_matheo[,f] <- new_splines[,k]
          }
        }
      }
      # Multiply the spline values by the covariates when there is an interaction
      for (k in 1:ncol(new_data_cov_int)) {
        for (f in 1:ncol(new_matheo)) {
          if(grepl(colnames(new_data_cov_int)[k], colnames(new_matheo)[f])){
            new_matheo[,f] <- new_matheo[,f]*new_data_cov_int[,k]
          }
        }
      }


      # Calculate V(w(t))
      effect <- as.data.frame(effect)
      effect$var_eff <- 0
      for(z in 1:nrow(effect)){
        effect$var_eff[z] <- t(new_matheo[z,]) %*% eff_varCov_tot %*% new_matheo[z,]
      }

    }



    # Calculate the confidence interval of the effects

    effect$var_eff <- sqrt(effect$var_eff)
    colnames(effect) <- c("Time","Effect","sd")
    effect <- as.data.frame(effect)

    # Calculate the confidence intervals (CI)
    effect$conf.low <- effect$Effect - qnorm(0.975)  * (effect$sd)
    effect$conf.high <- effect$Effect + qnorm(0.975)  * (effect$sd)


    # Calculate the average effect from 0 to -T
    # Compute the mean effect w_bar = (1 / (T + 1)) * sum(w(u))
    mean_effect <- mean(effect$Effect)
    #mean_effect <- mean(effect[[2]], na.rm = T)

    #################################################################################
    ############## Calculate the variance of the average effect #####################
    #################################################################################

    # if no interactions
    if (length(var)==0) {

    means_splines <- apply(new_splines,2,FUN = mean)  #1/T+1(B(t))

    # Calculate the intra-bootstrap variance of the average effect using the Delta method
    var_intra_boot <- numeric(length(boot_results))
    for (h in 1:length(boot_results)) {
      mat_vcov <- boot_results[[h]][2]$V  # Vcov matrix for each bootstrap sample
      wcie_names <- grep("^WCIE", colnames(mat_vcov), value = TRUE) #keep the WCIE var
      d <- mat_vcov[wcie_names,wcie_names] # variance parameters of wcie
      var_intra_boot[h] <- t(means_splines) %*% d %*% (means_splines)
    }
    var_intra_mean_effect <- mean(var_intra_boot)


    # Calculate the inter-bootstrap variance of the average effect using Rubin's formula
    #(M+1)/(M(M-1))sum(Xbar-mean(Xbar)²)

    data_splines <- seq(from=time.frame[1],to=time.frame[2],by=time.frame[3]) # sequence de mesure dans la fenêtre choisis

    ## splines recompile with the same parameters than put in the WCEland function
    new_splines <- as.matrix(ns(unlist(data_splines),knots = WCIE$splines.quantiles,
                                Boundary.knots = WCIE$boundary.quantiles,
                                intercept = TRUE))
    # rename WCIE1:WCIEk
    colnames(new_splines) <- paste0("WCIE", 1:ncol(new_splines))

    effect_mean <- data.frame(time=data_splines,eff_expo=0)
    effect_list <- replicate(length(boot_results), effect_mean, simplify = FALSE)
    # Perform a first loop over all bootstrap estimates
    for (h in 1:length(boot_results)) {

      # Calculate the total effect using a loop
      effect_mean <- data.frame(time=data_splines,eff_expo=0)
      effect_mean <- effect_mean %>% arrange(effect_mean[1])
      # If there are no interactions, just loop over the WCIE variables
      for (s in 1:dim(new_splines)[2]){
        effect_mean$eff_expo <- effect_mean$eff_expo + boot_results[[h]]$coef[paste0("WCIE",s)]*new_splines[,s]
      }
    effect_list[[h]]$eff_expo <- effect_mean$eff_expo
    }

    mean_effect_boot <- sapply(effect_list, function(df) mean(df$eff_expo))

    somme_var_var_effect <- Reduce("+", lapply(1:ncol(boot_est), function(i) {
      (mean_effect_boot[i] - mean(mean_effect_boot)) %*% t(mean_effect_boot[i] - mean(mean_effect_boot))# do the sum of wbar-mean(wbar)²
    }))

    var_inter_mean_effect <- ((n_boot + 1)/(n_boot*(n_boot-1)))*somme_var_var_effect #applicate the (M+1)/(M(M-1))
    real_mean_var_effect <- var_intra_mean_effect + as.numeric(var_inter_mean_effect)


    ###### Simple version to calculate the variance of the mean effect using the delta method
    # real_mean_var_effect <- t(means_splines) %*% eff_varCov_tot %*% means_splines  # v(w̄)
    # This estimate is almost similar to the one obtained by the delta method combined with Rubin's formula

    mean_variable_effect <- sqrt(real_mean_var_effect)
    }

    # Add a condition to check if there are interactions


    ###########################################################
    ############ do the graph if is ok ########################
    ###########################################################


    ### BOX_PLOT

    # Effect - Intercept
    graph_effect <- ggplot(effect, aes(x=as.factor(Time), y=Effect, group = 1)) +
      geom_line(color="darkgray", linewidth=1)+  # Replace the boxplot with a line plot
      geom_line(aes(y=conf.low), linetype="dashed", color="black") +  # Add a dashed lower boundary line
      geom_line(aes(y=conf.high), linetype="dashed", color="black") +  # Add a dashed upper boundary line
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high), fill="gray", alpha=0.3, linetype="dashed") +  # ad IC in dash
      theme(legend.position = "none") +
      xlab("Time preceding landmark") +
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
                call=WCIE$call,
                knots.quantile=WCIE$splines.quantiles,V=var_tot,var.time=var.time,AIC=mean_AIC
                ,loglike=mean_loglike,n=n,nb.subj.del=nb_subject_delete,
                time.processing=cost)

  class(result) <- "WCEland"
  return(result)
}




