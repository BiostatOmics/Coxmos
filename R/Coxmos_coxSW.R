#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' coxSW
#' @description The `coxSW` function conducts a stepwise Cox regression analysis on survival data,
#' leveraging the capabilities of the `My.stepwise` R package. The primary objective of this function
#' is to identify the most significant predictors for survival data by iteratively adding or removing
#' predictors based on their statistical significance in the model. The resulting model is of class
#' "Coxmos" with an attribute model labeled as "coxSW".
#'
#' @details
#' The `coxSW` function employs a stepwise regression technique tailored for survival data. This
#' method is particularly beneficial when dealing with a plethora of predictors, and there's a
#' necessity to distill the model to its most impactful variables. The stepwise procedure can be
#' configured to operate in forward, backward, or a hybrid mode, contingent on the parameters
#' specified by the user.
#'
#' During the iterative process, variables are evaluated for inclusion or exclusion based on
#' predefined significance levels (`alpha_ENT` for entry and `alpha_OUT` for removal). This ensures
#' that the model retains only those predictors that meet the significance criteria, thereby
#' enhancing the model's interpretability and predictive power.
#'
#' Additionally, the function offers several preprocessing options, such as centering and scaling of
#' the predictor matrix, removal of variables with near-zero or zero variance, and the ability to
#' enforce the inclusion of specific variables in the model. These preprocessing steps are crucial
#' for ensuring the robustness and stability of the resulting Cox regression model.
#'
#' It's worth noting that the function is equipped to handle both numeric and binary categorical
#' predictors. However, it's imperative that categorical variables are appropriately transformed
#' into binary format before analysis. The outcome or response variable should comprise two columns:
#' "time" representing the survival time and "event" indicating the occurrence of the event of interest.
#'
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If
#' NULL, the number of columns of X matrix is selected. When MIN_EPV is not meet, the value will be
#' change automatically (default: NULL).
#' @param BACKWARDS Logical. If BACKWARDS = TRUE, backward strategy is performed (default: TRUE).
#' @param alpha_ENT Numeric. Maximum P-Value threshold for an ANOVA test when comparing a more complex model to a simpler model that includes a new variable. If the p-value is less than or equal to this threshold, the new variable is considered significantly important and will be added to the model (default: 0.10).
#' @param alpha_OUT Numeric. Minimum P-Value threshold for an ANOVA test when comparing a simpler model to a more complex model that excludes an existing variable. If the p-value is greater than or equal to this threshold, the existing variable is considered not significantly important and will be removed from the model (default: 0.15).
#' @param toKeep.sw Character vector. Name of variables in X to not be deleted by Step-wise
#' selection (default: NULL).
#' @param initialModel Character vector. Name of variables in X to include in the initial model
#' (default: NULL).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance
#' variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will
#' be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance
#' filtering (default: NULL).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant
#' variables/components in final cox model will be removed until all variables are significant by
#' forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "Coxmos" and model "coxSW". The class contains the following elements:
#'
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(data)}: normalized Y matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix
#'  }
#' \code{survival_model}: List of survival model information
#' \itemize{
#'  \item \code{fit}: coxph object.
#'  \item \code{AIC}: AIC of cox model.
#'  \item \code{BIC}: BIC of cox model.
#'  \item \code{lp}: linear predictors for train data.
#'  \item \code{coef}: Coefficients for cox model.
#'  \item \code{YChapeau}: Y Chapeau residuals.
#'  \item \code{Yresidus}: Y residuals.
#' }
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{nsv}: Variables removed by remove_non_significant if any.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{nz_coeffvar}: Variables removed by coefficient variation near zero.
#'
#' \code{removed_variables_correlation}: Variables removed by being high correlated with other
#' variables.
#'
#' \code{class}: Model class.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @references
#' \insertRef{Efroymson_SW}{Coxmos}
#' \insertRef{coxSW_CRAN}{Coxmos}
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' X <- X_proteomic[,1:10]
#' Y <- Y_proteomic
#' coxSW(X, Y, x.center = TRUE, x.scale = TRUE)

coxSW <- function(X, Y,
                  max.variables = NULL, BACKWARDS = TRUE,
                  alpha_ENT = 0.1, alpha_OUT = 0.15, toKeep.sw = NULL,
                  initialModel = NULL,
                  x.center = TRUE, x.scale = FALSE,
                  remove_near_zero_variance = TRUE, remove_zero_variance = FALSE, toKeep.zv = NULL,
                  remove_non_significant = FALSE, alpha = 0.05,
                  MIN_EPV = 5, returnData = TRUE, verbose = FALSE){

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  if(is.null(initialModel)){
    initialModel = "NULL"
  }

  #### Check values classes and ranges
  params_with_limits <- list("alpha_ENT" = alpha_ENT, "alpha_OUT" = alpha_OUT, "alpha" = alpha)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("MIN_EPV" = MIN_EPV)
  if(!is.null(max.variables)){
    numeric_params$max.variables <- max.variables
  }

  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = x.center, "x.scale" = x.scale,
                         #"y.center" = y.center, "y.scale" = y.scale,
                          "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                          "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose,
                          "BACKWARDS" = BACKWARDS)
  check_class(logical_params, class = "logical")

  #### Check rownames
  lst_check <- checkXY.rownames(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Check colnames in X for Illegal Chars (affect cox formulas)
  X <- checkColnamesIllegalChars(X)

  #### REQUIREMENTS
  checkX.colnames(X)
  checkY.colnames(Y)
  lst_check <- checkXY.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### ZERO VARIANCE - ALWAYS
  if(remove_near_zero_variance || remove_zero_variance){
    lst_dnz <- deleteZeroOrNearZeroVariance(X = X,
                                            remove_near_zero_variance = remove_near_zero_variance,
                                            remove_zero_variance = remove_zero_variance,
                                            toKeep.zv = toKeep.zv,
                                            freqCut = FREQ_CUT)
    X <- lst_dnz$X
    variablesDeleted <- lst_dnz$variablesDeleted
  }else{
    variablesDeleted <- NULL
  }

  #### COEF VARIATION
  lst_dnzc <- deleteNearZeroCoefficientOfVariation(X = X)
  X <- lst_dnzc$X
  variablesDeleted_cvar <- lst_dnzc$variablesDeleted

  #### MAX PREDICTORS
  if(is.null(max.variables)){
    max.variables <- ncol(X)
  }

  max.variables <- check.ncomp(X, max.variables, verbose = verbose)
  max.variables <- check.maxPredictors(X, Y, MIN_EPV, max.variables, verbose = verbose)

  #### SCALING
  lst_scale <- XY.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh
  data <- cbind(Xh, Yh)

  #### INITIALISING VARIABLES
  data.all <- NULL
  oneToDelete <- c("")
  lstMeetAssumption = NULL
  problem_flag = FALSE
  removed_variables = NULL
  removed_variables_cor = NULL

  in.variable = initialModel

  #### ###
  # MODEL #
  #### ###
  model <- tryCatch(
    # Specifying expression
    expr ={
      stepwise.coxph(Time = "time", Status = "event", variable.list = NULL,
                     sle = alpha_ENT, sls = alpha_OUT, startBackwards = BACKWARDS,
                     in.variable = in.variable, data = data, max.variables = max.variables,
                     verbose = verbose)
    },
    # Specifying error message
    error = function(e){
      message(paste0("coxSW: ", e))
      # invisible(gc())
      return(NA)
    }
  )

  ## STEPWISE MODEL IS A NOT-ROBUST COX MODEL
  ## make it robust to have the same P-VAL
  if(!all(is.na(model))){
    model <- cox(X = Xh[,names(model$coefficients), drop = FALSE], Y = Yh,
                 x.center = x.center, x.scale = x.scale,
                 alpha = alpha, MIN_EPV = MIN_EPV,
                 remove_near_zero_variance = FALSE, remove_zero_variance = FALSE, remove_non_significant = FALSE,
                 returnData = returnData, verbose = verbose)
  }

  # if all NA, returna NULL model
  if(all(is.na(model)) | all(is.null(model$survival_model))){
    problem_flag = TRUE
    survival_model = NULL
    func_call <- match.call()
    return(coxSW_class(list(X = list("data" = if(returnData) data.all else NA, "x.mean" = xmeans, "x.sd" = xsds),
                            Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                            survival_model = survival_model,
                            call = func_call,
                            X_input = if(returnData) X_original else NA,
                            Y_input = if(returnData) Y_original else NA,
                            alpha = alpha,
                            nsv = NULL,
                            nzv = variablesDeleted,
                            removed_variables_correlation = NULL,
                            class = pkg.env$coxSW,
                            time = time)))
  }

  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  lst_model <- removeNAorINFcoxmodel(model, data)
  coxph.sw <- lst_model$model$survival_model$fit
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  #### ### ###
  # CHECK SIG #
  #### ### ###

  # already check for the SW algorithm. If you want to be more restrictive,
  # change the enter/leave parameters.

  #### ### ### ### ### ### #
  # NUMBER OF VARIABLES >> #
  #### ### ### ### ### ### #

  # already check for the SW algorithm.

  #### ### ### ### #
  # if NO PROBLEMS #
  #### ### ### ### #

  data.all <- data #must be returned
  message("Final Model Reached!")

  removed_variables <- NULL
  if(isa(coxph.sw,"coxph")){

    #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
    if(remove_non_significant){
      if(all(c("time", "event") %in% colnames(data))){
        lst_rnsc <- removeNonSignificativeCox(cox = coxph.sw, alpha = alpha, cox_input = data[,colnames(data) %in% c(names(coxph.sw$coefficients),"time", "event")], time.value = NULL, event.value = NULL)
      }else{
        lst_rnsc <- removeNonSignificativeCox(cox = coxph.sw, alpha = alpha, cox_input = cbind(data[,colnames(data) %in% names(coxph.sw$coefficients)], Yh), time.value = NULL, event.value = NULL)
      }

      coxph.sw <- lst_rnsc$cox
      removed_variables <- lst_rnsc$removed_variables
    }

    survival_model <- getInfoCoxModel(coxph.sw)
  }else{
    survival_model <- NULL
    message("No variables have been selected for coxSW model. Returning a NULL Survival Model.")
  }

  func_call <- match.call()

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(coxSW_class(list(X = list("data" = if(returnData) data.all else NA, "x.mean" = xmeans, "x.sd" = xsds),
                        Y = list("data" = Yh, "y.mean" = ymeans, "y.sd" = ysds),
                        survival_model = survival_model,
                        call = func_call,
                        X_input = if(returnData) X_original else NA,
                        Y_input = if(returnData) Y_original else NA,
                        alpha = alpha,
                        nsv = removed_variables,
                        removed_variables_correlation = removed_variables_cor,
                        nzv = variablesDeleted,
                        nz_coeffvar = variablesDeleted_cvar,
                        class = pkg.env$coxSW,
                        time = time)))
}

#### ### ##
# METHODS #
#### ### ##

computeFirstModelSW <- function(Time = NULL, Status = NULL, variable.list,
                                in.variable = "NULL", data, max.variables = NULL,
                                startBackwards = FALSE, verbose = FALSE){
  #### ###
  # BACKWARDS CODE
  #### ###

  if(startBackwards){
    if(!all(in.variable == "NULL")){
      in.variable <- colnames(data)[colnames(data) %in% in.variable]
      if(length(in.variable)==0){
        stop("No variables were selected. Check that 'in.variable' parameter must have variables that belong to X data.")
      }
    }else{
      in.variable <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
    }

    # If X matrix has more variables than EPV and BACKWARDS...
    if(!is.null(max.variables) & length(in.variable) > max.variables){

      icox <- getIndividualCox(data = data[,colnames(data) %in% c(in.variable, "time", "event", "status")], time_var = "time", event_var = "event")
      in.variable <- rownames(icox[1:min(max.variables, length(rownames(icox))),])

      # we cannot compute a standard cox bc the p.val and coefficients are not well computed for HD models
      # we are using a FORCE = TRUE model
      # we should run a Individual COX model for starting

      aux_data <- as.data.frame(data[,colnames(data) %in% c(in.variable, "time", "event", "status"),drop = FALSE])

      f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(in.variable, collapse = " + ")))
      initial.model <- survival::coxph(formula = f, data = aux_data,
                                       method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)

      # REMOVE NA-PVAL VARIABLES
      # p_val could be NA for some variables (if NA, then change it to P-VAL=1)
      # DO IT ALWAYS, we do not want problems in COX models
      if(all(c("time", "event") %in% colnames(aux_data))){
        lst_model <- removeNAorINFcoxmodel(model = initial.model, data = aux_data, time.value = NULL, event.value = NULL)
      }else{
        lst_model <- removeNAorINFcoxmodel(model = initial.model, data = cbind(cbind(aux_data, Time), Status), time.value = NULL, event.value = NULL)
      }
      initial.model <- lst_model$model

    }else{
      # we CAN compute a standard cox because the data has less variables than observations
      aux_data <- as.data.frame(data[,colnames(data) %in% c(in.variable, "time", "event", "status"),drop = FALSE])

      f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(in.variable, collapse = " + ")))
      initial.model <- survival::coxph(formula = f, data = aux_data,
                                       method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)

      lst_model <- removeNAorINFcoxmodel(initial.model, aux_data)
      initial.model <- lst_model$model
    }

  #### ###
  # FORWARDS CODE
  #### ###
  }else{

    # forward selection - one variable
    # we also need to compute individual cox for the first selection

    if(!all(in.variable == "NULL")){
      in.variable <- colnames(data)[colnames(data) %in% in.variable]
      if(length(in.variable)==0){
        stop("Any variable was selected. Check that 'initialModel' variables belong to X data.")
      }
    }else{
      in.variable <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
    }

    # we cannot compute a standard cox bc the p.val and coefficients are not well computed for HD models
    # we are using a FORCE = TRUE model
    # we should run a Individual COX model for starting

    icox <- getIndividualCox(data[,colnames(data) %in% c(in.variable, "time", "event", "status")])
    in.variable <- rownames(icox[1,])

    aux_data <- as.data.frame(data[,colnames(data) %in% c(in.variable, "time", "event", "status"),drop = FALSE])
    initial.model <- survival::coxph(as.formula(paste("Surv(", Time, ", ", Status, ") ~ .")), data = aux_data,
                                     method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)
    lst_model <- removeNAorINFcoxmodel(initial.model, aux_data)
    initial.model <- lst_model$model
  }

  return(initial.model)
}

stepwise.coxph <- function(Time = NULL, Status = NULL, variable.list,
                           in.variable = "NULL", data, max.variables = NULL,
                           sle = 0.15, sls = 0.15,
                           startBackwards = FALSE, verbose = FALSE){

  if(is.null(variable.list)){
    variable.list <- colnames(data)[!colnames(data) %in% c("time", "event", "status")]
  }

  initial.model <- computeFirstModelSW(Time = Time, Status = Status, variable.list,
                                       in.variable = in.variable, data, max.variables = max.variables,
                                       startBackwards = startBackwards, verbose = verbose)

  # save model
  temp.model <- initial.model

  #Stepwise algorithm model based on My.stepwise
  i <- 0
  break.rule <- TRUE
  broke = FALSE
  iter = 0

  last_removed = ""
  historic = NULL
  last_AIC <- getInfoCoxModel(temp.model)$AIC #minimum, best

  # WHILE after we got the best model
  while(break.rule){

    if(verbose){
      message(paste0("Stepwise loop: ", iter))
      #print(temp.model)
      print(head(temp.model$coefficients))
      print(paste0("AIC: ", last_AIC))
    }

    iter = iter + 1
    i <- i + 1

    # choose the list of variable than can enter to the model in this iteration
    if(i == 1){
      variable.list2 <- setdiff(variable.list, all.vars(temp.model$formula))
    }else{
      # after first iteration, cannot enter out.x from last iteration
      variable.list2 <- setdiff(variable.list, c(all.vars(temp.model$formula), out.x))
      out.x <- NULL
    }

    #### ### ### ##
    # CHECK ENTER #
    #### ### ### ##
    if(length(variable.list2) != 0 & length(names(temp.model$coefficients)) < max.variables){
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      AIC.value <- NULL
      for(k in 1:length(variable.list2)){

        if(verbose){
          message(paste0("Testing Entering: ",k, "/", length(variable.list2), " - iteration ", iter))
        }
        utils::flush.console()

        model <- tryCatch({
          # generate new model
          coeff <- names(coefficients(temp.model))
          f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(c(coeff, variable.list2[k]), collapse = " + ")))
          survival::coxph(formula = f, data = as.data.frame(data),
                          method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)
          },
          # save the error
          error=function(e){
            # Choose a return value in case of error
            message(paste0("COXSW: ", e$message))
            return(NA)
          })

        if(!isa(model,"coxph")){# an error happens, next variable
          next
        }

        if(length(model$coefficients) > 1){
          if(sum(is.na(model$coefficients)) != 0){
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
            AIC.value[k] <- 1000 #max value
          }else{
            # compare temporal model with current iteration model
            # different between models is significant if p.val lesser than sig
            anova.pvalue[k] <- anova(temp.model, model)[2,"Pr(>|Chi|)"] #model sig
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients),"Pr(>|z|)"] #last variable cox significance
            AIC.value[k] <- getInfoCoxModel(model)$AIC #aic model
          }
        }
      }
      # list of pvalues complete
      # it uses ANOVA P-VALUE for selecting the variables are going to enter
      # select indexes of variables whose cox coeff is less than 0.9 or not NA
      variable.list2.1 <- variable.list2[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      anova.pvalue2 <- anova.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      mv.pvalue2 <- mv.pvalue[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]
      AIC.value2 <- AIC.value[mv.pvalue <= 0.9 & !is.na(mv.pvalue)]

      # MinAnova Value have to be lesser than SLE
      enter.x <- variable.list2.1[anova.pvalue2 == min(anova.pvalue2, na.rm = TRUE) & anova.pvalue2 <= sle]
      enter.x.AIC <- AIC.value2[anova.pvalue2 == min(anova.pvalue2, na.rm = TRUE) & anova.pvalue2 <= sle]

      #select which variable enters if any or AIC is lesser than last_AIC
      if(length(enter.x)==0 || last_AIC < enter.x.AIC){
        #no variables able to enter into the model
        enter.x <- NULL

      }else{

        # if last removed want to enter again, do not allow it
        if(last_removed %in% enter.x){
          enter.x <- enter.x[-which(last_removed==enter.x)]
        }

        # But, if there is a tie, select that variable with the less p-val
        wald.p <- mv.pvalue2[anova.pvalue2 == min(anova.pvalue2, na.rm = TRUE) & anova.pvalue2 <= sle]
        if(length(setdiff(enter.x, NA)) != 0){

          # if another tie, select the first one
          if(length(enter.x) > 1){
            enter.x <- enter.x[which.min(wald.p)]
          }

          if(verbose){
            message(paste0("ENTER: ", paste0(enter.x, collapse = ", "), " AIC: ", enter.x.AIC, " LAST_AIC: ", last_AIC))
            historic <- c(historic, paste0("ENTER: ", paste0(enter.x, collapse = ", "), " AIC: ", enter.x.AIC, " LAST_AIC: ", last_AIC))
          }

          coeff <- names(coefficients(temp.model))
          f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(c(coeff, enter.x), collapse = " + ")))
          temp.model <- survival::coxph(formula = f, data = as.data.frame(data),
                                        method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)

          # update bc a new model is used
          last_AIC <- getInfoCoxModel(temp.model)$AIC #minimum, best
        }
      }

    }else{
      enter.x <- NULL
    }

    #### ### ### ##
    # CHECK LEAVE #
    #### ### ### ##

    out.x  <- NULL
    # at least one variable can leave the model
    if(length(names(temp.model$coefficients)) > 1){
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients),
                                c(enter.x))

      if(length(variable.list3) != 0){
        anova.pvalue <- NULL
        mv.pvalue <- NULL
        AIC.value <- NULL
        for(k in 1:length(variable.list3)){

          if(verbose){
            message(paste0("Testing Leaving: ",k, "/", length(variable.list3), " - iteration ", iter))
          }

          model <- tryCatch({

            # generate new model
            coeff <- names(coefficients(temp.model))
            coeff <- coeff[-which(coeff == variable.list3[k])]
            f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(coeff, collapse = " + ")))
            survival::coxph(formula = f, data = as.data.frame(data),
                            method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)
            },

            # save the error
            error=function(e){
              # Choose a return value in case of error
              message(paste0("COXSW: ", e$message))
              return(NA)
            })

          if(!isa(model,"coxph")){# an error happend, next variable
            next
          }

          if(length(model$coefficients) >= 1){
            if(sum(is.na(model$coefficients)) != 0){
              anova.pvalue[k] <- 1
              mv.pvalue[k] <- getPvalFromCox(model)[variable.list3[k]]
              AIC.value[k] <- 1000
            }else{
              # compare temporal model with current iteration model
              # different between models is significant if p.val lesser than sig
              anova.pvalue[k] <- anova(temp.model, model)[2,"Pr(>|Chi|)"] #model sig
              mv.pvalue[k] <- getPvalFromCox(temp.model)[variable.list3[k]] #significant coeff from TEMPORAL model bc variable exist only there
              AIC.value[k] <- getInfoCoxModel(model)$AIC #aic model
            }
          }
        }

        if(is.null(anova.pvalue)){
          message("LEAVE: anova is null")
          broke = TRUE
          break.rule = FALSE
          next
        }

        # it uses ANOVA P-VALUE for selecting the variables are going to leave
        # in this case, max P-val bc is the lesser important variable
        # Valor mínimo de Anova y Anova menor que SLE, pero creo que debería ser el coeficiente del modelo de COX, no del modelo comparativo.
        out.x <- variable.list3[anova.pvalue == max(anova.pvalue, na.rm = TRUE) & anova.pvalue > sls]
        out.x.AIC <- AIC.value[anova.pvalue == max(anova.pvalue, na.rm = TRUE) & anova.pvalue > sls]
        out.x <- setdiff(out.x, NA)
        out.x.AIC <- setdiff(out.x.AIC, NA)
        mv.pvalue <- mv.pvalue[which(variable.list3 %in% out.x)] # select just pval selected

        # mod: some out.x and out.x.AIC lesser than current value
        if(length(out.x) != 0 && last_AIC >= out.x.AIC){

          # if tie, then uses the pval
          if(length(out.x) > 1){
            # if another tie, select the first one
            if(length(out.x) > 1){
              out.x <- out.x[which.max(mv.pvalue)]
              mv.pvalue <- mv.pvalue[which.max(mv.pvalue)]
            }
          }else{
            out.x <- out.x
            mv.pvalue <- mv.pvalue
          }

          #check enter==out
          equal = FALSE
          if(length(enter.x)>0){
            equal <- enter.x == out.x
          }
          if(!equal){ #if the P-Val is grater than sls then remove, else, althouth the anova wants to remove, dont do it

            if(verbose){
              message(paste0("OUT: ", paste0(out.x, collapse = ", "), " AIC: ", out.x.AIC, " LAST_AIC: ", last_AIC))
              historic <- c(historic, paste0("OUT: ", paste0(out.x, collapse = ", "), " AIC: ", out.x.AIC, " LAST_AIC: ", last_AIC))
            }

            #temp.model <- update(temp.model, as.formula(paste(". ~ . - ", out.x.column, sep = "")))
            coeff <- names(coefficients(temp.model))
            coeff <- coeff[-which(coeff == out.x)]
            f <- as.formula(paste0("survival::Surv(", Time, ", ", Status, ") ~ ", paste0(coeff, collapse = " + ")))
            temp.model <- survival::coxph(formula = f, data = as.data.frame(data),
                            method = "efron", model = TRUE, singular.ok = TRUE, x = TRUE)

          }else{
            # enter.x is the same as out.x (no-sense)?
            message(paste0("Variable selected to enter the model and to leave in the same iteration is the same (", out.x, ")"))
            message("Performing one more iteration.")
          }

        # any p-val to leave
        }else{
          out.x <- NULL
        }

      # any variable to check for leaving
      }else{
        out.x <- NULL
      }
    }

    last_AIC <- getInfoCoxModel(temp.model)$AIC #minimum, best

    if((length(enter.x) + length(out.x)) == 0){
      final.model <- temp.model
      break.rule <- FALSE
    }
    enter.x <- NULL
  }

  if(broke){
    final.model <- temp.model #the model has not be updated because any variable enter to the model
  }

  return(final.model)
}

#### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Get the real name for the cox model variable that it is used in the matrix/data.frame
# because for categorical variables, cox add the level at the end and we can not do any match.
# We are not paying attention to different order, but yes for multiple interactions +2
#### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# data.frame/matrix  -  var_cox  -  interaction symbol #

getRealNameVarFromCoxModel <- function(x_sw_train, var, symbol = ":"){
  real.symbol <- symbol

  if(symbol=="."){
    symbol = "\\."
  }

  if(!var %in% colnames(x_sw_train)){

    if(isInteraction(var)){#check ifeach individual_variable is a factor
      n <- strsplit(var, symbol)[[1]]
      in_vector <- sapply(n, function(x){x %in% colnames(x_sw_train)})
      names_updated <- sapply(names(in_vector[in_vector==FALSE]), function(x){x <- substr(x, 0, nchar(x)-1)})
      names(in_vector)[in_vector==FALSE] <- names_updated

      res_var <- paste(names(in_vector), collapse = real.symbol)

    }else{#not interaction
      res_var <- substr(var, 0, nchar(var)-1)
      if(res_var %in% colnames(x_sw_train)){
        return(res_var)
      }

      res_var <- paste0("`",var,"`")
      if(res_var %in% colnames(x_sw_train)){
        return(res_var)
      }
    }
  }else{#var in colnames
    res_var = var
  }

  if(res_var %in% colnames(x_sw_train)){
    return(res_var)
  }else{
    return(NULL)
  }
}

isInteraction <- function(cn, symbol = ":"){
  if(symbol==".")
    symbol = "\\."
  if(length(grep(symbol,cn))>0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

deleteVariablesCox <- function(x_sw_train, x_sw_test=NULL, td, interactions = FALSE, symbol = ":"){
  real.symbol = symbol

  if(symbol==".")
    symbol = "\\."

  if(length(td)>0){
    index <- match(td,colnames(x_sw_train)) #perfect matchs
    index <- index[!is.na(index)]
    #ifit is NOT an interaction, delete all the interactions where the variable appears
    if(interactions){
      subtd <- unlist(lapply(td, function(x){if(!isInteraction(x))return(x)})) #for individual variables...
      inter <- NA
      splited <- lapply(colnames(x_sw_train), function(x){
        if(isInteraction(x))
          strsplit(x, symbol)})
      for(cn in splited){
        if(is.null(cn))
          next
        if(cn[[1]][1] %in% subtd | cn[[1]][2] %in% subtd){
          inter <- c(inter, paste0(cn[[1]][1], real.symbol, cn[[1]][2]))
        }
      }
      inter <- inter[!is.na(inter)]
      index <- unique(c(index, match(inter, colnames(x_sw_train))))
    }
    if(length(index)>0){
      if(class(x_sw_train)[1] %in% "data.frame"){
        x_sw_train[,index] <- NULL
        if(!is.null(x_sw_test))
          x_sw_test[,index] <- NULL
      }else if(class(x_sw_train)[1] %in% "matrix"){
        x_sw_train <- x_sw_train[,-index]
        if(!is.null(x_sw_test))
          x_sw_test <- x_sw_test[,-index]
      }
    }
  }
  return(list(cox.train = x_sw_train, cox.test = x_sw_test))
}

#my own fit_zph function
plotZPH <- function(fit_zph, oneToDelete, df=3){

  #DFCALLS
  x <- y <- NULL

  #https://rpubs.com/Cristina_Gil/Regr_no_lineal
  aux <- as.data.frame(cbind(fit_zph$time, fit_zph$y[,oneToDelete]))
  colnames(aux) <- c("x","y")
  ggp <- ggplot(aux, aes(x=x, y=y)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ splines::ns(x, df = df, intercept = TRUE), color = "red", #should be a base R package
                se = TRUE, level = 0.95) +
    xlab(label = "Time") + ylab(label = paste0("Beta(t) for ", oneToDelete))
  return(ggp)
}

### ## ##
# CLASS #
### ## ##

coxSW_class = function(cox_model, ...){
  model = structure(cox_model, class = pkg.env$model_class,
                    model = pkg.env$coxSW)
  return(model)
}
