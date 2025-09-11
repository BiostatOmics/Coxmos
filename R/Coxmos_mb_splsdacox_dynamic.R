#### ### ##
# METHODS #
#### ### ##

#' MB.sPLS-DACOX
#' @description
#' The MB.sPLS-DACOX function conducts a multi-block sparse partial least squares discriminant
#' analysis Cox (MB.sPLS-DACOX) using a dynamic variable selection approach. This analysis is
#' particularly suited for high-dimensional datasets where the goal is to identify the relationship
#' between explanatory variables and survival outcomes. The function outputs a model of class
#' "Coxmos" with an attribute labeled "MB.sPLS-DACOX".
#'
#' @details
#' The MB.sPLS-DACOX methodology is designed to handle multi-block datasets, where each block
#' represents a set of related variables. By employing a sparse partial least squares approach, the
#' function efficiently selects relevant variables from each block, ensuring that the final model is
#' both interpretable and predictive. The Cox proportional hazards model is then applied to the
#' selected variables to assess their association with survival outcomes.
#'
#' The function offers flexibility in terms of parameter tuning. For instance, users can specify the
#' number of latent components to compute, the range of variables to consider for optimal selection,
#' and the evaluation metric (either AUC or c-index). Additionally, data preprocessing options are
#' available, such as centering and scaling of the explanatory variables, and removal of variables
#' with near-zero or zero variance.
#'
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event
#' observations.
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 4).
#' @param vector Numeric vector or list. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL). If
#' vector is a list, must be named as the names of X param followed by the number of variables to select.
#' @param design Numeric matrix. Matrix of size (number of blocks in X) x (number of blocks in X) with
#' values between 0 and 1. Each value indicates the strength of the relationship to be modeled between
#' two blocks; a value of 0 indicates no relationship, 1 is the maximum value. If NULL, auto-design is computed (default: NULL).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 1).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use (default: NULL).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
#' as n.cut_points^n.blocks models will be computed as minimum (default: 5).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: "AUC").
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
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01).
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @return Instance of class "Coxmos" and model "MB.sPLS-DACOX". The class contains the following
#' elements:
#' \code{X}: List of normalized X data information.
#' \itemize{
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(weightings)}: PLS weights
#'  \item \code{(weightings_norm)}: PLS normalize weights
#'  \item \code{(W.star)}: PLS W* vector
#'  \item \code{(scores)}: PLS scores/variates
#'  \item \code{(E)}: error matrices
#'  \item \code{(x.mean)}: mean values for X matrix
#'  \item \code{(x.sd)}: standard deviation for X matrix
#'  }
#' \code{Y}: List of normalized Y data information.
#' \itemize{
#'  \item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
#'  \item \code{(dr.mean)}: mean values for deviance residuals Y matrix
#'  \item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
#'  \item \code{(data)}: normalized X matrix
#'  \item \code{(y.mean)}: mean values for Y matrix
#'  \item \code{(y.sd)}: standard deviation for Y matrix'
#'  }
#' \code{survival_model}: List of survival model information.
#' \itemize{
#'  \item \code{fit}: coxph object.
#'  \item \code{AIC}: AIC of cox model.
#'  \item \code{BIC}: BIC of cox model.
#'  \item \code{lp}: linear predictors for train data.
#'  \item \code{coef}: Coefficients for cox model.
#'  \item \code{YChapeau}: Y Chapeau residuals.
#'  \item \code{Yresidus}: Y residuals.
#' }
#'
#' \code{mb.model}: List of sPLS models computed for each block.
#'
#' \code{n.comp}: Number of components selected.
#'
#' \code{n.varX}: Number of variables selected for each block.
#'
#' \code{call}: call function
#'
#' \code{X_input}: X input matrix
#'
#' \code{Y_input}: Y input matrix
#'
#' \code{design}: Design used for the strength of the relationship to be modeled between blocks.
#'
#' \code{B.hat}: PLS beta matrix
#'
#' \code{R2}: sPLS acumulate R2
#'
#' \code{alpha}: Significance threshold used.
#'
#' \code{nsv}: Variables removed due to non-significance.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{nz_coeffvar}: Variables removed by coefficient variation near zero.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' \code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.
#'
#' \code{nz_coeffvar}: Variables removed by coefficient variation near zero.
#'
#' \code{time}: time consumed for running the cox analysis.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @references
#' \insertRef{MixOmics}{Coxmos}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("X_multiomic")
#' data("Y_multiomic")
#' X <- X_multiomic
#' X$mirna <- X$mirna[,1:50]
#' X$proteomic <- X$proteomic[,1:50]
#' Y <- Y_multiomic
#' mb.splsdacox(X, Y, n.comp = 2, vector = NULL, x.center = TRUE, x.scale = TRUE)
#' }

mb.splsdacox <- function (X, Y,
                          n.comp = 4, vector = NULL, design = NULL,
                          MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5, EVAL_METHOD = "AUC",
                          x.center = TRUE, x.scale = FALSE,
                          remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                          remove_non_significant = FALSE, alpha = 0.05,
                          MIN_AUC_INCREASE = 0.01, pred.method = "cenROC", max.iter = 200,
                          times = NULL, max_time_points = 15,
                          MIN_EPV = 5, returnData = TRUE, verbose = FALSE){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### Check values classes and ranges
  params_with_limits <- list("alpha" = alpha, "MIN_AUC_INCREASE" = MIN_AUC_INCREASE)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("n.comp" = n.comp, "MIN_NVAR" = MIN_NVAR, "n.cut_points" = n.cut_points,
                         "max_time_points" = max_time_points,
                         "MIN_EPV" = MIN_EPV, "tol" = tol, "max.iter" = max.iter)

  if(!is.null(MAX_NVAR)){
    numeric_params$MAX_NVAR <- MAX_NVAR
  }

  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                         #"y.center" = y.center, "y.scale" = y.scale,
                         "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                         "remove_non_significant" = remove_non_significant, "returnData" = returnData, "verbose" = verbose)
  check_class(logical_params, class = "logical")

  character_params <- list("EVAL_METHOD" = EVAL_METHOD, "pred.method" = pred.method)
  check_class(character_params, class = "character")

  #### Check rownames
  lst_check <- checkXY.rownames.mb(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Check colnames
  X <- checkColnamesIllegalChars.mb(X)

  #### REQUIREMENTS
  checkX.colnames.mb(X)
  checkY.colnames(Y)
  lst_check <- checkXY.mb.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### ZERO VARIANCE - ALWAYS
  lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
                                             remove_near_zero_variance = remove_near_zero_variance,
                                             remove_zero_variance = remove_zero_variance,
                                             toKeep.zv = toKeep.zv,
                                             freqCut = FREQ_CUT)
  X <- lst_dnz$X
  variablesDeleted <- lst_dnz$variablesDeleted

  #### COEF VARIATION
  lst_dnzc <- deleteNearZeroCoefficientOfVariation.mb(X = X, LIMIT = 0.1)
  X <- lst_dnzc$X
  variablesDeleted_cvar <- lst_dnzc$variablesDeleted

  #### SCALING
  lst_scale <- XY.mb.scale(X, Y, x.center, x.scale, y.center, y.scale)
  Xh <- lst_scale$Xh
  Yh <- lst_scale$Yh
  xmeans <- lst_scale$xmeans
  xsds <- lst_scale$xsds
  ymeans <- lst_scale$ymeans
  ysds <- lst_scale$ysds

  X_norm <- Xh

  #### MAX PREDICTORS
  n.comp <- check.mb.maxPredictors(X, Y, MIN_EPV, n.comp, verbose = verbose)
  max_comps <- min(unlist(purrr::map(X, ~ncol(.))))
  n.comp <- min(n.comp, max_comps)

  E <- list()
  R2 <- list()
  SCR <- list()
  SCT <- list()

  XXNA <- purrr::map(Xh, ~is.na(.)) #TRUE is NA
  YNA <- is.na(Y) #TRUE is NA

  #### ### ### ### ### ### ### ### ### ### ### ##
  ### ##         MB.sPLS-DACOX             ### ##
  #### ### ### ### ### ### ### ### ### ### ### ##

  # AUTO DESIGN - https://mixomicsteam.github.io/mixOmics-Vignette/id_06.html#id_06:diablo-design
  if(is.null(design)){
    design <- getDesign.MB(Xh)
  }

  #### ### ### ### ### ### ### ### ### ###
  # DIVIDE Y VENCERAS - BEST VECTOR SIZE #
  #### ### ### ### ### ### ### ### ### ###

  if(is.null(times)){
    times <- getTimesVector(Yh, max_time_points)
  }

  DR_coxph = NULL #not used in plsda

  if(is.null(vector)){
    lst_BV <- getBestVectorMB(Xh = Xh, DR_coxph = DR_coxph, Yh = Yh, n.comp = n.comp, max.iter = max.iter, vector = vector,
                              MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                              EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = FALSE, n_cores = NULL, mode = "splsda", times = times,
                              max_time_points = max_time_points, verbose = verbose)
    keepX <- lst_BV$best.keepX
    plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
  }else{
    if(isa(vector, "list")){
      keepX <- vector
      #if list, but not n.comp length... and just one value in each block
      if(!all(unlist(purrr::map(keepX, ~length(.)==n.comp))) & all(unlist(purrr::map(keepX, ~length(.)==1)))){
        keepX <- purrr::map(keepX, ~rep(., n.comp))
      }else if(!all(unlist(purrr::map(keepX, ~length(.)==1)))){
        #more than one value... just take the first one
        keepX <- purrr::map(keepX, ~rep(.[[1]], n.comp))
      }
    }else{
      #vector is the same length of blocks in X (so each value correspond to each vector)
      if(length(vector)==length(X)){
        keepX <- list()
        for(e in 1:length(vector)){
          keepX[[e]] <- rep(vector[[e]], n.comp)
        }
        names(keepX) <- names(X)
      }else{
        message("Vector does not has the proper structure. Optimizing best n.variables by using your vector as start vector.")
        lst_BV <- getBestVectorMB(Xh = Xh, DR_coxph = DR_coxph, Yh = Yh, n.comp = n.comp, max.iter = max.iter, vector = vector,
                                  MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, cut_points = n.cut_points,
                                  EVAL_METHOD = EVAL_METHOD, EVAL_EVALUATOR = pred.method, PARALLEL = FALSE, n_cores = NULL, mode = "splsda", times = times,
                                  max_time_points = max_time_points, verbose = verbose)
        keepX <- lst_BV$best.keepX
        plotVAR <- plot_VAR_eval(lst_BV, EVAL_METHOD = EVAL_METHOD)
      }
    }
  }

  mb.splsda <- mixOmics::block.splsda(Xh, Yh[,"event"], scale = FALSE, ncomp = n.comp, keepX = keepX, design = design,
                                      max.iter = max.iter, near.zero.var = FALSE, all.outputs = TRUE)

  #PREDICTION
  predplsfit <- predict_mixOmics.mb.pls(mb.spls = mb.splsda, Xh, n.comp)

  for(block in names(predplsfit$predict)){
    # E[[block]] <- list()
    # SCR[[block]] <- list()
    # SCT[[block]] <- list()
    R2[[block]] <- list()
    for(h in 1:n.comp){
      # E[[block]][[h]] <- Yh[,"event",drop=F] - predplsfit$predict[[block]][,h]
      #
      # SCR[[block]][[h]] = sum(apply(E[[block]][[h]],2,function(x) sum(x**2)))
      # SCT[[block]][[h]] = sum(apply(as.matrix(Yh[,"event"]),2,function(x) sum(x**2))) #equivalent sum((Yh[,"event"] - mean(Yh[,"event"]))**2)
      #
      # R2[[block]][[h]] = 1 - (SCR[[block]][[h]]/SCT[[block]][[h]]) #deviance residuals explanation
      R2[[block]][[h]] = mb.splsda$prop_expl_var[[block]][[h]] #deviance residuals explanation
    }
    R2[[block]] = mb.splsda$prop_expl_var[[block]]
  }

  #### ### ### ### ### ### ### ### ### ### ### #
  #                                            #
  #      Computation of the coefficients       #
  #      of the model with kk components       #
  #                                            #
  #### ### ### ### ### ### ### ### ### ### ### #

  #### ### ### ### ### ### ### ### ### ### ### ###
  ### ###         MB:sPLS-DACOX            ### ###
  #### ### ### ### ### ### ### ### ### ### ### ###
  n.comp_used <- ncol(mb.splsda$variates$Y) #can be lesser than expected because we have lesser variables to select because penalization

  n.varX_used <- list()
  for(i in names(Xh)){
    aux <- list()
    for(j in 1:n.comp){
      aux[[j]] <- rownames(mb.splsda$loadings[[i]][which(mb.splsda$loadings[[i]][,j]!=0),j,drop = FALSE])
    }
    names(aux) <- colnames(mb.splsda$loadings[[i]])
    n.varX_used[[i]] <- aux
  }

  data <- as.data.frame(mb.splsda$variates[[1]][,,drop = FALSE])
  for(b in names(Xh)[2:length(Xh)]){
    data <- cbind(data, as.data.frame(mb.splsda$variates[[b]][,,drop = FALSE]))
  }

  update_colnames <- paste0("comp_", 1:ncol(mb.splsda$variates[[1]]))
  colnames(data) <- apply(expand.grid(update_colnames, names(Xh)), 1, paste, collapse="_")
  cox_model <- cox(X = data, Y = Yh,
                   x.center = FALSE, x.scale = FALSE,
                   #y.center = FALSE, y.scale = FALSE,
                   alpha = alpha, remove_non_significant = FALSE, FORCE = TRUE)

  # RETURN a MODEL with ALL significant Variables from complete, deleting one by one
  removed_variables <- NULL
  removed_variables_cor <- NULL

  # REMOVE NA-PVAL VARIABLES
  # p_val could be NA for some variables (if NA change to P-VAL=1)
  # DO IT ALWAYS, we do not want problems in COX models
  if(all(c("time", "event") %in% colnames(data))){
    lst_model <- removeNAorINFcoxmodel(model = cox_model$survival_model$fit, data = data, time.value = NULL, event.value = NULL)
  }else{
    lst_model <- removeNAorINFcoxmodel(model = cox_model$survival_model$fit, data = cbind(data, Yh), time.value = NULL, event.value = NULL)
  }
  cox_model$survival_model$fit <- lst_model$model
  removed_variables_cor <- c(removed_variables_cor, lst_model$removed_variables)

  #RETURN a MODEL with ALL significant Variables from complete, deleting one by one in backward method
  if(remove_non_significant){
    if(all(c("time", "event") %in% colnames(data))){
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$survival_model$fit, alpha = alpha, cox_input = data, time.value = NULL, event.value = NULL)
    }else{
      lst_rnsc <- removeNonSignificativeCox(cox = cox_model$survival_model$fit, alpha = alpha, cox_input = cbind(data, Yh), time.value = NULL, event.value = NULL)
    }

    cox_model$survival_model$fit <- lst_rnsc$cox
    removed_variables <- lst_rnsc$removed_variables
  }

  survival_model <- cox_model$survival_model

  if(isa(survival_model$fit,"coxph")){
    survival_model <- getInfoCoxModel(survival_model$fit)
  }else{
    survival_model <- NULL
  }

  #get W.star
  Tmat <- Pmat <- Cmat <- Wmat <- W.star <- B.hat <- list()
  for(i in 1:length(Xh)){
    #select just features != 0 (selected features)
    names <- purrr::map(1:n.comp_used, ~rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,.,drop = FALSE]!=0)])
    all_names <- unique(unlist(names))

    aux_Pmat = matrix(data = 0, nrow = ncol(Xh[[i]]), ncol = n.comp_used)
    rownames(aux_Pmat) <- colnames(Xh[[i]])
    colnames(aux_Pmat) <- colnames(mb.splsda$loadings[[i]])

    for(c in 1:n.comp_used){
      names <- rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,c,drop = FALSE]!=0)]
      aux <- crossprod(Xh[[i]][,names,drop = FALSE], mb.splsda$variates[[i]][,c])
      aux_Pmat[names,c] = aux
    }

    Pmat[[i]] = aux_Pmat
    Cmat[[i]] = crossprod(Yh[,"event"], mb.splsda$variates[[i]])
    Wmat[[i]] = mb.splsda$loadings[[i]]
    Tmat[[i]] = mb.splsda$variates[[i]]

    colnames(Wmat[[i]]) <- paste0("comp_", 1:ncol(Wmat[[i]]))
    colnames(Pmat[[i]]) <- paste0("comp_", 1:ncol(Pmat[[i]]))
    colnames(Tmat[[i]]) <- paste0("comp_", 1:ncol(Tmat[[i]]))

    # W.star[[i]] <- lapply(1:n.comp, function(x){Wmat[[i]][,1:x,drop = FALSE] %*% solve(t(Pmat[[i]][,1:x,drop = FALSE]) %*% Wmat[[i]][, 1:x,drop = FALSE])})
    # B.hat[[i]] <- lapply(1:n.comp, function(x){W.star[[i]][[x]][,1:x,drop = FALSE] %*% t(Cmat[[i]][,1:x,drop = FALSE])})

    aux_W.star = matrix(data = 0, nrow = ncol(Xh[[i]]), ncol = n.comp_used)
    rownames(aux_W.star) <- colnames(Xh[[i]])
    colnames(aux_W.star) <- colnames(mb.splsda$loadings[[i]])

    for(c in 1:n.comp_used){
      names <- rownames(mb.splsda$loadings[[i]])[which(mb.splsda$loadings[[i]][,c,drop = FALSE]!=0)]

      if(is.null(Pmat[[i]][names,c,drop = FALSE]) | is.null(Wmat[[i]][names,c,drop = FALSE])){
        message(paste0(pkg.env$mb.splsdacox, " model cannot be computed because P or W vectors are NULL. Returning NA."))
        # invisible(gc())
        return(NA)
      }

      #aux <- Wmat[[i]][names,c,drop = FALSE] %*% solve(t(Pmat[[i]][names,c,drop = FALSE]) %*% Wmat[[i]][names,c,drop = FALSE])
      #W.star
      #sometimes solve(t(P) %*% W)
      #system is computationally singular: reciprocal condition number = 6.24697e-18
      # PW <- tryCatch(expr = {solve(t(Pmat[[i]][names,c,drop = FALSE]) %*% Wmat[[i]][names,c,drop = FALSE], tol = tol)},
      #                error = function(e){
      #                  if(verbose){
      #                    message(e$message)
      #                  }
      #                  NA
      #                })
      PW <- tryCatch(expr = {MASS::ginv(t(Pmat[[i]][names,c,drop = FALSE]) %*% Wmat[[i]][names,c,drop = FALSE])},
                     error = function(e){
                       if(verbose){
                         message(e$message)
                       }
                       NA
                     })

      if(all(is.na(PW))){
        message(paste0(pkg.env$mb.splsdacox," model cannot be computed due to ginv(t(P) %*% W). Multicollineality could be present in your data. Returning NA."))
        # invisible(gc())
        return(NA)
      }

      # What happen when you cannot compute W.star but you have P and W?
      aux <- Wmat[[i]][names,c,drop = FALSE] %*% PW
      aux_W.star[names,c] = aux
    }

    W.star[[i]] <- aux_W.star
    B.hat[[i]] <- W.star[[i]] %*% t(Cmat[[i]][,1:n.comp_used,drop = FALSE])

    colnames(W.star[[i]]) <- paste0("comp_", 1:ncol(W.star[[i]]))
  }

  # #get W.star
  # Tmat <- Pmat <- Cmat <- Wmat <- W.star <- B.hat <- list()
  # for(i in 1:length(Xh)){
  #   Pmat[[i]] = crossprod(Xh[[i]], mb.splsda$variates[[i]])
  #   Cmat[[i]] = crossprod(Yh[,"event"], mb.splsda$variates[[i]])
  #   Wmat[[i]] = mb.splsda$loadings[[i]]
  #   Tmat[[i]] = mb.splsda$variates[[i]]
  #
  #   colnames(Wmat[[i]]) <- paste0("comp_", 1:ncol(Wmat[[i]]))
  #   colnames(Pmat[[i]]) <- paste0("comp_", 1:ncol(Pmat[[i]]))
  #   colnames(Tmat[[i]]) <- paste0("comp_", 1:ncol(Tmat[[i]]))
  #
  #   # W.star[[i]] <- lapply(1:n.comp, function(x){Wmat[[i]][,1:x,drop = FALSE] %*% solve(t(Pmat[[i]][,1:x,drop = FALSE]) %*% Wmat[[i]][, 1:x,drop = FALSE])})
  #   # B.hat[[i]] <- lapply(1:n.comp, function(x){W.star[[i]][[x]][,1:x,drop = FALSE] %*% t(Cmat[[i]][,1:x,drop = FALSE])})
  #
  #   W.star[[i]] <- Wmat[[i]][,1:n.comp_used,drop = FALSE] %*% solve(t(Pmat[[i]][,1:n.comp_used,drop = FALSE]) %*% Wmat[[i]][, 1:n.comp_used,drop = FALSE])
  #   B.hat[[i]] <- W.star[[i]] %*% t(Cmat[[i]][,1:n.comp_used,drop = FALSE])
  # }

  names(Pmat) <- names(Xh)
  names(Cmat) <- names(Xh)
  names(Wmat) <- names(Xh)
  names(Tmat) <- names(Xh)
  names(W.star) <- names(Xh)
  names(B.hat) <- names(Xh)

  #MIX Omics, a la hora de generar los nuevos scores para nuevas X (o las mismas de entrenamiento),
  #a parte de realizar la multiplicacion X*W.STAR, realiza luego una normalización de los scores en base a la norma de la propia X usada,
  #de esa manera, en el multiblock de SPLS los resultados no coinciden con los de la funcion predict de MIXOMICS. La siguiente linea es
  #la que se ejecuta una vez realizado el calculo de los nuevos SCORES.

  # head(predplsfit$variates$genes)
  # head(mb.splsda$X$genes %*% W.star[[1]][[n.comp]])
  # head(mb.splsda$X$genes %*% W.star[[1]][[n.comp]])
  # head(mb.splsda$X$genes %*% Wmat[[1]] %*% solve(t(Pmat[[1]]) %*% Wmat[[1]]))
  # #
  # Pmat[[1]] = crossprod(Xh$genes, tt_mbsplsDR[[1]])
  # Wmat[[1]] = mb.splsda$AVE$AVE_inner
  # head(mb.splsda$X$genes %*% Wmat[[1]] %*% solve(t(Pmat[[1]]) %*% Wmat[[1]]))
  #
  # new_t <- mb.splsda$X$genes %*% W.star$genes[[n.comp_used]]
  # new_t2 <- matrix(data = sapply(1:ncol(new_t),
  #                                function(x) {new_t[, x] * apply(mb.splsda$variates$genes, 2,
  #                                                                function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(Xh$genes), ncol = ncol(new_t))
  # head(new_t2)

  # Si lo aplicamos a SPLS normal, también falla el cáclulo de la W*. Puede ser que sea debido a que los cálculos de
  # los loadings de X se estén realizando con la normalización de la C y por tanto la corrección de la norma soluciona el problema.
  # Sin embargo, hubiera sido más sencillo trabajar directamente con una metodología correcta. En mi caso, si utilizo mixomics, debo usar
  # su función siempre para predecir los scores de las nuevas X y NO LO ESTOY HACIENDO!

  #get W.star
  W <- Wmat
  P <- Pmat
  W.star <- W.star
  B.hat <- B.hat
  Ts <- Tmat

  func_call <- match.call()

  if(!returnData){
    survival_model <- removeInfoSurvivalModel(survival_model)
  }

  all_scores <- NULL
  for(b in names(Ts)){
    aux_scores <- Ts[[b]]
    colnames(aux_scores) <- paste0(colnames(aux_scores), "_", b)
    all_scores <- cbind(all_scores, aux_scores)
  }

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  return(mb.splsdacox_class(list(X = list("data" = if(returnData) X_norm else NA,
                                          "loadings" = P,
                                          "weightings" = if(returnData) W else NA,
                                          "W.star" = W.star,
                                          "scores" = Ts,
                                          "scores_all" = all_scores,
                                          "E" = if(returnData) E else NA,
                                          "x.mean" = xmeans, "x.sd" = xsds),
                                 Y = list("data" = Yh,
                                          "y.mean" = ymeans, "y.sd" = ysds),
                                 survival_model = survival_model,
                                 mb.model = mb.splsda,
                                 n.comp = n.comp_used, #number of components
                                 n.varX = n.varX_used,
                                 call = if(returnData) func_call else NA,
                                 X_input = if(returnData) X_original else NA,
                                 Y_input = if(returnData) Y_original else NA,
                                 design = design,
                                 B.hat = B.hat,
                                 R2 = R2,
                                 # SCR = SCR,
                                 # SCT = SCT,
                                 alpha = alpha,
                                 nsv = removed_variables,
                                 nzv = variablesDeleted,
                                 nz_coeffvar = variablesDeleted_cvar,
                                 class = pkg.env$mb.splsdacox,
                                 time = time)))
}

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' MB.sPLS-DACOX Cross-Validation
#' @description The cv.mb.splsdacox function performs cross-validation for the MB.sPLS-DACOX model,
#' a specialized model tailored for survival analysis with high-dimensional data. This function
#' systematically evaluates the performance of the model across different hyperparameters and
#' configurations to determine the optimal settings for the given data.
#'
#' @details The function operates by partitioning the data into multiple subsets (folds) and
#' iteratively holding out one subset for validation while training on the remaining subsets. The
#' cross-validation process is repeated for a specified number of runs, ensuring a robust assessment
#' of the model's performance. The function offers flexibility in terms of the number of PLS components,
#' the range of variables considered, and the evaluation metrics used.
#'
#' The function provides an option to center and scale the explanatory variables, which can be crucial
#' for ensuring consistent performance, especially when the variables are measured on different scales.
#' Additionally, the function incorporates features to handle near-zero and zero variance variables,
#' which can be problematic in high-dimensional datasets.
#'
#' For model evaluation, users can choose between various metrics, including AUC, c-index, and Brier
#' Score. The function also allows for the specification of weights for these metrics, enabling users
#' to prioritize certain metrics over others based on the research context.
#'
#' The function's design also emphasizes computational efficiency. It offers a parallel processing
#' option to expedite the cross-validation process, especially beneficial for large datasets. However,
#' users should be cautious about potential high RAM consumption when using this option.
#'
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation
#' (default: 8).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL). If
#' vector is a list, must be named as the names of X param followed by the number of variables to select.
#' @param design Numeric matrix. Matrix of size (number of blocks in X) x (number of blocks in X) with
#' values between 0 and 1. Each value indicates the strength of the relationship to be modeled between
#' two blocks; a value of 0 indicates no relationship, 1 is the maximum value. If NULL, auto-design is computed (default: NULL).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 1).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use (default: NULL).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used (default: 5).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: "AUC").
#' @param n_run Numeric. Number of runs for cross validation (default: 3).
#' @param k_folds Numeric. Number of folds for cross validation (default: 10).
#' @param x.center Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If remove_near_zero_variance = TRUE, near zero variance
#' variables will be removed (default: TRUE).
#' @param remove_zero_variance Logical. If remove_zero_variance = TRUE, zero variance variables will
#' be removed (default: TRUE).
#' @param toKeep.zv Character vector. Name of variables in X to not be deleted by (near) zero variance
#' filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If remove_variance_at_fold_level = TRUE, (near) zero
#' variance will be removed at fold level. Not recommended. (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE,
#' non-significant models are removed before computing the evaluation. A non-significant model is a
#' model with at least one component/variable with a P-Value higher than the alpha cutoff.
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant
#' variables/components in final cox model will be removed until all variables are significant by
#' forward selection (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_C.Index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_I.BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is
#' reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the
#' 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC
#' improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the
#' evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following:
#' "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200).
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated
#' simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test
#' observations. Once all have their linear predictors, the evaluation is perform across all the
#' observations together (default: FALSE).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param PARALLEL Logical. Run the cross validation with multicore option. As many cores as your
#' total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).
#' @param n_cores Numeric. Number of cores to use for parallel processing. This parameter is only
#' used if `PARALLEL` is `TRUE`. If `NULL`, it will use all available cores minus one. Otherwise,
#' it will use the minimum between the value specified and the total number of cores - 1. The fewer
#' cores used, the less RAM memory will be used.(default: NULL).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return Instance of class "Coxmos" and model "cv.MB.sPLS-DACOX".
#' \code{best_model_info}: A data.frame with the information for the best model.
#' \code{df_results_folds}: A data.frame with fold-level information.
#' \code{df_results_runs}: A data.frame with run-level information.
#' \code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha
#' information).
#'
#' \code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
#' \code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.
#'
#' \code{opt.comp}: Optimal component selected by the best_model.
#' \code{opt.nvar}: Optimal number of variables selected by the best_model.
#' \code{design}: Design used for the strength of the relationship to be modeled between blocks.
#'
#' \code{plot_AIC}: AIC plot by each hyper-parameter.
#' \code{plot_C.Index}: C-Index plot by each hyper-parameter.
#' \code{plot_I.BRIER}: Integrative Brier Score plot by each hyper-parameter.
#' \code{plot_AUC}: AUC plot by each hyper-parameter.
#'
#' \code{class}: Cross-Validated model class.
#'
#' \code{lst_train_indexes}: List (of lists) of indexes for the observations used in each run/fold
#' for train the models.
#' \code{lst_test_indexes}: List (of lists) of indexes for the observations used in each run/fold
#' for test the models.
#'
#' \code{time}: time consumed for running the cross-validated function.
#'
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("X_multiomic")
#' data("Y_multiomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_multiomic$event, p = .5, list = FALSE, times = 1)
#' X_train <- X_multiomic
#' X_train$mirna <- X_train$mirna[index_train,1:50]
#' X_train$proteomic <- X_train$proteomic[index_train,1:50]
#' Y_train <- Y_multiomic[index_train,]
#' cv.mb.splsdacox_model <- cv.mb.splsdacox(X_train, Y_train, max.ncomp = 2, vector = NULL,
#' n_run = 1, k_folds = 2, x.center = TRUE, x.scale = TRUE)
#' }

cv.mb.splsdacox <- function(X, Y,
                            max.ncomp = 8, vector = NULL, design = NULL,
                            MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5, EVAL_METHOD = "AUC",
                            n_run = 3, k_folds = 10,
                            x.center = TRUE, x.scale = FALSE,
                            remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                            remove_variance_at_fold_level = FALSE,
                            remove_non_significant_models = FALSE, remove_non_significant = FALSE, alpha = 0.05,
                            w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL,
                            max_time_points = 15,
                            MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                            pred.attr = "mean", pred.method = "cenROC", max.iter= 200, fast_mode = FALSE,
                            MIN_EPV = 5, return_models = FALSE, returnData = FALSE,
                            PARALLEL = FALSE, n_cores = NULL, verbose = FALSE, seed = 123){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

  t1 <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

  #### ### ###
  # WARNINGS #
  #### ### ###

  #### Check evaluator installed:
  checkLibraryEvaluator(pred.method)

  #### Check values classes and ranges
  params_with_limits <- list("MIN_AUC_INCREASE" = MIN_AUC_INCREASE, "MIN_AUC" = MIN_AUC, "alpha" = alpha,
                             "w_AIC" = w_AIC, "w_C.Index" = w_C.Index, "w_AUC" = w_AUC, "w_I.BRIER" = w_I.BRIER)
  check_min0_max1_variables(params_with_limits)

  numeric_params <- list("max.ncomp" = max.ncomp, "MIN_NVAR" = MIN_NVAR, "n.cut_points" = n.cut_points,
                         "n_run" = n_run, "k_folds" = k_folds, "max_time_points" = max_time_points,
                         "MIN_COMP_TO_CHECK" = MIN_COMP_TO_CHECK, "MIN_EPV" = MIN_EPV, "seed" = seed, "tol" = tol)

  if(!is.null(MAX_NVAR)){
    numeric_params$MAX_NVAR <- MAX_NVAR
  }

  check_class(numeric_params, class = "numeric")

  logical_params <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                         #"y.center" = y.center, "y.scale" = y.scale,
                         "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                         "remove_variance_at_fold_level" = remove_variance_at_fold_level,
                         "remove_non_significant_models" = remove_non_significant_models,
                         "remove_non_significant" = remove_non_significant,
                         "return_models" = return_models,"returnData" = returnData, "verbose" = verbose, "PARALLEL" = PARALLEL)
  check_class(logical_params, class = "logical")

  character_params <- list("EVAL_METHOD" = EVAL_METHOD, "pred.attr" = pred.attr, "pred.method" = pred.method)
  check_class(character_params, class = "character")

  #### Check cv-folds
  lst_checkFR <- checkFoldRuns(Y, n_run, k_folds, fast_mode)
  n_run <- lst_checkFR$n_run
  fast_mode <- lst_checkFR$fast_mode

  #### Check rownames
  lst_check <- checkXY.rownames.mb(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  #### Illegal chars in colnames
  X <- checkColnamesIllegalChars.mb(X)

  #### REQUIREMENTS
  checkX.colnames.mb(X)
  checkY.colnames(Y)
  lst_check <- checkXY.mb.class(X, Y, verbose = verbose)

  check.cv.weights(c(w_AIC, w_C.Index, w_I.BRIER, w_AUC))
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### ZERO VARIANCE - ALWAYS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnz <- deleteZeroOrNearZeroVariance.mb(X = X,
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
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnzc <- deleteNearZeroCoefficientOfVariation.mb(X = X)
    X <- lst_dnzc$X
    variablesDeleted_cvar <- lst_dnzc$variablesDeleted
  }else{
    variablesDeleted_cvar <- NULL
  }

  #### MAX PREDICTORS
  max.ncomp <- check.mb.ncomp(X, max.ncomp)
  max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)
  if(MIN_COMP_TO_CHECK >= max.ncomp){
    MIN_COMP_TO_CHECK = max(max.ncomp-1, 1)
  }

  # AUTO DESIGN - https://mixomicsteam.github.io/mixOmics-Vignette/id_06.html#id_06:diablo-design
  if(is.null(design)){
    #### SCALING
    lst_scale <- XY.mb.scale(X, Y, x.center, x.scale, y.center, y.scale)
    Xh <- lst_scale$Xh
    design <- getDesign.MB(Xh)
  }

  #### #
  # CV #
  #### #
  # lst_data <- splitData_Iterations_Folds.mb(X, Y, n_run = n_run, k_folds = k_folds, seed = seed) #FOR TEST
  # lst_X_train <- lst_data$lst_X_train
  # lst_Y_train <- lst_data$lst_Y_train
  # lst_X_test <- lst_data$lst_X_test
  # lst_Y_test <- lst_data$lst_Y_test
  #
  # lst_train_indexes <- lst_data$lst_train_index
  # lst_test_indexes <- lst_data$lst_test_index

  lst_data <- splitData_Iterations_Folds_indexes(Y, n_run = n_run, k_folds = k_folds, seed = seed) #FOR TEST

  lst_train_indexes <- lst_data$lst_train_index
  lst_test_indexes <- lst_data$lst_test_index

  #### ### ### ###
  # TRAIN MODELS #
  #### ### ### ###
  total_models <- 1 * k_folds * n_run
  #total_models <- max.ncomp * k_folds * n_run

  comp_model_lst <- get_Coxmos_models2.0(method = pkg.env$mb.splsdacox,
                                         X_train = X, Y_train = Y,
                                         lst_X_train = lst_train_indexes, lst_Y_train = lst_train_indexes,
                                         max.ncomp = max.ncomp, penalty.list = NULL, EN.alpha.list = NULL, max.variables = NULL, vector = vector, design = design,
                                         n_run = n_run, k_folds = k_folds,
                                         MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, MIN_AUC_INCREASE = MIN_AUC_INCREASE, EVAL_METHOD = EVAL_METHOD,
                                         n.cut_points = n.cut_points,
                                         x.center = x.center, x.scale = x.scale,
                                         y.center = y.center, y.scale = y.scale,
                                         remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = FALSE, toKeep.zv = NULL,
                                         alpha = alpha, MIN_EPV = MIN_EPV,
                                         remove_non_significant = remove_non_significant, tol = tol,
                                         max.iter = max.iter, times = times, pred.method = pred.method, max_time_points = max_time_points,
                                         returnData = returnData, total_models = total_models,
                                         PARALLEL = PARALLEL, n_cores = n_cores, verbose = verbose)

  # already check in Coxmos_models
  # if(all(is.na(unlist(lst_model)))){
  #   message(paste0("Best model could NOT be obtained. All models computed present problems. Try to remove variance at fold level. If problem persists, try to delete manually some problematic variables."))
  #
  #   t2 <- Sys.time()
  #   time <- difftime(t2,t1,units = "mins")
  #   if(return_models){
  #     return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = lst_model, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_C.Index = NULL, plot_I.BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  #   }else{
  #     return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_C.Index = NULL, plot_I.BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  #   }
  # }

  #### ### ### ### ### ### #
  # BEST MODEL FOR CV DATA #
  #### ### ### ### ### ### #
  total_models <- max.ncomp * k_folds * n_run
  df_results_evals <- get_COX_evaluation_AIC_CINDEX(comp_model_lst = comp_model_lst, alpha = alpha,
                                                    max.ncomp = max.ncomp, penalty.list = NULL, n_run = n_run, k_folds = k_folds,
                                                    total_models = total_models, remove_non_significant_models = remove_non_significant_models, verbose = verbose)

  if(all(is.null(df_results_evals))){
    message(paste0("Best model could NOT be obtained. All models computed present problems."))

    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(return_models){
      return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_C.Index = NULL, plot_I.BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }else{
      return(cv.mb.splsdacox_class(list(best_model_info = NULL, df_results_folds = NULL, df_results_runs = NULL, df_results_comps = NULL, lst_models = NULL, pred.method = pred.method, opt.comp = NULL, opt.nvar = NULL, plot_AIC = NULL, plot_C.Index = NULL, plot_I.BRIER = NULL, plot_AUC = NULL, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
    }
  }

  #### ### ### ### ### ### #
  # EVALUATING BRIER SCORE #
  #### ### ### ### ### ### #
  df_results_evals_comp <- NULL
  df_results_evals_run <- NULL
  df_results_evals_fold <- NULL
  optimal_comp_index <- NULL
  optimal_comp_flag <- FALSE
  optimal_eta_index <- NULL
  optimal_eta <- NULL

  if(TRUE){ #compute always BRIER SCORE
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    #As we are measuring just one evaluator and one method - PARALLEL = FALSE, n_cores = NULL
    lst_df <- get_COX_evaluation_BRIER(comp_model_lst = comp_model_lst,
                                       fast_mode = fast_mode,
                                       X_test = X, Y_test = Y,
                                       lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                       df_results_evals = df_results_evals, times = times,
                                       pred.method = pred.method, pred.attr = pred.attr,
                                       max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                       MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                       w_I.BRIER = w_I.BRIER, method.train = pkg.env$mb.splsdacox, PARALLEL = FALSE, n_cores = NULL, verbose = verbose)

    df_results_evals_comp <- lst_df$df_results_evals_comp
    df_results_evals_run <- lst_df$df_results_evals_run
    df_results_evals_fold <- lst_df$df_results_evals_fold
  }

  #### ### ### ### #
  # EVALUATING AUC #
  #### ### ### ### #

  if(w_AUC!=0){
    #total_models <- ifelse(!fast_mode, n_run * max.ncomp, k_folds * n_run * max.ncomp)#inside get_COX_evaluation_AUC

    #times should be the same for all folds
    #calculate time vector if still NULL
    if(is.null(times)){
      times <- getTimesVector(Y, max_time_points = max_time_points)
    }

    lst_df <- get_COX_evaluation_AUC(comp_model_lst = comp_model_lst,
                                     X_test = X, Y_test = Y,
                                     lst_X_test = lst_test_indexes, lst_Y_test = lst_test_indexes,
                                     df_results_evals = df_results_evals, times = times,
                                     fast_mode = fast_mode, pred.method = pred.method, pred.attr = pred.attr,
                                     max.ncomp = max.ncomp, n_run = n_run, k_folds = k_folds,
                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE, MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                     w_AUC = w_AUC, method.train = pkg.env$mb.splsdacox, PARALLEL = FALSE, n_cores = NULL, verbose = verbose)

    if(is.null(df_results_evals_comp)){
      df_results_evals_comp <- lst_df$df_results_evals_comp
    }else{
      df_results_evals_comp$AUC <- lst_df$df_results_evals_comp$AUC
    }

    if(is.null(df_results_evals_run)){
      df_results_evals_run <- lst_df$df_results_evals_run
    }else{
      df_results_evals_run$AUC <- lst_df$df_results_evals_run$AUC
    }

    if(is.null(df_results_evals_fold)){
      df_results_evals_fold <- lst_df$df_results_evals_fold
    }else{
      df_results_evals_fold$AUC <- lst_df$df_results_evals_fold$AUC
    }

    optimal_comp_index <- lst_df$optimal_comp_index
    optimal_comp_flag <- lst_df$optimal_comp_flag
    optimal_eta <- lst_df$optimal_eta
    optimal_eta_index <- lst_df$optimal_eta_index
  }

  #### ### ### #
  # BEST MODEL #
  #### ### ### #

  df_results_evals_comp <- cv.getScoreFromWeight(df_results_evals_comp, w_AIC, w_C.Index, w_I.BRIER, w_AUC,
                                                 colname_AIC = "AIC", colname_c_index = "C.Index", colname_AUC = "AUC", colname_BRIER = "IBS")

  if(optimal_comp_flag){
    best_model_info <- df_results_evals_comp[df_results_evals_comp[,"n.comps"]==optimal_comp_index,, drop = FALSE][1,]
    best_model_info <- as.data.frame(best_model_info)
  }else{
    best_model_info <- df_results_evals_comp[which(df_results_evals_comp[,"score"] == max(df_results_evals_comp[,"score"], na.rm = TRUE)),, drop = FALSE][1,]
    best_model_info <- as.data.frame(best_model_info)
  }

  best_n_var <- list()
  aux_n_var <- as.numeric(strsplit(as.character(best_model_info$n.var), "_")[[1]])
  for(e in 1:length(aux_n_var)){
    best_n_var[[e]] <- aux_n_var[[e]]
  }
  names(best_n_var) <- names(X)

  #### ###
  # PLOT #
  #### ###
  class = pkg.env$mb.splsdacox
  lst_EVAL_PLOTS <- get_EVAL_PLOTS(fast_mode = fast_mode, best_model_info = best_model_info, w_AUC = w_AUC, w_I.BRIER = w_I.BRIER, max.ncomp = max.ncomp, penalty.list = NULL,
                                   df_results_evals_fold = df_results_evals_fold, df_results_evals_run = df_results_evals_run, df_results_evals_comp = df_results_evals_comp,
                                   colname_AIC = "AIC", colname_c_index = "C.Index", colname_AUC = "AUC", colname_BRIER = "IBS", x.text = "Component",
                                   class = class)

  ggp_AUC <- lst_EVAL_PLOTS$ggp_AUC
  ggp_IBS <- lst_EVAL_PLOTS$ggp_IBS
  ggp_C.Index <- lst_EVAL_PLOTS$ggp_C.Index
  ggp_AIC <- lst_EVAL_PLOTS$ggp_AIC

  df_results_evals_comp <- lst_EVAL_PLOTS$df_results_evals_comp

  #### ### #
  # RETURN #
  #### ### #

  message(paste0("Best model obtained."))

  t2 <- Sys.time()
  time <- difftime(t2,t1,units = "mins")

  # invisible(gc())
  if(return_models){
    return(cv.mb.splsdacox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = comp_model_lst, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_n_var, design = design, plot_AIC = ggp_AIC, plot_C.Index = ggp_C.Index, plot_I.BRIER = ggp_IBS, plot_AUC = ggp_AUC, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }else{
    return(cv.mb.splsdacox_class(list(best_model_info = best_model_info, df_results_folds = df_results_evals_fold, df_results_runs = df_results_evals_run, df_results_comps = df_results_evals_comp, lst_models = NULL, pred.method = pred.method, opt.comp = best_model_info$n.comps, opt.nvar = best_n_var, design = design, plot_AIC = ggp_AIC, plot_C.Index = ggp_C.Index, plot_I.BRIER = ggp_IBS, plot_AUC = ggp_AUC, class = pkg.env$cv.mb.splsdacox, lst_train_indexes = lst_train_indexes, lst_test_indexes = lst_test_indexes, time = time)))
  }
}

### ## ##
# CLASS #
### ## ##

mb.splsdacox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$mb.splsdacox)
  return(model)
}

cv.mb.splsdacox_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.mb.splsdacox)
  return(model)
}
