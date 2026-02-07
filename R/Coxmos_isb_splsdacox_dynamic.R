#### ### ##
# METHODS #
#### ### ##

#' Iterative single-block sPLS-DACOX Dynamic
#' @description This function performs a single-block sparse partial least squares deviance residual
#' Cox analysis (sPLS-DACOX) using the optimal components and variables identified in a previous
#' cross-validation process. It builds the final model based on the selected hyperparameters.
#'
#' @details
#' The `isb.splsdacox_dynamic` function fits a single-block sPLS-DACOX model using the input data
#' and the optimal components and variables determined from cross-validation. The function allows
#' for centering and scaling of the data, and it offers the option to remove variables with near-zero
#' variance, zero variance, or those that are non-significant based on a specified alpha level.
#'
#' This method is particularly suited for high-dimensional data where there are many more variables
#' than observations. The function can handle multiple blocks of data, and integrates them into a
#' single model for Cox proportional hazards analysis.
#'
#' @param X List of numeric matrices or data.frames. Explanatory variables. If the variables are qualitative,
#' they must be transformed into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables with two columns: "time" and "event".
#' Accepted values for the event column are 0/1 or FALSE/TRUE for censored and event observations, respectively.
#' @param cv.isb Instance of class "Coxmos" and model "cv.iSB.sPLS-DACOX-Dynamic". Used to retrieve the
#' optimal components and variables for the sPLS Cox model.
#' @param x.center Logical. If TRUE, the X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If TRUE, the X matrix is scaled to unit variance (default: FALSE).
#' @param remove_near_zero_variance Logical. If TRUE, near-zero variability variables are removed (default: TRUE).
#' @param remove_zero_variance Logical. If TRUE, zero-variance variables are removed (default: TRUE).
#' @param toKeep.zv Character vector. Names of variables in X to retain despite near-zero variability filtering (default: NULL).
#' @param remove_non_significant Logical. If TRUE, non-significant variables/components in the final Cox model
#' are removed through forward selection (default: FALSE).
#' @param alpha Numeric. Significance threshold (default: 0.05).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) for the final Cox model. Limits the
#' number of variables/components allowed (default: 5).
#' @param returnData Logical. If TRUE, returns the original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If TRUE, extra messages will be displayed (default: FALSE).
#'
#' @return An object of class "Coxmos" and model "isb.splsdacox_dynamic", containing:
#' \itemize{
#'   \item \code{X}: List with normalized X data:
#'     \itemize{
#'       \item \code{data}: Normalized X matrix (or NA if not returned).
#'       \item \code{x.mean}: Mean values of the X matrix.
#'       \item \code{x.sd}: Standard deviations of the X matrix.
#'     }
#'   \item \code{Y}: List with normalized Y data:
#'     \itemize{
#'       \item \code{data}: Normalized Y matrix.
#'       \item \code{y.mean}: Mean values of the Y matrix.
#'       \item \code{y.sd}: Standard deviations of the Y matrix.
#'     }
#'   \item \code{survival_model}: Fitted survival model (Cox proportional hazards model).
#'   \item \code{list_spls_models}: List of sPLS models computed for each block.
#'   \item \code{n.comp}: Number of components selected.
#'   \item \code{n.varX}: Number of variables selected per block.
#'   \item \code{call}: Function call.
#'   \item \code{X_input}: Original X matrix (or NA if not returned).
#'   \item \code{Y_input}: Original Y matrix (or NA if not returned).
#'   \item \code{alpha}: Significance threshold used.
#'   \item \code{nsv}: Variables removed due to non-significance.
#'   \item \code{nzv}: Variables removed due to near-zero variability.
#'   \item \code{nz_coeffvar}: Variables removed due to near-zero coefficient of variation.
#'   \item \code{class}: Model class.
#'   \item \code{time}: Time taken to run the analysis.
#' }
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' data("Y_multiomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_multiomic$event, p = .25, list = FALSE, times = 1)
#' X_train <- X_multiomic
#' X_train$mirna <- X_train$mirna[index_train,1:15]
#' X_train$proteomic <- X_train$proteomic[index_train,1:15]
#' Y_train <- Y_multiomic[index_train,]
#' vector <- list()
#' vector$mirna <- c(10, 15)
#' vector$proteomic <- c(10, 15)
#' cv <- cv.isb.splsdacox(X_train, Y_train, max.ncomp = 1, vector = vector, n_run = 1, k_folds = 3)
#' model <- isb.splsdacox(X_train, Y_train, cv)

isb.splsdacox <- function(X, Y,
                                  cv.isb,
                                  x.center = TRUE, x.scale = FALSE,
                                  remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                                  remove_non_significant = FALSE, alpha = 0.05,
                                  MIN_EPV = 5, returnData = TRUE, verbose = FALSE){

        # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
        tol = 1e-10

        t1 <- Sys.time()
        y.center = y.scale = FALSE
        FREQ_CUT <- 95/5

        # Check class of CV
        if(!class(cv.isb) %in% pkg.env$model_class || !attr(cv.isb, "model") %in% pkg.env$cv.isb.splsdacox_dynamic){
          stop(paste0("Cross validation object 'cv.isb' must belong to '", pkg.env$cv.isb.splsdacox_dynamic, "' class."))
        }

        # Compute max.ncomp from cv.isb
        max.ncomp <- max(unlist(purrr::map(names(cv.isb$list_cv_spls_models), ~cv.isb$list_cv_spls_models[[.]]$opt.comp)))
        max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)

        #### Check values classes and ranges
        params_with_limits <- list("alpha" = alpha)
        check_min0_max1_variables(params_with_limits)

        numeric_params <- list("tol" = tol)
        check_class(numeric_params, class = "numeric")

        logical_params <- list("x.center" = unlist(x.center), "x.scale" = unlist(x.scale),
                               #"y.center" = y.center, "y.scale" = y.scale,
                               "remove_near_zero_variance" = remove_near_zero_variance, "remove_zero_variance" = remove_zero_variance,
                               "remove_non_significant" = remove_non_significant,
                               "returnData" = returnData, "verbose" = verbose)
        check_class(logical_params, class = "logical")

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

        max.ncomp <- check.mb.ncomp(X, max.ncomp)

        #### Original data
        X_original <- X
        Y_original <- Y

        time <- Y[,"time"]
        event <- Y[,"event"]

        #### SCALE
        if(length(x.center)==1){
          x.center <- rep(x.center, length(names(X)))
          names(x.center) <- names(X)
        }
        if(length(x.scale)==1){
          x.scale <- rep(x.scale, length(names(X)))
          names(x.scale) <- names(X)
        }

        #### ZERO VARIANCE - ALWAYS
        if((remove_near_zero_variance | remove_zero_variance)){
          lst_dnz <- deleteZeroOrNearZeroVariability.mb(X = X,
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
        if((remove_near_zero_variance | remove_zero_variance)){
          lst_dnzc <- deleteNearZeroCoefficientOfVariation.mb(X = X)
          X <- lst_dnzc$X
          variablesDeleted_cvar <- lst_dnzc$variablesDeleted
        }else{
          variablesDeleted_cvar <- NULL
        }

        #### SCALING
        lst_scale <- XY.mb.scale(X, Y, x.center, x.scale, y.center, y.scale)
        Xh <- lst_scale$Xh
        Yh <- lst_scale$Yh
        xmeans <- lst_scale$xmeans
        xsds <- lst_scale$xsds
        ymeans <- lst_scale$ymeans
        ysds <- lst_scale$ysds

        X_norm <- Xh

        # CREATE INDIVIDUAL MODELS
        lst_sb.spls <- list()

        message(paste0("\n"))

        for(b in names(Xh)){

          message(paste0("Creating ", pkg.env$sb.splsdacox_dynamic ," model for block: ", b, "..."))
          t1 <- Sys.time()

          cv.splsdacox_dynamic_res <- cv.isb$list_cv_spls_models[[b]]

          lst_sb.spls[[b]] <- splsdacox(X = Xh[[b]],
                                                Y = Yh,
                                                n.comp = cv.splsdacox_dynamic_res$opt.comp,
                                                vector = cv.splsdacox_dynamic_res$opt.nvar,
                                                remove_near_zero_variance = FALSE, remove_zero_variance = FALSE, toKeep.zv = NULL,
                                                remove_non_significant = remove_non_significant, alpha = alpha,
                                                returnData = FALSE,
                                                x.center = x.center[[b]], x.scale = x.scale[[b]],
                                                #y.scale = y.scale, y.center = y.center,
                                                verbose = verbose)

          t2 <- Sys.time()
          time <- difftime(t2,t1,units = "mins")
          if(verbose){
            message(paste0("\nTime for block ", b , ": ", as.character(round(time, 2)), "\n"))
          }
        }

        # CHECK ALL MODELS SAME COMPONENTS
        aux_ncomp <- purrr::map(lst_sb.spls, ~.$n.comp)
        aux_nvar <- purrr::map(lst_sb.spls, ~.$opt.nvar)

        # CREATE COMBINE MODEL
        data <- NULL
        cn.merge <- NULL
        for(b in names(Xh)){
          if(!is.null(lst_sb.spls[[b]]$survival_model)){
            data <- cbind(data, lst_sb.spls[[b]]$X$scores)
            cn.merge <- c(cn.merge, paste0(colnames(lst_sb.spls[[b]]$X$scores), "_", b))
          }else{
            next
          }
        }

        #colnames(data) <- apply(expand.grid(colnames(lst_sb.spls[[1]]$X$scores), names(Xh)), 1, paste, collapse="_")
        colnames(data) <- cn.merge
        cox_model <- cox(X = data, Y = Yh,
                         x.center = FALSE, x.scale = FALSE,
                         #y.center = FALSE, y.scale = FALSE,
                         remove_non_significant = remove_non_significant, FORCE = TRUE)

        if(remove_non_significant){
          removed_variables <- cox_model$nsv
        }else{
          removed_variables <- NULL
        }

        #### ### #
        # RETURN #
        #### ### #
        func_call <- match.call()

        if(!returnData){
          survival_model <- removeInfoSurvivalModel(cox_model$survival_model)
        }else{
          survival_model <- cox_model$survival_model
        }

        all_scores <- NULL
        for(b in names(lst_sb.spls)){
          aux_scores <- lst_sb.spls[[b]]$X$scores
          colnames(aux_scores) <- paste0(colnames(aux_scores), "_", b)
          all_scores <- cbind(all_scores, aux_scores)
        }

        t2 <- Sys.time()
        time <- difftime(t2,t1,units = "mins")

        # invisible(gc())
        return(isb.splsdacox_dynamic_class(list(X = list("data" = if(returnData) X_norm else NA,
                                                         "scores" = all_scores,
                                                         "x.mean" = xmeans, "x.sd" = xsds),
                                      Y = list("data" = Yh,
                                               "y.mean" = ymeans, "y.sd" = ysds),
                                      survival_model = survival_model,
                                      list_spls_models = lst_sb.spls,
                                      n.comp = aux_ncomp, #number of components used, but could be lesser than expected because not computed models
                                      n.varX = aux_nvar,
                                      call = func_call,
                                      X_input = if(returnData) X_original else NA,
                                      Y_input = if(returnData) Y_original else NA,
                                      alpha = alpha,
                                      nsv = removed_variables,
                                      nzv = variablesDeleted,
                                      nz_coeffvar = variablesDeleted_cvar,
                                      class = pkg.env$isb.splsdacox_dynamic,
                                      time = time)))
      }

#### ### ### ### ###
# CROSS-EVALUATION #
#### ### ### ### ###

#' Iterative SB.sPLS-DACOX-Dynamic Cross-Validation
#' @description This function performs cross-validated sparse partial least squares single-block for
#' sPLS-DACOX-Dynamic. It returns the optimal number of components and the optimal sparsity penalty value based
#' on cross-validation. Performance can be evaluated using multiple metrics, such as Area Under the Curve
#' (AUC), I. Brier Score, or C-Index. Users can also specify more than one metric simultaneously.
#'
#' @details
#' The `cv.isb.splsdacox_dynamic` function performs cross-validation for the single-block sparse partial least
#' squares deviance residual Cox analysis (sPLS-DACOX). Cross-validation evaluates different hyperparameter
#' combinations, including the number of components (`max.ncomp`) and the number of variables selected (`vector`).
#' The function systematically evaluates models across multiple runs and folds to determine the best configuration.
#' It allows flexibility in metrics, preprocessing steps (centering, scaling, variance filtering), and stopping criteria.
#'
#' For each run, the dataset is divided into training and test sets for the specified number of folds (`k_folds`).
#' Various metrics, such as AIC, C-Index, I. Brier Score, and AUC, are computed to assess model performance. The
#' function identifies the optimal hyperparameters that yield the best performance based on the selected evaluation metrics.
#'
#' Additionally, it offers options to control the evaluation algorithm method (`pred.method`), whether to return
#' all models, and parallel processing (`PARALLEL`). The function also allows the user to control the verbosity of
#' output messages and set the minimum threshold for Events Per Variable (`MIN_EPV`).
#'
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be transformed into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Must contain two columns: "time" and "event".
#' For the event column, accepted values are 0/1 or FALSE/TRUE for censored and event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute during cross-validation (default: 8).
#' @param vector Numeric vector. A vector indicating the number of variables to select for each block and component (default: NULL).
#' @param MIN_NVAR Numeric. Minimum number of variables to select in the model (default: 10).
#' @param MAX_NVAR Numeric. Maximum number of variables to select in the model (default: NULL).
#' @param n.cut_points Numeric. Number of cut points to evaluate the number of variables (default: 5).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement in AUC required between models to continue evaluation (default: 0.01).
#' @param EVAL_METHOD Character. Method for evaluating performance. Must be one of "AUC", "C-Index", etc. (default: "AUC").
#' @param n_run Numeric. Number of runs for cross-validation (default: 3).
#' @param k_folds Numeric. Number of folds for cross-validation (default: 10).
#' @param x.center Logical. If TRUE, the X matrix is centered to zero means (default: TRUE).
#' @param x.scale Logical. If TRUE, the X matrix is scaled to unit variances (default: FALSE).
#' @param remove_near_zero_variance Logical. If TRUE, near-zero variability variables are removed (default: TRUE).
#' @param remove_zero_variance Logical. If TRUE, zero-variance variables are removed (default: TRUE).
#' @param toKeep.zv Character vector. Names of variables in X to retain despite variance filtering (default: NULL).
#' @param remove_variance_at_fold_level Logical. If TRUE, variance filtering is applied at the fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If TRUE, models with non-significant components are removed before evaluation (default: FALSE).
#' @param remove_non_significant Logical. If TRUE, non-significant components in the final Cox model are removed (default: FALSE).
#' @param alpha Numeric. Significance threshold for selecting variables/components (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC in the evaluation. All weights must sum to 1 (default: 0).
#' @param w_C.Index Numeric. Weight for C-Index in the evaluation. All weights must sum to 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC in the evaluation. All weights must sum to 1 (default: 1).
#' @param w_I.BRIER Numeric. Weight for Integrative Brier Score in the evaluation. All weights must sum to 1 (default: 0).
#' @param times Numeric vector. Time points for AUC evaluation (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points for AUC evaluation (default: 15).
#' @param MIN_AUC Numeric. Minimum AUC to achieve during cross-validation (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of components to evaluate before stopping if no improvement is observed (default: 3).
#' @param pred.attr Character. Method for evaluating performance. Must be one of "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation method. Must be one of: "risksetROC", "survivalROC", "cenROC", etc. (default: "cenROC").
#' @param fast_mode Logical. If TRUE, only one fold is evaluated per run; otherwise, all folds are evaluated simultaneously (default: FALSE).
#' @param max.iter Numeric. Maximum number of iterations for convergence (default: 200).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable for the final Cox model (default: 5).
#' @param return_models Logical. If TRUE, returns all models computed during cross-validation (default: FALSE).
#' @param returnData Logical. If TRUE, returns original and normalized X and Y matrices (default: FALSE).
#' @param PARALLEL Logical. If TRUE, runs cross-validation in parallel using multiple cores (default: FALSE).
#' @param n_cores Numeric. Number of cores to use for parallel processing. This parameter is only
#' used if `PARALLEL` is `TRUE`. If `NULL`, it will use all available cores minus one. Otherwise,
#' it will use the minimum between the value specified and the total number of cores - 1. The fewer
#' cores used, the less RAM memory will be used.(default: NULL).
#' @param verbose Logical. If TRUE, extra messages are displayed during execution (default: FALSE).
#' @param seed Numeric. Seed for reproducibility (default: 123).
#'
#' @return An instance of class "Coxmos" and model "cv.SB.sPLS-DACOX-Dynamic", containing:
#' \itemize{
#'   \item \code{best_model_info}: Data frame with the best model's information.
#'   \item \code{df_results_folds}: Data frame with fold-level results.
#'   \item \code{df_results_runs}: Data frame with run-level results.
#'   \item \code{df_results_comps}: Data frame with component-level results.
#'   \item \code{list_cv_spls_models}: List of cross-validated models for each block.
#'   \item \code{opt.comp}: Optimal number of components.
#'   \item \code{opt.nvar}: Optimal number of variables selected.
#'   \item \code{class}: Model class.
#'   \item \code{time}: Time taken to run the cross-validation.
#' }
#' @author Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' data("Y_multiomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_multiomic$event, p = .25, list = FALSE, times = 1)
#' X_train <- X_multiomic
#' X_train$mirna <- X_train$mirna[index_train,1:20]
#' X_train$proteomic <- X_train$proteomic[index_train,1:20]
#' Y_train <- Y_multiomic[index_train,]
#' vector <- list()
#' vector$mirna <- c(10)
#' vector$proteomic <- c(10)
#' cv.isb.splsdacox_model <- cv.isb.splsdacox(X_train, Y_train, max.ncomp = 1, vector = vector,
#' n_run = 1, k_folds = 3, x.center = TRUE, x.scale = TRUE)

cv.isb.splsdacox <- function(X, Y,
                                     max.ncomp = 8, vector = NULL,
                                     MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5,
                                     MIN_AUC_INCREASE = 0.01,
                                     EVAL_METHOD = "AUC",
                                     n_run = 3, k_folds = 10,
                                     x.center = TRUE, x.scale = FALSE,
                                     remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                                     remove_variance_at_fold_level = FALSE,
                                     remove_non_significant_models = FALSE, remove_non_significant = FALSE,
                                     alpha = 0.05,
                                     w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL,
                                     max_time_points = 15,
                                     MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                                     pred.attr = "mean", pred.method = "cenROC", fast_mode = FALSE,
                                     max.iter = 200,
                                     MIN_EPV = 5, return_models = FALSE, returnData = FALSE,
                                     PARALLEL = FALSE, n_cores = NULL, verbose = FALSE, seed = 123){
  # tol Numeric. Tolerance for solving: solve(t(P) %*% W) (default: 1e-15).
  tol = 1e-10

  t1_true <- Sys.time()
  y.center = y.scale = FALSE
  FREQ_CUT <- 95/5

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

  #### Check colnames
  X <- checkColnamesIllegalChars.mb(X)

  #### REQUIREMENTS
  checkX.colnames.mb(X)
  checkY.colnames(Y)
  lst_check <- checkXY.mb.class(X, Y, verbose = verbose)
  X <- lst_check$X
  Y <- lst_check$Y

  check.cv.weights(c(w_AIC, w_C.Index, w_I.BRIER, w_AUC))
  max.ncomp <- check.mb.ncomp(X, max.ncomp)
  # if(!pred.method %in% c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I")){
  #   stop_quietly(paste0("pred.method must be one of the following: ", paste0(c("risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C", "smoothROCtime_I"), collapse = ", ")))
  # }
  if(!pred.method %in% pkg.env$AUC_evaluators){
    stop_quietly(paste0("pred.method must be one of the following: ", paste0(pkg.env$AUC_evaluators, collapse = ", ")))
  }

  #### Original data
  X_original <- X
  Y_original <- Y

  time <- Y[,"time"]
  event <- Y[,"event"]

  #### SCALE
  if(length(x.center)==1){
    x.center <- rep(x.center, length(names(X)))
    names(x.center) <- names(X)
  }
  if(length(x.scale)==1){
    x.scale <- rep(x.scale, length(names(X)))
    names(x.scale) <- names(X)
  }

  #### ZERO VARIANCE - ALWAYS
  if(!remove_variance_at_fold_level & (remove_near_zero_variance | remove_zero_variance)){
    lst_dnz <- deleteZeroOrNearZeroVariability.mb(X = X,
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
  max.ncomp <- check.mb.maxPredictors(X, Y, MIN_EPV, max.ncomp, verbose = verbose)
  if(MIN_COMP_TO_CHECK >= max.ncomp){
    MIN_COMP_TO_CHECK = max(max.ncomp-1, 1)
  }

  # CREATE INDIVIDUAL MODELS
  lst_cv.sb.spls <- list()
  for(b in names(Xh)){

    message(paste0("\nRunning cross validation ", pkg.env$sb.splsdacox_dynamic ," for block: ", b, "\n"))
    t1 <- Sys.time()

    # vector as list in SB
    aux_vector <- NULL
    if(b %in% names(vector)){
      aux_vector <- vector[[b]]
    }else{
      aux_vector <- vector
    }

    lst_cv.sb.spls[[b]] <- cv.splsdacox(X = Xh[[b]], Y = Yh,
                                                     max.ncomp = max.ncomp,
                                                     vector = aux_vector,
                                                     MIN_NVAR = MIN_NVAR, MAX_NVAR = MAX_NVAR, n.cut_points = n.cut_points,
                                                     MIN_AUC_INCREASE = MIN_AUC_INCREASE,
                                                     EVAL_METHOD = EVAL_METHOD,
                                                     n_run = n_run, k_folds = k_folds, alpha = alpha, remove_non_significant_models = remove_non_significant_models,
                                                     remove_non_significant = remove_non_significant,
                                                     w_AIC = w_AIC, w_C.Index = w_C.Index, w_I.BRIER = w_I.BRIER, w_AUC = w_AUC, times = times, max_time_points = max_time_points,
                                                     MIN_AUC = MIN_AUC, MIN_COMP_TO_CHECK = MIN_COMP_TO_CHECK,
                                                     x.scale = x.scale[[b]], x.center = x.center[[b]],
                                                     #y.scale = y.scale, y.center = y.center,
                                                     remove_near_zero_variance = remove_variance_at_fold_level, remove_zero_variance = FALSE, toKeep.zv = NULL,
                                                     remove_variance_at_fold_level = remove_variance_at_fold_level,
                                                     fast_mode = fast_mode, return_models = return_models,
                                                     MIN_EPV = MIN_EPV, verbose = verbose,
                                                     pred.attr = pred.attr, pred.method = pred.method, seed = seed, PARALLEL = PARALLEL, n_cores = n_cores, returnData = FALSE)
    t2 <- Sys.time()
    time <- difftime(t2,t1,units = "mins")
    if(verbose){
      message(paste0("\nTime for block ", b , ": ", as.character(round(time, 2)), "\n"))
    }
  }

  #### ### #
  # RETURN #
  #### ### #
  func_call <- match.call()

  t2_true <- Sys.time()
  time <- difftime(t2_true,t1_true,units = "mins")

  # invisible(gc())
  return(cv.isb.splsdacox_dynamic_class(list(X = list("data" = if(returnData) X_norm else NA,
                                                   "x.mean" = xmeans, "x.sd" = xsds),
                                          Y = list("data" = Yh,
                                                   "y.mean" = ymeans, "y.sd" = ysds),
                                          list_cv_spls_models = lst_cv.sb.spls,
                                          call = func_call,
                                          X_input = if(returnData) X_original else NA,
                                          Y_input = if(returnData) Y_original else NA,
                                          alpha = alpha,
                                          nzv = variablesDeleted,
                                          nz_coeffvar = variablesDeleted_cvar,
                                          class = pkg.env$cv.isb.splsdacox_dynamic,
                                          time = time)))
}

### ## ##
# CLASS #
### ## ##

isb.splsdacox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$isb.splsdacox_dynamic)
  return(model)
}

cv.isb.splsdacox_dynamic_class = function(pls_model, ...) {
  model = structure(pls_model, class = pkg.env$model_class,
                    model = pkg.env$cv.isb.splsdacox_dynamic)
  return(model)
}
