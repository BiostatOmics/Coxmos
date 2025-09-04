#' Multiblock COX Modeling Function
#'
#' @description This function provides a unified interface for all multiblock HD-COX modeling methods in the package.
#'
#' @param method Modeling method to use: "sb.splsicox", "sb.splsdrcox", "sb.splsdrcox_penalty", "sb.splsdacox",
#' "isb.splsicox", "isb.splsdrcox", "isb.splsdrcox_penalty", "isb.splsdacox",
#' "mb.splsdrcox", or "mb.splsdacox".
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param cv.isb Instance of class "Coxmos" and isb model. Used to retrieve the
#' optimal components and variables for the sPLS Cox model
#' (isb.splsicox, isb.splsdrcox, isb.splsdrcox_penalty and isb.splsdacox).
#' @param design Numeric matrix. Matrix of size (number of blocks in X) x (number of blocks in X) with
#' values between 0 and 1. Each value indicates the strength of the relationship to be modeled between
#' two blocks; a value of 0 indicates no relationship, 1 is the maximum value. If NULL, auto-design is
#' computed (default: NULL) (mb.splsdrcox and mb.splsdacox).
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
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the threshold
#' (default: 0.05).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model. Used to restrict the number of variables/components can be computed in final cox models.
#' If the minimum is not meet, the model cannot be computed (default: 5).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#'
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 4).
#' @param penalty Numeric. Penalty for variable selection for the individual cox models. Variables
#' with a lower P-Value than 1 - "penalty" in the individual cox analysis will be keep for the approach
#' (default: 0.5) (sb.splsicox and sb.splsdrcox_penalty).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL)
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 1) (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use. If NULL, the number of variables is selected (default: NULL)
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
#' as n.cut_points^n.blocks models will be computed as minimum (default: 5)
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01) (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: "AUC")
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC")  (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200)
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL)
#' (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15) (sb.splsdrcox, sb.splsdacox, mb.splsdrcox and mb.splsdacox).
#'
#' @return A Coxmos model of the specified multiblock type.
#'
#' @seealso
#' \code{\link{sb.splsicox}} for Single-Block SPLS-ICOX,
#' \code{\link{sb.splsdrcox_penalty}} for Single-Block SPLS-DRCOX with penalty,
#' \code{\link{sb.splsdrcox}} for Single-Block SPLS-DRCOX,
#' \code{\link{sb.splsdacox}} for Single-Block SPLS-DACOX,
#' \code{\link{isb.splsicox}} for Integrated Single-Block SPLS-ICOX,
#' \code{\link{isb.splsdrcox_penalty}} for Integrated Single-Block SPLS-DRCOX with penalty,
#' \code{\link{isb.splsdrcox}} for Integrated Single-Block SPLS-DRCOX,
#' \code{\link{isb.splsdacox}} for Integrated Single-Block SPLS-DACOX,
#' \code{\link{mb.splsdrcox}} for Multi-Block SPLS-DRCOX,
#' \code{\link{mb.splsdacox}} for Multi-Block SPLS-DACOX
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' data("Y_multiomic")
#' X <- X_multiomic
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_multiomic$event, p = .25, list = FALSE, times = 1)
#' X_train <- X_multiomic
#' X_train$mirna <- X_train$mirna[index_train,1:30]
#' X_train$proteomic <- X_train$proteomic[index_train,1:30]
#' Y_train <- Y_multiomic[index_train,]
#' mb.coxmos(method = "sb.splsicox", X_train, Y_train, n.comp = 2, x.center = TRUE, x.scale = FALSE)

mb.coxmos <- function(method = c("sb.splsicox", "sb.splsdrcox", "sb.splsdrcox_penalty", "sb.splsdacox",
                                 "isb.splsicox", "isb.splsdrcox", "isb.splsdrcox_penalty", "isb.splsdacox",
                                 "mb.splsdrcox", "mb.splsdacox"),
                              X, Y,
                              # Argumentos específicos de familia
                              cv.isb = NULL, design = NULL,
                              # Argumentos comunes
                              x.center = TRUE, x.scale = FALSE,
                              remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                              remove_non_significant = FALSE, alpha = 0.05,
                              MIN_EPV = 5, returnData = TRUE, verbose = FALSE,
                              # Argumentos específicos de método
                              n.comp = 4, penalty = 0.5, vector = NULL,
                              MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5,
                              MIN_AUC_INCREASE = 0.01, EVAL_METHOD = "AUC",
                              pred.method = "cenROC", max.iter = 200,
                              times = NULL, max_time_points = 15) {

  method <- match.arg(method)

  # Objective function
  target_func <- get(method)

  # Formal arguments for objective function
  valid_args <- formalArgs(target_func)

  # Get all current arguments
  all_args <- as.list(match.call())[-1]

  # Filter current arguments
  filtered_args <- all_args[names(all_args) %in% valid_args]

  # Run the function
  model <- do.call(target_func, filtered_args)

  return(model)
}

#' Multiblock COX Cross-Validation Function
#'
#' @description This function provides a unified interface for all multiblock HD-COX cross-validation methods in the package.
#'
#' @param method Cross-validation method to use: "sb.splsicox", "sb.splsdrcox", "sb.splsdrcox_penalty", "sb.splsdacox",
#' "isb.splsicox", "isb.splsdrcox", "isb.splsdrcox_penalty", "isb.splsdacox", "mb.splsdrcox", or "mb.splsdacox".
#' @param X List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param max.ncomp Numeric. Maximum number of PLS components to compute for the cross validation
#' (default: 8).
#' @param penalty.list Numeric vector. Penalty for variable selection for the individual cox
#' models. Variables with a lower P-Value than 1 - "penalty" in the individual cox analysis will
#' be keep for the sPLS-ICOX approach (default: seq(0.1,0.9,0.2)) (sb.splsicox, sb.splsdrcox_penalty, isb.splsicox, isb.splsdrcox_penalty).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 1) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use (default: NULL) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used (default: 5) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
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
#' variance will be removed at fold level (default: FALSE).
#' @param remove_non_significant_models Logical. If remove_non_significant_models = TRUE,
#' non-significant models are removed before computing the evaluation (default: FALSE).
#' @param remove_non_significant Logical. If remove_non_significant = TRUE, non-significant
#' variables/components in final cox model will be removed (default: FALSE).
#' @param alpha Numeric. Numerical values are regarded as significant if they fall below the
#' threshold (default: 0.05).
#' @param w_AIC Numeric. Weight for AIC evaluator (default: 0).
#' @param w_C.Index Numeric. Weight for C-Index evaluator (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator (default: 1).
#' @param w_I.BRIER Numeric. Weight for BRIER SCORE evaluator (default: 0).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values (default: 0.01).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables (default: "AUC") (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param pred.method Character. AUC evaluation algorithm method (default: "cenROC") (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param pred.attr Character. Way to evaluate the metric selected (default: "mean").
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC
#' improves (default: 3).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param times Numeric vector. Time points where the AUC will be evaluated (default: NULL) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15) (sb.splsdrcox, sb.splsdacox, isb.splsdrcox, isb.splsdacox, mb.splsdrcox, mb.splsdacox).
#' @param design Numeric matrix. Matrix of size (number of blocks in X) x (number of blocks in X) with
#' values between 0 and 1 (default: NULL) (mb.splsdrcox and mb.splsdacox).
#' @param fast_mode Logical. If fast_mode = TRUE, for each run, only one fold is evaluated
#' simultaneously (default: FALSE).
#' @param MIN_EPV Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
#' cox model (default: 5).
#' @param return_models Logical. Return all models computed in cross validation (default: FALSE).
#' @param returnData Logical. Return original and normalized X and Y matrices (default: TRUE).
#' @param PARALLEL Logical. Run the cross validation with multicore option (default: FALSE).
#' @param n_cores Numeric. Number of cores to use for parallel processing. This parameter is only
#' used if `PARALLEL` is `TRUE`. If `NULL`, it will use all available cores minus one. Otherwise,
#' it will use the minimum between the value specified and the total number of cores - 1. The fewer
#' cores used, the less RAM memory will be used.(default: NULL).
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @return A cross-validation object of the specified multiblock type.
#'
#' @seealso
#' \code{\link{cv.sb.splsicox}} for Single-Block SPLS-ICOX cross-validation,
#' \code{\link{cv.sb.splsdrcox_penalty}} for Single-Block SPLS-DRCOX with penalty cross-validation,
#' \code{\link{cv.sb.splsdrcox}} for Single-Block SPLS-DRCOX cross-validation,
#' \code{\link{cv.sb.splsdacox}} for Single-Block SPLS-DACOX cross-validation,
#' \code{\link{cv.isb.splsicox}} for Integrated Single-Block SPLS-ICOX cross-validation,
#' \code{\link{cv.isb.splsdrcox_penalty}} for Integrated Single-Block SPLS-DRCOX with penalty cross-validation,
#' \code{\link{cv.isb.splsdrcox}} for Integrated Single-Block SPLS-DRCOX cross-validation,
#' \code{\link{cv.isb.splsdacox}} for Integrated Single-Block SPLS-DACOX cross-validation,
#' \code{\link{cv.mb.splsdrcox}} for Multi-Block SPLS-DRCOX cross-validation,
#' \code{\link{cv.mb.splsdacox}} for Multi-Block SPLS-DACOX cross-validation
#'
#' @export
#'
#' @examples
#' data("X_multiomic")
#' data("Y_multiomic")
#' set.seed(123)
#' X_train <- X_multiomic
#' X_train$mirna <- X_train$mirna[1:20,1:30]
#' X_train$proteomic <- X_train$proteomic[1:20,1:30]
#' Y_train <- Y_multiomic[1:20,]
#' cv_mb <- cv.mb.coxmos(method = "sb.splsicox", X = X_train, Y = Y_train,
#'                      max.ncomp = 1, n_run = 1, k_folds = 2)
cv.mb.coxmos <- function(method = c("sb.splsicox", "sb.splsdrcox", "sb.splsdrcox_penalty", "sb.splsdacox",
                                    "isb.splsicox", "isb.splsdrcox", "isb.splsdrcox_penalty", "isb.splsdacox",
                                    "mb.splsdrcox", "mb.splsdacox"),
                         X, Y,
                         max.ncomp = 8,
                         penalty.list = seq(0.1, 0.9, 0.2),
                         vector = NULL,
                         MIN_NVAR = 1,
                         MAX_NVAR = NULL,
                         n.cut_points = 5,
                         n_run = 3,
                         k_folds = 10,
                         x.center = TRUE,
                         x.scale = FALSE,
                         remove_near_zero_variance = TRUE,
                         remove_zero_variance = TRUE,
                         toKeep.zv = NULL,
                         remove_variance_at_fold_level = FALSE,
                         remove_non_significant_models = FALSE,
                         remove_non_significant = FALSE,
                         alpha = 0.05,
                         w_AIC = 0,
                         w_C.Index = 0,
                         w_AUC = 1,
                         w_I.BRIER = 0,
                         MIN_AUC_INCREASE = 0.01,
                         EVAL_METHOD = "AUC",
                         pred.method = "cenROC",
                         pred.attr = "mean",
                         MIN_AUC = 0.8,
                         MIN_COMP_TO_CHECK = 3,
                         max.iter = 200,
                         times = NULL,
                         max_time_points = 15,
                         design = NULL,
                         fast_mode = FALSE,
                         MIN_EPV = 5,
                         return_models = FALSE,
                         returnData = FALSE,
                         PARALLEL = FALSE, n_cores = NULL,
                         verbose = FALSE,
                         seed = 123) {

  method <- match.arg(method)

  # Objective function
  target_func <- get(paste0("cv.", method))

  # Formal arguments for objective function
  valid_args <- formalArgs(target_func)

  # Get all current arguments
  all_args <- as.list(match.call())[-1]

  # Filter current arguments
  filtered_args <- all_args[names(all_args) %in% valid_args]

  # Run the function
  cv.model <- do.call(target_func, filtered_args)

  return(cv.model)
}
