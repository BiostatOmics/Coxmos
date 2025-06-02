#' Coxmos Modeling Function
#'
#' @description This function provides a unified interface for all HD-COX modeling methods in the package.
#'
#' @param method Modeling method to use: "cox", "coxSW", "coxEN", "splsicox", "splsdrcox", "splsdrcox_penalty", or "splsdacox" (default: cox).
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
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
#' @param max.variables Numeric. Maximum number of variables you want to keep in the cox model. If
#' MIN_EPV is not meet, the value will be change automatically (default: 20) (coxSW, coxEN).
#' @param BACKWARDS Logical. If BACKWARDS = TRUE, backward strategy is performed (default: TRUE) (coxSW).
#' @param alpha_ENT Numeric. Maximum P-Value threshold for an ANOVA test when comparing a more complex
#' model to a simpler model that includes a new variable. If the p-value is less than or equal to this
#' threshold, the new variable is considered significantly important and will be added to the model
#' (default: 0.10) (coxSW).
#' @param alpha_OUT Numeric. Minimum P-Value threshold for an ANOVA test when comparing a simpler model
#' to a more complex model that excludes an existing variable. If the p-value is greater than or equal
#' to this threshold, the existing variable is considered not significantly important and will be removed
#' from the model (default: 0.15) (coxSW).
#' @param toKeep.sw Character vector. Name of variables in X to not be deleted by Step-wise
#' selection (default: NULL) (coxSW).
#' @param initialModel Character vector. Name of variables in X to include in the initial model
#' (default: NULL) (coxSW).
#'
#' @param EN.alpha Numeric. Elastic net mixing parameter. If EN.alpha = 1 is the lasso penalty, and
#' EN.alpha = 0 the ridge penalty (default: 0.5). NOTE: When ridge penalty is used, EVP and
#' max.variables will not be used (coxEN).
#'
#' @param n.comp Numeric. Number of latent components to compute for the (s)PLS model (default: 4) (splsicox, splsdrcox, splsdrcox_penalty, splsdacox).
#' @param penalty Numeric. Penalty for variable selection for the individual cox models. Variables
#' with a lower P-Value than 1 - "penalty" in the individual cox analysis will be keep for the approach
#' (default: 0.5) (splsicox, splsdrcox_penalty).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL) (splsdrcox, splsdacox).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: 1) (splsdrcox, splsdacox).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use. If NULL, the number of variables is selected (default: NULL) (splsdrcox, splsdacox).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
#' as n.cut_points^n.blocks models will be computed as minimum (default: 5) (splsdrcox, splsdacox).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models to
#' continue evaluating higher values in the multiple tested parameters. If it is not reached for next
#' 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01) (splsdrcox, splsdacox).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: "AUC") (splsdrcox, splsdacox).
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC") (splsdrcox, splsdacox).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: 200) (splsdrcox, splsdacox).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL) (splsdrcox, splsdacox).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15) (splsdrcox, splsdacox).
#'
#' @param FORCE Logical. In case the MIN_EPV is not meet, it allows to compute the model (default: FALSE) (cox).
#'
#' @return A Coxmos model specified by the method argument.
#'
#' @seealso
#' \code{\link{cox}} for standard Cox model,
#' \code{\link{coxSW}} for stepwise selection,
#' \code{\link{coxEN}} for elastic net,
#' \code{\link{splsicox}} for SPLS-COX,
#' \code{\link{splsdrcox}} for SPLS-DRCOX,
#' \code{\link{splsdrcox_penalty}} for penalized SPLS-DRCOX,
#' \code{\link{splsdacox}} for SPLS-DACOX
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_proteomic$event, p = .5, list = FALSE, times = 1)
#' X_train <- X_proteomic[index_train,1:50]
#' Y_train <- Y_proteomic[index_train,]
#'
#' # Elastic Net Cox
#' model_en <- coxmos(method = "coxEN", X = X_train, Y = Y_train, EN.alpha = 0.7)

coxmos <- function(method = c("cox", "coxSW", "coxEN", "splsicox", "splsdrcox", "splsdrcox_penalty", "splsdacox"),
                   X, Y,
                   # Common arguments
                   x.center = TRUE, x.scale = FALSE,
                   remove_near_zero_variance = TRUE, remove_zero_variance = FALSE, toKeep.zv = NULL,
                   remove_non_significant = FALSE, alpha = 0.05,
                   MIN_EPV = 5, returnData = TRUE, verbose = FALSE,
                   # Method-specific arguments
                   max.variables = NULL, BACKWARDS = TRUE,
                   alpha_ENT = 0.1, alpha_OUT = 0.15, toKeep.sw = NULL,
                   initialModel = NULL, EN.alpha = 0.5,
                   n.comp = 4, penalty = 0.5, vector = NULL,
                   MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5,
                   MIN_AUC_INCREASE = 0.01, EVAL_METHOD = "AUC",
                   pred.method = "cenROC", max.iter = 200,
                   times = NULL, max_time_points = 15,
                   FORCE = FALSE) {

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

#' Cross-Validation for COX Models
#'
#' @description This function provides a unified interface for cross-validation of all COX modeling methods in the package.
#'
#' @param method Cross-validation method to use: "coxEN", "splsicox", "splsdrcox", "splsdrcox_penalty", or "splsdacox" (default: coxEN).
#' @param X Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
#' transform into binary variables.
#' @param Y Numeric matrix or data.frame. Response variables. Object must have two columns named as
#' "time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
#' event observations.
#' @param n_run Numeric. Number of runs for cross validation (default: 3).
#' @param k_folds Numeric. Number of folds for cross-validation (default: 10).
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
#' @param w_AIC Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).
#' @param w_C.Index Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).
#' @param w_AUC Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).
#' @param w_I.BRIER Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).
#' @param times Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
#' 'max_time_points' points will be selected equally distributed (default: NULL).
#' @param max_time_points Numeric. Maximum number of time points to use for evaluating the model
#' (default: 15).
#' @param MIN_AUC_INCREASE Numeric. Minimum improvement between different cross validation models
#' to continue evaluating higher values in the multiple tested parameters. If it is not reached for
#' next 'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
#' (default: 0.01).
#' @param MIN_AUC Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is
#' reached, the evaluation could stop if the improvement does not reach an AUC higher than adding
#' the 'MIN_AUC_INCREASE' value (default: 0.8).
#' @param MIN_COMP_TO_CHECK Numeric. Number of penalties/components to evaluate to check if the AUC
#' improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet,
#' the evaluation could stop (default: 3).
#' @param pred.attr Character. Way to evaluate the metric selected. Must be one of the following:
#' "mean" or "median" (default: "mean").
#' @param pred.method Character. AUC evaluation algorithm method for evaluate the model performance.
#' Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
#' "smoothROCtime_I" (default: "cenROC").
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
#' @param verbose Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).
#' @param seed Number. Seed value for performing runs/folds divisions (default: 123).
#'
#' @param EN.alpha.list Numeric vector. Elastic net mixing parameter values for cross-validation.
#' If NULL, default seq(0,1,0.1) will be used (default: NULL) (coxEN).
#' @param max.ncomp Numeric. Maximum number of latent components to compute for the (s)PLS model
#' (default: NULL) (splsicox, splsdrcox, splsdrcox_penalty, splsdacox).
#' @param penalty.list Numeric vector. Penalty values for variable selection.
#' If NULL, default seq(0,0.9,0.1) will be used (default: NULL) (splsicox, splsdrcox_penalty).
#' @param vector Numeric vector. Used for computing best number of variables. As many values as
#' components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL) (splsdrcox, splsdacox).
#' @param MIN_NVAR Numeric. Minimum range size for computing cut points to select the best number of
#' variables to use (default: NULL) (splsdrcox, splsdacox).
#' @param MAX_NVAR Numeric. Maximum range size for computing cut points to select the best number of
#' variables to use. If NULL, the number of variables is selected (default: NULL) (splsdrcox, splsdacox).
#' @param n.cut_points Numeric. Number of cut points for searching the optimal number of variables.
#' If only two cut points are selected, minimum and maximum size are used (default: NULL) (splsdrcox, splsdacox).
#' @param EVAL_METHOD Character. The selected metric will be use to compute the best
#' number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: NULL) (splsdrcox, splsdacox).
#' @param max.iter Numeric. Maximum number of iterations for PLS convergence (default: NULL) (splsdrcox, splsdacox).
#' @param max.variables Numeric. Maximum number of variables to keep in the cox model (default: NULL) (coxEN).
#'
#' @return A cross-validation object of the specified type.
#'
#'
#' @seealso
#' \code{\link{cv.coxEN}} for Elastic Net cross-validation,
#' \code{\link{cv.splsicox}} for SPLS-COX cross-validation,
#' \code{\link{cv.splsdrcox}} for SPLS-DRCOX cross-validation,
#' \code{\link{cv.splsdrcox_penalty}} for penalized SPLS-DRCOX cross-validation,
#' \code{\link{cv.splsdacox}} for SPLS-DACOX cross-validation
#'
#' @export
#'
#' @examples
#' data("X_proteomic")
#' data("Y_proteomic")
#' set.seed(123)
#' index_train <- caret::createDataPartition(Y_proteomic$event, p = .5, list = FALSE, times = 1)
#' X_train <- X_proteomic[index_train,1:50]
#' Y_train <- Y_proteomic[index_train,]
#'
#' # Elastic Net Cross-Validation
#' # cv_model <- cv_coxmos(method = "coxEN", X = X_train, Y = Y_train, EN.alpha.list = seq(0, 1, 0.2))

cv.coxmos <- function(method = c("coxEN", "splsicox", "splsdrcox", "splsdrcox_penalty", "splsdacox"),
                      X, Y,
                      # Argumentos comunes
                      n_run = 3, k_folds = 10,
                      x.center = TRUE, x.scale = FALSE,
                      remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
                      remove_variance_at_fold_level = FALSE,
                      remove_non_significant_models = FALSE, remove_non_significant = FALSE, alpha = 0.05,
                      w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL,
                      max_time_points = 15,
                      MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
                      pred.attr = "mean", pred.method = "cenROC", fast_mode = FALSE,
                      MIN_EPV = 5, return_models = FALSE, returnData = FALSE,
                      PARALLEL = FALSE, verbose = FALSE, seed = 123,
                      # Argumentos específicos de método
                      EN.alpha.list = seq(0,1,0.2), max.ncomp = 8, penalty.list = seq(0.1,0.9,0.2),
                      vector = NULL, MIN_NVAR = 1, MAX_NVAR = NULL, n.cut_points = 5,
                      EVAL_METHOD = "AUC", max.iter = 200, max.variables = NULL) {

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
