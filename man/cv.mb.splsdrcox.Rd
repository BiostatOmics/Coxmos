% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_mb_splsdrcox_dynamic.R
\name{cv.mb.splsdrcox}
\alias{cv.mb.splsdrcox}
\title{MB.sPLS-DRCOX Cross-Validation}
\usage{
cv.mb.splsdrcox(
  X,
  Y,
  max.ncomp = 8,
  vector = NULL,
  design = NULL,
  MIN_NVAR = 1,
  MAX_NVAR = NULL,
  n.cut_points = 5,
  EVAL_METHOD = "AUC",
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
  times = NULL,
  max_time_points = 15,
  MIN_AUC_INCREASE = 0.01,
  MIN_AUC = 0.8,
  MIN_COMP_TO_CHECK = 3,
  pred.attr = "mean",
  pred.method = "cenROC",
  max.iter = 200,
  fast_mode = FALSE,
  MIN_EPV = 5,
  return_models = FALSE,
  returnData = FALSE,
  PARALLEL = FALSE,
  verbose = FALSE,
  seed = 123
)
}
\arguments{
\item{X}{List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be
transform into binary variables.}

\item{Y}{Numeric matrix or data.frame. Response variables. Object must have two columns named as
"time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and
event observations.}

\item{max.ncomp}{Numeric. Maximum number of PLS components to compute for the cross validation
(default: 8).}

\item{vector}{Numeric vector. Used for computing best number of variables. As many values as
components have to be provided. If vector = NULL, an automatic detection is perform (default: NULL). If
vector is a list, must be named as the names of X param followed by the number of variables to select.}

\item{design}{Numeric matrix. Matrix of size (number of blocks in X) x (number of blocks in X) with
values between 0 and 1. Each value indicates the strength of the relationship to be modeled between
two blocks; a value of 0 indicates no relationship, 1 is the maximum value. If NULL, auto-design is computed (default: NULL).}

\item{MIN_NVAR}{Numeric. Minimum range size for computing cut points to select the best number of
variables to use (default: 1).}

\item{MAX_NVAR}{Numeric. Maximum range size for computing cut points to select the best number of
variables to use (default: NULL).}

\item{n.cut_points}{Numeric. Number of cut points for searching the optimal number of variables.
If only two cut points are selected, minimum and maximum size are used. For MB approaches as many
as n.cut_points^n.blocks models will be computed as minimum (default: 5).}

\item{EVAL_METHOD}{Character. The selected metric will be use to compute the best
number of variables. Must be one of the following: "AUC", "IBS" or "C.Index" (default: "AUC").}

\item{n_run}{Numeric. Number of runs for cross validation (default: 3).}

\item{k_folds}{Numeric. Number of folds for cross validation (default: 10).}

\item{x.center}{Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).}

\item{x.scale}{Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).}

\item{remove_near_zero_variance}{Logical. If remove_near_zero_variance = TRUE, near zero variance
variables will be removed (default: TRUE).}

\item{remove_zero_variance}{Logical. If remove_zero_variance = TRUE, zero variance variables will
be removed (default: TRUE).}

\item{toKeep.zv}{Character vector. Name of variables in X to not be deleted by (near) zero variance
filtering (default: NULL).}

\item{remove_variance_at_fold_level}{Logical. If remove_variance_at_fold_level = TRUE, (near)
zero variance will be removed at fold level. Not recommended. (default: FALSE).}

\item{remove_non_significant_models}{Logical. If remove_non_significant_models = TRUE,
non-significant models are removed before computing the evaluation. A non-significant model is a
model with at least one component/variable with a P-Value higher than the alpha cutoff.}

\item{remove_non_significant}{Logical. If remove_non_significant = TRUE, non-significant
variables/components in final cox model will be removed until all variables are significant by
forward selection (default: FALSE).}

\item{alpha}{Numeric. Numerical values are regarded as significant if they fall below the
threshold (default: 0.05).}

\item{w_AIC}{Numeric. Weight for AIC evaluator. All weights must sum 1 (default: 0).}

\item{w_C.Index}{Numeric. Weight for C-Index evaluator. All weights must sum 1 (default: 0).}

\item{w_AUC}{Numeric. Weight for AUC evaluator. All weights must sum 1 (default: 1).}

\item{w_I.BRIER}{Numeric. Weight for BRIER SCORE evaluator. All weights must sum 1 (default: 0).}

\item{times}{Numeric vector. Time points where the AUC will be evaluated. If NULL, a maximum of
'max_time_points' points will be selected equally distributed (default: NULL).}

\item{max_time_points}{Numeric. Maximum number of time points to use for evaluating the model
(default: 15).}

\item{MIN_AUC_INCREASE}{Numeric. Minimum improvement between different cross validation models to
continue evaluating higher values in the multiple tested parameters. If it is not reached for next
'MIN_COMP_TO_CHECK' models and the minimum 'MIN_AUC' value is reached, the evaluation stops
(default: 0.01).}

\item{MIN_AUC}{Numeric. Minimum AUC desire to reach cross-validation models. If the minimum is
reached, the evaluation could stop if the improvement does not reach an AUC higher than adding the
'MIN_AUC_INCREASE' value (default: 0.8).}

\item{MIN_COMP_TO_CHECK}{Numeric. Number of penalties/components to evaluate to check if the AUC
improves. If for the next 'MIN_COMP_TO_CHECK' the AUC is not better and the 'MIN_AUC' is meet, the
evaluation could stop (default: 3).}

\item{pred.attr}{Character. Way to evaluate the metric selected. Must be one of the following:
"mean" or "median" (default: "mean").}

\item{pred.method}{Character. AUC evaluation algorithm method for evaluate the model performance.
Must be one of the following: "risksetROC", "survivalROC", "cenROC", "nsROC", "smoothROCtime_C",
"smoothROCtime_I" (default: "cenROC").}

\item{max.iter}{Numeric. Maximum number of iterations for PLS convergence (default: 200).}

\item{fast_mode}{Logical. If fast_mode = TRUE, for each run, only one fold is evaluated
simultaneously. If fast_mode = FALSE, for each run, all linear predictors are computed for test
observations. Once all have their linear predictors, the evaluation is perform across all the
observations together (default: FALSE).}

\item{MIN_EPV}{Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
cox model. Used to restrict the number of variables/components can be computed in final cox models.
If the minimum is not meet, the model cannot be computed (default: 5).}

\item{return_models}{Logical. Return all models computed in cross validation (default: FALSE).}

\item{returnData}{Logical. Return original and normalized X and Y matrices (default: TRUE).}

\item{PARALLEL}{Logical. Run the cross validation with multicore option. As many cores as your
total cores - 1 will be used. It could lead to higher RAM consumption (default: FALSE).}

\item{verbose}{Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).}

\item{seed}{Number. Seed value for performing runs/folds divisions (default: 123).}
}
\value{
Instance of class "Coxmos" and model "cv.MB.sPLS-DRCOX".
\code{best_model_info}: A data.frame with the information for the best model.
\code{df_results_folds}: A data.frame with fold-level information.
\code{df_results_runs}: A data.frame with run-level information.
\code{df_results_comps}: A data.frame with component-level information (for cv.coxEN, EN.alpha
information).

\code{lst_models}: If return_models = TRUE, return a the list of all cross-validated models.
\code{pred.method}: AUC evaluation algorithm method for evaluate the model performance.

\code{opt.comp}: Optimal component selected by the best_model.
\code{opt.nvar}: Optimal number of variables selected by the best_model.
\code{design}: Design matrix used for computing the MultiBlocks models.

\code{plot_AIC}: AIC plot by each hyper-parameter.
\code{plot_C.Index}: C-Index plot by each hyper-parameter.
\code{plot_I.BRIER}: Integrative Brier Score plot by each hyper-parameter.
\code{plot_AUC}: AUC plot by each hyper-parameter.

\code{class}: Cross-Validated model class.

\code{lst_train_indexes}: List (of lists) of indexes for the observations used in each run/fold
for train the models.
\code{lst_test_indexes}: List (of lists) of indexes for the observations used in each run/fold
for test the models.

\code{time}: time consumed for running the cross-validated function.
}
\description{
The cv.mb.splsdrcox function performs cross-validation for the MB.sPLS-DRCOX model,
a specialized model for survival analysis with high-dimensional data. This function
systematically evaluates the performance of the model across different hyperparameters and
configurations to determine the optimal settings for the given data.
}
\details{
The function operates by partitioning the data into multiple subsets (folds) and
iteratively holding out one subset for validation while training on the remaining subsets. The
cross-validation process is repeated for a specified number of runs, ensuring a robust assessment
of the model's performance. The function offers flexibility in terms of the number of PLS components,
the range of variables considered, and the evaluation metrics used.

The function provides an option to center and scale the explanatory variables, which can be crucial
for ensuring consistent performance, especially when the variables are measured on different scales.
Additionally, the function incorporates features to handle near-zero and zero variance variables,
which can be problematic in high-dimensional datasets.

For model evaluation, users can choose between various metrics, including AUC, c-index, and Brier
Score. The function also allows for the specification of weights for these metrics, enabling users
to prioritize certain metrics over others based on the research context.

The function's design also emphasizes computational efficiency. It offers a parallel processing
option to expedite the cross-validation process, especially beneficial for large datasets. However,
users should be cautious about potential high RAM consumption when using this option.
}
\examples{
\donttest{
data("X_multiomic")
data("Y_multiomic")
set.seed(123)
index_train <- caret::createDataPartition(Y_multiomic$event, p = .5, list = FALSE, times = 1)
X_train <- X_multiomic
X_train$mirna <- X_train$mirna[index_train,1:50]
X_train$proteomic <- X_train$proteomic[index_train,1:50]
Y_train <- Y_multiomic[index_train,]
cv.mb.splsdrcox_model <- cv.mb.splsdrcox(X_train, Y_train, max.ncomp = 2, vector = NULL,
n_run = 1, k_folds = 2, x.center = TRUE, x.scale = TRUE)
}
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
