## ----include = FALSE----------------------------------------------------------
dpi = 125

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=dpi, 
  fig.retina=1, 
  fig.width=1440/dpi, #4:3 FHD
  fig.height=1080/dpi, 
  out.width="100%",
  crop = NULL,
  warning = T, 
  error = T
)

rm(dpi)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("Coxmos")

## ----eval = FALSE-------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("BiostatOmics/Coxmos", build_vignettes = TRUE)

## ----setup, results = "hide"--------------------------------------------------
library(Coxmos)

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("RColorConesa")
# library(RColorConesa)

## ----load data, warning=FALSE-------------------------------------------------
### load dataset
data("X_proteomic")
data("Y_proteomic")

X <- X_proteomic
Y <- Y_proteomic

rm(X_proteomic, Y_proteomic)

## ----data dimensions, echo = FALSE, warning=FALSE-----------------------------
knitr::kable(X[1:5,1:5])
knitr::kable(Y[1:5,])

## ----echo = FALSE, warning=FALSE----------------------------------------------
knitr::kable(dim(X), col.names = "X")
knitr::kable(dim(Y), col.names = "Y")

## -----------------------------------------------------------------------------
ggp_density.event <- plot_events(Y = Y, 
                                 categories = c("Censored","Event"),
                                 y.text = "Number of observations", 
                                 roundTo = 0.5, 
                                 max.breaks = 15, 
                                 txt.x.angle = 90)

## ----fig.small = T------------------------------------------------------------
ggp_density.event$plot

## -----------------------------------------------------------------------------
EPV <- getEPV(X, Y)

## -----------------------------------------------------------------------------
EPV

## ----eval=FALSE---------------------------------------------------------------
# factorToBinary(X = X, all = TRUE, sep = "_")

## -----------------------------------------------------------------------------
split_data <- getTrainTest(X, Y, p = 0.7)

X_train <- split_data$X_train #106x369
Y_train <- split_data$Y_train

X_test <- split_data$X_test #44x369
Y_test <- split_data$Y_test

## ----eval=FALSE, message=T----------------------------------------------------
# cox_model <- coxmos(method = "cox",
#                     X = X_train, Y = Y_train,
#                     x.center = T, x.scale = F,
#                     MIN_EPV = 5)

## ----message=T----------------------------------------------------------------
coxSW_model <- coxmos(method = "coxSW",
                    X = X_train, Y = Y_train, 
                    x.center = T, x.scale = F,
                    MIN_EPV = 5)

## -----------------------------------------------------------------------------
coxSW_model

## ----warning=F, eval=FALSE----------------------------------------------------
# cv.coxen_res <- cv.coxmos(method = "coxEN",
#                           X = X_train, Y = Y_train,
#                           EN.alpha.list = c(0.1, 0.5, 0.9),
#                           n_run = 2, k_folds = 5)

## ----eval=FALSE---------------------------------------------------------------
# cv.coxen_res

## ----warning=F----------------------------------------------------------------
coxen_model <- coxmos(method = "coxEN",
                      X = X_train, Y = Y_train, 
                      EN.alpha = 0.5, #cv.coxen_res$opt.EN.alpha,
                      max.variables = 8) #cv.coxen_res$opt.nvar

## -----------------------------------------------------------------------------
coxen_model

## ----warning=F----------------------------------------------------------------
coxen_model <- coxmos(method = "coxEN",
                      X = X_train, Y = Y_train, 
                      EN.alpha = 0.5, #cv.coxen_res$opt.EN.alpha
                      max.variables = 8, #cv.coxen_res$opt.nvar
                      remove_non_significant = T, alpha = 0.05)

## -----------------------------------------------------------------------------
coxen_model

## -----------------------------------------------------------------------------
coxen_model$nsv

## ----warning=F, eval=FALSE----------------------------------------------------
# cv.splsicox_res <- cv.coxmos(method = "splsicox",
#                              X = X_train, Y = Y_train,
#                              max.ncomp = 2, penalty.list = c(0.5, 0.9),
#                              n_run = 2, k_folds = 5)

## ----eval=FALSE---------------------------------------------------------------
# cv.splsicox_res

## ----fig.small=T, warning=F, eval=FALSE---------------------------------------
# ### plot cv.plsicox
# cv.splsicox_res$plot_AUC

## -----------------------------------------------------------------------------
splsicox_model <- coxmos(method = "splsicox",
                         X = X_train, Y = Y_train, 
                         n.comp = 2, #cv.splsicox_res$opt.comp, 
                         penalty  = 0.9) #cv.splsicox_res$opt.penalty

splsicox_model

## ----warning=FALSE, eval=FALSE------------------------------------------------
# cv.splsdrcox_dynamic_res <- cv.coxmos(method = "splsdrcox",
#                                       X = X_train, Y = Y_train,
#                                       max.ncomp = 2, vector = NULL,
#                                       MIN_NVAR = 1, MAX_NVAR = 50,
#                                       n.cut_points = 10,
#                                       n_run = 2, k_folds = 5)

## ----eval=FALSE---------------------------------------------------------------
# cv.splsdrcox_dynamic_res

## -----------------------------------------------------------------------------
splsdrcox_dynamic_model <- coxmos(method = "splsdrcox",
                                  X = X_train, Y = Y_train, 
                                  n.comp = 2, #cv.splsdrcox_dynamic_res$opt.comp,
                                  vector = 50, #cv.splsdrcox_dynamic_res$opt.nvar,
                                  n.cut_points = 10)

splsdrcox_dynamic_model

## ----warning=FALSE, eval=FALSE------------------------------------------------
# cv.splsdacox_dynamic_res <- cv.coxmos(method = "splsdacox",
#                                       X = X_train, Y = Y_train,
#                                       max.ncomp = 2, vector = NULL,
#                                       MIN_NVAR = 1, MAX_NVAR = 50,
#                                       n.cut_points = 10,
#                                       n_run = 2, k_folds = 5)

## ----eval=FALSE, eval=FALSE---------------------------------------------------
# cv.splsdacox_dynamic_res

## -----------------------------------------------------------------------------
splsdacox_dynamic_model <- coxmos(method = "splsdacox",
                                  X = X_train, Y = Y_train, 
                                  n.comp = 2, #cv.splsdacox_dynamic_res$opt.comp, 
                                  vector = 50) #cv.splsdacox_dynamic_res$opt.nvar,
                                  
splsdacox_dynamic_model

## -----------------------------------------------------------------------------
lst_models <- list("COX-EN" = coxen_model,
                   "sPLS-ICOX" = splsicox_model,
                   "sPLS-DRCOX-Dynamic" = splsdrcox_dynamic_model,
                   "sPLS-DACOX-Dynamic" = splsdacox_dynamic_model)

eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                   X_test = X_test, Y_test = Y_test, 
                                   times = NULL)

## -----------------------------------------------------------------------------
eval_results

## ----warning=F----------------------------------------------------------------
lst_eval_results_auc <- plot_evaluation(eval_results, evaluation = "AUC", pred.attr = "mean")
lst_eval_results_brier <- plot_evaluation(eval_results, evaluation = "IBS", pred.attr = "mean")

## ----fig.small=T, warning=F---------------------------------------------------
lst_eval_results_auc$lst_plots$lineplot.mean
lst_eval_results_auc$lst_plot_comparisons$anova

## -----------------------------------------------------------------------------
lst_models_time <- list(#cv.coxen_res,
                        coxen_model,
                        #cv.splsicox_res,
                        splsicox_model,
                        #cv.splsdrcox_dynamic_res,
                        splsdrcox_dynamic_model,
                        #cv.splsdacox_dynamic_res,
                        splsdacox_dynamic_model, 
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.list(lst_models_time, txt.x.angle = 90)

## ----fig.small=T, warning=F---------------------------------------------------
ggp_time

## ----results='hide'-----------------------------------------------------------
lst_ph_ggplot <- plot_proportionalHazard(lst_models$`sPLS-DRCOX-Dynamic`)

## ----fig.small=T, warning=F---------------------------------------------------
lst_ph_ggplot

## -----------------------------------------------------------------------------
lst_forest_plot <- plot_forest(lst_models$`sPLS-DRCOX-Dynamic`)

## ----fig.small=T, warning=F---------------------------------------------------
lst_forest_plot

## -----------------------------------------------------------------------------
density.plots.lp <- plot_cox.event(lst_models$`sPLS-DRCOX-Dynamic`, type = "lp")

## ----fig.small=T, warning=F---------------------------------------------------
density.plots.lp$plot.density
density.plots.lp$plot.histogram

## ----warning=F----------------------------------------------------------------
variable_auc_results <- eval_Coxmos_model_per_variable(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                                                       X_test = lst_models$`sPLS-DRCOX-Dynamic`$X_input, 
                                                       Y_test = lst_models$`sPLS-DRCOX-Dynamic`$Y_input)

variable_auc_plot_train <- plot_evaluation(variable_auc_results, evaluation = "AUC")

## ----fig.small=T, warning=F---------------------------------------------------
variable_auc_plot_train$lst_plots$lineplot.mean

## -----------------------------------------------------------------------------
ggp_scores <- plot_sPLS_Coxmos(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                             comp = c(1,2), mode = "scores")

## ----fig.small=T, warning=F---------------------------------------------------
ggp_scores

## -----------------------------------------------------------------------------
ggp_loadings <- plot_sPLS_Coxmos(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10)

## ----fig.small=T, warning=F---------------------------------------------------
ggp_loadings

## -----------------------------------------------------------------------------
ggp_biplot <- plot_sPLS_Coxmos(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                               comp = c(1,2), mode = "biplot",
                               top = 15,
                               only_top = T,
                               overlaps = 20)

## ----fig.small=T, warning=F---------------------------------------------------
ggp_biplot

## -----------------------------------------------------------------------------
ggp.simulated_beta <- plot_pseudobeta(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                                      error.bar = T, onlySig = T, alpha = 0.05, 
                                      zero.rm = T, auto.limits = T, top = 20,
                                      show_percentage = T, size_percentage = 2)

## ----fig.small=T, warning=F---------------------------------------------------
ggp.simulated_beta$plot

## ----warning=FALSE------------------------------------------------------------
LST_KM_RES_LP <- getAutoKM(type = "LP",
                           model = lst_models$`sPLS-DRCOX-Dynamic`)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_RES_LP$LST_PLOTS$LP

## ----warning=FALSE------------------------------------------------------------
LST_KM_RES_COMP <- getAutoKM(type = "COMP",
                             model = lst_models$`sPLS-DRCOX-Dynamic`,
                             comp = 1:2)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_RES_COMP$LST_PLOTS$comp_1
LST_KM_RES_COMP$LST_PLOTS$comp_2

## ----warning=FALSE, warning=FALSE---------------------------------------------
LST_KM_RES_VAR <- getAutoKM(type = "VAR",
                            model = lst_models$`sPLS-DRCOX-Dynamic`,
                            top = 10,
                            ori_data = T, #original data selected
                            only_sig = T, alpha = 0.05)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_RES_VAR$LST_PLOTS$var_840

## -----------------------------------------------------------------------------
new_pat <- X_test[1,,drop=FALSE]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat),])

## -----------------------------------------------------------------------------
pat_density <- plot_observation.eventDensity(observation = new_pat,
                                             model = lst_models$`sPLS-DRCOX-Dynamic`,
                                             type = "lp")

## ----fig.small=T, warning=F---------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_observation.eventHistogram(observation = new_pat, 
                                                 model = lst_models$`sPLS-DRCOX-Dynamic`, 
                                                 type = "lp")

## ----fig.small=T, warning=F---------------------------------------------------
pat_histogram

## ----warning=F----------------------------------------------------------------
ggp.simulated_beta_newObs <- plot_observation.pseudobeta(model = lst_models$`sPLS-DRCOX-Dynamic`,
                                                         observation = new_pat,
                                                         error.bar = T, onlySig = T, alpha = 0.05,
                                                         zero.rm = T, auto.limits = T, show.betas = T, top = 20,
                                                         txt.x.angle = 90)

## ----fig.small=T, warning=F---------------------------------------------------
ggp.simulated_beta_newObs$plot

## ----warning=F----------------------------------------------------------------
variable_auc_results_test <- eval_Coxmos_model_per_variable(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                                                       X_test = X_test, 
                                                       Y_test = Y_test)

variable_auc_plot_test <- plot_evaluation(variable_auc_results_test, evaluation = "AUC")

## ----fig.small=T, warning=F---------------------------------------------------
variable_auc_plot_test$lst_plots$lineplot.mean

## -----------------------------------------------------------------------------
lst_observations <- X_test[1:5,]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[1:5,])

## ----warning=F----------------------------------------------------------------
lst_cox.comparison <- plot_multipleObservations.LP(model = lst_models$`sPLS-DRCOX-Dynamic`,
                                                   observations = lst_observations,
                                                   top = 10)

## ----fig.small=T, warning=F---------------------------------------------------
lst_cox.comparison$plot

## ----warning=FALSE------------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                            X_test = X_test, Y_test = Y_test, 
                            type = "LP",
                            BREAKTIME = NULL, n.breaks = 20,
                            cutoff = lst_cutoff)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_TEST_LP

## ----warning=FALSE------------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                              X_test = X_test, Y_test = Y_test, 
                              type = "COMP",
                              BREAKTIME = NULL, n.breaks = 20,
                              cutoff = lst_cutoff)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_TEST_COMP$comp_1
LST_KM_TEST_COMP$comp_2

## ----warning=FALSE------------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM(model = lst_models$`sPLS-DRCOX-Dynamic`, 
                             X_test = X_test, Y_test = Y_test, 
                             type = "VAR", ori_data = T,
                             BREAKTIME = NULL, n.breaks = 20,
                             cutoff = lst_cutoff)

## ----fig.small=T, warning=F---------------------------------------------------
LST_KM_TEST_VAR$var_840

