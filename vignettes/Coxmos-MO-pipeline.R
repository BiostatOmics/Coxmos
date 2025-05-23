## ---- include = FALSE---------------------------------------------------------
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

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("Coxmos")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("BiostatOmics/Coxmos", build_vignettes = TRUE)

## ----setup, results = "hide"--------------------------------------------------
library(Coxmos)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("RColorConesa")
#  library(RColorConesa)

## -----------------------------------------------------------------------------
# load dataset
data("X_multiomic", package = "Coxmos")
data("Y_multiomic", package = "Coxmos")

X <- X_multiomic
Y <- Y_multiomic

rm(X_multiomic, Y_multiomic)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(X$mirna[1:5,1:5]);knitr::kable(X$proteomic[1:5,1:5])
knitr::kable(Y[1:5,])

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
EPV <- getEPV.mb(X, Y)

## -----------------------------------------------------------------------------
EPV

## ---- eval=FALSE--------------------------------------------------------------
#  X$mirna <- factorToBinary(X = X$mirna, all = TRUE, sep = "_")

## -----------------------------------------------------------------------------
split_data <- getTrainTest(X, Y, p = 0.7)

X_train <- split_data$X_train #106x642 and 106x369
Y_train <- split_data$Y_train

X_test <- split_data$X_test #44x642 and 44x369
Y_test <- split_data$Y_test

## ---- warning=T, eval=F-------------------------------------------------------
#  cv.sb.splsicox_res <- cv.mb.coxmos(method = "sb.splsicox",
#                                     X = X_train, Y = Y_train,
#                                     max.ncomp = 2, penalty.list = c(0.5,0.9),
#                                     n_run = 1, k_folds = 5)

## ---- eval=FALSE--------------------------------------------------------------
#  cv.sb.splsicox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.sb.splsicox_res$plot_AUC

## -----------------------------------------------------------------------------
sb.splsicox_model <- mb.coxmos(method = "sb.splsicox",
                               X = X_train, Y = Y_train,
                               n.comp = 2, #cv.sb.splsicox_res$opt.comp
                               penalty = 0.9) #cv.sb.splsicox_res$opt.penalty

sb.splsicox_model

## -----------------------------------------------------------------------------
sb.splsicox_model <- mb.coxmos(method = "sb.splsicox",
                               X = X_train, Y = Y_train,
                               n.comp = 2, #cv.sb.splsicox_res$opt.comp
                               penalty = 0.9, #cv.sb.splsicox_res$opt.penalty
                               remove_non_significant = TRUE)

sb.splsicox_model

## ----warning=T, eval=F--------------------------------------------------------
#  cv.isb.splsicox_res <- cv.mb.coxmos(method = "isb.splsicox",
#                                      X = X_train, Y = Y_train,
#                                      max.ncomp = 2, penalty.list = c(0.5,0.9),
#                                      n_run = 1, k_folds = 5,
#                                      remove_non_significant = TRUE)
#  
#  cv.isb.splsicox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.isb.splsicox_res$list_cv_spls_models$mirna$plot_AUC
#  cv.isb.splsicox_res$list_cv_spls_models$proteomic$plot_AUC

## ---- warning=F, eval=F-------------------------------------------------------
#  isb.splsicox_model <- mb.coxmos(method = "isb.splsicox",
#                                  X = X_train, Y = Y_train, cv.isb = cv.isb.splsicox_res,
#                                  remove_non_significant = FALSE)
#  
#  isb.splsicox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.sb.splsdrcox_res <- cv.mb.coxmos(method = "sb.splsdrcox",
#                                      X = X_train, Y = Y_train,
#                                      max.ncomp = 2, vector = NULL,
#                                      n_run = 1, k_folds = 5)
#  
#  cv.sb.splsdrcox_res

## -----------------------------------------------------------------------------
sb.splsdrcox_model <- mb.coxmos(method = "sb.splsdrcox", 
                                X = X_train, Y = Y_train, 
                                n.comp = 2, #cv.sb.splsdrcox_res$opt.comp, 
                                vector = list("mirna" = 161, "proteomic" = 185), #cv.sb.splsdrcox_res$opt.nvar,
                                remove_non_significant = T)

sb.splsdrcox_model

## ----warning=T, eval=F--------------------------------------------------------
#  cv.isb.splsdrcox_res <- cv.mb.coxmos(method = "isb.splsdrcox",
#                                       X = X_train, Y = Y_train,
#                                       max.ncomp = 2, vector = NULL,
#                                       n_run = 1, k_folds = 5)
#  
#  cv.isb.splsdrcox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.isb.splsdrcox_res$list_cv_spls_models$mirna$plot_AUC
#  cv.isb.splsdrcox_res$list_cv_spls_models$proteomic$plot_AUC

## ---- warning=F, eval=F-------------------------------------------------------
#  isb.splsdrcox_model <- mb.coxmos(method = "isb.splsdrcox",
#                                   X = X_train, Y = Y_train, cv.isb = cv.isb.splsdrcox_res,
#                                   remove_non_significant = TRUE)
#  
#  isb.splsdrcox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.sb.splsdacox_res <- cv.mb.coxmos(method = "sb.splsdacox",
#                                      X = X_train, Y = Y_train,
#                                      max.ncomp = 2, vector = NULL,
#                                      n_run = 1, k_folds = 5)
#  
#  cv.sb.splsdacox_res

## -----------------------------------------------------------------------------
sb.splsdacox_model <- mb.coxmos(method = "sb.splsdacox", 
                                X = X_train, Y = Y_train, 
                                n.comp = 1, #cv.sb.splsdacox_res$opt.comp, 
                                vector = list("mirna" = 321, "proteomic" = 369), #cv.sb.splsdacox_res$opt.nvar,
                                remove_non_significant = T)

sb.splsdacox_model

## ----warning=T, eval=F--------------------------------------------------------
#  cv.isb.splsdacox_res <- cv.mb.coxmos(method = "isb.splsdacox",
#                                       X = X_train, Y = Y_train,
#                                       max.ncomp = 2, vector = NULL,
#                                       n_run = 1, k_folds = 5)
#  
#  cv.isb.splsdacox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.isb.splsdacox_res$list_cv_spls_models$mirna$plot_AUC
#  cv.isb.splsdacox_res$list_cv_spls_models$proteomic$plot_AUC

## ---- warning=F, eval=F-------------------------------------------------------
#  isb.splsdacox_model <- mb.coxmos(method = "isb.splsdacox",
#                                   X = X_train, Y = Y_train, cv.isb = cv.isb.splsdacox_res,
#                                   remove_non_significant = TRUE)
#  
#  isb.splsdacox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.mb.splsdrcox_res <- cv.mb.coxmos(method = "mb.splsdrcox",
#                                      X = X_train, Y = Y_train,
#                                      max.ncomp = 2, vector = NULL,
#                                      MIN_NVAR = 10, MAX_NVAR = NULL, n.cut_points = 10, EVAL_METHOD = "AUC",
#                                      n_run = 1, k_folds = 5)
#  
#  cv.mb.splsdrcox_res

## -----------------------------------------------------------------------------
mb.splsdrcox_model <- mb.coxmos(method = "mb.splsdrcox", 
                                X = X_train, Y = Y_train, 
                                n.comp = 1, #cv.mb.splsdrcox_res$opt.comp,
                                vector = list("mirna" = 10, "proteomic" = 369)) #cv.mb.splsdrcox_res$opt.nvar
mb.splsdrcox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.mb.splsdacox_res <- cv.mb.coxmos(method = "mb.splsdacox",
#                                      X = X_train, Y = Y_train,
#                                      max.ncomp = 2, vector = NULL,
#                                      n_run = 1, k_folds = 5)
#  
#  cv.mb.splsdacox_res

## -----------------------------------------------------------------------------
mb.splsdacox_model <- mb.coxmos(method = "mb.splsdacox", 
                                X = X_train, Y = Y_train, 
                                n.comp = 2, #cv.mb.splsdacox_res$opt.comp,
                                vector = list("mirna" = 1, "proteomic" = 93)) #cv.mb.splsdacox_res$opt.nvar

mb.splsdacox_model

## -----------------------------------------------------------------------------
lst_models <- list("SB.sPLS-ICOX" = sb.splsicox_model,
                   #"iSB.sPLS-ICOX" = isb.splsicox_model,
                   "SB.sPLS-DRCOX-Dynamic" = sb.splsdrcox_model,
                   #"iSB.sPLS-DRCOX-Dynamic" = isb.splsdrcox_model,
                   "SB.sPLS-DACOX-Dynamic" = sb.splsdacox_model,
                   #"iSB.sPLS-DACOX-Dynamic" = isb.splsdacox_model,
                   "MB.sPLS-DRCOX" = mb.splsdrcox_model,
                   "MB.sPLS-DACOX" = mb.splsdacox_model)

eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                  X_test = X_test, Y_test = Y_test, 
                                  times = NULL)

## -----------------------------------------------------------------------------
eval_results

## ---- warning=F---------------------------------------------------------------
lst_eval_results_auc <- plot_evaluation(eval_results, evaluation = "AUC", pred.attr = "mean")
lst_eval_results_brier <- plot_evaluation(eval_results, evaluation = "IBS", pred.attr = "mean")

## ---- fig.small=T, warning=F--------------------------------------------------
lst_eval_results_auc$lst_plots$lineplot.mean
lst_eval_results_auc$lst_plot_comparisons$anova

## -----------------------------------------------------------------------------
lst_models_time <- list(#cv.sb.splsicox_res,
                        sb.splsicox_model,
                        #isb.splsicox_model,
                        #cv.sb.splsdrcox_res,
                        sb.splsdrcox_model,
                        #isb.splsdrcox_model,
                        #cv.sb.splsdacox_res,
                        sb.splsdacox_model,
                        #isb.splsdacox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdrcox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdacox_model,
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.list(lst_models_time, txt.x.angle = 90)

## ---- fig.small=T-------------------------------------------------------------
ggp_time

## ---- results='hide'----------------------------------------------------------
lst_ph_ggplot <- plot_proportionalHazard(lst_models$`SB.sPLS-ICOX`)

## ---- fig.small=T, warning=F--------------------------------------------------
lst_ph_ggplot

## -----------------------------------------------------------------------------
lst_forest_plot <- plot_forest(lst_models$`SB.sPLS-ICOX`)

## ---- fig.small=T-------------------------------------------------------------
lst_forest_plot

## -----------------------------------------------------------------------------
density.plots.lp <- plot_cox.event(lst_models$`SB.sPLS-ICOX`, type = "lp")

## ---- fig.small=T-------------------------------------------------------------
density.plots.lp$plot.density
density.plots.lp$plot.histogram

## ---- warning=F---------------------------------------------------------------
variable_auc_results <- eval_Coxmos_model_per_variable(model = lst_models$`SB.sPLS-ICOX`, 
                                                      X_test = lst_models$`SB.sPLS-ICOX`$X_input, 
                                                      Y_test = lst_models$`SB.sPLS-ICOX`$Y_input)

variable_auc_plot_train <- plot_evaluation(variable_auc_results, evaluation = "AUC")

## ---- fig.small=T, warning=F--------------------------------------------------
variable_auc_plot_train$lst_plots$lineplot.mean

## -----------------------------------------------------------------------------
ggp_scores <- plot_sPLS_Coxmos(model = lst_models$`SB.sPLS-ICOX`,
                               comp = c(1,2), mode = "scores")

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_scores

## -----------------------------------------------------------------------------
ggp_loadings <- plot_sPLS_Coxmos(model = lst_models$`SB.sPLS-ICOX`, 
                                 comp = c(1,2), mode = "loadings",
                                 top = 10)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_loadings

## -----------------------------------------------------------------------------
ggp_biplot <- plot_sPLS_Coxmos(model = lst_models$`SB.sPLS-ICOX`, 
                               comp = c(1,2), mode = "biplot",
                               top = 15,
                               only_top = T,
                               overlaps = 20)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_biplot

## ---- warning=FALSE-----------------------------------------------------------
ggp.simulated_beta <- plot_pseudobeta(model = lst_models$`SB.sPLS-ICOX`, 
                                      error.bar = T, onlySig = F, alpha = 0.05, 
                                      zero.rm = T, auto.limits = T, top = 20,
                                      show_percentage = T, size_percentage = 2)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp.simulated_beta$plot

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp.simulated_beta$mb_plot$plot

## ---- warning=F---------------------------------------------------------------
LST_KM_RES_LP <- getAutoKM(type = "LP",
                           model = lst_models$`SB.sPLS-ICOX`)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_LP$LST_PLOTS$LP

## ---- warning=F---------------------------------------------------------------
LST_KM_RES_COMP <- getAutoKM(type = "COMP",
                             model = lst_models$`SB.sPLS-ICOX`,
                             comp = 1:2)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_COMP$LST_PLOTS$mirna$comp_1
LST_KM_RES_COMP$LST_PLOTS$proteomic$comp_2

## ---- warning=F---------------------------------------------------------------
LST_KM_RES_VAR <- getAutoKM(type = "VAR",
                            model = lst_models$`SB.sPLS-ICOX`,
                            top = 10,
                            ori_data = T, #original data selected
                            only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_VAR$LST_PLOTS$mirna$hsa.minus.miR.minus.491.minus.3p
LST_KM_RES_VAR$LST_PLOTS$proteomic$var_5080

## -----------------------------------------------------------------------------
new_pat <- X_test
new_pat$mirna <- new_pat$mirna[1,,drop=F]
new_pat$proteomic <- new_pat$proteomic[1,,drop=F]

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat$mirna),])

## -----------------------------------------------------------------------------
pat_density <- plot_observation.eventDensity(observation = new_pat,
                                             model = lst_models$`SB.sPLS-ICOX`,
                                             type = "lp")

## ---- fig.small=T, warning=F--------------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_observation.eventHistogram(observation = new_pat, 
                                                 model = lst_models$`SB.sPLS-ICOX`, 
                                                 type = "lp")

## ---- fig.small=T, warning=F--------------------------------------------------
pat_histogram

## ---- warning=F---------------------------------------------------------------
ggp.simulated_beta_newObs <- plot_observation.pseudobeta(model = lst_models$`SB.sPLS-ICOX`,
                                                         observation = new_pat,
                                                         error.bar = T, onlySig = T, alpha = 0.05,
                                                         zero.rm = T, auto.limits = T, show.betas = T, top = 20,
                                                         txt.x.angle = 90)

## ---- fig.small=T, warning=F--------------------------------------------------
ggp.simulated_beta_newObs$plot

## ---- warning=FALSE-----------------------------------------------------------
variable_auc_results_test <- eval_Coxmos_model_per_variable(model = lst_models$`SB.sPLS-ICOX`, 
                                                       X_test = X_test, 
                                                       Y_test = Y_test)

variable_auc_plot_test <- plot_evaluation(variable_auc_results_test, evaluation = "AUC")

## ---- fig.small=T, warning=F--------------------------------------------------
variable_auc_plot_test$lst_plots$lineplot.mean

## -----------------------------------------------------------------------------
lst_observations <- list()
for(b in names(X_test)){
  lst_observations[[b]] <- X_test[[b]][1:5,]
}

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(lst_observations$mirna),])

## ---- warning=F---------------------------------------------------------------
lst_cox.comparison <- plot_multipleObservations.LP(model = lst_models$`SB.sPLS-ICOX`,
                                                   observations = lst_observations,
                                                   top = 10)

## ---- fig.small=T, warning=F--------------------------------------------------
lst_cox.comparison$plot

## ---- warning=FALSE-----------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM(model = lst_models$`SB.sPLS-ICOX`, 
                            X_test = X_test, Y_test = Y_test, 
                            type = "LP",
                            BREAKTIME = NULL, n.breaks = 20,
                            cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
LST_KM_TEST_LP

## ---- warning=FALSE-----------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM(model = lst_models$`SB.sPLS-ICOX`, 
                              X_test = X_test, Y_test = Y_test, 
                              type = "COMP",
                              BREAKTIME = NULL, n.breaks = 20,
                              cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
LST_KM_TEST_COMP$comp_1_mirna
LST_KM_TEST_COMP$comp_2_proteomic

## ---- warning=FALSE-----------------------------------------------------------
lst_cutoff <- getCutoffAutoKM(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM(model = lst_models$`SB.sPLS-ICOX`, 
                             X_test = X_test, Y_test = Y_test, 
                             type = "VAR", ori_data = T,
                             BREAKTIME = NULL, n.breaks = 20,
                             cutoff = lst_cutoff)

## ---- fig.small=T, warning=F--------------------------------------------------
LST_KM_TEST_VAR$mirna$hsa.minus.miR.minus.491.minus.3p

