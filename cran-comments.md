# CRAN Submission Comments (2025-06-02)
This resubmission incorporates the following adjustments:

- New 'coxmos()' and 'cv.coxmos()' functions.
- New 'getTrainTest()' function for splitting data into train and test for single-omic and multi-omic X data.
- Updated 'coxSW()'. Now the max.variables parameter admits NULL values and it is the default value.
- Updated default values for all dynamic functions MIN_NVAR = 1 and MAX_NVAR = NULL.
- Updated EN.alpha.list default value.
- Updated 'plot_PLS_Coxmos()' function name to 'plot_sPLS_Coxmos()'. 
- Fixed 'plot_Coxmos.PLS.model()' and 'plot_Coxmos.MB.PLS.modelFurthermore()'. Outlier result was delete as not proportionate outlier information.
- Fixed 'plot_events()' description.
- Fixed 'plot_sPLS_Coxmos()' for biplot mode when selecting custom components.
- Fixed 'plot_evaluation()' for IBS evaluation. Now Y axes show "Brier Score" instead "IBS".
- Fixed 'plot_multipleObservations.LP()' for MB approaches.
- Updated 'plot_sPLS_Coxmos()' to add the ellipses parameter, allowing the user to show or hide the ellipses when coloring by factor.
- Miscellaneous changes in MB.sPLS default parameters.
- Miscellaneous changes in sPLS plot axis.

# CRAN Submission Comments (2025-03-05)
This resubmission incorporates the following adjustments:

- Changed 'plot_BRIER' for 'plot_I.BRIER' in all CV methods and fix their descriptions.
- Changed 'plot_c_index' for 'plot_C.Index' in all CV methods.
- Changed 'w_c.index' for 'w_C.Index' in all CV methods.
- Changed 'w_BRIER' for 'w_I.BRIER' in all CV methods.
- Changed 'brier.cox' for 'i.brier.cox' in 'evaluation_list_Coxmos' function.
- Changed 'c_index.cox' for 'c.index.cox' in 'evaluation_list_Coxmos' function.
- Changed labels and titles for CV plots in 'get_EVAL_PLOTS' function and their respective CV methods.
- Changed X label for 'comboplot.performance2.0' function.
- Changed ggpubr::ggarrange by patchwork R package.
- Fixed R2 values for all sPLS plots and description updated for methods that compute this value. 
- Fixed 'getCutoffAutoKM' function for COMP mode.
- Fixed KM plots. Now, if a variable was affected by transformIllegalChars(), now is reverted to plot it.
- Fixed 'plot_pseudobeta' for MB models.
- Fixed 'coxweightplot.fromVector.Coxmos'. Now, if a variable was affected by transformIllegalChars(), now is reverted to plot it.
- Fixed 'plot_pseudobeta.newObservation' to work with observation with more variables than the model selection.
- Fixed 'plot_LP.multipleObservations' to work with observation with more variables than the model selection.
- Fixed 'pseudobetas' plot for y axis.
- Updated 'plot_events' function to manage different number of decimals and angle of x.axis. Now return 4 plots and 4 data.frames instead 1. Percentage relative to observations, to group and to time added.
- Updated 'getAutoKM', 'getAutoKM.list', 'getTestKM' and 'getTestKM.list' functions to manage subtitles. The 'getAutoKM.list' and and 'getTestKM.list' functions now use the model name as title when title is NULL.
- Updated 'plot_pseudobeta' for work with a selection of cox variables instead of working only with full model.
- Updated 'plot_pseudobeta' for MB models with a new plot which mix all omics together.
- Updated 'plot_observation.pseudobeta' for better legend, arrange and text.
- Updated 'plot_LP.multipleObservations' to manage more parameters for plot design.
- Updated 'getEPV' to show error message when using a MB data instead of HD.
- Renamed multiple parameters for 'plot_LP.multipleObservations' and 'plot_pseudobeta_newObservation' functions.
- Miscellaneous changes in variables names.
- Miscellaneous changes in all descriptions.
- Miscellaneous changes in plot titles and axis.
- Miscellaneous changes in function names.
- README updated.

# CRAN Submission Comments (2025-02-06)
This resubmission incorporates the following adjustments:

- Function 'predict.Coxmos' fixed for iSB and SB methods ('mean' and 'sd', and when not all block are present in final model and for NULL data).
- Function 'isb.splsicox' now return scores for all blocks.
- Function 'isb.splsdrcox' now return scores for all blocks.
- Function 'isb.splsdrcox_dynamic' now return scores for all blocks.
- Function 'isb.splsdacox_dynamic' now return scores for all blocks.
- Functions 'cv.isb' now compute the consumed time correctly.

# CRAN Submission Comments (2025-01-07)
This resubmission incorporates the following adjustments:

- Functions 'save_ggplot' and 'save_ggplot_lst' updated to support multiple formats.
- Functions 'save_ggplot.svg' and 'save_ggplot_lst.svg' removed, as their functionality is now included in 'save_ggplot'.
- Function 'plot_time.list' updated to accept additional parameters, enabling greater customization.
- Function 'checkTestTimesVSTrainTimes' updated to handle NA cases.
- Function 'eval_Coxmos_models' updated to exclude models with NA values.
- Function 'evaluation_list_Coxmos' fixed to ensure compatibility with 'checkColnamesIllegalChars()'.
- Function 'transformIllegalChars' updated.
- Function 'plot_pseudobeta.list' updated to handle non-PLS models.
- Functions 'getAUC_vector' and 'getAUC_from_LP_2.0' fixed for 'smoothROCtime_C' and 'smoothROCtime_I'.
- Functions 'getTestKM.list' and 'getTestKM' fixed to handle X_test column names using 'checkColnamesIllegalChars()'.
- Function 'coxSW' has been improved:
  - 'AIC' metric is now considered when selecting variables to add or remove from the model.
  - Updated order of Cox model checks to handle 'Infinity' and 'NA' values more effectively.
  - Verbose messages and documentation have been updated.
- Function 'splsdacox_dynamic' updated to handle 'NA' and 'NULL' values in the final Cox model.
- Function 'removeNAorINFcoxmodel' fixed.
- Function 'getBestVector' fixed and messages updated.
- Function 'getAUC_from_LP_2.0' fixed.
- Function 'eval_Coxmos_models' fixed to handle cases where 'I. Brier' cannot be computed for the first model.
  - Issues with repeated times in 'I. Brier' Score results have been resolved.
- Function 'getCIndex_AUC_CoxModel_spls' fixed to handle non-converging models.
- Function 'get_Coxmos_models2.0' updated to manage the 'EVAL_METHOD' parameter in MB functions.
- Descriptions of 'Dynamic' functions updated to include all possible values for 'EVAL_METHOD'.
- Function 'boxplot.performance' updated to separate different types of evaluators.
- New functions 'checkX.colnames' and 'checkX.colnames.mb' added to all relevant methods.
- Parameter 'EVAL_METHOD' updated across all functions that include it.
- Function 'plot_evaluation' fixed to handle missing values at certain time points.
- Functions 'deleteZeroOrNearZeroVariance' and 'deleteZeroOrNearZeroVariance.mb' updated to return the percentage of unique values (overall and per block, respectively).
- Function 'deleteNearZeroCoefficientOfVariation.mb' fixed to return deleted variables.
- Updated 'MIN_COMP_TO_CHECK' value in all methods.
- Function 'plot_Coxmos.MB.PLS.model' fixed for MB models.
  - Biplot method now normalizes loadings to align with score values.
- Functions in the 'MB' family now generate an optimized design matrix when 'design = NULL'.
  - Design matrix usage has been fixed.
- New function 'getDesign.MB' added to compute the design matrix automatically.
- Function 'deleteNearZeroCoefficientOfVariation' updated to remove NA values, if present.
- Function 'splsdrcox' updated to align with [doi:10.1093/bioinformatics/btu660], except for the penalty.
- Function 'coxEN' updated to handle illegal characters in column names and to handle problems with one variable models.
- Function 'predict.Coxmos' fixed and verified for all methods.
- Minor changes made to 'sb.splsicox'.
- All algorithms now transform illegal characters in column names to manageable alternatives.
- Function 'getAUC_RUN_AND_COMP' fixed.
- Functions 'getAUC_RUN_AND_COMP' and 'getAUC_RUN_AND_COMP_sPLS' updated to work with specific method groups.
- Function 'print.Coxmos' updated to support iSB models.
- Miscellaneous changes made to plotting functions related to 'deleteIllegalChars'.
- Function renaming:
  - Original 'drcox' functions now use the '_penalty' suffix.
  - '_dynamic' functions renamed to reflect only the algorithm name.
- Function 'coxEN' updated to select the maximum number of variables by default.
- Functions 'getBestVector' and 'getBestVector2' fixed.
- Functions 'getCompKM' and 'getLPVarKM' fixed for MB algorithms.
- New algorithms added:
  - 'isb.splsicox_dynamic', 'cv.isb.splsdacox_dynamic'
  - 'isb.splsdrcox', 'cv.isb.splsdrcox'
  - 'sb.splsdrcox_dynamic', 'isb.splsdrcox_dynamic', 'cv.sb.splsdrcox_dynamic', 'cv.isb.splsdrcox_dynamic'
  - 'sb.splsdacox_dynamic', 'isb.splsdacox_dynamic', 'cv.sb.splsdacox_dynamic', 'cv.isb.splsdacox_dynamic'
  - New annotations added in common functions.
- Function 'get_Coxmos_models2.0' updated to support all new methods.
- Functions 'get_COX_evaluation_AIC_CINDEX' and 'get_COX_evaluation_BRIER_sPLS' updated to manage all new methods.
- Functions 'get_COX_evaluation_BRIER' and 'get_COX_evaluation_BRIER_sPLS' fixed to work correctly with BRIER metric.
- Function 'getTimesVector' updated to select an optimal time interval for evaluation.
- New parameter added to Kaplan-Meier functions. Now, minimum proportion of observations can be included for optimal cutoff in numerical variables.
- Function 'get_COX_evaluation_AUC' fixed.
- New 'eval_Coxmos_model_per_variable.list' function.
- New 'is.binaryMatrix' function.
- Functions 'getCutoffAutoKM' and 'getCutoffAutoKM.list' now deal qualitative KM.
- Functions 'getCutoffAutoKM' and 'getCutoffAutoKM.list' now uses a CV with 5 folds to obtain the cutoff.
- Functions 'getTestKM' and 'getTestKM.list' now deal qualitative KM.
- Function 'getLPKM' and 'getCompKM' fixed when obtaining LP for specific variabes.
- Function 'Csurv' in 'cenROC' updated when 'bw.SJ(M, method = "dpi")' cannot work.
- Function 'eval_Coxmos_model_per_variable' and 'eval_Coxmos_model_per_variable.list' updated to work by default with TRAIN times when 'times' vector is NULL.
- Function 'timesAsumption_AUC_Eval' fixed to not manage eval. times greater than the maximum of Y.
- Function 'splitData_Iterations_Folds', 'splitData_Iterations_Folds_indexes', 'splitData_Iterations_Folds.mb' and 'getLogRank_NumVariables' fixed when all observations have the same event.

# CRAN Submission Comments (2024-03-20 & 2024-03-22)
This resubmission incorporates the following adjustments:

- Vignettes have been changed in order to use the correct parameters and compute all the results.

# CRAN Submission Comments (2024-02-29)

## R CMD check results
- 0 errors.
- 0 warnings.
- 1 note: Received a note regarding "checking for detritus in the temp directory." 

I think this is due to the development environment on Windows 10, and it should not affect the package's functionality or integrity.

## CRAN Submission Comments  (2024-03-05)
This resubmission incorporates the following adjustments:

- **R/Coxmos_plot_functions.R - getwd() Issue**: The folder argument is now mandatory. Examples use `tempdir()`.

- **Authors@R**: All authors are listed in the DESCRIPTION file, including thesis advisers and contributors following their license agreements.

- **Seed Control in Functions**: Introduced a parameter in all functions to specify the seed.

- **References in DESCRIPTION**: Key references are now cited in the DESCRIPTION file. Additional references are noted to be available within each function's documentation.

- **Examples using dontrun{}**: Converted all applicable `dontrun{}` instances to `donttest{}`.

- **Dependency Checks**: Revised to ensure dependencies are appropriately managed without the explicit need for `requireNamespace()` checks in code.

- **Console Messages**: Transitioned non-object `print()` and `cat()` calls to `message()`, and where appropriate, to `warning()`.

- **Acronyms**: Defined all acronyms in the DESCRIPTION file and corrected "COX" to "Cox".
