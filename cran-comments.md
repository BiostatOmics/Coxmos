# CRAN Submission Comments (2024-11-26)
This resubmission incorporates the following adjustments:

- Functions 'save_ggplot' and 'save_ggplot_lst' updated to support multiple formats.
- Functions 'save_ggplot.svg' and 'save_ggplot_lst.svg' removed, as their functionality is now included in 'save_ggplot'.
- Function 'plot_time.list' updated to accept additional parameters, enabling greater customization.
- Function 'checkTestTimesVSTrainTimes' updated to handle NA cases.
- Function 'eval_Coxmos_models' updated to exclude models with NA values.
- Function 'evaluation_list_Coxmos' fixed to ensure compatibility with 'checkColnamesIllegalChars()'.
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
- Function 'eval_Coxmos_models' fixed to handle cases where 'I.Brier' cannot be computed for the first model.
  - Issues with repeated times in 'I.Brier' Score results have been resolved.
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
