# CRAN Submission Comments (2024-11-11)
This resubmission incorporates the following adjustments:

- Functions 'save_ggplot' and 'save_ggplot_lst' updated to work with multiples formats.
- Functions 'save_ggplot.svg' and 'save_ggplot_lst.svg' deleted due to they are included in 'save_ggplot'.
- Function 'plot_time.list' have been updated to work with extra parameters to allow more customization.
- Function 'checkTestTimesVSTrainTimes' updated to consider NA cases.
- Function 'eval_Coxmos_models' updated to remove NA models cases.
- Function 'evaluation_list_Coxmos' fixed for all cases with  'checkColnamesIllegalChars()'.
- Function 'plot_pseudobeta.list' updated to work no PLS model case.
- Functions 'getAUC_vector' and 'getAUC_from_LP_2.0' fixed for smoothROCtime_C and smoothROCtime_I.
- Function 'getTestKM.list' and 'getTestKM' fixed for all cases: X_test colnames with 'checkColnamesIllegalChars()'.
- Function 'coxSW' has been fixed in several aspects:
  - Now, AIC metric is also considered when selecting a new variable to enter or leave the old model.
  - Updated order of checking Cox models to avoid problems removing Infinity or NA values.
  - Verbose messages and documentation have been updated.
- Function 'splsdacox_dynamic' has been update for manage problems with NA or NUL values in final Cox model.
- Function 'removeNAorINFcoxmodel' fixed.
- Function 'getBestVector' fixed.
- Function 'getAUC_from_LP_2.0' fixed.
- Function 'eval_Coxmos_models' fixed when I.Brier cannot be computed for the first model.
  - Fixed problems with repeated times in I.Brier Score results.
- Function 'getCIndex_AUC_CoxModel_spls' fixed for managing models that do not converge.
- Function 'get_Coxmos_models2.0' updated for MB functions to manage 'EVAL_METHOD' parameter.
- Description of 'Dynamic' Functions updated to include all values for 'EVAL_METHOD'.
- Function 'boxplot.performance' updated to split different type of evaluators.
- New function 'checkX.colnames' and 'checkX.colnames.mb' added in all methods.
- Parameter 'EVAL_METHOD' updated in all functions that contains this parameter.
- Function 'plot_evaluation' fixed to manage missin values for some times.
- Function 'deleteZeroOrNearZeroVariance' updated to return percentage of unique values.
- Function 'deleteZeroOrNearZeroVariance.mb' updated to return percentage of unique values per block.
- Function 'deleteNearZeroCoefficientOfVariation.mb' fixed to return the deleted variables.
- Updated MIN_COMP_TO_CHECK value in all methods.
- Function 'plot_Coxmos.MB.PLS.model' fixed for MB models.
- The function 'plot_Coxmos.MB.PLS.model' has been fixed in the Biplot method. It now normalizes loadings to align with score values.
- Functions 'MB' now generates an optimized DESIGN matrix if design is NULL.
  - Fixed the use of design matrix.
- New function 'getDesign.MB' for computing design automatically.
- Function 'deleteNearZeroCoefficientOfVariation' updated to remove NAs if any.
- Function 'splsdrcox' updated to work as doi:10.1093/bioinformatics/btu660 except for the penalty.
- Function 'coxEN' updated to manage illegal chars in column names.
- Function 'predict.Coxmos' fixed and checked for all methods.
- Minor changes in 'sb.splsicox'.
- Now, all algorithms transform Illegal Chars in columns to something manageable.
- Functions 'getAUC_RUN_AND_COMP' and 'getAUC_RUN_AND_COMP_sPLS' updated to work with specific group of methods.
- Function 'print.Coxmos' updated to manage iSB.
- Miscellaneous changes in Plot functiones related to 'deleteIllegalChars'.
- Function Rename: 'drcox' original functions were renamed to use the suffix '_penalty'.
- Function Rename: '_dynamic' functions were renamed to use only the algorithm name.
- Function 'coxEN' updated to select, as default, the maximum number of variables.
- Function 'getBestVector' and 'getBestVector2' fixed.
- Functions 'getCompKM' and 'getLPVarKM' fixed for MB algorithms.
- New algorithms:
  - New 'isb.splsicox_dynamic' and 'cv.isb.splsdacox_dynamic'.
  - New 'isb.splsdrcox' and 'cv.isb.splsdrcox' algorithms.
  - New 'sb.splsdrcox_dynamic', 'isb.splsdrcox_dynamic', 'cv.sb.splsdrcox_dynamic' and 'cv.isb.splsdrcox_dynamic' algorithms.
  - New 'sb.splsdacox_dynamic', 'isb.splsdacox_dynamic', 'cv.sb.splsdacox_dynamic' and 'cv.isb.splsdacox_dynamic' algorithms.
  - New annotations in common functions.
- Function 'get_Coxmos_models2.0' updated to manage all new methods.
- Function 'get_COX_evaluation_AIC_CINDEX' updated to manage all new methods.
- Function 'get_COX_evaluation_BRIER_sPLS' updated to manage all new methods.

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
