# CRAN Submission Comments (2024-02-29)

## R CMD check results
- 0 errors.
- 0 warnings.
- 1 note: Received a note regarding "checking for detritus in the temp directory." 

I think this is due to the development environment on Windows 10, and it should not affect the package's functionality or integrity.

## Resubmission (2024-03-05)
This resubmission incorporates the following adjustments:

- **R/Coxmos_plot_functions.R - getwd() Issue**: The folder argument is now mandatory. Examples use `tempdir()`.

- **Authors@R**: All authors are listed in the DESCRIPTION file, including thesis advisers and contributors following their license agreements.

- **Seed Control in Functions**: Introduced a parameter in all functions to specify the seed.

- **References in DESCRIPTION**: Key references are now cited in the DESCRIPTION file. Additional references are noted to be available within each function's documentation.

- **Examples using dontrun{}**: Converted all applicable `dontrun{}` instances to `donttest{}`.

- **Dependency Checks**: Revised to ensure dependencies are appropriately managed without the explicit need for `requireNamespace()` checks in code.

- **Console Messages**: Transitioned non-object `print()` and `cat()` calls to `message()`, and where appropriate, to `warning()`.

- **Acronyms**: Defined all acronyms in the DESCRIPTION file and corrected "COX" to "Cox".
