% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{plot_pseudobeta.list}
\alias{plot_pseudobeta.list}
\title{plot_pseudobeta.list}
\usage{
plot_pseudobeta.list(
  lst_models,
  error.bar = TRUE,
  onlySig = FALSE,
  alpha = 0.05,
  zero.rm = TRUE,
  top = NULL,
  auto.limits = TRUE,
  selected_variables = NULL,
  show_percentage = TRUE,
  size_percentage = 3,
  title = NULL,
  title_size_text = 15,
  subtitle = NULL,
  subtitle_size_text = 12,
  legend.position = "right",
  legend_title = "Method",
  legend_size_text = 12,
  x_axis_size_text = 10,
  y_axis_size_text = 10,
  label_x_axis_size = 10,
  label_y_axis_size = 10,
  verbose = FALSE
)
}
\arguments{
\item{lst_models}{List of Coxmos models.}

\item{error.bar}{Logical. Show error bar (default: TRUE).}

\item{onlySig}{Logical. Compute pseudobetas using only significant components (default: FALSE).}

\item{alpha}{Numeric. Numerical values are regarded as significant if they fall below the
threshold (default: 0.05).}

\item{zero.rm}{Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).}

\item{top}{Numeric. Show "top" first variables with the higher pseudobetas in absolute value.
If top = NULL, all variables are shown (default: NULL).}

\item{auto.limits}{Logical. If "auto.limits" = TRUE, limits are detected automatically
(default: TRUE).}

\item{selected_variables}{Character. Name of survival model variables to performed a custom selection (default: NULL).}

\item{show_percentage}{Logical. If show_percentage = TRUE, it shows the contribution percentage
for each variable to the full model (default: TRUE).}

\item{size_percentage}{Numeric. Size of percentage text (default: 3).}

\item{title}{Character. Plot title (default: NULL).}

\item{title_size_text}{Numeric. Text size for title (default: 15).}

\item{subtitle}{Character. Plot subtitle (default: NULL).}

\item{subtitle_size_text}{Numeric. Text size for subtitle (default: 12).}

\item{legend.position}{Character. Legend position. Must be one of the following: "top", "bottom", "right" or "left (default: "right").}

\item{legend_title}{Character. Legend title (default: "Method").}

\item{legend_size_text}{Numeric. Text size for legend title (default: 12).}

\item{x_axis_size_text}{Numeric. Text size for x axis (default: 10).}

\item{y_axis_size_text}{Numeric. Text size for y axis (default: 10).}

\item{label_x_axis_size}{Numeric. Text size for x label axis (default: 10).}

\item{label_y_axis_size}{Numeric. Text size for y label axis (default: 10).}

\item{verbose}{Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).}
}
\value{
A list containing the following elements per model:
\code{plot}: Depending on the model type, this can either be a single ggplot object visualizing the pseudo-beta coefficients for the original variables in a single block PLS-Cox model, or a list of ggplot objects for each block in a multiblock PLS-Cox model. Each plot provides a comprehensive visualization of the pseudo-beta coefficients, potentially including error bars, significance filtering, and variable contribution percentages.
\code{beta}: A matrix or list of matrices (for multiblock models) containing the computed pseudo-beta coefficients for the original variables. These coefficients represent the influence of each original variable on the survival prediction.
\code{sd.min}: A matrix or list of matrices (for multiblock models) representing the lower bounds of the error bars for the pseudo-beta coefficients.
\code{sd.max}: A matrix or list of matrices (for multiblock models) representing the upper bounds of the error bars for the pseudo-beta coefficients.
}
\description{
Run the function "plot_pseudobeta" for a list of models. More information in
"?plot_pseudobeta".
}
\examples{
data("X_proteomic")
data("Y_proteomic")
X <- X_proteomic[,1:50]
Y <- Y_proteomic
splsicox.model <- splsicox(X, Y, n.comp = 2, penalty = 0.5, x.center = TRUE, x.scale = TRUE)
splsdrcox.model <- splsdrcox_penalty(X, Y, n.comp = 2, penalty = 0.5,
x.center = TRUE, x.scale = TRUE)
lst_models = list("sPLSICOX" = splsicox.model, "sPLSDRCOX" = splsdrcox.model)
plot_pseudobeta.list(lst_models = lst_models)
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
