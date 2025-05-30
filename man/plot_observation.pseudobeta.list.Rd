% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{plot_observation.pseudobeta.list}
\alias{plot_observation.pseudobeta.list}
\title{plot_observation.pseudobeta.list}
\usage{
plot_observation.pseudobeta.list(
  lst_models,
  observation,
  error.bar = TRUE,
  onlySig = TRUE,
  alpha = 0.05,
  zero.rm = TRUE,
  txt.x.angle = 0,
  top = NULL,
  auto.limits = TRUE,
  show.betas = FALSE,
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

\item{observation}{Numeric matrix or data.frame. New explanatory variables (raw data) for one
observation. Qualitative variables must be transform into binary variables.}

\item{error.bar}{Logical. Show error bar (default: TRUE).}

\item{onlySig}{Logical. Compute pseudobetas using only significant components (default: TRUE).}

\item{alpha}{Numeric. Numerical values are regarded as significant if they fall below the
threshold (default: 0.05).}

\item{zero.rm}{Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).}

\item{txt.x.angle}{Numeric. Angle of X text (default: 0).}

\item{top}{Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If
top = NULL, all variables are shown (default: NULL).}

\item{auto.limits}{Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).}

\item{show.betas}{Logical. Show original betas (default: FALSE).}

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
A list of lst_models length with a list of four elements per each model:
\code{plot}: Linear prediction per variable.
\code{lp.var}: Value of each linear prediction per variable.
\code{norm_observation}: Observation normalized using the model information.
\code{observation}: Observation used.
}
\description{
Run the function "plot_observation.pseudobeta" for a list of models. More information
in "?plot_observation.pseudobeta".
}
\examples{
data("X_proteomic")
data("Y_proteomic")
set.seed(123)
index_train <- caret::createDataPartition(Y_proteomic$event, p = .5, list = FALSE, times = 1)
X_train <- X_proteomic[index_train,1:50]
Y_train <- Y_proteomic[index_train,]
X_test <- X_proteomic[-index_train,1:50]
Y_test <- Y_proteomic[-index_train,]
splsicox.model <- splsicox(X_train, Y_train, n.comp = 2, penalty = 0.5, x.center = TRUE,
x.scale = TRUE)
splsdrcox.model <- splsdrcox_penalty(X_train, Y_train, n.comp = 2, penalty = 0.5, x.center = TRUE,
x.scale = TRUE)
lst_models = list("sPLSICOX" = splsicox.model, "sPLSDRCOX" = splsdrcox.model)
plot_observation.pseudobeta.list(lst_models, observation = X_test[1,,drop=FALSE])
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
