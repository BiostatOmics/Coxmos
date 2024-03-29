% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_plot_functions.R
\name{plot_pseudobeta_newObservation}
\alias{plot_pseudobeta_newObservation}
\title{plot_pseudobeta.newObservation}
\usage{
plot_pseudobeta_newObservation(
  model,
  new_observation,
  error.bar = TRUE,
  onlySig = TRUE,
  alpha = 0.05,
  zero.rm = TRUE,
  top = NULL,
  auto.limits = TRUE,
  show.betas = FALSE
)
}
\arguments{
\item{model}{Coxmos model.}

\item{new_observation}{Numeric matrix or data.frame. New explanatory variables (raw data) for one
observation. Qualitative variables must be transform into binary variables.}

\item{error.bar}{Logical. Show error bar (default: TRUE).}

\item{onlySig}{Logical. Compute pseudobetas using only significant components (default: TRUE).}

\item{alpha}{Numeric. Numerical values are regarded as significant if they fall below the
threshold (default: 0.05).}

\item{zero.rm}{Logical. Remove variables with a pseudobeta equal to 0 (default: TRUE).}

\item{top}{Numeric. Show "top" first variables with the higher pseudobetas in absolute value. If
top = NULL, all variables are shown (default: NULL).}

\item{auto.limits}{Logical. If "auto.limits" = TRUE, limits are detected automatically (default: TRUE).}

\item{show.betas}{Logical. Show original betas (default: FALSE).}
}
\value{
A list of four elements:
\code{plot}: Linear prediction per variable.
\code{lp.var}: Value of each linear prediction per variable.
\code{norm_observation}: Observation normalized using the model information.
\code{observation}: Observation used.
}
\description{
Generates a visual representation comparing the pseudobeta values derived from the Coxmos model
with the values of a new observation. This function provides insights into how the new observation
aligns with the established model, offering a graphical comparison of the pseudobeta directions.
}
\details{
The function \code{plot_pseudobeta.newObservation} is designed to visually compare the pseudobeta values
from the Coxmos model with those of a new observation. The generated plot is based on the ggplot2
framework and offers a comprehensive view of the relationship between the model's pseudobeta values
and the new observation's values.

The function first checks the validity of the provided model and ensures that it belongs to the
appropriate class. Depending on the type of the model (either PLS or MB Coxmos methods).

For the actual plotting, the function computes the linear predictor values for the new observation
and juxtaposes them with the pseudobeta values from the model. If the \code{show.betas} parameter is
set to TRUE, the original beta values are also displayed on the plot. Error bars can be included
to represent the variability in the pseudobeta values, providing a more comprehensive view of the
data's distribution.

The resulting plot serves as a valuable tool for researchers and statisticians to visually assess
the alignment of a new observation with an established Coxmos model, facilitating better
interpretation and understanding of the data in the context of the model.
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
plot_pseudobeta_newObservation(model = splsicox.model, new_observation = X_test[1,,drop=FALSE])
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
