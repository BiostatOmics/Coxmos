% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_splsdrcox.R, R/Coxmos_splsdrcox_median.R
\name{splsdrcox_penalty}
\alias{splsdrcox_penalty}
\title{sPLS-DRCOX}
\usage{
splsdrcox_penalty(
  X,
  Y,
  n.comp = 4,
  penalty = 0.5,
  x.center = TRUE,
  x.scale = FALSE,
  remove_near_zero_variance = TRUE,
  remove_zero_variance = FALSE,
  toKeep.zv = NULL,
  remove_non_significant = FALSE,
  alpha = 0.05,
  MIN_EPV = 5,
  returnData = TRUE,
  verbose = FALSE
)

splsdrcox_penalty(
  X,
  Y,
  n.comp = 4,
  penalty = 0.5,
  x.center = TRUE,
  x.scale = FALSE,
  remove_near_zero_variance = TRUE,
  remove_zero_variance = FALSE,
  toKeep.zv = NULL,
  remove_non_significant = FALSE,
  alpha = 0.05,
  MIN_EPV = 5,
  returnData = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{X}{Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be
transform into binary variables.}

\item{Y}{Numeric matrix or data.frame. Response variables. Object must have two columns named as
"time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event
observations.}

\item{n.comp}{Numeric. Number of latent components to compute for the (s)PLS model (default: 4).}

\item{penalty}{Numeric. Penalty for sPLS-DRCOX. If penalty = 0 no penalty is applied, when
penalty = 1 maximum penalty (no variables are selected) based on 'plsRcox' penalty. Equal or greater
than 1 cannot be selected (default: 0.5).}

\item{x.center}{Logical. If x.center = TRUE, X matrix is centered to zero means (default: TRUE).}

\item{x.scale}{Logical. If x.scale = TRUE, X matrix is scaled to unit variances (default: FALSE).}

\item{remove_near_zero_variance}{Logical. If remove_near_zero_variance = TRUE, near zero variance
variables will be removed (default: TRUE).}

\item{remove_zero_variance}{Logical. If remove_zero_variance = TRUE, zero variance variables will
be removed (default: TRUE).}

\item{toKeep.zv}{Character vector. Name of variables in X to not be deleted by (near) zero
variance filtering (default: NULL).}

\item{remove_non_significant}{Logical. If remove_non_significant = TRUE, non-significant
variables/components in final cox model will be removed until all variables are significant by
forward selection (default: FALSE).}

\item{alpha}{Numeric. Numerical values are regarded as significant if they fall below the
threshold (default: 0.05).}

\item{MIN_EPV}{Numeric. Minimum number of Events Per Variable (EPV) you want reach for the final
cox model. Used to restrict the number of variables/components can be computed in final cox models.
If the minimum is not meet, the model cannot be computed (default: 5).}

\item{returnData}{Logical. Return original and normalized X and Y matrices (default: TRUE).}

\item{verbose}{Logical. If verbose = TRUE, extra messages could be displayed (default: FALSE).}
}
\value{
Instance of class "Coxmos" and model "sPLS-DRCOX". The class contains the following
elements:
\code{X}: List of normalized X data information.
\itemize{
\item \code{(data)}: normalized X matrix
\item \code{(weightings)}: sPLS weights
\item \code{(weightings_norm)}: sPLS normalize weights
\item \code{(W.star)}: sPLS W* vector
\item \code{(loadings)}: sPLS loadings
\item \code{(scores)}: sPLS scores/variates
\item \code{(E)}: error matrices
\item \code{(x.mean)}: mean values for X matrix
\item \code{(x.sd)}: standard deviation for X matrix
}
\code{Y}: List of normalized Y data information.
\itemize{
\item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
\item \code{(dr.mean)}: mean values for deviance residuals Y matrix
\item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
\item \code{(data)}: normalized X matrix
\item \code{(weightings)}: sPLS weights
\item \code{(loadings)}: sPLS loadings
\item \code{(scores)}: sPLS scores/variates
\item \code{(ratio)}: r value for the sPLS model (used to perform predictions)
\item \code{(y.mean)}: mean values for Y matrix
\item \code{(y.sd)}: standard deviation for Y matrix'
}

\code{survival_model}: List of survival model information.
\itemize{
\item \code{fit}: coxph object.
\item \code{AIC}: AIC of cox model.
\item \code{BIC}: BIC of cox model.
\item \code{lp}: linear predictors for train data.
\item \code{coef}: Coefficients for cox model.
\item \code{YChapeau}: Y Chapeau residuals.
\item \code{Yresidus}: Y residuals.
}

\code{penalty}: Penalty value selected.

\code{n.comp}: Number of components selected.

\code{var_by_component}: Variables selected in each PLS component.

\code{call}: call function

\code{X_input}: X input matrix

\code{Y_input}: Y input matrix

\code{B.hat}: sPLS beta matrix

\code{R2}: sPLS acumulate R2

\code{SCR}: sPLS SCR

\code{SCT}: sPLS SCT

\code{alpha}: alpha value selected

\code{nsv}: Variables removed by cox alpha cutoff.

\code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.

\code{nz_coeffvar}: Variables removed by coefficient variation near zero.

\code{class}: Model class.

\code{time}: time consumed for running the cox analysis.

Instance of class "Coxmos" and model "sPLS-DRCOX". The class contains the following
elements:
\code{X}: List of normalized X data information.
\itemize{
\item \code{(data)}: normalized X matrix
\item \code{(weightings)}: sPLS weights
\item \code{(weightings_norm)}: sPLS normalize weights
\item \code{(W.star)}: sPLS W* vector
\item \code{(loadings)}: sPLS loadings
\item \code{(scores)}: sPLS scores/variates
\item \code{(E)}: error matrices
\item \code{(x.mean)}: mean values for X matrix
\item \code{(x.sd)}: standard deviation for X matrix
}
\code{Y}: List of normalized Y data information.
\itemize{
\item \code{(deviance_residuals)}: deviance residual vector used as Y matrix in the sPLS.
\item \code{(dr.mean)}: mean values for deviance residuals Y matrix
\item \code{(dr.sd)}: standard deviation for deviance residuals Y matrix'
\item \code{(data)}: normalized X matrix
\item \code{(weightings)}: sPLS weights
\item \code{(loadings)}: sPLS loadings
\item \code{(scores)}: sPLS scores/variates
\item \code{(ratio)}: r value for the sPLS model (used to perform predictions)
\item \code{(y.mean)}: mean values for Y matrix
\item \code{(y.sd)}: standard deviation for Y matrix'
}

\code{survival_model}: List of survival model information.
\itemize{
\item \code{fit}: coxph object.
\item \code{AIC}: AIC of cox model.
\item \code{BIC}: BIC of cox model.
\item \code{lp}: linear predictors for train data.
\item \code{coef}: Coefficients for cox model.
\item \code{YChapeau}: Y Chapeau residuals.
\item \code{Yresidus}: Y residuals.
}

\code{penalty}: Penalty value selected.

\code{n.comp}: Number of components selected.

\code{var_by_component}: Variables selected in each PLS component.

\code{call}: call function

\code{X_input}: X input matrix

\code{Y_input}: Y input matrix

\code{B.hat}: sPLS beta matrix

\code{R2}: sPLS acumulate R2

\code{SCR}: sPLS SCR

\code{SCT}: sPLS SCT

\code{alpha}: alpha value selected

\code{nsv}: Variables removed by cox alpha cutoff.

\code{nzv}: Variables removed by remove_near_zero_variance or remove_zero_variance.

\code{nz_coeffvar}: Variables removed by coefficient variation near zero.

\code{class}: Model class.

\code{time}: time consumed for running the cox analysis.
}
\description{
This function performs a sparse partial least squares deviance residual Cox (sPLS-DRCOX)
(based on plsRcox R package). The function returns a Coxmos model with the attribute model as
"sPLS-DRCOX".

This function performs a sparse partial least squares deviance residual Cox (sPLS-DRCOX)
(based on plsRcox R package). The function returns a Coxmos model with the attribute model as
"sPLS-DRCOX".
}
\details{
The \code{sPLS-DRCOX} function implements the sparse partial least squares deviance residual Cox
(sPLS-DRCOX) model, a specialized approach tailored for survival analysis. This method integrates
the strengths of the sparse partial least squares (sPLS) technique with the Cox proportional hazards
model, leveraging deviance residuals as a bridge.

The function's core lies in its ability to handle high-dimensional data, often encountered in
genomics or other omics studies. By incorporating the \code{penalty} parameter, which governs the sparsity
level, the function offers a fine-grained control over variable selection. This ensures that only
the most informative predictors contribute to the model, enhancing interpretability and reducing
overfitting.

Data preprocessing is seamlessly integrated, with options to center and scale the predictors, and
to remove variables exhibiting near-zero or zero variance. The function also provides a mechanism
to retain specific variables, regardless of their variance, ensuring that domain-specific knowledge
can be incorporated.

The output is comprehensive, detailing both the sPLS and Cox model components. It provides insights
into the selected variables, their contributions across latent components, and the overall fit of
the survival model. This rich output aids in understanding the underlying relationships between
predictors and survival outcomes.

The \code{sPLS-DRCOX} function is grounded in established methodologies and is a valuable tool for
researchers aiming to unravel complex survival associations in high-dimensional datasets.

The \code{sPLS-DRCOX} function implements the sparse partial least squares deviance residual Cox
(sPLS-DRCOX) model, a specialized approach tailored for survival analysis. This method integrates
the strengths of the sparse partial least squares (sPLS) technique with the Cox proportional hazards
model, leveraging deviance residuals as a bridge.

The function's core lies in its ability to handle high-dimensional data, often encountered in
genomics or other omics studies. By incorporating the \code{penalty} parameter, which governs the sparsity
level, the function offers a fine-grained control over variable selection. This ensures that only
the most informative predictors contribute to the model, enhancing interpretability and reducing
overfitting.

Data preprocessing is seamlessly integrated, with options to center and scale the predictors, and
to remove variables exhibiting near-zero or zero variance. The function also provides a mechanism
to retain specific variables, regardless of their variance, ensuring that domain-specific knowledge
can be incorporated.

The output is comprehensive, detailing both the sPLS and Cox model components. It provides insights
into the selected variables, their contributions across latent components, and the overall fit of
the survival model. This rich output aids in understanding the underlying relationships between
predictors and survival outcomes.

The \code{sPLS-DRCOX} function is grounded in established methodologies and is a valuable tool for
researchers aiming to unravel complex survival associations in high-dimensional datasets.
}
\examples{
data("X_proteomic")
data("Y_proteomic")
X <- X_proteomic[,1:50]
Y <- Y_proteomic
splsdrcox_penalty(X, Y, n.comp = 3, penalty = 0.25, x.center = TRUE, x.scale = TRUE)
data("X_proteomic")
data("Y_proteomic")
X <- X_proteomic[,1:50]
Y <- Y_proteomic
splsdrcox_penalty(X, Y, n.comp = 3, penalty = 0.25, x.center = TRUE, x.scale = TRUE)
}
\references{
\insertRef{Bastien_2008}{Coxmos}
\insertRef{Bastien_2015}{Coxmos}

\insertRef{Bastien_2008}{Coxmos}
\insertRef{Bastien_2015}{Coxmos}
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
