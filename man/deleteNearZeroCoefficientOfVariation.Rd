% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_common_functions.R
\name{deleteNearZeroCoefficientOfVariation}
\alias{deleteNearZeroCoefficientOfVariation}
\title{deleteNearZeroCoefficientOfVariation}
\usage{
deleteNearZeroCoefficientOfVariation(X, LIMIT = 0.1)
}
\arguments{
\item{X}{Numeric matrix or data.frame. Explanatory variables. Qualitative variables must be transform
into binary variables.}

\item{LIMIT}{Numeric. Cutoff for minimum variation. If coefficient is lesser than the limit, the
variables are removed because not vary enough (default: 0.1).}
}
\value{
Return a list of two objects:
\code{X}: The new data.frame X filtered.
\code{variablesDeleted}: The variables that have been removed by the filter.
\code{coeff_variation}: The coefficient variables per each variable tested.
}
\description{
Filters out variables from a dataset that exhibit a coefficient of variation below a
specified threshold, ensuring the retention of variables with meaningful variability.
}
\details{
The \code{deleteNearZeroCoefficientOfVariation} function is a pivotal tool in data preprocessing,
especially when dealing with high-dimensional datasets. The coefficient of variation (CoV) is a
normalized measure of data dispersion, calculated as the ratio of the standard deviation to the mean.
In many scientific investigations, variables with a low CoV might be considered as offering
limited discriminative information, potentially leading to noise in subsequent statistical analyses.
By setting a threshold through the \code{LIMIT} parameter, this function provides a systematic approach
to identify and exclude variables that do not meet the desired variability criteria. The underlying
rationale is that variables with a CoV below the set threshold might not contribute significantly
to the variability of the dataset and could be redundant or even detrimental for certain analyses.
The function returns a modified dataset, a list of deleted variables, and the computed coefficients
of variation for each variable. This comprehensive output ensures that researchers are well-informed
about the preprocessing steps and can make subsequent analytical decisions with confidence.
}
\examples{
data("X_proteomic")
X <- X_proteomic
filter <- deleteNearZeroCoefficientOfVariation(X, LIMIT = 0.1)
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
