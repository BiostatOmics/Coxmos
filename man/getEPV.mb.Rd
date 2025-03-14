% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coxmos_multiblock_functions.R
\name{getEPV.mb}
\alias{getEPV.mb}
\title{getEPV.mb}
\usage{
getEPV.mb(X, Y)
}
\arguments{
\item{X}{List of numeric matrices or data.frames. Explanatory variables. Qualitative variables must be transform
into binary variables.}

\item{Y}{Numeric matrix or data.frame. Response variables. Object must have two columns named as
"time" and "event". For event column, accepted values are: 0/1 or FALSE/TRUE for censored and event
observations.}
}
\value{
Return the EPV value for a specific X (explanatory variables) and Y (time and censored variables) data.
}
\description{
Provides a quantitative assessment of the dataset by computing the Events per Variable
(EPV) metric for multi-block data, which gauges the proportionality between observed events and the
number of explanatory variables.
}
\details{
In the realm of survival analysis, the balance between observed events and explanatory
variables is paramount. The \code{getEPV} function serves as a tool for researchers to ascertain this
balance, which can be pivotal in determining the robustness and interpretability of subsequent
statistical models. By evaluating the ratio of events in the \code{Y} matrix to the variables in the \code{X}
matrix, the function yields the EPV metric. It is of utmost importance that the \code{Y} matrix encompasses
two distinct columns, namely "time" and "event". The latter, "event", should strictly encapsulate
binary values, delineating censored (either 0 or FALSE) and event (either 1 or TRUE) observations.
To ensure the integrity of the data and the precision of the computation, the function is equipped
with an error mechanism that activates if the "event" column remains undetected.
}
\examples{
data("X_multiomic")
data("Y_multiomic")
X <- X_multiomic
Y <- Y_multiomic
getEPV.mb(X,Y)
}
\author{
Pedro Salguero Garcia. Maintainer: pedsalga@upv.edu.es
}
