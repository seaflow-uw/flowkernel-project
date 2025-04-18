% Generated by roxygen2: do not edit by hand
% Please edit documentation in _main.Rmd.
\name{compute_log_like}
\alias{compute_log_like}
\title{Compute the Average Log-Likelihood for Predicted Parameters}
\usage{
compute_log_like(y, pred)
}
\arguments{
\item{y}{A list of length T, where each element is a matrix of observations at a given time point. Each matrix should have dimensions corresponding to the number of observations at that time point by the data dimension.}

\item{pred}{A list containing the predicted parameters from the model. It should include:
\describe{
\item{mu}{An array of predicted means with dimensions \link{T, K, d}, where T is the number of time points, K is the number of clusters, and d is the data dimension.}
\item{Sigma}{An array of predicted covariance matrices with dimensions \link{T, K, d, d}.}
\item{pi}{A matrix of predicted mixing proportions with dimensions \link{T, K}.}
}}
}
\value{
A numeric value representing the average log-likelihood over all time points.
}
\description{
This function computes the average log-likelihood of the observed data over all time points using the predicted parameters.
}