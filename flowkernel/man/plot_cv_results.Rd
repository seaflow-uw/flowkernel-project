% Generated by roxygen2: do not edit by hand
% Please edit documentation in _main.Rmd.
\name{plot_cv_results}
\alias{plot_cv_results}
\title{Plot Cross-Validation Results}
\usage{
plot_cv_results(results_list)
}
\arguments{
\item{results_list}{A list returned by \code{cv_grid_search}. This list should include:
\describe{
\item{results}{A data frame of cross-validation scores, where the row names correspond to evaluated \code{hmu} values and the column names to evaluated \code{hpi} values.}
\item{hmu_vals}{A numeric vector of the \code{hmu} bandwidth values that were evaluated.}
\item{hpi_vals}{A numeric vector of the \code{hpi} bandwidth values that were evaluated.}
}}
}
\value{
A Plotly object representing a 3D scatter plot with:
\describe{
\item{x-axis}{Evaluated \code{hmu} bandwidth values.}
\item{y-axis}{Evaluated \code{hpi} bandwidth values.}
\item{z-axis}{Corresponding cross-validation scores.}
}
}
\description{
This function generates an interactive 3D scatter plot of cross-validation results using Plotly.
It accepts a results list returned by \code{cv_grid_search} and visualizes the grid of cross-validation
scores along with the evaluated bandwidth values (\code{hmu} and \code{hpi}).
}