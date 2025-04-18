% Generated by roxygen2: do not edit by hand
% Please edit documentation in _main.Rmd.
\name{generate_smooth_gauss_mix}
\alias{generate_smooth_gauss_mix}
\title{Generate data from smoothly-varying mixture of Gaussians model}
\usage{
generate_smooth_gauss_mix(
  mu_function,
  Sigma_function,
  pi_function,
  num_points,
  start_date = "2017-05-28 21:00:00 UTC"
)
}
\arguments{
\item{mu_function}{a function that maps a vector of times to a T-by-K-by-d
array of means}

\item{Sigma_function}{a function that maps a vector of times to a
T-K-by-d-by-d array of covariance matrices}

\item{pi_function}{a function that maps a vector of times to a T-by-K vector
of probabilities}

\item{num_points}{a T vector of integers giving the number of points n_t to
generate at each time point t.}

\item{start_date}{the starting date for the generated data. This is created for when running a test requires dates.}
}
\description{
The smoothly-varying mixture of Gaussians model is defined as follows:
}
\details{
At time t there are n_t points generated as follows:

Y_{it}|\{Z_{it}=k\} ~ N_d(mu_{kt},Sigma_{kt})
where
P(Z_{it}=k)=pi_{kt}
and the parameters (mu_{kt},Sigma_{kt}, pi_{kt}) are all slowly varying in time.

This function generates Y and Z.
}