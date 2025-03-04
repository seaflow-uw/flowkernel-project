# Generated from _main.Rmd: do not edit by hand

#' Compute the Average Log-Likelihood for Predicted Parameters
#'
#' This function computes the average log-likelihood of the observed data over all time points using the predicted parameters.
#'
#' @param y A list of length T, where each element is a matrix of observations at a given time point. Each matrix should have dimensions corresponding to the number of observations at that time point by the data dimension.
#' @param pred A list containing the predicted parameters from the model. It should include:
#'   \describe{
#'     \item{mu}{An array of predicted means with dimensions [T, K, d], where T is the number of time points, K is the number of clusters, and d is the data dimension.}
#'     \item{Sigma}{An array of predicted covariance matrices with dimensions [T, K, d, d].}
#'     \item{pi}{A matrix of predicted mixing proportions with dimensions [T, K].}
#'   }
#'
#' @return A numeric value representing the average log-likelihood over all time points.
#' @export
compute_log_like <- function(y, pred) {
  d <- ncol(y[[1]])  # Dimension of the data
  ntimes <- length(y)  # Number of time points
  K <- dim(pred$pi)[2]  # Number of clusters
  log_like <- numeric(ntimes)
  for (tt in seq(ntimes)) {
    likelihood_tt <- 0
    for (k in seq(K)) {
      if (d == 1) {
        # Univariate case
        likelihood_k <- stats::dnorm(y[[tt]], mean = pred$mu[tt, k, 1], 
                                     sd = sqrt(pred$Sigma[tt, k, 1, 1]))
      } else {
        # Multivariate case
        likelihood_k <- mvtnorm::dmvnorm(y[[tt]], mean = pred$mu[tt, k, ], 
                                         sigma = pred$Sigma[tt, k, , ])
      }
      # Weighted likelihood with mixing proportions
      likelihood_tt <- likelihood_tt + pred$pi[tt, k] * likelihood_k
    }
    log_like[tt] <- sum(log(likelihood_tt))
  }
  
  return(mean(log_like))
}
