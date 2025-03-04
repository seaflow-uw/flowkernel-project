# Generated from _main.Rmd: do not edit by hand

#' Cross-Validated Log-Likelihood for Kernel EM Model
#'
#' This function performs K-fold cross-validation (by default, 5-fold) for selecting the bandwidth parameters in a kernel-based EM model. The data is partitioned into folds by leaving out every \code{leave_out_every}-th time point, where the \eqn{\ell}-th fold is given by 
#' \deqn{\mathcal{I}_\ell = \{\ell, \ell+ \text{leave_out_every}, \ell+2 \times \text{leave_out_every}, \dots\}.}
#' For each fold, the model is fit on the remaining time points using \code{kernel_em} and predictions for the left-out points are generated using \code{kernel_em_predict}. The log-likelihood of the left-out data is computed using \code{compute_log_like}, and the final cross-validation score is the average log-likelihood over all folds.
#'
#' @param y A list of length T, where each element is a matrix of observations at a given time point.
#' @param K Number of mixture components (clusters) in the model.
#' @param dates A vector of date-time objects (e.g., of class \code{POSIXct}) of length T corresponding to the time points in \code{y}.
#' @param hmu Bandwidth parameter for smoothing the \code{mu} (mean) estimates.
#' @param hSigma Bandwidth parameter for smoothing the \code{Sigma} (covariance) estimates.
#' @param hpi Bandwidth parameter for smoothing the \code{pi} (mixing proportions) estimates.
#' @param biomass A list of length T, where each element is a numeric vector of length \eqn{n_t} containing the biomass (or count) of particles for the observations at that time point.
#' @param leave_out_every The number of folds for cross-validation (default is 5). The \eqn{\ell}-th fold is defined as the set of time points \eqn{\{\ell, \ell + \text{leave_out_every}, \ell + 2 \times \text{leave_out_every}, \dots\}}.
#' @param scale_dates Logical flag indicating whether to rescale the date-time objects. If \code{TRUE}, the dates are converted to numeric values (by subtracting the minimum date and dividing by 3600) before smoothing; if \code{FALSE}, the dates are used directly.
#'
#' @return A numeric value representing the average log-likelihood over all cross-validation folds.
#'
#' @details For each fold, the model is fit using the training data (all time points not in the current fold) and then predictions are made for the left-out time points. The log-likelihood for each fold is computed by summing the log-likelihoods over all observations and then averaging over the time points in the fold. The overall cross-validation score is obtained by averaging the log-likelihoods over all folds.
#' @export
kernel_em_cv <- function(y, 
                         K, 
                         dates, 
                         hmu, 
                         hSigma, 
                         hpi, 
                         biomass, 
                         leave_out_every = 5,
                         scale_dates = TRUE) 
{
  n <- length(y)
  if (leave_out_every <= 0 || leave_out_every >= n) {
    stop("leave_out_every must be a positive integer smaller than the number of data points")
  }
  
  pred <- vector("list", leave_out_every)
  log_like <- vector("list", leave_out_every)
  
  for (i in seq_len(leave_out_every)) {
    test_indices <- seq(i, n, by = leave_out_every)
    train_data   <- lapply(y[-test_indices], ensure_matrix)
    train_date   <- dates[-test_indices]
    train_biomass <- biomass[-test_indices]
    test_data    <- lapply(y[test_indices], ensure_matrix)
    test_date    <- dates[test_indices]
    test_biomass <- biomass[test_indices]
    
    if (scale_dates){
      fit_train <- kernel_em(y      = train_data,
                             K      = K,
                             hmu    = hmu,
                             hSigma = hSigma,
                             hpi    = hpi,
                             dates  = train_date,
                             biomass = train_biomass)
    } else {
      fit_train <- kernel_em(y      = train_data,
                             K      = K,
                             hmu    = hmu,
                             hSigma = hSigma,
                             hpi    = hpi,
                             biomass = train_biomass)
    }
    
    pred[[i]] <- kernel_em_predict(fit           = fit_train,
                                   test_dates    = test_date,
                                   train_dates   = train_date,
                                   train_data    = train_data,
                                   train_biomass = train_biomass,
                                   hmu           = hmu,
                                   hSigma        = hSigma,
                                   hpi           = hpi,
                                   scale_dates   = scale_dates)
    
    log_like[[i]] <- compute_log_like(test_data, pred[[i]])
  }
  
  avg_log_likelihood <- mean(unlist(log_like))
  return(avg_log_likelihood)
}
