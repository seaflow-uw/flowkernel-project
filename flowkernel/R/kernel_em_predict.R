# Generated from _main.Rmd: do not edit by hand

#' @param fit A list returned by `kernel_em` containing the fitted parameters (mu, Sigma, pi, and responsibilities) from the training data.
#' @param test_dates A vector of POSIXct (or POSIXt) date-time objects of length T_test representing the time points at which predictions are required.
#' @param train_dates A vector of POSIXct (or POSIXt) date-time objects of length T_train representing the time points corresponding to the training data.
#' @param train_data A list of length T_train, where each element is a matrix containing the observations for the training set at the corresponding time point.
#' @param train_biomass A list of length T_train, where each element is a numeric vector of length n_t containing the biomass (or count) for the observations in the training set.
#' @param hmu Bandwidth parameter for smoothing the mu estimates.
#' @param hSigma Bandwidth parameter for smoothing the Sigma estimates.
#' @param hpi Bandwidth parameter for smoothing the pi (mixing proportions) estimates.
#' @param scale_dates Logical flag indicating whether to rescale the date-time objects. If TRUE, the dates are converted to numeric values and rescaled (by subtracting the minimum training date and dividing by 3600); if FALSE, the dates are used directly. Defaults to TRUE.
#' @export
kernel_em_predict <- function(fit, 
                                test_dates, 
                                train_dates, 
                                train_data, 
                                train_biomass, 
                                hmu, 
                                hSigma, 
                                hpi,
                                scale_dates = TRUE) {
  num_test <- length(test_dates)
  num_train <- length(train_dates)
  K <- ncol(fit$pi)
  d <- dim(fit$mu)[3]
  mu <- array(NA, c(num_test, K, d))
  Sigma <- array(NA, c(num_test, K, d, d))
  pi <- matrix(NA, num_test, K)
  
  # Convert to numeric
  numeric_train_dates <- as.numeric(train_dates)
  numeric_test_dates  <- as.numeric(test_dates)
  
  if (scale_dates) {
    # If we are given real dates, then we do the original rescaling
    min_date <- min(numeric_train_dates)
    rescaled_train_dates <- (numeric_train_dates - min_date) / 3600
    rescaled_test_dates  <- (numeric_test_dates  - min_date) / 3600
  } else {
    # If dates are already in hours or integer indices, treat them directly
    rescaled_train_dates <- numeric_train_dates
    rescaled_test_dates  <- numeric_test_dates
  }
  
  # 1) Predict pi (mixing proportions) ----------------------------------------
  # Weighted responsibilities (weights = biomass)
  resp_weighted <- purrr::map2(train_biomass, fit$resp, ~ .y * .x)
  # Sum over each mixture component across all observations
  resp_sum <- purrr::map(resp_weighted, ~ colSums(.x)) %>%
    unlist() %>%
    matrix(ncol = K, byrow = TRUE)
  
  # Smooth those sums across time
  resp_sum_smooth <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hpi,
                     x.points = rescaled_test_dates)$y
  )
  
  # Normalize to get pi
  pi <- resp_sum_smooth / rowSums(resp_sum_smooth)
  
  # 2) M-step for mu ---------------------------------------------------------
  # Weighted sum of y's for each mixture component
  y_sum <- purrr::map2(resp_weighted, train_data, ~ crossprod(.x, .y)) %>% 
    unlist() %>% 
    array(c(K, d, num_train)) %>% 
    aperm(c(3,1,2))
  
  # Smooth each dimension across time
  y_sum_smoothed <- apply(
    y_sum, 2:3, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hmu,
                     x.points = rescaled_test_dates)$y
  )
  
  # Smooth the sum of responsibilities (again, but with bandwidth hmu)
  resp_sum_smooth_mu <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hmu,
                     x.points = rescaled_test_dates)$y
  )
  
  # Combine smoothed sums to get mu
  for (j in seq(d)) {
    mu[, , j] <- y_sum_smoothed[, , j] / resp_sum_smooth_mu
  }
  
  # 3) M-step for Sigma ------------------------------------------------------
  mat_sum <- array(NA, c(num_train, K, d, d))
  for (tt in seq(num_train)) {
    # Prepare a matrix for each observation's difference from mu
    yy <- matrix(NA, nrow(train_data[[tt]]), d)
    for (k_idx in seq(K)) {
      for (dd in seq(d)) {
        yy[, dd] <- train_data[[tt]][, dd] - fit$mu[tt, k_idx, dd]
      }
      mat_sum[tt, k_idx, , ] <- crossprod(yy, yy * resp_weighted[[tt]][, k_idx])
    }
  }
  
  # Smooth Sigma the same way
  mat_sum_smoothed <- apply(
    mat_sum, 2:4, function(x)
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hSigma,
                     x.points = rescaled_test_dates)$y
  )
  
  resp_sum_smooth_Sigma <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hSigma,
                     x.points = rescaled_test_dates)$y
  )
  
  for (j in seq(d)) {
    for (l in seq(d)) {
      Sigma[, , j, l] <- mat_sum_smoothed[, , j, l] / resp_sum_smooth_Sigma
    }
  }
  
  list(mu = mu, Sigma = Sigma, pi = pi)
}
