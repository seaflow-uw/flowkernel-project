# Generated from _main.Rmd: do not edit by hand

#' Initialize the Kernel EM-algorithm using constant parameters
#' 
#' @param y length T list with `y[[t]]` being a n_t-by-d matrix
#' @param K number of components
#' @param times_to_sample number of time points to sample
#' @param points_to_sample number of bins to sample from each sampled time point
#' @export
init_const <- function (y, K, times_to_sample = 50, points_to_sample = 50){
  num_times <- length(y)
  d <- ncol(y[[1]])
  mu <- array(NA, c(num_times, K, d))
  Sigma <- array(NA, c(num_times, K, d, d))
  pi <- matrix(NA, num_times, K)
  # Repeatedly call Mclust until it gives a non-NULL fit:
  init_fit <- NULL
  while (is.null(init_fit)) {
    # subsample data:
    sampled_times <- sample(num_times, times_to_sample, replace=TRUE)
    if (d == 1) {
      sample_data <- y[sampled_times] %>%
        purrr::map(~ .x[sample(nrow(.x), points_to_sample, replace=TRUE)]) %>% 
        unlist()
    }
    else {
      sample_data <- y[sampled_times] %>%
        purrr::map(~ t(.x[sample(nrow(.x), points_to_sample, replace=TRUE), ])) %>% 
        unlist() %>% 
        matrix(ncol = d, byrow = TRUE)
    }
    
    
    if (d == 1) {
      init_fit <- mclust::Mclust(sample_data, G = K, modelNames = "V")
      for (tt in seq(num_times)) {
        mu[tt, , 1] <- init_fit$parameters$mean
        Sigma[tt, , 1, 1] <- init_fit$parameters$variance$sigmasq
        pi[tt, ] <- init_fit$parameters$pro
      }
    } else if (d > 1) {
      init_fit <- mclust::Mclust(sample_data, G = K, modelNames = "VVV")
      if (is.matrix(init_fit$parameters$mean)){
        for (tt in seq(num_times)) {
          mu[tt, ,] <- t(init_fit$parameters$mean)
          pi[tt, ] <- init_fit$parameters$pro
          Sigma[tt, , , ] <- aperm(init_fit$parameters$variance$sigma, c(3,1,2))
        }
      }
    }
  }
  
  #calculate responsibilities
  resp <- calculate_responsibilities(y, mu, Sigma, pi)
  zest <- resp %>% purrr::map(~ max.col(.x))
  list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
}
