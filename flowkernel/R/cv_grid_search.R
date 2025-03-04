# Generated from _main.Rmd: do not edit by hand

#' Grid Search Cross-Validation for Kernel EM Model
#'
#' This function performs a grid search over bandwidth parameters for the kernel-based EM model via cross-validation.
#' The grid search is performed over the \code{hmu} (for the mean parameter) and \code{hpi} (for the mixing proportions)
#' while keeping \code{hSigma} fixed. Cross-validation is implemented by partitioning the data into folds by leaving out every
#' \code{leave_out_every}-th time point. For each parameter combination, the model is fit on the training data and evaluated
#' on the left-out data using the log-likelihood. The best parameters are selected as those that maximize the cross-validation score.
#'
#' @param y A list of length T, where each element is a matrix of observations at a given time point.
#' @param K Number of mixture components (clusters) in the model.
#' @param hSigma Bandwidth parameter for smoothing the \code{Sigma} (covariance) estimates.
#' @param dates A vector of date-time objects (e.g., of class \code{POSIXct}) of length T corresponding to the time points in \code{y}.
#'   If \code{NULL}, the function uses the sequence of indices.
#' @param biomass A list of length T, where each element is a numeric vector of length \eqn{n_t} containing the biomass (or count)
#'   for the observations at that time point. Defaults to \code{default_biomass(y)}.
#' @param leave_out_every The number of folds for cross-validation (default is 5). The \eqn{\ell}-th fold is defined as the set of
#'   time points \eqn{\{\ell, \ell + \text{leave_out_every}, \ell + 2 \times \text{leave_out_every}, \dots\}}.
#' @param grid_size The number of grid points to evaluate for the bandwidth parameters.
#' @param ncores Number of cores to use for parallel processing. If \code{NULL} or less than 2, the grid search is performed sequentially.
#' @param log_grid Logical flag indicating whether the grid for the bandwidth parameters should be logarithmically spaced (if \code{TRUE})
#'   or linearly spaced (if \code{FALSE}). Defaults to \code{TRUE}.
#' @param grid_length_half Logical flag indicating whether to use half the total time length (or median time) for computing the end of the grid.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{A data frame (converted from a matrix) of cross-validation scores, with row names corresponding to evaluated \code{hmu} values
#'     and column names to evaluated \code{hpi} values.}
#'     \item{best_hmu}{The best \code{hmu} value (bandwidth for smoothing \code{mu}) according to the highest CV score.}
#'     \item{best_hpi}{The best \code{hpi} value (bandwidth for smoothing \code{pi}) according to the highest CV score.}
#'     \item{max_cv_score}{The highest cross-validation score observed.}
#'     \item{hmu_vals}{The vector of \code{hmu} values that were evaluated.}
#'     \item{hpi_vals}{The vector of \code{hpi} values that were evaluated.}
#'     \item{time_elapsed}{The total time elapsed (in seconds) during the grid search.}
#'   }
#'
#' @details The function first determines whether the dates should be scaled based on whether \code{dates} is provided.
#'   A grid of bandwidth values is then constructed (either logarithmically or linearly spaced) from the calculated range.
#'   Cross-validation is performed for each parameter combination via the \code{kernel_em_cv} function. The grid search can be run
#'   either sequentially or in parallel, depending on the value of \code{ncores}.
#'
#' @export
cv_grid_search <- function(y, 
                                        K, 
                                        hSigma = 10, 
                                        dates = NULL, 
                                        biomass = default_biomass(y), 
                                        leave_out_every = 5, 
                                        grid_size = 10,
                                        ncores = NULL,
                                        log_grid = TRUE,
                                        grid_length_half = TRUE) {
  start_time <- Sys.time()
  
  # Decide if we are scaling dates or not
  if (is.null(dates)) {
    numeric_dates <- seq_along(y)
    scale_dates <- FALSE
    starting_date <- min(numeric_dates)
    middle_date   <- stats::median(numeric_dates)
    
    if (grid_length_half) {
      end_cv <- middle_date - starting_date
    } else {
      end_cv <- max(numeric_dates) - starting_date
    }
  } else {
    numeric_dates <- as.numeric(dates)
    scale_dates   <- TRUE
    starting_date <- min(numeric_dates)
    middle_date   <- stats::median(numeric_dates)
    
    if (grid_length_half) {
      end_cv <- (middle_date - starting_date) / 3600
    } else {
      end_cv <- (max(numeric_dates) - starting_date) / 3600
    }
  }

  if (log_grid) {
    # LOG-SPACED GRID
    h_values <- exp(seq(log(end_cv), log(2), length = grid_size))
  } else {
    # LINEARLY-SPACED GRID
    h_values <- seq(end_cv, 2, length.out = grid_size)
  }
  
  h_values <- round(h_values)
  h_values <- unique(h_values)
  h_values <- h_values[h_values > 0]
  
  hmu_vals <- h_values
  hpi_vals <- h_values
  
  # Create a data frame of all parameter combinations
  param_grid <- expand.grid(hmu = hmu_vals, hpi = hpi_vals)
  
  # Decide whether to run in parallel or sequentially
  if (is.null(ncores) || as.numeric(ncores) < 2) {
    message("Running sequentially...")
    
    # Prepare an empty data frame to store results
    cv_scores <- data.frame(hmu = numeric(0), hpi = numeric(0), cv_score = numeric(0))
    total <- nrow(param_grid)
    
    for (i in seq_len(nrow(param_grid))) {
      hmu <- param_grid$hmu[i]
      hpi <- param_grid$hpi[i]
      processed <- i
      remaining <- total - i
      
      # Print detailed progress information
      message("Evaluating parameter combination ", processed, " of ", total, 
              ": hmu = ", hmu, ", hpi = ", hpi, ". Remaining: ", remaining)
      
      cv_score <- kernel_em_cv(
        y             = y,
        K             = K,
        dates         = numeric_dates,
        hmu           = hmu,
        hSigma        = hSigma,
        hpi           = hpi,
        biomass       = biomass,
        leave_out_every = leave_out_every,
        scale_dates   = scale_dates
      )
      
      cv_scores <- rbind(cv_scores, data.frame(hmu = hmu, hpi = hpi, cv_score = cv_score))
    }
    
  } else {
    message("Running in parallel on ", ncores, " cores...")
    
    cl <- parallel::makeCluster(as.integer(ncores))
    doParallel::registerDoParallel(cl)
    
    # Export needed variables/functions
    parallel::clusterExport(
      cl, 
      varlist = c("kernel_em_cv", "ensure_matrix", "kernel_em", "kernel_em_predict",
                  "kernel_em_predict", "calculate_responsibilities", 
                  "compute_log_like", "y", "K", "hSigma", "biomass", 
                  "leave_out_every", "numeric_dates", "scale_dates"),
      envir = environment()
    )
    
    # Load required packages on the workers
    parallel::clusterEvalQ(cl, {
      library(flowkernel) 
    })
    
    cv_scores <- foreach::foreach(i = seq_len(nrow(param_grid)), .combine = rbind) %dopar% {
      hmu <- param_grid$hmu[i]
      hpi <- param_grid$hpi[i]
      
      cv_score <- kernel_em_cv(
        y             = y,
        K             = K,
        dates         = numeric_dates,
        hmu           = hmu,
        hSigma        = hSigma,
        hpi           = hpi,
        biomass       = biomass,
        leave_out_every = leave_out_every,
        scale_dates   = scale_dates
      )
      
      data.frame(hmu = hmu, hpi = hpi, cv_score = cv_score)
    }
    
    parallel::stopCluster(cl)
  }
  
  # Convert cv_scores to matrix form
  if (nrow(cv_scores) == 0) {
    stop("No cross-validation scores were computed.")
  }
  
  results <- reshape2::acast(cv_scores, hmu ~ hpi, value.var = "cv_score")
  results_df <- as.data.frame(results)
  
  # Re-label columns and rows for clarity
  colnames(results_df) <- paste0("hpi_", colnames(results_df))
  rownames(results_df) <- paste0("hmu_", rownames(results_df))
  
  # Identify best hmu/hpi
  max_pos  <- which(results == max(results, na.rm = TRUE), arr.ind = TRUE)
  best_hmu <- as.numeric(rownames(results)[max_pos[1]])
  best_hpi <- as.numeric(colnames(results)[max_pos[2]])
  
  # Timing
  end_time  <- Sys.time()
  run_time  <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours     <- floor(run_time / 3600)
  minutes   <- floor((run_time %% 3600) / 60)
  seconds   <- round(run_time %% 60)
  
  cat("Best hmu:", best_hmu, 
      "\nBest hpi:", best_hpi, 
      "\nHighest CV score:", max(results, na.rm = TRUE), 
      "\nTime Elapsed:", run_time, "seconds =>",
      hours, "hours", minutes, "minutes", seconds, "seconds\n")
  
  results_list <- list(
    results      = results_df,
    best_hmu     = best_hmu,
    best_hpi     = best_hpi,
    max_cv_score = max(results, na.rm = TRUE),
    hmu_vals = hmu_vals,
    hpi_vals = hpi_vals,
    time_elapsed = run_time
  )
  return(results_list)
}
