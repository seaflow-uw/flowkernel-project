# Generated from _main.Rmd: do not edit by hand

#' Plot Cross-Validation Results
#'
#' This function generates an interactive 3D scatter plot of cross-validation results using Plotly.
#' It accepts a results list returned by \code{cv_grid_search} and visualizes the grid of cross-validation
#' scores along with the evaluated bandwidth values (\code{hmu} and \code{hpi}).
#'
#' @param results_list A list returned by \code{cv_grid_search}. This list should include:
#'   \describe{
#'     \item{results}{A data frame of cross-validation scores, where the row names correspond to evaluated \code{hmu} values and the column names to evaluated \code{hpi} values.}
#'     \item{hmu_vals}{A numeric vector of the \code{hmu} bandwidth values that were evaluated.}
#'     \item{hpi_vals}{A numeric vector of the \code{hpi} bandwidth values that were evaluated.}
#'   }
#'
#' @return A Plotly object representing a 3D scatter plot with:
#'   \describe{
#'     \item{x-axis}{Evaluated \code{hmu} bandwidth values.}
#'     \item{y-axis}{Evaluated \code{hpi} bandwidth values.}
#'     \item{z-axis}{Corresponding cross-validation scores.}
#'   }
#' @export
plot_cv_results <- function(results_list) {
  # Pull out the results
  results_df <- results_list$results
  
  # Convert row/column names from strings ("hmu_99") to numeric (99)
  row_nums <- as.numeric(sub("hmu_", "", rownames(results_df)))
  col_nums <- as.numeric(sub("hpi_", "", colnames(results_df)))
  
  # Sort the numeric row and column names
  row_sorted <- sort(row_nums)
  col_sorted <- sort(col_nums)
  
  # Create a new, reordered matrix
  # match() finds the indices of row_nums that correspond to row_sorted, etc.
  row_order <- match(row_sorted, row_nums)
  col_order <- match(col_sorted, col_nums)
  
  # Reorder results_df so row i truly corresponds to row_sorted[i], etc.
  results_mat <- as.matrix(results_df[row_order, col_order])
  
  # Flatten in column-major order
  cv_scores <- c(results_mat)
  
  # Build a data frame matching the sorted bandwidths to their CV scores
  df <- expand.grid(hmu = row_sorted, hpi = col_sorted)
  df$cv_score <- cv_scores
  
  # Plot in 3D
  plotly::plot_ly(
    data = df,
    x = ~hmu,
    y = ~hpi,
    z = ~cv_score,
    type = "scatter3d",
    mode = "markers"
  ) %>%
    plotly::layout(
      title = "Cross-Validation Results",
      scene = list(
        xaxis = list(title = "hmu"),
        yaxis = list(title = "hpi"),
        zaxis = list(title = "CV Score")
      )
    )
}
