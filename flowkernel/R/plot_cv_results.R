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
  # Extract the matrix of results
  results_df <- results_list$results
  hmu_vals   <- results_list$hmu_vals
  hpi_vals   <- results_list$hpi_vals
  
  # Flatten in row-major order by first transposing, then c(...)
  cv_scores <- c(t(as.matrix(results_df)))
  
  # Build a data frame that enumerates (hmu, hpi) in row-major order
  df <- expand.grid(hmu = hmu_vals, hpi = hpi_vals)
  df$cv_score <- cv_scores
  
  # Create a 3D scatter plot
  fig <- plotly::plot_ly(
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
      ),
      margin = list(l = 0, r = 0, b = 0, t = 50)
    )
  
  return(fig)
}
