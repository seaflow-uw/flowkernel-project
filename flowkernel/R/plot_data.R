# Generated from _main.Rmd: do not edit by hand

#' Plot raw data
#' 
#' @param y length T list with `y[[t]]` being a n_t-by-d matrix
#' 
#' @export
plot_data <- function(y) {
  d <- ncol(y[[1]])
  if (d == 1){
    y_label <- ifelse(is.null(colnames(y[[1]])), "V1", colnames(y[[1]]))
    fig <- purrr::map_dfr(y, ~ tibble::tibble(y = .x), .id = "time") %>% 
    dplyr::mutate(time = as.numeric(.data$time)) %>% 
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$y)) +
      ggplot2::geom_point(alpha = 0.2) + 
      ggplot2::labs(x = "Time", y = y_label, title = "Raw Data")
  }
  else if (d == 3){
      d <- ncol(y[[1]])
      max_val <- list()
      max_val_time <- list()
      min_val = list()
      min_val_time = list()
      for (dd in seq(d)) {
        max_val[[dd]] <- sapply(y, function(mat) max(mat[, dd]))
        max_val_time[[dd]] <- max(max_val[[dd]])
        min_val[[dd]] <- sapply(y, function(mat) min(mat[, dd]))
        min_val_time[[dd]] <- min(min_val[[dd]])
      }
      y = unname(y)
      
      # Determine axis labels
      if (is.null(colnames(y[[1]]))) {
        x_label <- "V1"
        y_label <- "V2"
        z_label <- "V3"
      } else {
        x_label <- colnames(y[[1]])[1]
        y_label <- colnames(y[[1]])[2]
        z_label <- colnames(y[[1]])[3]
      }
      y <- y %>% 
        purrr::map_dfr(~ tibble::tibble(x = .x[, 1], y = .x[, 2], z = .x[, 3]), .id = "time")
      y$time = as.integer(y$time)
      
      fig <- plotly::plot_ly(
        data = y,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d", frame = ~time, mode = "markers", size = 80,
        colors = colorRamp(c("blue", "red", "purple", "cyan", "magenta", "brown", "gray", "darkgreen", "darkblue", "darkred", "darkorange"))) %>%
        plotly::layout(title = 'Raw Data', scene = list(
          xaxis = list(range = c(1.1 * min_val_time[[1]], 1.1 * max_val_time[[1]]), title = x_label), 
          yaxis = list(range = c(1.1 * min_val_time[[2]], 1.1 * max_val_time[[2]]), title = y_label),
          zaxis = list(range = c(1.1 * min_val_time[[3]], 1.1 * max_val_time[[3]]), title = z_label),  
          aspectmode = "manual", 
          aspectratio = list(x = 1, y = 1, z = 1)  # Specify the fixed aspect ratio
        ))
  }
  return(fig)
}
