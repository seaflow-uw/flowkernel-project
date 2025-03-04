# Generated from _main.Rmd: do not edit by hand

#' Ensure an Object is a Matrix
#'
#' This helper function checks whether the provided object is a matrix and, if it is not, converts it to a matrix using \code{as.matrix()}.
#'
#' @param x An object that is expected to be a matrix. If \code{x} is not already a matrix, it will be coerced into one.
#'
#' @return A matrix. If the input was not a matrix, it is converted using \code{as.matrix()}.
#'
#' @examples
#' # Example 1: x is a vector
#' vec <- c(1, 2, 3)
#' mat <- ensure_matrix(vec)
#' is.matrix(mat)  # Should return TRUE
#'
#' # Example 2: x is already a matrix
#' mat2 <- matrix(1:4, nrow = 2)
#' identical(ensure_matrix(mat2), mat2)  # Should return TRUE
#'
#' @export
ensure_matrix <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  return(x)
}
