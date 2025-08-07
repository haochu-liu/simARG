#' The square of correlation coefficient for LD
#'
#' Compute the square of Pearson’s correlation coefficient for linkage disequilibrium.
#'
#' @param mat An incidence matrix with two columns.
#' @return r square.
#' @export
#'
#' @examples
#' r_square <- LD(mat)
LD_r <- function(mat) {
  if (ncol(mat) != 2) {
    cli::cli_abort("`mat` must have two columns.")
  }

  n <- nrow(mat)

  n_A <- sum(mat[, 1])
  n_a <- n - n_A

  n_B <- sum(mat[, 2])
  n_b <- n - n_B

  n_AB <- sum(mat[, 1] & mat[, 2])

  D_AB <- n_AB / n - (n_A * n_B) / (n^2)
  r_square <- D_AB^2 * n^4 / (n_A * n_a * n_B * n_b)

  return(r_square)
}
