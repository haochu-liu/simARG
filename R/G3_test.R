#' Three-gamete test
#'
#' Compute the number of combinations that have patterns 01, 10, 11.
#'
#' @param mat An incidence matrix with two columns.
#' @return Number of combinations.
#' @export
#'
#' @examples
#' r_square <- G3_test(mat)
G3_test <- function(mat) {
  if (ncol(mat) != 2) {
    cli::cli_abort("`mat` must have two columns.")
  } else if (nrow(mat) < 3) {
    cli::cli_abort("`mat` must have at least 3 rows.")
  }

  num01 <- sum(!mat[, 1] & mat[, 2])
  num10 <- sum(mat[, 1] & !mat[, 2])
  num11 <- sum(mat[, 1] & mat[, 2])

  num_logical <- as.logical(c(num01, num10, num11))
  if (all(num_logical)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
