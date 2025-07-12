#' ClonalOrigin Recombination Nodes
#'
#' Find the row index for recombination nodes on a given edge.
#'
#' @param matrix_data A recomb_edge matrix.
#' @param n The edge index.
#' @return vector of where are the nodes
ClonalOrigin_nodes <- function(matrix_data, n) {
  rows_type_a <- which(matrix_data[, 3] == n)
  rows_type_b <- which(matrix_data[, 1] == n)
  found_nodes <- matrix(NA, nrow=(length(rows_type_a)+length(rows_type_b)), ncol=2)

  if (length(rows_type_a) > 0) {
    found_nodes[1:length(rows_type_a), 1] <- rows_type_a
    found_nodes[1:length(rows_type_a), 2] <- matrix_data[rows_type_a, 4]
  }

  if (length(rows_type_b) > 0) {
    found_nodes[(length(rows_type_a)+1):(length(rows_type_a)+length(rows_type_b)), 1] <- -rows_type_b
    found_nodes[(length(rows_type_a)+1):(length(rows_type_a)+length(rows_type_b)), 2] <- matrix_data[rows_type_b, 2]
  }

  if (nrow(found_nodes) == 0) {
    return(numeric(0))
  }

  if (nrow(found_nodes) > 1) {
    found_nodes <- found_nodes[order(found_nodes[, 2]), ]
  }

  return(found_nodes[, 1])
}
