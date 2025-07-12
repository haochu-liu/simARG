#' Local tree height for `ClonalOrigin` object
#'
#' Find the height of a locat site at give site for a `ClonalOrigin` object.
#'
#' @param ARG A `ClonalOrigin` object.
#' @param location The site location for local tree.
#' @return Height of the local tree
#' @export
local_height_ClonalOrigin <- function(ARG, location) {
  if (!inherits(ARG, "ClonalOrigin")) {
    cli::cli_abort("Object must be of class 'ClonalOrigin'!")
  }

  if (is.null(ARG$recomb_edge)) {
    return(ARG$clonal_time)
  }
  root_node_index <- which(ARG$node_info[, 2]>=ARG$clonal_time)
  if (length(root_node_index)==1) {
    return(ARG$clonal_time)
  }
  root_node <- ARG$node_info[root_node_index, 1]
  root_node_mat <- ARG$node_mat[root_node, ]
  tar_node_index <- root_node_index[which(root_node_mat[, location])]

  return(min(ARG$node_info[tar_node_index, 2]))
}
