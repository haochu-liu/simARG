#' Select edges for local tree graph
#'
#' Construct a local tree from a `FSM_ARG` object with given site.
#'
#' @param ARG A `FSM_ARG` object.
#' @param location The site location for local tree.
#' @return `localtree` object; Local tree at a chosen site.
local_tree_FSM <- function(ARG, location) {
  edge_index <- 1:nrow(ARG$edge)
  keep_edge <- c()
  for (i in 1:nrow(ARG$edge_mat)) {
    if (ARG$edge_mat[i, location]) {
      keep_edge <- c(keep_edge, i)
    }
  }

  drop_index <- c()
  duplicated_nodes <- duplicated(ARG$edge[keep_edge, 1]) |
    duplicated(ARG$edge[keep_edge, 1], fromLast=TRUE)
  for (i in length(duplicated_nodes):1) {
    if (duplicated_nodes[i]) {
      break
    } else {
      drop_index <- c(drop_index, i)
    }
  }
  if (length(drop_index)) {keep_edge <- keep_edge[-drop_index]}

  ARG$edge <- ARG$edge[keep_edge, ]
  ARG$edge_mat <- ARG$edge_mat[keep_edge, ]
  ARG$edge_index <- edge_index[keep_edge]

  class(ARG) <- "localtree"
  return(ARG)
}
