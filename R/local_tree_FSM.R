#' Select edges for local tree graph
#'
#' Construct a local tree from a `FSM_ARG` object with given site.
#'
#' @param ARG A `FSM_ARG` object.
#' @param location The site location for local tree.
#' @return `localtree` object; Local tree at a chosen site.
local_tree_FSM <- function(ARG, location) {
  keep_edge <- which(as.logical(ARG$edge_mat[, location]))
  edge_index <- 1:nrow(ARG$edge)

  ARG$edge <- ARG$edge[keep_edge, ]
  ARG$edge_mat <- ARG$edge_mat[keep_edge, ]
  ARG$edge_index <- edge_index[keep_edge]

  ARG$edge <- ARG$edge[order(ARG$edge[, 1]), ]
  duplicated_edge <- duplicated(ARG$edge[, 1]) | duplicated(ARG$edge[, 1], fromLast = T)
  last_duplicated <- tail(which(duplicated_edge), 1)
  ARG$edge <- ARG$edge[1:last_duplicated, ]

  class(ARG) <- "localtree"
  return(ARG)
}
