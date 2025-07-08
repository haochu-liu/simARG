#' Construct clonal tree from ARG
#'
#' Select clonal lineages to create a clonal tree.
#'
#' @param ARG An `FSM_ARG` object.
#' @return `localtree` object; the clonal tree from ARG.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' clonal_tree <- clonal_tree_FSM(ARG)
clonal_tree_FSM <- function(ARG) {
  if (!inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'FSM_ARG'!")
  } else if (is.null(ARG$node_clonal)) {
    cli::cli_abort("`ARG` must have clonal information!")
  }

  edge_clonal <- which(ARG$node_clonal[ARG$edge[, 1]] &
                       ARG$node_clonal[ARG$edge[, 2]])
  edge_index <- 1:nrow(ARG$edge)

  ARG$edge <- ARG$edge[edge_clonal, ]
  ARG$edge_mat <- ARG$edge_mat[edge_clonal, ]
  ARG$edge_index <- edge_index[edge_clonal]

  ARG$edge <- ARG$edge[order(ARG$edge[, 1]), ]
  duplicated_edge <- duplicated(ARG$edge[, 1]) | duplicated(ARG$edge[, 1], fromLast = T)
  last_duplicated <- tail(which(duplicated_edge), 1)
  ARG$edge <- ARG$edge[1:last_duplicated, ]

  class(ARG) <- "localtree"

  ARG$waiting_time <- NULL
  ARG$k <- NULL
  ARG$sum_time <- tree_height(ARG)
  if (!is.null(ARG$mutation)) {
    ARG$mutation[ARG$mutation[, 1] %in% ARG$edge_index, ]
  }

  return(ARG)
}
