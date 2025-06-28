#' Construct a layout matrix for tree-like ARG
#'
#' Provide the layout coordinates for `igraph::plot.igraph`.
#'
#' @param ARG An `ISM_ARG` or `FSM_ARG` object.
#' @return A layout matrix for `igraph::plot.igraph`.
ARG_igraph_layout <- function(ARG) {
  if (!inherits(ARG, "ISM_ARG") & !inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'ISM_ARG' or 'FSM_ARG'")
  } else if (inherits(ARG, "sim_ISM_ARG")) {
    n <- nrow(ARG$node)
  } else {
    n <- length(ARG$node_height)
  }

  # create layout matrix and leaf nodes
  edge_matrix <- as.matrix(ARG$edge[, 1:3])
  layout_coord <- matrix(data=NA, nrow=n, ncol=2)
  leaf_nodes <- unique(edge_matrix[edge_matrix[, 2] %in% 1:ARG$n, 2])
  layout_coord[leaf_nodes, 1] <- 1:ARG$n * 3
  layout_coord[leaf_nodes, 2] <- 0
  current_level <- 1

  # store the recombination nodes
  recomb_nodes <- edge_matrix[duplicated(edge_matrix[, 2]) | duplicated(edge_matrix[, 2],
                                                                        fromLast=TRUE), 2]
  recomb_nodes <- unique(recomb_nodes)
  recomb_parents <- edge_matrix[edge_matrix[, 2] %in% recomb_nodes, 1]

  # add node coordinate to layout matrix
  for (i in (ARG$n+1):n) {
    target_node <- i
    if (target_node %in% recomb_parents) {
      # if the target node is from recombination
      if (is.na(layout_coord[i, 1])) {
        base_node <- edge_matrix[which(target_node == edge_matrix[, 1]), 2]
        target_node <- edge_matrix[edge_matrix[, 2] %in% base_node, 1] # get all two nodes
        layout_coord[target_node, 1] <- c(-1, 1) + layout_coord[base_node, 1]
        layout_coord[target_node, 2] <- current_level
        current_level <- current_level + 1
      }
    } else {
      # if the target node is from coalescent
      base_node <- edge_matrix[which(target_node == edge_matrix[, 1]), 2]
      layout_coord[target_node, 1] <- mean(layout_coord[base_node, 1])
      layout_coord[target_node, 2] <- current_level
      current_level <- current_level + 1
    }
  }
  return(layout_coord)
}
