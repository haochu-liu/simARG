#' Local tree height
#'
#' Compute the length from leaf node to root node in the local tree.
#'
#' @param tree An `localtree` object.
#' @return Tree height.
tree_height <- function(tree) {
  if (!inherits(tree, "localtree")) {
    cli::cli_abort("Object must be of class `localtree`")
  }

  # provide the height (time to MRCA) of the tree
  traj <- tree_traj(1, tree)
  return(sum(tree$edge[traj$edge_index, 3]))
}
