#' Input: tree (localtree), label
#' Convert the local tree to a phylo object that can be plotted by ape::plot.phylo
#' Output: a phylo object

#' Convert a `localtree` object to a `phylo` object
#'
#' Change a local tree to a `phylo` object in `ape` package.
#'
#' @param tree A `localtree` object.
#' @param label logical; If TRUE, add string labels on leaf lineages.
#' @return `phylo` object in `ape` package.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' local_tree <- local_tree(ARG_mutation, 2L)
#' phylo_tree1 <- localtree_to_phylo(local_tree, label=TRUE)
localtree_to_phylo <- function(tree, label=FALSE) {
  if (!inherits(tree, "localtree")) {
    cli::cli_abort("Object must be of class 'localtree'")
  }

  # find the node1 that only occurs once
  edge_df <- as.matrix(tree$edge[, 1:3])
  counts_node1 <- table(edge_df[, 1])
  single_node1 <- as.numeric(names(counts_node1[counts_node1 == 1]))

  if (length(single_node1)) {
    # for loop to modify the edges
    for (i in single_node1) {
      delete_edge <- which(edge_df[, 1] %in% i)
      target_edge <- which(edge_df[, 2] %in% i)
      edge_df[target_edge, 2] <- edge_df[delete_edge, 2]
      edge_df[target_edge, 3] <- edge_df[target_edge, 3] +
        edge_df[delete_edge, 3]
      edge_df[delete_edge, ] <- NA
    }
    edge_df <- na.omit(edge_df)
  }

  # update the node index for phylo object
  unique_node1 <- sort(unique(edge_df[, 1]))
  edge_matrix <- edge_df[, 1:2]
  for (i in 1:length(unique_node1)) {
    edge_matrix[edge_df[, 1:2] == unique_node1[i]] <- 2*tree$n - i
  }

  if (label & !is.null(tree$node$gene_str)) {
    leaf_labels <- c()
    for (i in 1:tree$n) {
      node_str <- paste(i, tree$node$gene_str[i], sep=":")
      leaf_labels <- c(leaf_labels, node_str)
    }
  } else {
    leaf_labels <- as.character(1:tree$n)
  }

  # convert the local tree object to phylo object for ape::plot.phylo
  tree_phylo <- list(edge=as.matrix(edge_matrix),
                     edge.length=edge_df[, 3],
                     tip.label=leaf_labels,
                     Nnode=as.integer(tree$n - 1))
  class(tree_phylo) <- "phylo"

  return(tree_phylo)
}
