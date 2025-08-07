#' Clonal Genealogy simulation
#'
#' Simulate coalescent backwards in time to construct a clonal tree.
#'
#' @param n An integer for the number of leaf lineages.
#' @return A list containing clonal tree.
#' @export
#'
#' @examples
#' tree1 <- clonal_genealogy(10L)
#' tree2 <- clonal_genealogy(15L)
clonal_genealogy <- function(n) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  }

  k = n
  t_sum <- 0

  # Initialize varables for clonal tree
  clonal_edge <- matrix(NA, nrow=2*(n-1), ncol=3) # root and leaf nodes, length
  colnames(clonal_edge) <- c("node1", "node2", "length")
  clonal_node_height <- rep(NA, 2*n-1)            # node height to recent time
  clonal_node_height[1:n] <- 0                    # initialize first n nodes

  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)

  # clonal tree by coalescent only
  while (k > 1) {
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1)/2)
    t_sum <- t_sum + event_time
    # coalescent event
    leaf_node <- sample(pool, size=2, replace=FALSE)

    # append edges
    clonal_edge[c(edge_index, edge_index+1), 1] <- node_index
    clonal_edge[c(edge_index, edge_index+1), 2] <- leaf_node
    clonal_edge[c(edge_index, edge_index+1), 3] <- t_sum-clonal_node_height[leaf_node]

    # append root node
    clonal_node_height[node_index] <- t_sum

    # updates for iteration
    pool <- c(setdiff(pool, leaf_node), node_index)
    edge_index <- edge_index + 2L
    node_index <- node_index + 1L
    k <- k - 1
  }

  tree = list(edge=clonal_edge,
              node=clonal_node_height,
              n=n)
  class(tree) <- "clonal_tree"
  return(tree)
}
