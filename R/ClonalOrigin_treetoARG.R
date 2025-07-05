#' An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm
#'
#' Simulate coalescent an recombination events by ClonalOrigin algorithm.
#' Recombinations are added to a clonal tree to construct an ARG.
#'
#' @param n An integer for the number of leaf lineages.
#' @param rho The recombination parameter.
#' @param L An integer for the number of sites.
#' @param delta numeric; If bacteria = TRUE, delta is the mean of recombinant segment length.
#' @param node_max numeric; Initial maximal node size (default = 1000).
#' @param optimise_recomb numeric; If TRUE, skip the recombinations that containing no effective segment.
#' @param edgemat numeric; If TRUE, return full edge material matrix.
#' @return A list containing matrices of edges and nodes, and other information about ARG.
#' @export
#'
#' @examples
#' ARG1 <- ClonalOrigin_treetoARG(100L, 5, 100L, 5)
#' ARG2 <- ClonalOrigin_treetoARG(5L, 1, 10L, 1, optimise_recomb=TRUE)
ClonalOrigin_treetoARG <- function(n, rho, L, delta, node_max=1000,
                                   optimise_recomb=FALSE, edgemat=TRUE) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_integer(L, n=1)) {
    cli::cli_abort("`L` must be a single integer!")
  } else if (n >= node_max) {
    cli::cli_abort("Maximal node size must greater than the number of leaf lineages!")
  }

  k = n
  t_sum <- 0

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
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
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

  # number of recombination edges
  l <- sum(clonal_edge[, 3])
  n_recomb <- rpois(1, rho*l/2) # num of recombs | l ~ Poisson(rho*l/2)

  if (edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta)
  } else {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat_index=edge_mat_index[1:(edge_index-1)],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               delta=delta)
  }
  class(ARG) <- "FSM_ARG"
  return(ARG)
}
