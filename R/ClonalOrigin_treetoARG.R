#' An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm
#'
#' Simulate coalescent an recombination events by ClonalOrigin algorithm.
#' Recombinations are added to a clonal tree to construct an ARG.
#'
#' @param n An integer for the number of leaf lineages.
#' @param rho The recombination parameter.
#' @param L An integer for the number of sites.
#' @param delta numeric; If bacteria = TRUE, delta is the mean of recombinant segment length.
#' @param optimise_recomb numeric; If TRUE, skip the recombinations that containing no effective segment.
#' @param edgemat numeric; If TRUE, return full edge material matrix.
#' @return A list containing matrices of edges and nodes, and other information about ARG.
#' @export
#'
#' @examples
#' ARG1 <- ClonalOrigin_treetoARG(100L, 5, 100L, 5)
#' ARG2 <- ClonalOrigin_treetoARG(5L, 1, 10L, 1, optimise_recomb=TRUE)
ClonalOrigin_treetoARG <- function(n, rho, L, delta,
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

  # number of recombination edges
  tree_length <- sum(clonal_edge[, 3])
  # num of recombs | tree_length ~ Poisson(rho*l/2)
  n_recomb <- rpois(1, rho*tree_length/2)

  # stop when there is no recombination
  if (n_recomb == 0) {
    if (edgemat) {
      ARG = list(edge=clonal_edge,
                 edge_mat=matrix(1, nrow=nrow(clonal_edge), ncol=L),
                 node_height=clonal_node_height,
                 node_mat=matrix(1, nrow=length(clonal_node_height), ncol=L),
                 node_clonal=rep(TRUE, length(clonal_node_height)),
                 sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
    } else {
      ARG = list(edge=clonal_edge,
                 edge_mat_index=clonal_edge[1, ],
                 node_height=clonal_node_height,
                 node_mat=matrix(1, nrow=length(clonal_node_height), ncol=L),
                 node_clonal=rep(TRUE, length(clonal_node_height)),
                 sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
    }
    class(ARG) <- "FSM_ARG"
    return(ARG)
  }

  # Initialize recombination edges
  recomb_edge <- matrix(NA, nrow=n_recomb, ncol=6) # matrix for b, a, x, y
  colnames(recomb_edge) <- c("b_edge", "b_height",
                             "a_edge", "a_height",
                             "x", "y")
  probstartcum <- cumsum(rep(1/L, L))
  a_rexp <- rexp(n_recomb, rate=1)

  # simulate b_edge (similar to mutation)
  recomb_edge[, 1] <- sample(1:(2*(n-1)), n_recomb,
                             replace=TRUE, prob=clonal_edge[, 3])

  # simulate x and y
  recomb_edge[, 5] <- findInterval(runif(n_recomb), probstartcum) + 1
  recomb_edge[, 6] <- pmin(recomb_edge[, 5] + rgeom(n_recomb, 1/delta), L)
  for (i in 1:n_recomb) {
    # simulate b_height
    recomb_edge[i, 2] <- runif(1, max=clonal_edge[recomb_edge[i, 1], 3]) +
                         clonal_node_height[clonal_edge[recomb_edge[i, 1], 2]]
    # identify a_height
    t_above_b <- clonal_node_height[n:(2*n-1)] - recomb_edge[i, 2]
    i_above_b <- c(0, t_above_b[t_above_b >= 0])
    i_above_b <- i_above_b[2:length(i_above_b)] - i_above_b[1:(length(i_above_b)-1)]
    cuml_above_b <- cumsum(i_above_b * (1+length(i_above_b)):2)
    num_lineage <- (1+length(i_above_b)) - length(which(a_rexp[i] > cuml_above_b))
    if (num_lineage == (1+length(i_above_b))) {
      recomb_edge[i, 4] <- a_rexp[i] + recomb_edge[i, 2]
    } else {
      recomb_edge[i, 4] <- (a_rexp[i]-cuml_above_b[1+length(i_above_b)-num_lineage]) / num_lineage +
                           sum(i_above_b[1:(1+length(i_above_b)-num_lineage)]) +
                           recomb_edge[i, 2]
    }
    # simulate a_edge
    if (num_lineage > 1) {
      # print(which(recomb_edge[i, 4] < clonal_node_height)[1] == (2*n+1-num_lineage))
      pool_edge <- which((clonal_edge[, 1] >= (2*n+1-num_lineage)) &
                         (clonal_edge[, 2] < (2*n+1-num_lineage)))
      recomb_edge[i, 3] <- sample(pool_edge, 1, replace=TRUE)
    }
  }

  # Initialize output variables
  t_sum <- max(t_sum, recomb_edge[, 4])
  clonal_edge <- rbind(clonal_edge, matrix(NA, nrow=2*n_recomb, ncol=3))
  node_mat <- matrix(NA, nrow=(2*n_recomb+2*n-1), ncol=L)
  edge_mat_index <- c(clonal_edge[, 2], rep(NA, 2*n_recomb))
  node_clonal <- c(rep(TRUE, 2*n-1), rep(NA, 2*n_recomb))

  # rearrange and backwards in time
  for (i in 1:n_recomb){

  }



  if (edgemat) {
    ARG = list(edge=,
               edge_mat=,
               node_height=,
               node_mat=,
               node_clonal=,
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  } else {
    ARG = list(edge=,
               edge_mat_index=,
               node_height=,
               node_mat=,
               node_clonal=,
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  }
  class(ARG) <- "FSM_ARG"
  return(ARG)
}
