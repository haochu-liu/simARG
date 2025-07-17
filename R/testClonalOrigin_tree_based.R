#' An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm
#'
#' Simulate coalescent an recombination events by ClonalOrigin algorithm.
#' Recombinations are added to a clonal tree to construct an ARG.
#'
#' @param n An integer for the number of leaf lineages.
#' @param rho The recombination parameter.
#' @param L An integer for the number of sites.
#' @param delta numeric; If bacteria = TRUE, delta is the mean of recombinant segment length.
#' @return A list containing clonal tree and recombination edges.
#' @export
#'
#' @examples
#' ARG1 <- ClonalOrigin_tree_based(100L, 5, 100L, 5)
#' ARG2 <- ClonalOrigin_tree_based(5L, 1, 10L, 1)
testClonalOrigin_tree_based <- function(n, rho, L, delta) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_integer(L, n=1)) {
    cli::cli_abort("`L` must be a single integer!")
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
  # n_recomb <- rpois(1, rho*tree_length/2)
  n_recomb <- 1

  # stop when there is no recombination
  if (n_recomb == 0) {
    node_mat <- matrix(TRUE, nrow=(2*n-1), ncol=L)
    ARG = list(clonal_edge=clonal_edge,
               recomb_edge=NULL,
               clonal_node_height=clonal_node_height,
               node_mat=node_mat,
               sum_time=t_sum, n=n, rho=rho, L=L,
               delta=delta, clonal_time=t_sum)
    class(ARG) <- "ClonalOrigin"
    return(ARG)
  }

  # Initialize recombination edges
  recomb_edge <- matrix(NA, nrow=n_recomb, ncol=6) # matrix for b, a, x, y
  colnames(recomb_edge) <- c("b_edge", "b_height",
                             "a_edge", "a_height",
                             "x", "y")
  a_rexp <- rexp(n_recomb, rate=1)
  probstart <- rep(1, L)
  probstart[1] <- delta
  probstart <- probstart / sum(probstart)
  probstartcum <- cumsum(probstart)
  # simulate b_edge (similar to mutation)
  recomb_edge[, 1] <- sample(1:(2*(n-1)), n_recomb,
                             replace=TRUE, prob=clonal_edge[, 3])

  for (i in 1:n_recomb) {
    # simulate x and y
    recomb_edge[i, 5] <- which(runif(1) < probstartcum)[1]
    recomb_edge[i, 6] <- min(recomb_edge[i, 5] + rgeom(1, 1/delta), L)
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
      recomb_edge[i, 4] <- a_rexp[i] / num_lineage + recomb_edge[i, 2]
    } else {
      recomb_edge[i, 4] <- (a_rexp[i]-cuml_above_b[1+length(i_above_b)-num_lineage]) / num_lineage +
                            sum(i_above_b[1:(1+length(i_above_b)-num_lineage)]) +
                            recomb_edge[i, 2]
    }
    # simulate a_edge
    if (num_lineage > 1) {
      pool_edge <- which((clonal_node_height[clonal_edge[, 1]] >= recomb_edge[i, 4]) &
                         (clonal_node_height[clonal_edge[, 2]] < recomb_edge[i, 4]))
      recomb_edge[i, 3] <- sample(pool_edge, 1, replace=TRUE)
    } else {
      recomb_edge[i, 3] <- 2*n - 1
    }
  }

  # recombination segment and ancestral material
  node_mat <- matrix(NA, nrow=(2*n - 1 + 2*n_recomb), ncol=L)
  recomb_node_mat <- matrix(FALSE, nrow=n_recomb, ncol=L)
  node_info <- matrix(NA, nrow=(2*n - 1 + 2*n_recomb), ncol=3)
  colnames(node_info) <- c("index", "height", "recomb")
  node_mat[1:n, ] <- TRUE
  node_info[, 1] <- 1:(2*n - 1 + 2*n_recomb)
  node_info[1:(2*n-1), 2] <- clonal_node_height
  node_info[1:(2*n-1), 3] <- 0
  node_info[(2*n):(2*n-1+n_recomb), 2] <- recomb_edge[, 2]
  node_info[(2*n):(2*n-1+n_recomb), 3] <- -c(1:n_recomb)
  node_info[(2*n+n_recomb):(2*n-1+2*n_recomb), 2] <- recomb_edge[, 4]
  node_info[(2*n+n_recomb):(2*n-1+2*n_recomb), 3] <- 1:n_recomb
  node_info <- node_info[order(node_info[, 2]), ]
  # recombination nodes on every edge
  recomb_node <- lapply(1:(2*n - 1), function(n){
    ClonalOrigin_nodes(recomb_edge, n)
  })
  # Add ancestral material to every node
  for (i in (n+1):(2*n-1+2*n_recomb)) {
    node_index <- node_info[i, 1]
    if (node_info[i, 3]==0) {
      # clonal tree
      leaf_edge <- which(clonal_edge[, 1] == node_index)
      leaf_node <- rep(NA, 2)
      if (length(recomb_node[[leaf_edge[1]]])) {
        # target node is tail(recomb_node[[leaf_edge[1]]], 1)
        tar_node <- tail(recomb_node[[leaf_edge[1]]], 1)
        leaf_node[1] <- node_info[which(tar_node==node_info[, 3]), 1]
      } else {
        leaf_node[1] <- clonal_edge[leaf_edge[1], 2]
      }
      if (length(recomb_node[[leaf_edge[2]]])) {
        # target node is tail(recomb_node[[leaf_edge[2]]], 1)
        tar_node <- tail(recomb_node[[leaf_edge[2]]], 1)
        leaf_node[2] <- node_info[which(tar_node==node_info[, 3]), 1]
      } else {
        leaf_node[2] <- clonal_edge[leaf_edge[2], 2]
      }

      node_mat[node_index, ] <- node_mat[leaf_node[1], ] | node_mat[leaf_node[2], ]
    } else if (node_info[i, 3]<0) {
      # recombination edge out node
      node_index <- node_info[i, 1]
      leaf_edge <- recomb_edge[abs(node_info[i, 3]), 1]
      tar_node <- which(recomb_node[[leaf_edge]]==node_info[i, 3])
      if (tar_node==1) {
        leaf_node <- clonal_edge[leaf_edge, 2]
      } else {
        leaf_node <- node_info[which(recomb_node[[leaf_edge]][tar_node-1]==node_info[, 3]), 1]
      }

      x <- recomb_edge[abs(node_info[i, 3]), 5]
      y <- recomb_edge[abs(node_info[i, 3]), 6]
      node_mat[node_index, ] <- FALSE
      recomb_node_mat[abs(node_info[i, 3]), x:y] <- node_mat[leaf_node, x:y]
      node_mat[node_index, -(x:y)] <- node_mat[leaf_node, -(x:y)]
    } else if (node_info[i, 3]>0) {
      # recombination edge in node
      node_index <- node_info[i, 1]
      leaf_edge <- recomb_edge[node_info[i, 3], 3]
      tar_node <- which(recomb_node[[leaf_edge]]==node_info[i, 3])
      if (tar_node==1) {
        if (leaf_edge==(2*n - 1)) {
          leaf_node <- 2*n - 1
        } else {
         leaf_node <- clonal_edge[leaf_edge, 2]
        }
      } else {
        leaf_node <- node_info[which(recomb_node[[leaf_edge]][tar_node-1]==node_info[, 3]), 1]
      }

      node_mat[node_index, ] <- node_mat[leaf_node, ] | recomb_node_mat[node_info[i, 3], ]
    }
  }

  ARG = list(clonal_edge=clonal_edge,
             recomb_edge=recomb_edge,
             node_info=node_info,
             node_mat=node_mat,
             recomb_node_mat=recomb_node_mat,
             sum_time=node_info[(2*n - 1 + 2*n_recomb), 2], n=n, rho=rho, L=L,
             delta=delta, clonal_time=t_sum)
  class(ARG) <- "ClonalOrigin"
  return(ARG)
}
