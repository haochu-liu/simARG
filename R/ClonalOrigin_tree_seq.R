#' An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm
#'
#' Simulate coalescent an recombination events by ClonalOrigin algorithm.
#' Recombinations are added to a clonal tree to construct an ARG.
#'
#' @param n An integer for the number of leaf lineages.
#' @param rho_site The recombination parameter per site.
#' @param L An integer for the number of sites.
#' @param delta Numeric, delta is the mean of recombinant segment length.
#' @return A list containing clonal tree and recombination edges.
#' @export
#'
#' @examples
#' ARG1 <- ClonalOrigin_tree_seq(100L, 0.5, 100L, 5)
#' ARG2 <- ClonalOrigin_tree_seq(5L, 1, 10L, 1)
ClonalOrigin_tree_seq <- function(n, rho_site, L, delta) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_integer(L, n=1)) {
    cli::cli_abort("`L` must be a single integer!")
  }

  k = n
  t_sum <- 0
  rho <- L * rho_site

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

  # Initialize recombination edges
  nrow_max <- 1000
  recomb_edge <- matrix(NA, nrow=nrow_max, ncol=6) # matrix for b, a, x, y
  colnames(recomb_edge) <- c("b_edge", "b_height",
                             "a_edge", "a_height",
                             "x", "y")

  # Add recombination sequentially
  n_recomb <- 0
  for (i in 1:L) {
    if (i == 1) {
      R_new <- rpois(1, rho_site*delta*tree_length/2)
      R_old <- 0
    } else {
      survive_index <- which(recomb_edge[1:n_recomb, 6] == (i-1))
      R_new <- rpois(1, rho_site*tree_length/2)
      if (length(survive_index) >= 0) {
        R_old <- rbinom(1, length(survive_index), (1 - 1/delta))
        if (length(survive_index) == 1) {
          remain_index <- survive_index
        } else {
          remain_index <- sample(survive_index, R_old)
        }
      } else {
        R_old <- 0
      }
    }

    if (R_new > 0) {
      if (n_recomb+R_new >= nrow_max) {
        recomb_edge <- rbind(recomb_edge, matrix(NA, nrow=node_max, ncol=6))
        nrow_max <- 2 * nrow_max
      }

      recomb_edge[(n_recomb+1):(n_recomb+R_new), c(5, 6)] <- i
      a_rexp <- rexp(R_new, rate=1)
      # simulate b_edge (similar to mutation)
      recomb_edge[(n_recomb+1):(n_recomb+R_new), 1] <- sample(1:(2*(n-1)), R_new,
                                                              replace=TRUE, prob=clonal_edge[, 3])
      for (j in 1:R_new) {
        # simulate b_height
        recomb_edge[n_recomb+j, 2] <- runif(1, max=clonal_edge[recomb_edge[n_recomb+j, 1], 3]) +
          clonal_node_height[clonal_edge[recomb_edge[n_recomb+j, 1], 2]]
        # identify a_height
        t_above_b <- clonal_node_height[n:(2*n-1)] - recomb_edge[n_recomb+j, 2]
        i_above_b <- c(0, t_above_b[t_above_b >= 0])
        i_above_b <- i_above_b[2:length(i_above_b)] - i_above_b[1:(length(i_above_b)-1)]
        cuml_above_b <- cumsum(i_above_b * (1+length(i_above_b)):2)
        num_lineage <- (1+length(i_above_b)) - length(which(a_rexp[j] > cuml_above_b))
        if (num_lineage == (1+length(i_above_b))) {
          recomb_edge[n_recomb+j, 4] <- a_rexp[j] / num_lineage + recomb_edge[n_recomb+j, 2]
        } else {
          recomb_edge[n_recomb+j, 4] <- (a_rexp[j]-cuml_above_b[1+length(i_above_b)-num_lineage]) / num_lineage +
            sum(i_above_b[1:(1+length(i_above_b)-num_lineage)]) +
            recomb_edge[n_recomb+j, 2]
        }
        # simulate a_edge
        if (num_lineage > 1) {
          pool_edge <- which((clonal_node_height[clonal_edge[, 1]] >= recomb_edge[n_recomb+j, 4]) &
                               (clonal_node_height[clonal_edge[, 2]] < recomb_edge[n_recomb+j, 4]))
          recomb_edge[n_recomb+j, 3] <- sample(pool_edge, 1, replace=TRUE)
        } else {
          recomb_edge[n_recomb+j, 3] <- 2*n - 1
        }
      }
    }

    if (R_old > 0) {
      recomb_edge[remain_index, 6] <- i
    }

    n_recomb <- n_recomb + R_new
  }
  if (n_recomb == 0) {
    node_mat <- matrix(TRUE, nrow=(2*n-1), ncol=L)
    edge_mat <- matrix(TRUE, nrow=2*(n-1), ncol=L)
    ARG = list(edge=clonal_edge,
               edge_mat=edge_mat,
               node_height=clonal_node_height,
               node_mat=node_mat,
               node_clonal=rep(TRUE, (2*n-1)),
               sum_time=max(clonal_node_height),
               n=n, rho=rho, L=L, delta=delta)
    class(ARG) <- "FSM_ARG"
    return(ARG)
  } else {
    recomb_edge <- recomb_edge[1:n_recomb, , drop = FALSE]
  }

  # recombination segment and ancestral material
  node_max <- 2*n - 1 + 3*n_recomb
  edge_max <- 2*(n - 1) + 4*n_recomb
  edge_matrix <- matrix(NA, nrow=edge_max, ncol=3) # root and leaf nodes, length
  colnames(edge_matrix) <- c("node1", "node2", "length")
  edge_mat_index <- rep(NA, edge_max)              # edge material index
  node_mat <- matrix(NA, nrow=node_max, ncol=L)
  node_info <- matrix(NA, nrow=node_max, ncol=4)
  colnames(node_info) <- c("index", "height", "recomb", "clonal")
  node_mat[1:n, ] <- TRUE
  node_info[, 1] <- 1:node_max

  node_info[1:(2*n-1), 2] <- clonal_node_height
  node_info[1:(2*n-1), 3] <- 0
  node_info[1:(2*n-1), 4] <- TRUE

  node_info[(2*n):(2*n-1+2*n_recomb), 2] <- c(mapply(c, recomb_edge[, 2], recomb_edge[, 2]))
  node_info[(2*n):(2*n-1+2*n_recomb), 3] <- c(mapply(c, (-c(1:n_recomb)), rep(NA, n_recomb)))
  node_info[(2*n):(2*n-1+2*n_recomb), 4] <- c(mapply(c, rep(T, n_recomb), rep(F, n_recomb)))

  # node_info[(2*n+n_recomb):(2*n-1+2*n_recomb), 2] <- recomb_edge[, 2]
  # node_info[(2*n+n_recomb):(2*n-1+2*n_recomb), 3] <- NA
  # node_info[(2*n+n_recomb):(2*n-1+2*n_recomb), 4] <- FALSE

  node_info[(2*n+2*n_recomb):node_max, 2] <- recomb_edge[, 4]
  node_info[(2*n+2*n_recomb):node_max, 3] <- 1:n_recomb
  node_info[(2*n+2*n_recomb):node_max, 4] <- TRUE

  node_info <- node_info[order(node_info[, 2]), ]
  # recombination nodes on every edge
  recomb_node <- lapply(1:(2*n - 1), function(n){
    ClonalOrigin_nodes(recomb_edge, n)
  })
  # Add ancestral material to every node and construct full ARG edges
  i <- n + 1
  edge_index <- 1L
  repeat {
    if (node_info[i, 3]==0) {
      # clonal tree
      node_index <- node_info[i, 1]
      leaf_edge <- which(clonal_edge[, 1] == node_index)
      leaf_index <- rep(NA, 2)
      leaf_node <- rep(NA, 2)
      if (length(recomb_node[[leaf_edge[1]]])) {
        # target node is tail(recomb_node[[leaf_edge[1]]], 1)
        tar_node <- tail(recomb_node[[leaf_edge[1]]], 1)
        leaf_index[1] <- which(tar_node==node_info[, 3])
        leaf_node[1] <- node_info[leaf_index[1], 1]
      } else {
        leaf_node[1] <- clonal_edge[leaf_edge[1], 2]
        leaf_index[1] <- which(leaf_node[1]==node_info[, 1])
      }
      if (length(recomb_node[[leaf_edge[2]]])) {
        # target node is tail(recomb_node[[leaf_edge[2]]], 1)
        tar_node <- tail(recomb_node[[leaf_edge[2]]], 1)
        leaf_index[2] <- which(tar_node==node_info[, 3])
        leaf_node[2] <- node_info[leaf_index[2], 1]
      } else {
        leaf_node[2] <- clonal_edge[leaf_edge[2], 2]
        leaf_index[2] <- which(leaf_node[2]==node_info[, 1])
      }

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- i
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_index
      edge_matrix[c(edge_index, edge_index+1), 3] <- node_info[i, 2] - node_info[leaf_index, 2]
      edge_mat_index[c(edge_index, edge_index+1)] <- leaf_index

      # append root node
      node_mat[i, ] <- node_mat[leaf_index[1], ] | node_mat[leaf_index[2], ]

      edge_index <- edge_index + 2L
      i <- i + 1
    } else if (node_info[i, 3]<0) {
      # recombination edge out node
      node_index <- node_info[c(i, i+1L), 1]
      leaf_edge <- recomb_edge[abs(node_info[i, 3]), 1]
      tar_node <- which(recomb_node[[leaf_edge]]==node_info[i, 3])
      if (tar_node==1) {
        leaf_node <- clonal_edge[leaf_edge, 2]
      } else {
        leaf_node <- node_info[which(recomb_node[[leaf_edge]][tar_node-1]==node_info[, 3]), 1]
      }
      leaf_index <- which(leaf_node==node_info[, 1])

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- c(i, i+1)
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_index
      edge_matrix[c(edge_index, edge_index+1), 3] <- node_info[i, 2] - node_info[leaf_index, 2]
      edge_mat_index[c(edge_index, edge_index+1)] <- c(i, i+1)

      x <- recomb_edge[abs(node_info[i, 3]), 5]
      y <- recomb_edge[abs(node_info[i, 3]), 6]

      # append root node
      node_mat[c(i, i+1), ] <- FALSE
      node_mat[i+1, x:y] <- node_mat[leaf_index, x:y]
      node_mat[i, -(x:y)] <- node_mat[leaf_index, -(x:y)]

      edge_index <- edge_index + 2L
      i <- i + 2
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
      leaf_index <- rep(NA, 2)
      leaf_index[1] <- which(leaf_node==node_info[, 1])
      leaf_index[2] <- which(node_info[, 3]==(-node_info[i, 3])) + 1

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- i
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_index
      edge_matrix[c(edge_index, edge_index+1), 3] <- node_info[i, 2] - node_info[leaf_index, 2]
      edge_mat_index[c(edge_index, edge_index+1)] <- leaf_index

      # append root node
      node_mat[i, ] <- node_mat[leaf_index[1], ] | node_mat[leaf_index[2], ]

      edge_index <- edge_index + 2L
      i <- i + 1
    }

    if (i > node_max) {break}
  }

  ARG = list(edge=edge_matrix,
             edge_mat=node_mat[edge_mat_index, ],
             node_height=node_info[, 2],
             node_mat=node_mat,
             node_clonal=node_info[, 4],
             sum_time=node_info[node_max, 2],
             n=n, rho=rho, L=L, delta=delta)
  class(ARG) <- "FSM_ARG"
  return(ARG)
}
