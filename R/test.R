if (!inherits(tree, "clonal_tree")) {
  cli::cli_abort("Object must be of class 'clonal_tree'")
} else if (!rlang::is_integer(L, n=1)) {
  cli::cli_abort("`L` must be a single integer!")
} else if (!rlang::is_integer(k, n=1)) {
  cli::cli_abort("`k` must be a single integer!")
} else if (delta <= 0) {
  cli::cli_abort("`delta` must be greater than zero!")
}

rho <- L * rho_site
clonal_edge <- tree$edge
clonal_node_height <- tree$node
n <- tree$n

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
for (i in 1:2) {
  if (i == 1) {
    R_new <- rpois(1, rho_site*delta*tree_length/2)
    R_old <- 0
  } else { # i == 2
    survive_index <- which(recomb_edge[1:n_recomb, 6] == 1)
    delta2 <- sum((1 - 1/delta)^c(0:(k-1)))
    R_new <- rpois(1, rho_site*delta2*tree_length/2)
    R_old <- rbinom(1, length(survive_index), (1 - 1/delta)^k)
    if (length(survive_index) == 1) {
      remain_index <- survive_index
    } else {
      remain_index <- sample(survive_index, R_old)
    }
  }

  if (R_new > 0) {
    if (n_recomb+R_new >= nrow_max) {
      recomb_edge <- rbind(recomb_edge, matrix(NA, nrow=(nrow_max+n_recomb+R_new), ncol=6))
      nrow_max <- 2 * nrow_max + n_recomb + R_new
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
      num_lineage <- (1+length(i_above_b)) - sum(a_rexp[j] > cuml_above_b)
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
                             (clonal_node_height[clonal_edge[, 2]] <  recomb_edge[n_recomb+j, 4]))
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
  node_mat <- matrix(TRUE, nrow=(2*n-1), ncol=2)
  edge_mat <- matrix(TRUE, nrow=2*(n-1), ncol=2)
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
