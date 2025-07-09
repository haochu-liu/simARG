#' An approximated ancestral recombination graph (ARG) using ClonalOrigin algorithm
#'
#' Simulate coalescent an recombination events by ClonalOrigin algorithm.
#' The non-clonal lineages are not allowed to either recombine or coalesce with each other.
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
#' ARG1 <- ClonalOrigin_ARG(100L, 5, 100L, 5)
#' ARG2 <- ClonalOrigin_ARG(5L, 1, 10L, 1, optimise_recomb=TRUE)
ClonalOrigin_ARG <- function(n, rho, L, delta, node_max=1000,
                             optimise_recomb=FALSE, edgemat=TRUE) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_integer(L, n=1)) {
    cli::cli_abort("`L` must be a single integer!")
  } else if (n >= node_max) {
    cli::cli_abort("Maximal node size must greater than the number of leaf lineages!")
  }

  k = n
  k_vector <- c(k)
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0

  edge_matrix <- matrix(NA, nrow=node_max, ncol=3) # root and leaf nodes, length
  colnames(edge_matrix) <- c("node1", "node2", "length")
  edge_mat_index <- rep(NA, node_max)              # edge material index
  node_height <- rep(NA, node_max)                 # node height to recent time
  node_mat <- matrix(NA, nrow=node_max, ncol=L)    # node material
  node_clonal <- rep(NA, node_max)                 # node clonal
  node_height[1:n] <- 0                            # initialize first n nodes
  node_mat[1:n, ] <- TRUE
  node_clonal[1:n] <- TRUE

  # Probability of starting recombination at each site
  probstart <- rep(1/L, L)
  probstartcum <- cumsum(probstart)

  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)

  while (k > 1) {
    # sample a new event time
    k_clonal <- sum(node_clonal[pool])
    k_nonclonal <- k - k_clonal
    coale_rate <- k*(k-1)/2 - k_nonclonal*(k_nonclonal-1)/2
    recomb_rate <- k_clonal*rho/2
    event_rate <- coale_rate + recomb_rate
    coale_prob <- coale_rate / event_rate
    event_time <- rexp(1, rate=event_rate)
    t <- c(t, event_time)
    t_sum <- t_sum + event_time
    # sample whether the event is a coalescent
    p_coale <- rbinom(n=1, size=1, prob=coale_prob)
    if (p_coale == 1) {
      # coalescent event
      repeat {
        leaf_node <- sample(pool, size=2, replace=FALSE)
        if (any(node_clonal[leaf_node])) {break}
      }

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- node_index
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
      edge_mat_index[c(edge_index, edge_index+1)] <- leaf_node

      # append root node
      node_height[node_index] <- t_sum
      node_mat[node_index, ] <- node_mat[leaf_node[1], ] | node_mat[leaf_node[2], ]

      # update clonal lineage
      node_clonal[node_index] <- TRUE

      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), node_index)
      edge_index <- edge_index + 2L
      node_index <- node_index + 1L
      k <- k - 1
    } else {
      # recombination event
      repeat {
        leaf_node <- sample(pool, size=1, replace=FALSE)
        if (node_clonal[leaf_node]) {break}
      }

      x <- which(runif(1) < probstartcum)[1]
      y <- min(x + rgeom(1, 1/delta), L)

      if (optimise_recomb & !any(node_mat[leaf_node, x:y])) {next}

      edge_mat_index[c(edge_index, edge_index+1)] <- c(node_index, node_index+1L)

      node_mat[c(node_index, node_index+1), ] <- FALSE
      node_mat[node_index, x:y] <- node_mat[leaf_node, x:y]
      node_mat[node_index+1, -(x:y)] <- node_mat[leaf_node, -(x:y)]

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- c(node_index, node_index+1L)
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]

      # append root node
      node_height[c(node_index, node_index+1)] <- t_sum

      # update clonal lineage
      node_clonal[node_index] <- FALSE
      node_clonal[node_index+1] <- TRUE

      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), node_index, node_index+1L)
      edge_index <- edge_index + 2L
      node_index <- node_index + 2L
      k <- k + 1
    }
    k_vector <- c(k_vector, k)
    if (max(edge_index, node_index) >= node_max - 1) {
      # add empty rows or elements if more edges than expected
      edge_matrix <- rbind(edge_matrix, matrix(NA, nrow=node_max, ncol=3))
      edge_mat_index <- c(edge_mat_index, rep(NA, node_max))
      node_height <- c(node_height, rep(NA, node_max))
      node_mat <- rbind(node_mat, matrix(NA, nrow=node_max, ncol=L))
      node_clonal <- c(node_clonal, rep(NA, node_max))
      node_max <- 2 * node_max
    }
  }

  if (edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  } else {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat_index=edge_mat_index[1:(edge_index-1)],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               sum_time=t_sum, n=n, rho=rho, L=L, delta=delta)
  }
  class(ARG) <- "FSM_ARG"
  return(ARG)
}
