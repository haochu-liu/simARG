#' Input: number of leaf lineages, recombination parameter, number of sites
#' bacteria recombination or not,
#' if yes, delta is mean of the length of recombinant segment,
#' initial maximal node size (default = 1000),
#' optimise recombination edges or not, distinguish clonal lineages or not
#' Create a full ARG with coalescent and recombination
#' Edge dataframe: root node, leaf node, edge length, edge material interval
#' Node dataframe: node index, node height, node material interval
#' Output: edge dataframe, node dataframe, waiting time for each event,
#' total time, number of lineages at each event time, number of leaf alleles,
#' recombination parameter, bacteria recombination or not, and parameter delta
FSM_ARG <- function(n, rho, L, bacteria=FALSE, delta=NULL, node_max=1000,
                        optimise_recomb=FALSE, clonal=FALSE, edgemat=TRUE) {
  if (n!=as.integer(n)) {
    stop("Sample size must be an integer")
  } else if (L!=as.integer(L)) {
    stop("Number of sites must be an integer")
  } else if (bacteria & is.null(delta)) {
    stop("Must provide parameter for the length of recombinant segment")
  } else if (n >= node_max) {
    stop("Maximal node size must greater than the number of leaf lineages")
  } else if (clonal & (!bacteria)) {
    stop("Cannot consider clonal lineages for human gene")
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
  if (clonal) {                                    # node clonal
    node_clonal <- rep(NA, node_max)
  } else {
    node_clonal <- NULL
  }
  node_height[1:n] <- 0                            # initialize first n nodes
  node_mat[1:n, ] <- 1
  if (clonal) {node_clonal[1:n] <- TRUE}

  # Probability of starting recombination at each site
  probstart <- rep(1/L, L)
  probstartcum <- cumsum(probstart)

  # Initialize variables and vector
  edge_index <- 1L
  node_index <- as.integer(n + 1)
  pool <- as.integer(1:n)

  while (k > 1) {
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t <- c(t, event_time)
    t_sum <- t_sum + event_time
    # sample whether the event is a coalescent
    p_coale <- rbinom(n=1, size=1, prob=(k-1)/(k-1+rho))
    if (p_coale == 1) {
      # coalescent event
      leaf_node <- sample(pool, size=2, replace=FALSE)

      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- node_index
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]
      edge_mat_index[c(edge_index, edge_index+1)] <- leaf_node

      # append root node
      node_height[node_index] <- t_sum
      node_mat[node_index, ] <- as.integer(node_mat[leaf_node[1], ] |
                                             node_mat[leaf_node[2], ])

      # update clonal lineage
      if (clonal) {
        node_clonal[node_index] <- any(as.logical(node_clonal[leaf_node]))
      }

      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), node_index)
      edge_index <- edge_index + 2L
      node_index <- node_index + 1L
      k <- k - 1
    } else {
      # recombination event
      leaf_node <- sample(pool, size=1, replace=FALSE)

      if (bacteria) {
        x <- which(runif(1) < probstartcum)[1]
        y <- min(x + rgeom(1, 1/delta), L)

        if (clonal & optimise_recomb) {
          if ((sum(node_mat[leaf_node, x:y])==0 |
               sum(node_mat[leaf_node, -(x:y)])==0) &
              (!node_clonal[leaf_node])) {                   # not clonal
            next
          } else if (sum(node_mat[leaf_node, x:y])==0 &
                     node_clonal[leaf_node]) {               # clonal
            next
          }
        } else if (optimise_recomb) {
          if (sum(node_mat[leaf_node, x:y])==0 |
              sum(node_mat[leaf_node, -(x:y)])==0) {next}
        }

        edge_mat_index[c(edge_index, edge_index+1)] <- c(node_index, node_index+1)

        node_mat[c(node_index, node_index+1), ] <- 0
        node_mat[node_index, x:y] <- node_mat[leaf_node, x:y]
        node_mat[node_index+1, -(x:y)] <- node_mat[leaf_node, -(x:y)]
      } else {
        x <- which(runif(1) < probstartcum)[1]

        if ((sum(node_mat[leaf_node, 1:(x-1)])==0 |
             sum(node_mat[leaf_node, x:L])==0) & optimise_recomb) {
          next
        }

        edge_mat_index[c(edge_index, edge_index+1)] <- c(node_index, node_index+1)

        node_mat[c(node_index, node_index+1), ] <- 0
        node_mat[node_index, 1:(x-1)] <- node_mat[leaf_node, 1:(x-1)]
        node_mat[node_index+1, x:L] <- node_mat[leaf_node, x:L]
      }
      # append edges
      edge_matrix[c(edge_index, edge_index+1), 1] <- c(node_index, node_index+1L)
      edge_matrix[c(edge_index, edge_index+1), 2] <- leaf_node
      edge_matrix[c(edge_index, edge_index+1), 3] <- t_sum-node_height[leaf_node]

      # append root node
      node_height[c(node_index, node_index+1)] <- t_sum

      # update clonal lineage
      if (clonal) {
        node_clonal[node_index] <- FALSE
        node_clonal[node_index+1] <- node_clonal[leaf_node]
      }

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
      if (clonal) {node_clonal <- c(node_clonal, rep(NA, node_max))}
      node_max <- 2 * node_max
    }
  }

  if (edgemat) {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               edge_mat=node_mat[edge_mat_index[1:(edge_index-1)], ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               bacteria=bacteria, delta=delta)
  } else {
    ARG = list(edge=edge_matrix[1:(edge_index-1), ],
               node_height=node_height[1:(node_index-1)],
               node_mat=node_mat[1:(node_index-1), ],
               node_clonal=node_clonal[1:(node_index-1)],
               waiting_time=t, sum_time=t_sum, k=k_vector, n=n, rho=rho, L=L,
               bacteria=bacteria, delta=delta)
  }
  class(ARG) <- "FSM_ARG"
  return(ARG)
}
