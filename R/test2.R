if (n_recomb == 0) {
  node_mat <- matrix(TRUE, nrow=(2*n-1), ncol=k)
  edge_mat <- matrix(TRUE, nrow=2*(n-1), ncol=k)
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
node_mat <- matrix(NA, nrow=node_max, ncol=k)
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
