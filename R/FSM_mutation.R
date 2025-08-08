#' Add mutations to a `FSM_ARG` object
#'
#' Add mutations uniformly onto the edges of ARG with infinite site assumption.
#'
#' @param ARG A `FSM_ARG` object.
#' @param theta_site The mutation rate per site.
#' @param binary boolean; If TRUE, return the incidence matrix of nodes, instead of sequences.
#' @return Add a mutation matrix and a genotype dataframe to the `FSM_ARG` object.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' ARG_mutation <- FSM_mutation(ARG, 2)
FSM_mutation <- function(ARG, theta_site, binary=FALSE) {
  if (!inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'FSM_ARG'")
  }

  theta <- theta_site * ncol(ARG$node_mat)
  l <- sum(ARG$edge[, 3])
  n <- rpois(1, theta*l/2) # num of mutations | l ~ Poisson(theta*l/2)
  ARG$mutation <- matrix(NA, nrow=n, ncol=3)
  colnames(ARG$mutation) <- c("edge_index", "pos", "site")
  if (binary) {
    ARG$node_site <- matrix(FALSE, nrow=nrow(ARG$node_mat), ncol=ncol(ARG$node_mat))
  }
  ARG$node <- data.frame(gene_str=rep("[]", length(ARG$node_height)))
  ARG$node$gene <- list(numeric())

  # if there is no mutation
  if (n == 0) {return(ARG)}

  # if there are mutations
  mutate_edge <- sample(1:nrow(ARG$edge), n,
                        replace=TRUE, prob=ARG$edge[, 3])
  mutate_site <- sample(1:ncol(ARG$node_mat), n, replace=TRUE)
  keep_mutation <- c()

  # ignore mutations not in the edge material
  for (i in 1:n) {
    if (ARG$edge_mat[mutate_edge[i], mutate_site[i]]) {
      keep_mutation <- c(keep_mutation, i)
    }
  }
  mutate_edge <- mutate_edge[keep_mutation]
  mutate_site <- mutate_site[keep_mutation]
  ARG$mutation <- ARG$mutation[keep_mutation, ]
  if (!length(keep_mutation)) {return(ARG)}

  # dataframe mutations to store information of mutations
  ARG$mutation[, 1] <- mutate_edge
  ARG$mutation[, 3] <- mutate_site
  for (i in 1:length(mutate_edge)) {
    ARG$mutation[i, 2] <- runif(1, max=ARG$edge[mutate_edge[i], 3])
  }

  # simulate the mutations at every node
  for (i in nrow(ARG$edge):1) {
    edge_mutation <- as.integer(ARG$mutation[ARG$mutation[, 1]==i, 3])
    parent_seq <- ARG$node$gene[[ARG$edge[i, 1]]]
    parent_seq <- c(parent_seq, edge_mutation)
    for (j in 1:ncol(ARG$node_mat)) {
      if (ARG$edge_mat[i, j] == 0) {
        parent_seq <- parent_seq[parent_seq != j]
      }
    }
    ARG$node$gene[[ARG$edge[i, 2]]] <- sort(c(parent_seq,
                                              ARG$node$gene[[ARG$edge[i, 2]]]))
  }

  if (binary) {
    # convert to incidence matrix
    for (i in 1:nrow(ARG$node_mat)) {
      for (j in 1:ncol(ARG$node_mat)) {
        num_mutation <- sum(ARG$node$gene[[i]] == j)
        if (num_mutation %% 2) {
          ARG$node_site[i, j] <- TRUE
        }
      }
    }
  } else {
    # convert to string
    for (i in 1:nrow(ARG$node_mat)) {
      ARG$node$gene_str[i] <- paste0("[", paste(round(ARG$node$gene[[i]], 3),
                                                collapse = ", "), "]")
    }
  }

  return(ARG)
}
