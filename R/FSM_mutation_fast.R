#' Add mutations to a `FSM_ARG` object
#'
#' Add mutations uniformly onto the edges of ARG with infinite site assumption.
#' Return only the incidence matrix of nodes.
#'
#' No mutation matrix and only keep `mutate_edge` and `mutate_site`.
#' No `node` dataframe, simulate `node_site` through the tree.
#'
#' @param ARG A `FSM_ARG` object.
#' @param theta_site The mutation rate per site.
#' @return Add a site matrix for every node to the `FSM_ARG` object.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' ARG_mutation <- FSM_mutation_fast(ARG, 2)
FSM_mutation_fast <- function(ARG, theta_site) {
  if (!inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'FSM_ARG'")
  }

  theta <- theta_site * ncol(ARG$node_mat)
  l <- sum(ARG$edge[, 3])
  n <- rpois(1, theta*l/2) # num of mutations | l ~ Poisson(theta*l/2)

  ARG$node_site <- matrix(FALSE, nrow=nrow(ARG$node_mat), ncol=ncol(ARG$node_mat))

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
  if (!length(keep_mutation)) {return(ARG)}

  # simulate the mutations at every node
  for (i in nrow(ARG$edge):1) {
    edge_mutation <- as.integer(mutate_site[mutate_edge==i])
    parent_seq <- ARG$node_site[ARG$edge[i, 1], ]
    flip_counts <- tabulate(edge_mutation, nbins = length(parent_seq))
    parent_seq <- xor(parent_seq, flip_counts %% 2 == 1)
    material_range <- ARG$edge_mat[i, ]
    ARG$node_site[ARG$edge[i, 2], material_range] <- parent_seq[material_range]
  }

  return(ARG)
}
