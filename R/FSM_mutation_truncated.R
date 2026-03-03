#' Add mutations to a `FSM_ARG` object (truncated mutation distribution)
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
#' ARG_mutation <- FSM_mutation_truncated(ARG, 2)
FSM_mutation_truncated <- function(ARG, theta_site) {
  if (!inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'FSM_ARG'")
  }

  num_site <-  ncol(ARG$node_mat)
  n_vec <- rep(NA, num_site)
  mutate_edge <- c()
  mutate_site <- c()

  ARG$node_site <- matrix(FALSE, nrow=nrow(ARG$node_mat), ncol=ncol(ARG$node_mat))

  # Simulate mutations by sites
  for (i in 1:num_site) {
    # Length of local tree without reduction
    local_edge <- which(ARG$edge_mat[, i])
    local_length <- ARG$edge[local_edge, 3]
    # Truncated Poisson distribution
    local_n <- actuar::rztpois(n=1, lambda=theta_site*sum(local_length)/2)
    n_vec[i] <- local_n
    mutate_site <- c(mutate_site, rep(i, local_n))
    # Simulate edges
    mutate_edge <- c(mutate_edge, sample(local_edge, local_n,
                                         replace=TRUE, prob=local_length))
  }

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
