#' Homoplasy index
#'
#' Compute homoplasy index from ARG and sequence data.
#'
#' @param ARG A `FSM_ARG` object after mutation.
#' @return A numerical value, homoplasy index.
#' @export
#'
#' @examples
#' tree <- clonal_genealogy(10L)
#' ARG <- ClonalOrigin_pair_seq_fast(tree, 0.5, 100L, 5, 10L)
#' ARG_mutated <- FSM_mutation_fast(ARG, 0.3)
#' hi <- homoplasy_index(ARG_mutated)
homoplasy_index <- function(ARG) {
  if (!inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'FSM_ARG'")
  }

  n_node <- nrow(ARG$node_mat)
  n_leaf <- ARG$n
  n_site <- ncol(ARG$node_mat)
  m_vec <- rep(NA, n_site)
  s_vec <- rep(NA, n_site)

  for (site_loc in 1:n_site) {
    ARG_site <- local_tree_FSM(ARG, site_loc)
  }

}
