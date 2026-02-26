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
  m_vec <- rep(0, n_site)
  s_vec <- rep(NA, n_site)

  for (site_loc in 1:n_site) {
    ARG_site <- local_tree_FSM(ARG, site_loc)
    node_vec <- sort(unique(as.vector(ARG_site$edge[, 1:2])))
    node_site_vec <- ARG_site$node_site[node_vec, site_loc]

    # Compute minimum possible changes
    if (any(node_site_vec[1:n_leaf])) {m_vec[site_loc] <- 1}

    # Compute actual changes
    s_site <- 0
    site_list <- setNames(rep(list(NA), length(node_vec)), node_vec)
    site_list[1:n_leaf] <- as.integer(node_site_vec[1:n_leaf])
    for (i in (n_leaf+1):length(node_vec)) {
      parent_node <- node_vec[i]
      node_index <- which(ARG_site$edge[, 1] == parent_node)
      if (length(node_index) == 2) {
        # Coalescent structure
        children_node <- ARG_site$edge[node_index, 2]
        child_1 <- site_list[[as.character(children_node[1])]]
        child_2 <- site_list[[as.character(children_node[2])]]
        intersec <- intersect(child_1, child_2)
        if (length(intersec) == 0) {
          # intersection is empty -> mutation
          site_list[[as.character(parent_node)]] <- union(child_1, child_2)
          s_site <- s_site + 1
        } else {
          # intersection is not empty -> no mutation
          site_list[[as.character(parent_node)]] <- intersec
        }
      } else if (length(node_index) == 1) {
        # Recombination structure
        children_node <- ARG_site$edge[node_index, 2]
        site_list[[as.character(parent_node)]] <- site_list[[as.character(children_node)]]
      }
    }
    s_vec[site_loc] <- s_site
  }

  hi <- 1 - sum(m_vec) / sum(s_vec)
  if ((sum(m_vec) == 0) & (sum(s_vec) == 0)) {hi <- 0}
  return(hi)
}
