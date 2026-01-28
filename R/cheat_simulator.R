#' Simulate approximated ARG using ClonalOrigin with pair method
#'
#' Simulate pairs of sites in the ClonalOrigin approximation from a given clonal tree.
#' Provide ARG height, length and number of mutations as cheated summary statistics.
#'
#' @param tree The clonal genealogy.
#' @param rho_site The recombination parameter per site.
#' @param theta_site The mutation rate per site.
#' @param L An integer for the number of sites.
#' @param delta Numeric, delta is the mean of recombinant segment length.
#' @param N Number of pairs to simulate.
#' @param k_vec A vector for the distance values between two sites. Default is 50, 200, 2000.
#' @return A 9 dimension vector as the summary statistics of simulations.
#' @export
#'
#' @examples
#' tree <- clonal_genealogy(15L)
#' summary_stats <- cheat_simulator(tree, 0.5, 0.3, 100L, 5,
#'                                  2000, c(10L, 20L, 70L))
cheat_simulator <- function(tree, rho_site, theta_site, L, delta,
                            N, k_vec=c(50L, 200L, 2000L)) {
  if (max(k_vec) > L) {
    cli::cli_abort("Site distance cannot be greater than the number of sites!")
  } else if (!rlang::is_integer(k_vec, n=3)) {
    cli::cli_abort("`k_vec` must be a vector of three integer values!")
  }

  s_vec <- rep(NA, 9)
  tree_width <- tree$n
  for (j in 1:3) {
    v_hei <- rep(NA, N)
    v_len <- rep(NA, N)
    v_mut <- rep(NA, N)
    for (i in 1:N) {
      ARG <- ClonalOrigin_pair_seq(tree, rho_site, L, delta, k_vec[j])
      ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)

      v_hei[i] <- ARG$sum_time
      v_len[i] <- sum(ARG$edge[, 3])
      v_mut[i] <- nrow(ARG_mutated$mutation)
    }
    s_vec[j] <- mean(v_hei)
    s_vec[j+3] <- mean(v_len)
    s_vec[j+6] <- mean(v_mut)
  }

  return(s_vec)
}
