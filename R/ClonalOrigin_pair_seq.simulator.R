#' Simulate approximated ARG using ClonalOrigin with pair method
#'
#' Simulate pairs of sites in the ClonalOrigin approximation from a given clonal tree.
#' Provide summary statistics for pairs in three different distance.
#'
#' @param tree The clonal genealogy.
#' @param rho_site The recombination parameter per site.
#' @param theta_site The mutation rate per site.
#' @param L An integer for the number of sites.
#' @param delta Numeric, delta is the mean of recombinant segment length.
#' @param N Number of pairs to simulate.
#' @param k_vec A vector for the distance values between two sites. Default is 50, 200, 2000.
#' @return A 7 dimension vector as the summary statistics of simulations.
#' @export
#'
#' @examples
#' tree <- clonal_genealogy(15L)
#' summary_stats <- ClonalOrigin_pair_seq.simulator(tree, 0.5, 0.3, 100L, 5,
#'                                                  2000, c(10L, 20L, 70L))
ClonalOrigin_pair_seq.simulator <- function(tree, rho_site, theta_site,
                                            L, delta, N,
                                            k_vec=c(50L, 200L, 2000L)) {
  if (max(k_vec) > L) {
    cli::cli_abort("Site distance cannot be greater than the number of sites!")
  } else if (!rlang::is_integer(k_vec, n=3)) {
    cli::cli_abort("`k_vec` must be a vector of three integer values!")
  }

  s_vec <- rep(NA, 7)
  tree_width <- tree$n
  v_s <- rep(NA, N*3)
  for (j in 1:3) {
    v_r <- rep(NA, N)
    v_g3 <- rep(NA, N)
    for (i in 1:N) {
      ARG <- ClonalOrigin_pair_seq(tree, rho_site, L, delta, k_vec[j])
      ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
      mat <- ARG_mutated$node_site[1:tree_width, ]

      v_r[i] <- LD_r(mat)
      v_g3[i] <- G3_test(mat)
      v_s[i+(j-1)*N] <- any(as.logical(mat[, 1]))
    }
    s_vec[j] <- mean(v_r)
    s_vec[j+3] <- mean(v_g3)
  }
  s_vec[7] <- mean(v_s)

  return(s_vec)
}
