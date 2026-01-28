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

}

