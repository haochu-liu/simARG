#' Homoplasy index
#'
#' Compute homoplasy index from ARG and sequence data.
#'
#' @param ARG A `FSM_ARG` object after mutation.
#' @return A numerical value, homoplasy index.
#' @export
#'
#' @examples
#' tree <- clonal_genealogy(15L)
#' ARG <- ClonalOrigin_pair_seq_fast(tree, 0.5, 100L, 5, c(10L, 20L, 70L))
#' ARG_mutated <- FSM_mutation_fast(ARG, 0.3)
#' hi <- homoplasy_index(ARG_mutated)
homoplasy_index <- function(ARG) {

}
