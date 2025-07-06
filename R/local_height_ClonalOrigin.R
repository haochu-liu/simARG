#' Local tree height for `ClonalOrigin` object
#'
#' Find the height of a locat site at give site for a `ClonalOrigin` object.
#'
#' @param ARG A `ClonalOrigin` object.
#' @param location The site location for local tree.
#' @return Height of the local tree
local_height_ClonalOrigin <- function(ARG, location) {
  if (!inherits(ARG, "ClonalOrigin")) {
    cli::cli_abort("Object must be of class 'ClonalOrigin'!")
  }

  if (is.na(ARG$recomb_edge)) {
    return(ARG$sum_time)
  }
  edge_index <- which(recomb_edge[, 3] == -1 &
                      recomb_edge[, 5] <= location &
                      recomb_edge[, 6] >= location)
  if (length(edge_index)) {
    return(ARG$sum_time)
  } else {
    return(max(recomb_edge[edge_index, 4]))
  }
}
