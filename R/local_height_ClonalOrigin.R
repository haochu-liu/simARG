#' Local tree height for `ClonalOrigin` object
#'
#' Find the height of a locat site at give site for a `ClonalOrigin` object.
#'
#' @param ARG A `ClonalOrigin` object.
#' @param location The site location for local tree.
#' @return Height of the local tree
#' @export
local_height_ClonalOrigin <- function(ARG, location) {
  if (!inherits(ARG, "ClonalOrigin")) {
    cli::cli_abort("Object must be of class 'ClonalOrigin'!")
  }

  if (is.null(ARG$recomb_edge)) {
    return(ARG$clonal_time)
  }
  edge_index <- which(ARG$recomb_edge[, 3] == -1 &
                      ARG$recomb_edge[, 5] <= location &
                      ARG$recomb_edge[, 6] >= location)
  if (!length(edge_index)) {
    return(ARG$clonal_time)
  } else {
    return(max(c(ARG$recomb_edge[edge_index, 4], ARG$clonal_time)))
  }
}
