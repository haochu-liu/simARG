#' Select edges for local tree graph
#'
#' Construct a local tree from an `ISM_ARG` object with given site.
#'
#' @param ARG An `ISM_ARG` object.
#' @param location The site location for local tree.
#' @return `localtree` object; Local tree at a chosen site.
local_tree_ISM <- function(ARG, location) {
  interval <- iv(location, location+.Machine$double.eps)
  keep_edge <- c()
  for (i in 1:nrow(ARG$edge)) {
    if (iv_count_overlaps(interval, ARG$edge$material[[i]])) {
      keep_edge <- c(keep_edge, i)
    }
  }
  local_tree <- ARG
  local_tree$edge <- ARG$edge[keep_edge, ]

  class(local_tree) <- "localtree"
  return(local_tree)
}
