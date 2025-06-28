#' Plot a tree-like ARG
#'
#' `plot.ARG` is able to plot a tree-like network for an ARG.
#'
#' @param ARG An `ISM_ARG` or `FSM_ARG` object.
#' @return A plot from `igraph::plot.igraph`.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' plot.ARG(ARG)
plot.ARG <- function(ARG) {
  if (!inherits(ARG, "ISM_ARG") & !inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'ISM_ARG' or 'FSM_ARG'")
  }

  ARG_matrix <- as.matrix(ARG$edge[, c(1, 2)])
  g <- igraph::graph_from_edgelist(ARG_matrix, directed = FALSE)
  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
  layout_coord <- ARG_igraph_layout(ARG)
  igraph::plot.igraph(g, layout=layout_coord)
}
