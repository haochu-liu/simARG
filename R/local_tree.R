#' Construct local tree from ARG
#'
#' Call function `local_tree_ISM` or `local_tree_FSM`.
#' Select edges by a chosen site to create the local tree.
#'
#' @param ARG An `ISM_ARG` or `FSM_ARG` object.
#' @param location The site location for local tree.
#' @return `localtree` object; Local tree at a chosen site.
#' @export
#'
#' @examples
#' ARG <- FSM_ARG(20L, 1, 100L)
#' local_tree <- local_tree(ARG, 2L)
local_tree <- function(ARG, location) {
  if (!inherits(ARG, "ISM_ARG") & !inherits(ARG, "FSM_ARG")) {
    cli::cli_abort("Object must be of class 'ISM_ARG' or 'FSM_ARG'")
  } else if (inherits(ARG, "FSM_ARG") & (!rlang::is_integer(location, n=1))) {
    cli::cli_abort("`location` for a `FSM_ARG` object must be an integer!")
  } else if (inherits(ARG, "ISM_ARG") & (location < 0 | location >= 1)) {
    cli::cli_abort("`location` for a `ISM_ARG` object must be in [0, 1)!")
  }

  if (inherits(ARG, "ISM_ARG")) {
    local_tree <- local_tree_ISM(ARG, location)
  } else {
    local_tree <- local_tree_FSM(ARG, location)
  }

  local_tree$waiting_time <- NULL
  local_tree$k <- NULL
  local_tree$sum_time <- tree_height(local_tree)
  if (!is.null(local_tree$mutation)){
    local_tree$mutation[local_tree$mutation[, 1] %in% local_tree$edge_index, ]
  }

  return(local_tree)
}
