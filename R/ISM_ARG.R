#' Create a full ancestral recombination graph (ARG)
#'
#' Simulate coalescent an recombination events and construct an ARG.
#' Assume infinite site model.
#'
#' @param n An integer for the number of leaf lineages.
#' @param rho The recombination parameter.
#' @param bacteria logical; If TRUE, genes recombine by conversion.
#' @param delta numeric; If bacteria = TRUE, delta is the mean of recombinant segment length.
#' @return A list containing dataframes of edges and nodes, and other information about ARG.
#' @export
#'
#' @examples
#' ARG1 <- ISM_ARG(5L, 1)
#' ARG2 <- ISM_ARG(10L, 0.5, bacteria=TRUE, delta=0.1)
ISM_ARG <- function(n, rho, bacteria=FALSE, delta=NULL) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (bacteria & is.null(delta)) {
    cli::cli_abort("Must provide parameter delta for gene conversion!")
  }

  k = n
  k_vector <- c(k)
  t <- vector("numeric", length = 0) # vector of event times
  t_sum <- 0
  edge <- tibble::tibble(
    node1 = integer(), # root node
    node2 = integer(), # leaf node
    length = numeric(), # edge length
    material = list() # edge material interval
  )
  node <- tibble::tibble(
    height = rep(0, n), # node height to recent time
    material = list(ivs::iv(0, 1)) # node material interval
  )
  pool <- as.integer(1:n)
  next_node <- as.integer(n+1)

  while (k > 1) {
    # sample a new event time
    event_time <- rexp(1, rate=k*(k-1+rho)/2)
    t <- c(t, event_time)
    t_sum <- t_sum + event_time
    # sample whether the event is a coalescent
    p_coale <- rbinom(n=1, size=1, prob=(k-1)/(k-1+rho))
    if (p_coale == 1) {
      # coalescent event
      leaf_node <- sample(pool, size=2, replace=FALSE)
      # append edges
      append_edge <- tibble::tibble(
        node1 = rep(next_node, 2),
        node2 = leaf_node,
        length = c(t_sum-node$height[leaf_node[1]],
                   t_sum-node$height[leaf_node[2]]),
        material = list(node$material[[leaf_node[1]]],
                        node$material[[leaf_node[2]]])
      )
      edge <- dplyr::bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble::tibble(
        height = t_sum,
        material = list(ivs::iv_set_union(append_edge$material[[1]],
                                          append_edge$material[[2]]))
      )
      node <- dplyr::bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node)
      next_node <- next_node + 1L
      k <- k - 1
    } else {
      # recombination event
      leaf_node <- sample(pool, size=1, replace=FALSE)
      if (bacteria) {
        x <- runif(1, min=0, max=1)
        y <- min(1, x + rexp(1, rate=1/delta))
        iv1 <- ivs::iv(x, y)
        iv2 <- ivs::iv_set_complement(iv1, lower = 0, upper = 1)
        df_mat <- list(ivs::iv_set_intersect(iv1, node$material[[leaf_node]]),
                       ivs::iv_set_intersect(iv2, node$material[[leaf_node]]))
      } else {
        u <- runif(1, min=0, max=1)
        df_mat <- list(ivs::iv_set_intersect(ivs::iv(0, u), node$material[[leaf_node]]),
                       ivs::iv_set_intersect(ivs::iv(u, 1), node$material[[leaf_node]]))
      }
      # append edges
      append_edge <- tibble::tibble(
        node1 = c(next_node, next_node + 1L),
        node2 = rep(leaf_node, 2),
        length = c(t_sum-node$height[leaf_node],
                   t_sum-node$height[leaf_node]),
        material = df_mat
      )
      edge <- dplyr::bind_rows(edge, append_edge)
      # append root node
      append_node <- tibble::tibble(
        height = t_sum,
        material = df_mat
      )
      node <- dplyr::bind_rows(node, append_node)
      # updates for iteration
      pool <- c(setdiff(pool, leaf_node), next_node, next_node+1L)
      next_node <- next_node + 2L
      k <- k + 1
    }
    k_vector <- c(k_vector, k)
  }
  ARG = list(edge=edge, node=node, waiting_time=t, sum_time=t_sum, k=k_vector,
             n=n, rho=rho, bacteria=bacteria, delta=delta)
  class(ARG) <- "ISM_ARG"
  return(ARG)
}
