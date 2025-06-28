#' Compute the effective recombination rate for `simbac_ARG`
#'
#' Compute the effective recombination rate,
#' and provide the probability of recombination initiating sites.
#'
#' @param mat A vector of node ancestral material.
#' @param delta numeric; The mean of recombinant segment length.
#' @param rho The recombination parameter.
#' @param clonal logical; whether the lineage is clonal.
#' @param include_site A vector of optimised sites.
#' @return Effective recombination rate, cumulative sum of initiating sites probability.
effective_R <- function(mat, delta, rho, clonal, include_site=NULL) {
  if (!is.null(include_site)) {
    mat <- mat[include_site]
  }

  L <- length(mat)
  if (L < 2) {
    if (!is.null(include_site)) {
      return(list(R_eff = 0,
                  probstartcum=rep(0, length(include_site))))
    }
    return(list(R_eff = 0,
                probstartcum=rep(0, L)))
  }

  R <- rho / L
  v_s <- which(mat & (mat != c(0, mat[1:(L-1)]))) # compute s1, ..., sb
  v_e <- which(mat & (mat != c(mat[2:L], 0)))     # compute e1, ..., eb
  v_e <- c(v_e[length(v_e)] - L, v_e)             # add e0
  b <- length(v_s)                                # number of blocks

  # compute R_eff
  R_gap <- R * delta * (1 - (1 - 1/delta)^(v_s - v_e[1:b]))
  R_eff <- sum(R_gap) + R * (sum(mat) - b)
  if (!clonal) {
    R_eff <- R_eff - sum(R_gap * (1 - 1/delta)^(L - (v_s - v_e[1:b]))) -
      R * (1 - 1/delta)^(L - 1) * (sum(mat) - b)
  }

  # generate probstart
  probstart <- rep(0, L)
  if (clonal) {
    probstart[as.logical(mat)] <- R
    probstart[v_s] <- R_gap
  } else {
    probstart[as.logical(mat)] <- R * (1 - (1 - 1/delta)^(L - 1))
    probstart[v_s] <- R_gap * (1 - (1 - 1/delta)^(L - (v_s - v_e[1:b])))
  }

  if (sum(probstart) != 0) {probstart <- probstart / sum(probstart)}
  if (!is.null(include_site)) {
    full_probstart <- rep(0, length(include_site))
    full_probstart[include_site] <- probstart
    probstart <- full_probstart
  }

  return(list(R_eff = R_eff,
              probstartcum=cumsum(probstart)))
}
