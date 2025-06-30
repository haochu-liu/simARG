#' Input: n and rho
#' Simulate a birth death process until it hits the boundary
#' Output: the time when the process first hit k=1

#' Hitting time for a birth-death process
#'
#' Simulate a birth-death process and find its hitting time at boundary equals to 1.
#'
#' @param n A single integer for initial value.
#' @param rho The death rate.
#' @return Hitting time.
#' @export
#'
#' @examples
#' t <- birth_death_sim(20L, 1)
birth_death_sim <- function(n, rho) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  }

  t <- 0
  k <- n
  while (k > 1) {
    t <- t + rexp(1, rate=k*(k-1+rho)/2)
    if (runif(1) < (k - 1) / (k - 1 + rho)) {
      k <- k - 1
    } else {
      k <- k + 1
    }
  }
  return(t)
}
