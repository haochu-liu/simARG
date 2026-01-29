#' Hitting time for a birth-death process
#'
#' Simulate a birth-death process and find its hitting time at boundary equals to 1.
#'
#' @param n A single integer for initial value.
#' @param rho The birth rate.
#' @return Hitting time and number of jumps as a vector.
#' @export
#'
#' @examples
#' t <- birth_death_sim2(20L, 1)
birth_death_sim2 <- function(n, rho) {
  if (!rlang::is_integer(n, n=1)) {
    cli::cli_abort("`n` must be a single integer!")
  }

  t <- 0
  k <- n
  num <- 0
  while (k > 1) {
    t <- t + rexp(1, rate=k*(k-1+rho)/2)
    if (runif(1) < (k - 1) / (k - 1 + rho)) {
      k <- k - 1
      num <- num + 1
    } else {
      k <- k + 1
      num <- num + 2
    }
  }
  return(c(t, num))
}
