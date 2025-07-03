#' Convert decimal rows to binary rows
#'
#' Convert every decimal number to `n` binary elements.
#'
#' @param x A integer vector.
#' @param n An integer.
#' @return A boolean vector of binary numbers.
decimal2binary <- function(x, n) {
  if (!rlang::is_integer(n)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_integer(x)) {
    cli::cli_abort("`x` must be an integer vector!")
  } else if (n > 30) {
    cli::cli_abort("`n` must be smaller than 31.")
  } else if (any(x > sum(2^(0:(n-1))))) {
    cli::cli_abort("Cannot convert to `n` digits, try larger value.")
  }

  binary_x <- rep(NA, length(x) * n)

  for (i in 1:length(x)) {
    start_index <- (i - 1) * n + 1
    binary_x[start_index:(start_index+n-1)] <- rev(as.logical(intToBits(x[i]))[1:n])
  }

  return(binary_x)
}
