#' Convert binary rows to decimal rows
#'
#' For every `n` binary numbers, convert them into one decimal number.
#'
#' @param x A boolean vector.
#' @param n An integer that devides the length of `x`.
#' @return A vector of decimal numbers.
binary_to_decimal <- function(x, n) {
  if (!rlang::is_integer(n)) {
    cli::cli_abort("`n` must be a single integer!")
  } else if (!rlang::is_logical(x)) {
    cli::cli_abort("`x` must be a boolean vector!")
  } else if (n > 30) {
    cli::cli_abort("`n` must be smaller than 31.")
  } else if (length(x) %% n != 0) {
    cli::cli_abort("`n` must devide the length of x.")
  }

  num_groups <- length(x) / n

  decimal_x <- numeric(num_groups)

  for (i in 1:num_groups) {
    start_index <- (i - 1) * n + 1
    current_bits <- rev(x[start_index:(start_index+n-1)])

    packed_byte <- packBits(c(current_bits, rep(F, 32-n)), type = "integer")
    decimal_x[i] <- as.integer(packed_byte)
  }

  return(decimal_x)
}
