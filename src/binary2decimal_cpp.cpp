#include <Rcpp.h>
using namespace Rcpp;


//' Convert binary rows to decimal rows (Rcpp version)
//'
//' For every `n` binary numbers, convert them into one decimal number.
//'
//' @param x A boolean vector.
//' @param n An integer that devides the length of `x`.
//' @param n_last An integer as the number of binary elements for the last decimal value.
//' @return A vector of decimal numbers.
// [[Rcpp::export]]
IntegerVector binary2decimal_cpp(LogicalVector x, int n, int n_last) {
  if (n > 30) {
    Rcpp::stop("`n` must be smaller than 31.");
  } else if (n_last > n) {
    Rcpp::stop("`n_last` cannot be greater than `n`!");
  }

  int decimal_x_length = (x.length()-n_last) / n + 1;
  IntegerVector decimal_x(decimal_x_length);

  for (int i = 0; i < decimal_x_length; ++i) {
    int j_max = n;
    if (i == decimal_x_length-1) {
      j_max = n_last;
    }
    for (int j = j_max - 1; j >= 0; --j) {
      decimal_x[i] += x[i*n + j] * pow(2, j_max-1-j);
    }
  }

  return decimal_x;
}
