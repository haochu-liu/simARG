#include <Rcpp.h>
using namespace Rcpp;


//' Convert decimal rows to binary rows (Rcpp version)
//'
//' Convert every decimal number to `n` binary elements.
//'
//' @param x A integer vector.
//' @param n An integer.
//' @param n_last An integer as the number of binary elements for last value in `x`.
//' @return A boolean vector of binary numbers.
// [[Rcpp::export]]
LogicalVector decimal2binary_cpp(IntegerVector x, int n, int n_last) {
  if (n > 30) {
    Rcpp::stop("`n` must be smaller than 31.");
  }

  int binary_x_length = (x.length()-1) * n + n_last;
  LogicalVector binary_x(binary_x_length);

  for (int i = 0; i < x.length(); ++i) {
    int j_max = n;
    if (i == x.length()-1) {
      j_max = n_last;
    }
    int target_x = x[i];
    for (int j = j_max - 1; j >= 0; --j) {
      binary_x[i*n + j] = (target_x % 2 == 1) ? TRUE : FALSE;
      target_x /= 2;
    }
  }

  return binary_x;
}
