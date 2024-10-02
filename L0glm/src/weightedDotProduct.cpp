#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat weightedDotProduct(const arma::mat& Y, const arma::mat& weights) {
  arma::uword n = Y.n_cols;  // Change type to arma::uword
  arma::mat result(n, n, fill::zeros);
  
#pragma omp parallel for shared(Y, weights, result, n) default(none)
  for (arma::uword i = 0; i < n; ++i) {  // Change loop variable type to arma::uword
    for (arma::uword j = i; j < n; ++j) {
      for (arma::uword k = 0; k < Y.n_rows; ++k) {  // Change loop variable type to arma::uword
        result(i, j) += weights(k, i) * weights(k, j) * Y(k, i) * Y(k, j);
      }
      if (i != j) {
        result(j, i) = result(i, j); // Mirror the result for symmetry
      }
    }
  }
  
  return result;
}
