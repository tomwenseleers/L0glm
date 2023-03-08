// Fast Approximate Solution Trajectory nnls (FAST-nnls), code by Zach de Bruine
// https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
typedef unsigned int uint;

// [[Rcpp::export]]
arma::vec fastnnls_dense(arma::mat A, arma::vec b) {

  // initial x is the unbounded least squares solution
  arma::vec x = arma::solve(A, b, arma::solve_opts::likely_sympd + arma::solve_opts::fast);

  while (any(x < 0)) {

    // define the feasible set as all values greater than 0
    arma::uvec nz = find(x > 0);

    // reset x
    x.zeros();

    // solve the least squares solution for values in the feasible set
    x.elem(nz) = solve(A.submat(nz, nz), b.elem(nz), arma::solve_opts::likely_sympd + arma::solve_opts::fast);
  }
  return x;
}

