// Fast Approximate Solution Trajectory nnls (FAST-nnls), code by Zach de Bruine
// https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
typedef unsigned int uint;

// [[Rcpp::export]]
arma::vec fastnnls_dense(arma::mat &A, arma::vec b) {

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


// [[Rcpp::export]]
arma::vec fastnnls_sparse(arma::sp_mat &A, arma::vec b) {
  
  // initial x is the unbounded least squares solution
  arma::vec x = arma::spsolve(A, b, "superlu");

  while (any(x < 0)) {

    // define the feasible set as all values greater than 0
    arma::uvec nz = find(x > 0);

    // reset x
    x.zeros();
    
    // Compute the minimum and maximum indices from the 'nz' vector
    arma::uword row_start = nz.min();
    arma::uword row_end = nz.max();

    // Extract the submatrix using the computed indices and all columns
    arma::sp_mat A_sub = A.submat(row_start, 0, row_end, A.n_cols - 1);
    arma::sp_mat ATA_sub = A_sub.t() * A_sub;
    arma::vec ATb_sub = A_sub.t() * b.elem(nz);
    x.elem(nz) = spsolve(ATA_sub, ATb_sub);  // try , "superlu" to use spsolve with "superlu" option, but requires installation of superlu library & it being correctly linked
  }
  return x;
}
