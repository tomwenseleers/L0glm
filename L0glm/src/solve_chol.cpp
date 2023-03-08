// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// Solve Ax = b using Cholesky decomposition for sparse or dense covariate matrix A
// adapted from https://github.com/nk027/sanic/blob/master/src/solve_Cholesky.cpp

// Solve Ax = b using sparse LLT (pivot = 0) or LDLT (pivot = 1) Cholesky decomposition
// [[Rcpp::export]]
Rcpp::List solve_sparse_chol(
    const Eigen::MappedSparseMatrix<double> A,
    const Eigen::Map<Eigen::MatrixXd> b,
    unsigned int pivot = 1, unsigned int ordering = 0) {

  Eigen::SimplicialLDLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > solver;
  if(ordering == 1) { // use NaturalOrdering
    Eigen::SimplicialLDLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > solver;
  } else if(ordering > 1) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }

  if(pivot == 0) {
    Eigen::SimplicialLLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > solver;
    if(ordering == 1) { // use NaturalOrdering
      Eigen::SimplicialLLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > solver;
    } else if(ordering > 1) {
      Rcpp::warning("No valid ordering requested -- using default.");
    }
  } else if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  solver.compute(A);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}


// Solve Ax = b using dense LLT (pivot = 0) or LDLT (pivot = 1) Cholesky decomposition
// [[Rcpp::export]]
Rcpp::List solve_dense_chol(
    const Eigen::Map <Eigen::MatrixXd> A,
    const Eigen::Map <Eigen::MatrixXd> b,
    unsigned int pivot = 1) {

  Eigen::LDLT <Eigen::MatrixXd> solver;
  if(pivot == 0) {
    Eigen::LLT <Eigen::MatrixXd> solver;
  } else if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  solver.compute(A);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}
