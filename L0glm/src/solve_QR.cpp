// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SparseQR>

// Solve Ax = b using QR decomposition for sparse or dense covariate matrix A
// adapted from https://github.com/nk027/sanic/blob/master/src/solve_QR.cpp

// Solve Ax = b using sparse QR decomposition
// [[Rcpp::export]]
Rcpp::List solve_sparse_QR(const Eigen::MappedSparseMatrix<double> A,
                           const Eigen::Map<Eigen::MatrixXd> b,
                           unsigned int ordering = 1){

  Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
  if(ordering == 0) {
    Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
  } else if(ordering == 2) {
    Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::NaturalOrdering<int> > solver;
  } else if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }

  solver.compute(A);
  if(solver.info() != Eigen::Success) {
    // decomposition failed
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


// Solve Ax = b using dense QR decomposition(pivot=0: HouseHolderQR, pivot=1: ColPivHouseholderQR,
// pivot=2: FullPivHouseholderQR)
// [[Rcpp::export]]
Rcpp::List solve_dense_QR(const Eigen::Map<Eigen::MatrixXd> A,
                          const Eigen::Map<Eigen::MatrixXd> b,
                          unsigned int pivot = 1) {

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver;
  if(pivot == 0) {
      Eigen::HouseholderQR<Eigen::MatrixXd> solver;
    } else if(pivot == 2) {
      Eigen::FullPivHouseholderQR<Eigen::MatrixXd> solver;
    } else if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default ColPivHouseholderQR.");
  }

  solver.compute(A);

  if(!solver.isInvertible()) {
   // decomposition failed
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

