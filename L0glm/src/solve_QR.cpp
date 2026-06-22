// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SparseQR>

// Solve Ax = b using QR decomposition for sparse or dense covariate matrix A
// adapted from https://github.com/nk027/sanic/blob/master/src/solve_QR.cpp

template <typename SolverType, typename MatrixType, typename RhsType>
Rcpp::List run_sparse_qr_solver(const MatrixType &A, const RhsType &b) {
  SolverType solver;
  solver.compute(A);

  if(solver.info() != Eigen::Success) {
    return Rcpp::List::create(
      Rcpp::Named("status") = false,
      Rcpp::Named("info") = static_cast<int>(solver.info()));
  }

  Eigen::MatrixXd x = solver.solve(b);
  bool success = solver.info() == Eigen::Success;

  if(!x.allFinite()) {
    return Rcpp::List::create(
      Rcpp::Named("status") = false,
      Rcpp::Named("info") = static_cast<int>(solver.info()));
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = success,
    Rcpp::Named("info") = static_cast<int>(solver.info()),
    Rcpp::Named("coefficients") = x);
}

template <typename SolverType, typename MatrixType, typename RhsType>
Rcpp::List run_dense_qr_solver(const MatrixType &A, const RhsType &b) {
  SolverType solver;
  solver.compute(A);

  Eigen::MatrixXd x = solver.solve(b);

  if(!x.allFinite()) {
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = true,
    Rcpp::Named("coefficients") = x);
}

template <typename SolverType, typename MatrixType, typename RhsType>
Rcpp::List run_dense_ranked_qr_solver(const MatrixType &A, const RhsType &b) {
  SolverType solver;
  solver.compute(A);

  if(solver.rank() < A.cols()) {
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solve(b);

  if(!x.allFinite()) {
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = true,
    Rcpp::Named("coefficients") = x);
}

// Solve Ax = b using sparse QR decomposition
// [[Rcpp::export]]
Rcpp::List solve_sparse_QR(const Eigen::MappedSparseMatrix<double> &A,
                           const Eigen::Map<Eigen::MatrixXd> b,
                           unsigned int ordering = 1){

  if(ordering == 0) {
    Eigen::SparseMatrix<double> A_copy(A);
    typedef Eigen::SparseQR<Eigen::SparseMatrix<double>,
                            Eigen::AMDOrdering<int> > Solver;
    return run_sparse_qr_solver<Solver>(A_copy, b);
  }

  if(ordering == 2) {
    Eigen::SparseMatrix<double> A_copy(A);
    typedef Eigen::SparseQR<Eigen::SparseMatrix<double>,
                            Eigen::NaturalOrdering<int> > Solver;
    return run_sparse_qr_solver<Solver>(A_copy, b);
  }

  if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }

  Eigen::SparseMatrix<double> A_copy(A);
  typedef Eigen::SparseQR<Eigen::SparseMatrix<double>,
                          Eigen::COLAMDOrdering<int> > Solver;
  return run_sparse_qr_solver<Solver>(A_copy, b);
}


// Solve Ax = b using dense QR decomposition(pivot=0: HouseHolderQR, pivot=1: ColPivHouseholderQR,
// pivot=2: FullPivHouseholderQR)
// [[Rcpp::export]]
Rcpp::List solve_dense_QR(const Eigen::Map<Eigen::MatrixXd> &A,
                          const Eigen::Map<Eigen::MatrixXd> b,
                          unsigned int pivot = 1) {

  if(pivot == 0) {
    typedef Eigen::HouseholderQR<Eigen::MatrixXd> Solver;
    return run_dense_qr_solver<Solver>(A, b);
  }

  if(pivot == 2) {
    typedef Eigen::FullPivHouseholderQR<Eigen::MatrixXd> Solver;
    return run_dense_ranked_qr_solver<Solver>(A, b);
  }

  if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default ColPivHouseholderQR.");
  }

  typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXd> Solver;
  return run_dense_ranked_qr_solver<Solver>(A, b);
}
