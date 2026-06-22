// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// Solve Ax = b using Cholesky decomposition for sparse or dense covariate matrix A
// adapted from https://github.com/nk027/sanic/blob/master/src/solve_Cholesky.cpp

template <typename SolverType, typename MatrixType, typename RhsType>
Rcpp::List run_chol_solver(const MatrixType &A, const RhsType &b) {
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

// Solve Ax = b using sparse LLT (pivot = 0) or LDLT (pivot = 1) Cholesky decomposition
// [[Rcpp::export]]
Rcpp::List solve_sparse_chol(
    const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,
    unsigned int pivot = 1, unsigned int ordering = 0) {

  if(ordering > 1) {
    Rcpp::warning("No valid ordering requested -- using default.");
    ordering = 0;
  }

  if(pivot == 0) {
    if(ordering == 1) { // use NaturalOrdering
      typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,
                                   Eigen::Lower,
                                   Eigen::NaturalOrdering<int> > Solver;
      return run_chol_solver<Solver>(A, b);
    }

    typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,
                                 Eigen::Lower,
                                 Eigen::AMDOrdering<int> > Solver;
    return run_chol_solver<Solver>(A, b);
  }

  if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  if(ordering == 1) {
    typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,
                                  Eigen::Lower,
                                  Eigen::NaturalOrdering<int> > Solver;
    return run_chol_solver<Solver>(A, b);
  }

  typedef Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,
                                Eigen::Lower,
                                Eigen::AMDOrdering<int> > Solver;
  return run_chol_solver<Solver>(A, b);
}


// Solve Ax = b using dense LLT (pivot = 0) or LDLT (pivot = 1) Cholesky decomposition
// [[Rcpp::export]]
Rcpp::List solve_dense_chol(
    const Eigen::Map <Eigen::MatrixXd> &A,
    const Eigen::Map <Eigen::MatrixXd> b,
    unsigned int pivot = 1) {

  if(pivot == 0) {
    typedef Eigen::LLT<Eigen::MatrixXd> Solver;
    return run_chol_solver<Solver>(A, b);
  }

  if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  typedef Eigen::LDLT<Eigen::MatrixXd> Solver;
  return run_chol_solver<Solver>(A, b);
}
