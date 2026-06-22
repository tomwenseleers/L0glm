// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

// Solve Ax = b using dense SVD decompositions. Sparse inputs are converted to
// dense matrices because Eigen's core sparse solvers do not include sparse SVD.

template <typename SolverType, typename MatrixType, typename RhsType>
Rcpp::List run_svd_solver(const MatrixType &A,
                          const RhsType &b,
                          unsigned int compute_flags,
                          double threshold) {
  SolverType solver;
  if(threshold > 0) {
    solver.setThreshold(threshold);
  }

  solver.compute(A, compute_flags);

  Eigen::MatrixXd x = solver.solve(b);
  if(!x.allFinite()) {
    return Rcpp::List::create(
      Rcpp::Named("status") = false,
      Rcpp::Named("rank") = solver.rank());
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = true,
    Rcpp::Named("rank") = solver.rank(),
    Rcpp::Named("coefficients") = x);
}

template <typename MatrixType, typename RhsType>
Rcpp::List run_svd_dispatch(const MatrixType &A,
                            const RhsType &b,
                            unsigned int svd,
                            unsigned int precond,
                            double threshold) {
  if(svd == 0) {
    typedef Eigen::BDCSVD<Eigen::MatrixXd> Solver;
    return run_svd_solver<Solver>(A, b,
                                  Eigen::ComputeThinU | Eigen::ComputeThinV,
                                  threshold);
  }

  if(svd > 1) {
    Rcpp::warning("No valid SVD requested -- using JacobiSVD.");
  }

  if(precond == 1) {
    typedef Eigen::JacobiSVD<Eigen::MatrixXd,
                             Eigen::FullPivHouseholderQRPreconditioner> Solver;
    return run_svd_solver<Solver>(A, b,
                                  Eigen::ComputeFullU | Eigen::ComputeFullV,
                                  threshold);
  }

  if(precond == 2) {
    typedef Eigen::JacobiSVD<Eigen::MatrixXd,
                             Eigen::HouseholderQRPreconditioner> Solver;
    return run_svd_solver<Solver>(A, b,
                                  Eigen::ComputeThinU | Eigen::ComputeThinV,
                                  threshold);
  }

  if(precond == 3) {
    if(A.rows() != A.cols()) {
      return Rcpp::List::create(
        Rcpp::Named("status") = false,
        Rcpp::Named("message") =
          "JacobiSVD NoQRPreconditioner requires a square matrix.");
    }

    typedef Eigen::JacobiSVD<Eigen::MatrixXd,
                             Eigen::NoQRPreconditioner> Solver;
    return run_svd_solver<Solver>(A, b,
                                  Eigen::ComputeFullU | Eigen::ComputeFullV,
                                  threshold);
  }

  if(precond > 3) {
    Rcpp::warning("No valid JacobiSVD preconditioner requested -- using ColPivHouseholderQRPreconditioner.");
  }

  typedef Eigen::JacobiSVD<Eigen::MatrixXd,
                           Eigen::ColPivHouseholderQRPreconditioner> Solver;
  return run_svd_solver<Solver>(A, b,
                                Eigen::ComputeThinU | Eigen::ComputeThinV,
                                threshold);
}

// svd = 0: BDCSVD; svd = 1: JacobiSVD.
// For JacobiSVD, precond = 0: ColPivHouseholderQR, 1: FullPivHouseholderQR,
// 2: HouseholderQR, 3: NoQR (square matrices only).
// [[Rcpp::export]]
Rcpp::List solve_dense_SVD(const Eigen::Map<Eigen::MatrixXd> &A,
                           const Eigen::Map<Eigen::MatrixXd> b,
                           unsigned int svd = 0,
                           unsigned int precond = 0,
                           double threshold = 0) {
  return run_svd_dispatch(A, b, svd, precond, threshold);
}

// [[Rcpp::export]]
Rcpp::List solve_sparse_SVD(const Eigen::MappedSparseMatrix<double> &A,
                            const Eigen::Map<Eigen::MatrixXd> b,
                            unsigned int svd = 0,
                            unsigned int precond = 0,
                            double threshold = 0) {
  Eigen::MatrixXd A_dense(A);
  return run_svd_dispatch(A_dense, b, svd, precond, threshold);
}
