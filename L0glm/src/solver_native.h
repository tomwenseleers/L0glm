#ifndef L0GLM_SOLVER_NATIVE_H
#define L0GLM_SOLVER_NATIVE_H
#include <RcppEigen.h>

// Allocation-free solver interface shared by the native IRLS driver and the
// Rcpp entry points. R lists and warnings are assembled on the R thread only.
namespace L0glm {
using DenseRef = Eigen::Ref<const Eigen::MatrixXd>;
using SparseRef = Eigen::Ref<const Eigen::SparseMatrix<double>>;
struct SolverResult {
  Eigen::MatrixXd coefficients;
  bool success = false;
  int info = 0, iterations = 0, rank = 0;
  double error = 0;
  bool finite() const { return coefficients.size() && coefficients.allFinite(); }
  Rcpp::List list(int kind = 0) const {
    Rcpp::List value = Rcpp::List::create(Rcpp::_["status"] = success);
    if (kind != 2 && kind != 3) value["info"] = info;
    if (kind == 1) { value["iterations"] = iterations; value["error"] = error; }
    if (kind == 3) value["rank"] = rank;
    if (finite()) value["coefficients"] = coefficients;
    return value;
  }
};

template<class Solver, class Matrix> SolverResult factor_solve(const Matrix& A, DenseRef b) {
  SolverResult result; Solver solver; solver.compute(A);
  result.info = solver.info();
  if (solver.info() != Eigen::Success) return result;
  result.coefficients = solver.solve(b); result.info = solver.info();
  result.success = solver.info() == Eigen::Success && result.finite();
  return result;
}
template<class Solver, class Matrix> SolverResult iterative_solve(
    const Matrix& A, DenseRef b, DenseRef start, double tol, int maxit) {
  SolverResult result; Solver solver;
  if (tol > 0) solver.setTolerance(tol);
  if (maxit > 0) solver.setMaxIterations(maxit);
  solver.compute(A); result.info = solver.info();
  if (solver.info() != Eigen::Success) return result;
  result.coefficients = solver.solveWithGuess(b, start);
  result.info = solver.info(); result.iterations = solver.iterations(); result.error = solver.error();
  result.success = solver.info() == Eigen::Success && result.finite();
  return result;
}

SolverResult chol_solve(DenseRef, DenseRef, unsigned int, unsigned int = 0);
SolverResult chol_solve(SparseRef, DenseRef, unsigned int, unsigned int = 0);
SolverResult cg_solve(DenseRef, DenseRef, DenseRef, double, int, unsigned int);
SolverResult cg_solve(SparseRef, DenseRef, DenseRef, double, int, unsigned int);
SolverResult lscg_solve(DenseRef, DenseRef, DenseRef, double, int, unsigned int);
SolverResult lscg_solve(SparseRef, DenseRef, DenseRef, double, int, unsigned int);
SolverResult qr_solve(DenseRef, DenseRef, unsigned int, unsigned int = 1);
SolverResult qr_solve(SparseRef, DenseRef, unsigned int, unsigned int = 1);
SolverResult svd_solve(DenseRef, DenseRef, unsigned int, unsigned int, double);
SolverResult svd_solve(SparseRef, DenseRef, unsigned int, unsigned int, double);
}
#endif
