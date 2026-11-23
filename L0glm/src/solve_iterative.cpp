#include "solver_native.h"
#include <type_traits>

namespace L0glm {
namespace {
using SparseMatrix = Eigen::SparseMatrix<double>;
using DenseMatrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using IncompleteCholesky = Eigen::IncompleteCholesky<double>;

// Eigen's incomplete Cholesky takes sparse storage. Use it on the system
// for CG and on the normal system for LSCG, leaving the operators unchanged.
template<bool normal_system>
class SparseCholeskyPreconditioner : public IncompleteCholesky {
public:
  template<class Matrix>
  SparseCholeskyPreconditioner& compute(const Matrix& matrix) {
    SparseMatrix system;
    if constexpr (!normal_system) {
      system = matrix.sparseView(0., 0.);
    } else if constexpr (std::is_same_v<typename Matrix::StorageKind, Eigen::Sparse>) {
      system = matrix.transpose() * matrix;
      system.prune(0., 0.);
    } else {
      const DenseMatrix gram = matrix.transpose() * matrix;
      system = gram.sparseView(0., 0.);
    }
    IncompleteCholesky::compute(system);
    return *this;
  }
};

template<class Preconditioner, class Matrix>
SolverResult solve_cg(const Matrix& matrix, DenseRef rhs, DenseRef guess,
                      double tolerance, int max_iterations) {
  // The normal equations supplied by L0glm store both triangles.
  using CG = Eigen::ConjugateGradient<typename Matrix::PlainObject,
                                      Eigen::Lower | Eigen::Upper, Preconditioner>;
  return iterative_solve<CG>(matrix, rhs, guess, tolerance, max_iterations);
}

template<class Matrix>
SolverResult preconditioned_cg(const Matrix& matrix, DenseRef rhs, DenseRef guess,
                               double tolerance, int max_iterations, unsigned int precond) {
  if (precond == 0)
    return solve_cg<Eigen::IdentityPreconditioner>(matrix, rhs, guess, tolerance, max_iterations);
  if (precond == 2) {
    constexpr bool sparse = std::is_same_v<typename Matrix::StorageKind, Eigen::Sparse>;
    using IC = std::conditional_t<sparse, IncompleteCholesky, SparseCholeskyPreconditioner<false>>;
    return solve_cg<IC>(matrix, rhs, guess, tolerance, max_iterations);
  }
  return solve_cg<Eigen::DiagonalPreconditioner<double>>(matrix, rhs, guess, tolerance, max_iterations);
}

double normal_residual_error(double norm, double scale) {
  if (!std::isfinite(norm)) return std::numeric_limits<double>::infinity();
  if (scale > 0) return norm / scale;
  return norm == 0 ? 0 : std::numeric_limits<double>::infinity();
}

template<class Preconditioner, bool restart = false, class Matrix>
SolverResult solve_lscg(const Matrix& matrix, DenseRef rhs, DenseRef start,
                        double tolerance, int max_iterations) {
  Eigen::LeastSquaresConjugateGradient<typename Matrix::PlainObject, Preconditioner> solver;
  if (tolerance > 0) solver.setTolerance(tolerance);
  if (max_iterations > 0) solver.setMaxIterations(max_iterations);
  solver.compute(matrix);

  SolverResult result;
  result.info = solver.info();
  if (result.info != Eigen::Success) return result;
  const int budget = solver.maxIterations();
  const double target = solver.tolerance();
  result.coefficients.resize(start.rows(), rhs.cols());
  Vector residual(rhs.rows()), normal(matrix.cols());
  for (Eigen::Index col = 0; col < rhs.cols(); ++col) {
    const auto response = rhs.col(col);
    auto solution = result.coefficients.col(col);
    Vector guess = start.col(col);
    normal.noalias() = matrix.transpose() * response;
    const double scale = normal.norm();
    const auto normal_residual = [&](const Vector& beta) {
      residual.noalias() = matrix * beta;
      residual = response - residual;
      normal.noalias() = matrix.transpose() * residual;
      return normal.norm();
    };

    if constexpr (restart) {
      // A warm start from a different ridge system can be worse than zero.
      if (normal_residual(guess) > scale) guess.setZero();
    }
    solver.setMaxIterations(budget);
    // Evaluate the lazy solve into a plain vector so each RHS is a single
    // Eigen solve even when the destination is a column of a matrix.
    Vector estimate = solver.solveWithGuess(response, guess);
    const int first_iterations = solver.iterations();
    result.iterations += first_iterations;
    double error = normal_residual_error(normal_residual(estimate), scale);
    if constexpr (restart) {
      // A successful recursive residual may underestimate the true residual.
      // Count the first update even when Eigen reports zero iterations.
      if (error > target && solver.info() == Eigen::Success &&
          first_iterations < budget - 1 && estimate.allFinite()) {
        solver.setMaxIterations(budget - first_iterations - 1);
        guess = solver.solveWithGuess(response, estimate);
        result.iterations += solver.iterations();
        estimate.swap(guess);
        error = normal_residual_error(normal_residual(estimate), scale);
      }
    }
    solution = estimate;
    result.error = std::max(result.error, error);
  }
  result.success = result.finite() && result.error <= target;
  result.info = result.success ? Eigen::Success : Eigen::NoConvergence;
  return result;
}

template<class Matrix>
SolverResult preconditioned_lscg(const Matrix& matrix, DenseRef rhs, DenseRef guess,
                                 double tolerance, int max_iterations, unsigned int precond) {
  if (precond == 0)
    return solve_lscg<Eigen::IdentityPreconditioner>(matrix, rhs, guess, tolerance, max_iterations);
  if (precond == 2)
    return solve_lscg<SparseCholeskyPreconditioner<true>, true>(matrix, rhs, guess,
                                                              tolerance, max_iterations);
  return solve_lscg<Eigen::LeastSquareDiagonalPreconditioner<double>>(matrix, rhs, guess,
                                                                       tolerance, max_iterations);
}

unsigned int checked_preconditioner(unsigned int precond) {
  if (precond <= 2) return precond;
  Rcpp::warning("No valid preconditioner requested -- using incomplete Cholesky.");
  return 2;
}

Eigen::Map<const DenseMatrix> map_numeric_matrix(const Rcpp::NumericMatrix& matrix) {
  return {matrix.begin(), matrix.nrow(), matrix.ncol()};
}
}

SolverResult cg_solve(DenseRef A, DenseRef b, DenseRef guess, double tol,
                       int maxit, unsigned int precond) {
  return preconditioned_cg(A, b, guess, tol, maxit, precond);
}
SolverResult cg_solve(SparseRef A, DenseRef b, DenseRef guess, double tol,
                       int maxit, unsigned int precond) {
  return preconditioned_cg(A, b, guess, tol, maxit, precond);
}
SolverResult lscg_solve(DenseRef A, DenseRef b, DenseRef guess, double tol,
                         int maxit, unsigned int precond) {
  return preconditioned_lscg(A, b, guess, tol, maxit, precond);
}
SolverResult lscg_solve(SparseRef A, DenseRef b, DenseRef guess, double tol,
                         int maxit, unsigned int precond) {
  return preconditioned_lscg(A, b, guess, tol, maxit, precond);
}
}

// [[Rcpp::export]]
Rcpp::List solve_sparse_LSCG(const Eigen::MappedSparseMatrix<double> &A,
                             const Eigen::Map<Eigen::MatrixXd> b,
                             const Eigen::Map<Eigen::MatrixXd> x0,
                             double tol = 0, int maxit = 0, unsigned int precond = 2) {
  return L0glm::lscg_solve(L0glm::SparseRef(A), b, x0, tol, maxit,
                          L0glm::checked_preconditioner(precond)).list(1);
}

// [[Rcpp::export]]
Rcpp::List solve_dense_LSCG(const Rcpp::NumericMatrix& A,
                            const Rcpp::NumericMatrix& b,
                            const Rcpp::NumericMatrix& x0,
                            double tol = 0, int maxit = 0, unsigned int precond = 2) {
  return L0glm::lscg_solve(L0glm::map_numeric_matrix(A), L0glm::map_numeric_matrix(b),
                           L0glm::map_numeric_matrix(x0), tol, maxit,
                          L0glm::checked_preconditioner(precond)).list(1);
}

// [[Rcpp::export]]
Rcpp::List solve_sparse_CG(Eigen::MappedSparseMatrix<double> &A, Eigen::Map<Eigen::MatrixXd> b,
                           Eigen::Map<Eigen::MatrixXd> x0, double tol = 0, int maxit = 0,
                           unsigned int precond = 2) {
  return L0glm::cg_solve(L0glm::SparseRef(A), b, x0, tol, maxit,
                        L0glm::checked_preconditioner(precond)).list(1);
}

// [[Rcpp::export]]
Rcpp::List solve_dense_CG(Eigen::Map<Eigen::MatrixXd> &A, Eigen::Map<Eigen::MatrixXd> b,
                          Eigen::Map<Eigen::MatrixXd> x0, double tol = 0, int maxit = 0,
                          unsigned int precond = 2) {
  return L0glm::cg_solve(L0glm::DenseRef(A), b, x0, tol, maxit,
                        L0glm::checked_preconditioner(precond)).list(1);
}
