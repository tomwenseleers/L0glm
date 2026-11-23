#include "solver_native.h"

namespace L0glm {
namespace {
using SparseMatrix = Eigen::SparseMatrix<double>;

template<class Ordering>
SolverResult sparse_cholesky(SparseRef matrix, DenseRef rhs, bool llt) {
  using LLT = Eigen::SimplicialLLT<SparseMatrix, Eigen::Lower, Ordering>;
  using LDLT = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Ordering>;
  return llt ? factor_solve<LLT>(matrix, rhs)
             : factor_solve<LDLT>(matrix, rhs);
}

unsigned int cholesky_pivot(unsigned int pivot) {
  if (pivot <= 1) return pivot;
  Rcpp::warning("No valid pivoting scheme requested -- using default.");
  return 1;
}
}

SolverResult chol_solve(DenseRef matrix, DenseRef rhs, unsigned int pivot,
                        unsigned int) {
  using LLT = Eigen::LLT<Eigen::MatrixXd, Eigen::Lower>;
  using LDLT = Eigen::LDLT<Eigen::MatrixXd, Eigen::Lower>;
  return pivot == 0 ? factor_solve<LLT>(matrix, rhs)
                    : factor_solve<LDLT>(matrix, rhs);
}

SolverResult chol_solve(SparseRef matrix, DenseRef rhs, unsigned int pivot,
                        unsigned int ordering) {
  const bool llt = pivot == 0;
  if (ordering == 1)
    return sparse_cholesky<Eigen::NaturalOrdering<int>>(matrix, rhs, llt);
  return sparse_cholesky<Eigen::AMDOrdering<int>>(matrix, rhs, llt);
}
}

// [[Rcpp::export]]
Rcpp::List solve_sparse_chol(const Eigen::MappedSparseMatrix<double> &A,
                             const Eigen::Map<Eigen::MatrixXd> b,
                             unsigned int pivot = 1, unsigned int ordering = 0) {
  if (ordering > 1) {
    Rcpp::warning("No valid ordering requested -- using default.");
    ordering = 0;
  }
  return L0glm::chol_solve(L0glm::SparseRef(A), b, L0glm::cholesky_pivot(pivot),
                          ordering).list();
}

// [[Rcpp::export]]
Rcpp::List solve_dense_chol(const Eigen::Map<Eigen::MatrixXd> &A,
                            const Eigen::Map<Eigen::MatrixXd> b,
                            unsigned int pivot = 1) {
  return L0glm::chol_solve(L0glm::DenseRef(A), b, L0glm::cholesky_pivot(pivot)).list();
}
