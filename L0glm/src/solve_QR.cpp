#include "solver_native.h"
#include <Eigen/SparseQR>

namespace L0glm {
namespace {
using DenseMatrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

template<class Decomposition>
SolverResult pivoted_qr(DenseRef matrix, DenseRef rhs) {
  Decomposition qr(matrix);
  if (qr.rank() != matrix.cols()) return {};
  SolverResult result;
  result.coefficients = qr.solve(rhs);
  result.success = result.finite();
  return result;
}

template<class Ordering>
SolverResult sparse_qr(SparseRef matrix, DenseRef rhs) {
  using Decomposition = Eigen::SparseQR<SparseMatrix, Ordering>;
  return factor_solve<Decomposition>(matrix, rhs);
}
}

SolverResult qr_solve(DenseRef matrix, DenseRef rhs, unsigned int pivot,
                      unsigned int) {
  if (pivot == 0) {
    Eigen::HouseholderQR<DenseMatrix> qr(matrix);
    SolverResult result;
    result.coefficients = qr.solve(rhs);
    result.success = result.finite();
    return result;
  }
  if (pivot == 2)
    return pivoted_qr<Eigen::FullPivHouseholderQR<DenseMatrix>>(matrix, rhs);
  return pivoted_qr<Eigen::ColPivHouseholderQR<DenseMatrix>>(matrix, rhs);
}

SolverResult qr_solve(SparseRef matrix, DenseRef rhs, unsigned int,
                      unsigned int ordering) {
  // AMD symmetrizes the pattern, which requires a square matrix.
  if (ordering == 0 && matrix.rows() == matrix.cols())
    return sparse_qr<Eigen::AMDOrdering<int>>(matrix, rhs);
  if (ordering == 2)
    return sparse_qr<Eigen::NaturalOrdering<int>>(matrix, rhs);
  return sparse_qr<Eigen::COLAMDOrdering<int>>(matrix, rhs);
}
}

// [[Rcpp::export]]
Rcpp::List solve_sparse_QR(const Eigen::MappedSparseMatrix<double> &A,
                           const Eigen::Map<Eigen::MatrixXd> b,
                           unsigned int ordering = 1) {
  if (ordering > 2)
    Rcpp::warning("No valid ordering requested -- using default.");
  return L0glm::qr_solve(L0glm::SparseRef(A), b, 1, ordering).list();
}

// [[Rcpp::export]]
Rcpp::List solve_dense_QR(const Eigen::Map<Eigen::MatrixXd> &A,
                          const Eigen::Map<Eigen::MatrixXd> b,
                          unsigned int pivot = 1) {
  if (pivot > 2)
    Rcpp::warning("No valid pivoting scheme requested -- using default ColPivHouseholderQR.");
  return L0glm::qr_solve(L0glm::DenseRef(A), b, pivot).list(2);
}
