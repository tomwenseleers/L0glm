#ifndef L0GLM_RIDGE_PROFILE_H
#define L0GLM_RIDGE_PROFILE_H

// Observation-space ridge solves with optional unpenalized nuisance columns.
namespace L0glm {
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

inline MatrixXd row_gram(const MatrixXd& Z) {
  const int n = Z.rows(), p = Z.cols();
  MatrixXd A = MatrixXd::Zero(n, n);
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)("L", "N", &n, &p, &one, Z.data(), &n, &zero,
                  A.data(), &n FCONE FCONE);
  return A;
}

inline SparseMatrix<double> row_gram(const SparseMatrix<double>& Z) {
  return Z * Z.transpose();
}

inline void drop_exact_zeros(MatrixXd&) {}
inline void drop_exact_zeros(SparseMatrix<double>& Z) { Z.prune(0.0); }

inline bool solve_system(const MatrixXd& A, const MatrixXd& b,
                          MatrixXd& x, bool) {
  Eigen::LLT<MatrixXd> solver(A);
  if (solver.info() != Eigen::Success) return false;
  x = solver.solve(b);
  return solver.info() == Eigen::Success && x.allFinite();
}

inline bool solve_system(const SparseMatrix<double>& A, const MatrixXd& b,
                          MatrixXd& x, bool natural) {
  if (natural) {
    Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) return false;
    x = solver.solve(b);
    return solver.info() == Eigen::Success && x.allFinite();
  }
  Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) return false;
  x = solver.solve(b);
  return solver.info() == Eigen::Success && x.allFinite();
}

template <typename Design>
bool ridge_dual(const Design& X, const VectorXd& y, const VectorXd& penalty,
                 VectorXd center, VectorXd& beta, bool natural) {
  const int n = X.rows(), p = X.cols();
  const bool penalized = (penalty.array() > 0.0).any();
  std::vector<int> free;
  VectorXd invsqrt = VectorXd::Ones(p);
  for (int j = 0; j < p; ++j) {
    if (penalized) {
      invsqrt[j] = penalty[j] > 0.0 ? 1.0 / std::sqrt(penalty[j]) : 0.0;
      if (penalty[j] == 0.0) { free.push_back(j); center[j] = 0.0; }
    } else center[j] = 0.0;
  }
  const int q = free.size();
  Design Z = X * invsqrt.asDiagonal();
  // Zero-penalty nuisance columns have zero scale here. Retained sparse
  // entries would create structural fill in the observation-space Gram.
  drop_exact_zeros(Z);
  Design A = row_gram(Z);
  if (penalized) for (int i = 0; i < n; ++i) A.coeffRef(i, i) += 1.0;
  MatrixXd rhs(n, q + 1), U(n, q), solution;
  rhs.col(0) = y - X * center;
  for (int j = 0; j < q; ++j) {
    U.col(j) = X.col(free[j]);
    rhs.col(j + 1) = U.col(j);
  }
  if (!solve_system(A, rhs, solution, natural)) return false;
  VectorXd u = solution.col(0), alpha(q);
  if (q) {
    MatrixXd G = U.transpose() * solution.rightCols(q);
    alpha = G.colPivHouseholderQr().solve(U.transpose() * u);
    u.noalias() -= solution.rightCols(q) * alpha;
  }
  beta = center + invsqrt.asDiagonal() * (Z.transpose() * u);
  for (int j = 0; j < q; ++j) beta[free[j]] = alpha[j];
  return beta.allFinite();
}
} // namespace L0glm
#endif
