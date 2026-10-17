// [[Rcpp::depends(RcppEigen)]]
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <R_ext/BLAS.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace {

using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;

struct FitResult {
  VectorXd beta;
  VectorXd weights;
  VectorXd z;
  VectorXd lambdas;
};

inline bool finite_pos(double x) {
  return std::isfinite(x) && x > 0.0;
}

inline VectorXd constrain_vec(const VectorXd& x) {
  return x.cwiseMax(0.0);
}

inline std::vector<int> active_positions(const VectorXd& beta) {
  std::vector<int> active;
  active.reserve(beta.size());
  for (int i = 0; i < beta.size(); ++i) {
    if (beta[i] != 0.0) active.push_back(i);
  }
  return active;
}

inline void group_support_start(const MatrixXd& coef,
                                std::vector<int>& vars,
                                MatrixXd& start_current) {
  const int p = coef.rows();
  const int m = coef.cols();
  vars.clear();
  std::vector<double> starts;
  starts.reserve(p);

  for (int j = 0; j < p; ++j) {
    bool active = false;
    double signed_sum = 0.0;
    double abs_sum = 0.0;
    for (int col = 0; col < m; ++col) {
      const double value = coef(j, col);
      signed_sum += value;
      abs_sum += std::abs(value);
      active = active || value != 0.0;
    }
    if (active) {
      const double avg = signed_sum / static_cast<double>(m);
      const double mean_abs = abs_sum / static_cast<double>(m);
      vars.push_back(j);
      starts.push_back((avg == 0.0) ? mean_abs : avg);
    }
  }

  start_current.resize(static_cast<int>(vars.size()), m);
  for (int col = 0; col < m; ++col) {
    for (int a = 0; a < static_cast<int>(vars.size()); ++a) {
      start_current(a, col) = starts[a];
    }
  }
}

inline VectorXd compute_scaling_dense(const Map<MatrixXd>& X,
                                      const std::vector<int>& vars) {
  VectorXd scaling(vars.size());
  for (int j = 0; j < static_cast<int>(vars.size()); ++j) {
    double ss = 0.0;
    const int col = vars[j];
    for (int i = 0; i < X.rows(); ++i) {
      const double x = X(i, col);
      ss += x * x;
    }
    const double s = std::sqrt(ss);
    scaling[j] = (s == 0.0 || !std::isfinite(s)) ? 1.0 : s;
  }
  return scaling;
}

inline VectorXd compute_scaling_sparse(const Eigen::MappedSparseMatrix<double>& X,
                                       const std::vector<int>& vars) {
  VectorXd scaling(vars.size());
  for (int j = 0; j < static_cast<int>(vars.size()); ++j) {
    double ss = 0.0;
    const int col = vars[j];
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(X, col); it; ++it) {
      ss += it.value() * it.value();
    }
    const double s = std::sqrt(ss);
    scaling[j] = (s == 0.0 || !std::isfinite(s)) ? 1.0 : s;
  }
  return scaling;
}

inline VectorXd predict_dense(const Map<MatrixXd>& X, const VectorXd& offset,
                              const std::vector<int>& vars,
                              const VectorXd& beta) {
  VectorXd eta = offset;
  for (int j = 0; j < static_cast<int>(vars.size()); ++j) {
    const double bj = beta[j];
    if (bj == 0.0) continue;
    eta.noalias() += X.col(vars[j]) * bj;
  }
  return eta;
}

inline VectorXd predict_sparse(const Eigen::MappedSparseMatrix<double>& X,
                               const VectorXd& offset,
                               const std::vector<int>& vars,
                               const VectorXd& beta) {
  VectorXd eta = offset;
  for (int j = 0; j < static_cast<int>(vars.size()); ++j) {
    const double bj = beta[j];
    if (bj == 0.0) continue;
    const int col = vars[j];
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(X, col); it; ++it) {
      eta[it.row()] += it.value() * bj;
    }
  }
  return eta;
}

inline VectorXd init_start_dense(const Map<MatrixXd>& X, const VectorXd& y,
                                 const VectorXd& prior_weights,
                                 const VectorXd& offset,
                                 const VectorXd& prior_mean,
                                 const std::vector<int>& vars,
                                 double lambda) {
  const int n = X.rows();
  const int q = static_cast<int>(vars.size());
  VectorXd numer = VectorXd::Zero(q);
  VectorXd denom = VectorXd::Constant(q, lambda);
  for (int j = 0; j < q; ++j) numer[j] = lambda * prior_mean[j];

  for (int i = 0; i < n; ++i) {
    const double mu = y[i] + 0.1;
    const double w = prior_weights[i] / mu;
    if (!finite_pos(w)) continue;
    const double z = y[i] - offset[i];
    for (int j = 0; j < q; ++j) {
      const double x = X(i, vars[j]);
      numer[j] += w * x * z;
      denom[j] += w * x * x;
    }
  }

  VectorXd beta(q);
  for (int j = 0; j < q; ++j) {
    beta[j] = (denom[j] != 0.0 && std::isfinite(denom[j])) ?
      numer[j] / denom[j] : 0.0;
  }
  beta = constrain_vec(beta);
  return beta;
}

inline VectorXd init_start_sparse(const Eigen::MappedSparseMatrix<double>& X,
                                  const VectorXd& y,
                                  const VectorXd& prior_weights,
                                  const VectorXd& offset,
                                  const VectorXd& prior_mean,
                                  const std::vector<int>& vars,
                                  double lambda) {
  const int n = X.rows();
  const int q = static_cast<int>(vars.size());
  VectorXd w(n);
  VectorXd z(n);
  for (int i = 0; i < n; ++i) {
    const double wi = prior_weights[i] / (y[i] + 0.1);
    w[i] = finite_pos(wi) ? wi : 0.0;
    z[i] = y[i] - offset[i];
  }

  VectorXd numer(q);
  VectorXd denom(q);
  for (int j = 0; j < q; ++j) {
    numer[j] = lambda * prior_mean[j];
    denom[j] = lambda;
    const int col = vars[j];
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(X, col); it; ++it) {
      const int row = it.row();
      const double wi = w[row];
      if (wi == 0.0) continue;
      const double x = it.value();
      numer[j] += wi * x * z[row];
      denom[j] += wi * x * x;
    }
  }

  VectorXd beta(q);
  for (int j = 0; j < q; ++j) {
    beta[j] = (denom[j] != 0.0 && std::isfinite(denom[j])) ?
      numer[j] / denom[j] : 0.0;
  }
  beta = constrain_vec(beta);
  return beta;
}

inline VectorXd solve_ridge_dense(const Map<MatrixXd>& X, const VectorXd& z,
                                  const VectorXd& w,
                                  const std::vector<int>& vars,
                                  const std::vector<int>& active,
                                  const VectorXd& lambdas,
                                  const VectorXd& prior_mean,
                                  const VectorXd& scaling,
                                  const VectorXd& start) {
  const int n = X.rows();
  const int k = static_cast<int>(active.size());
  // Assemble weighted columns once, then use R's BLAS symmetric crossproduct.
  // The former a/i/b loop traversed X across columns in its innermost loop,
  // producing O(n*k^2) scalar operations with poor cache locality.
  const VectorXd sqrtw = w.array().max(0.0).sqrt();
  const VectorXd zw = z.array() * sqrtw.array();
  MatrixXd Xw(n, k);
  for (int a = 0; a < k; ++a) {
    const int pos = active[a];
    Xw.col(a) = X.col(vars[pos]).array() * sqrtw.array() / scaling[pos];
  }
  VectorXd penalty(k), center(k);
  for (int a = 0; a < k; ++a) {
    penalty[a] = lambdas[active[a]];
    center[a] = prior_mean[active[a]];
  }
  if (k > n && ((penalty.array() > 0.0).all() || (penalty.array() == 0.0).all())) {
    const bool penalized = (penalty.array() > 0.0).all();
    VectorXd invsqrt = VectorXd::Ones(k);
    if (penalized) invsqrt = penalty.array().sqrt().inverse();
    else center.setZero();
    MatrixXd Z = Xw * invsqrt.asDiagonal();
    MatrixXd dual = MatrixXd::Zero(n, n);
    const double one = 1.0, zero = 0.0;
    F77_CALL(dsyrk)("L", "N", &n, &k, &one, Z.data(), &n, &zero,
                    dual.data(), &n FCONE FCONE);
    if (penalized) dual.diagonal().array() += 1.0;
    Eigen::LLT<MatrixXd> solver(dual);
    if (solver.info() != Eigen::Success) return start;
    VectorXd u = solver.solve(zw - Xw * center);
    if (solver.info() != Eigen::Success || !u.allFinite()) return start;
    VectorXd sol = center + invsqrt.asDiagonal() * (Z.transpose() * u);
    sol = constrain_vec(sol);
    for (int a = 0; a < k; ++a) sol[a] /= scaling[active[a]];
    return sol;
  }
  MatrixXd A = MatrixXd::Zero(k, k);
  const double one = 1.0;
  const double zero = 0.0;
  // LLT reads the lower triangle only. BLAS computes it without the redundant
  // upper triangle, and uses the optimized BLAS selected by the R installation.
  F77_CALL(dsyrk)("L", "T", &k, &n, &one, Xw.data(), &n, &zero,
                  A.data(), &k FCONE FCONE);
  VectorXd b = Xw.transpose() * zw;
  for (int a = 0; a < k; ++a) {
    const int pos = active[a];
    A(a, a) += lambdas[pos];
    b[a] += lambdas[pos] * prior_mean[pos];
  }

  Eigen::LLT<MatrixXd> solver(A);
  if (solver.info() != Eigen::Success) return start;
  VectorXd sol = solver.solve(b);
  if (solver.info() != Eigen::Success) return start;

  sol = constrain_vec(sol);
  for (int a = 0; a < k; ++a) sol[a] /= scaling[active[a]];
  return sol;
}

inline VectorXd solve_ridge_sparse(const Eigen::MappedSparseMatrix<double>& X,
                                   const VectorXd& z,
                                   const VectorXd& w,
                                   const std::vector<int>& vars,
                                   const std::vector<int>& active,
                                   const VectorXd& lambdas,
                                   const VectorXd& prior_mean,
                                   const VectorXd& scaling,
                                   const VectorXd& start,
                                   bool natural_ordering) {
  const int n = X.rows();
  const int k = static_cast<int>(active.size());
  VectorXd sqrtw(n);
  VectorXd zw(n);
  for (int i = 0; i < n; ++i) {
    sqrtw[i] = (w[i] > 0.0) ? std::sqrt(w[i]) : 0.0;
    zw[i] = z[i] * sqrtw[i];
  }

  std::vector<Triplet<double> > triplets;
  for (int a = 0; a < k; ++a) {
    const int pos = active[a];
    const int col = vars[pos];
    const double scale = scaling[pos];
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(X, col); it; ++it) {
      const int row = it.row();
      const double sw = sqrtw[row];
      if (sw == 0.0) continue;
      triplets.emplace_back(row, a, it.value() * sw / scale);
    }
  }

  SparseMatrix<double> Xw(n, k);
  Xw.setFromTriplets(triplets.begin(), triplets.end());
  Xw.makeCompressed();

  VectorXd penalty(k), center(k);
  for (int a = 0; a < k; ++a) {
    penalty[a] = lambdas[active[a]];
    center[a] = prior_mean[active[a]];
  }
  if (k > n && ((penalty.array() > 0.0).all() || (penalty.array() == 0.0).all())) {
    const bool penalized = (penalty.array() > 0.0).all();
    VectorXd invsqrt = VectorXd::Ones(k);
    if (penalized) invsqrt = penalty.array().sqrt().inverse();
    else center.setZero();
    SparseMatrix<double> Z = Xw * invsqrt.asDiagonal();
    SparseMatrix<double> dual = Z * Z.transpose();
    if (penalized) for (int i = 0; i < n; ++i) dual.coeffRef(i, i) += 1.0;
    dual.makeCompressed();
    VectorXd rhs = zw - Xw * center;
    VectorXd u;
    if (natural_ordering) {
      Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::NaturalOrdering<int> > solver;
      solver.compute(dual);
      if (solver.info() != Eigen::Success) return start;
      u = solver.solve(rhs);
      if (solver.info() != Eigen::Success) return start;
    } else {
      Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::AMDOrdering<int> > solver;
      solver.compute(dual);
      if (solver.info() != Eigen::Success) return start;
      u = solver.solve(rhs);
      if (solver.info() != Eigen::Success) return start;
    }
    if (!u.allFinite()) return start;
    VectorXd sol = center + invsqrt.asDiagonal() * (Z.transpose() * u);
    sol = constrain_vec(sol);
    for (int a = 0; a < k; ++a) sol[a] /= scaling[active[a]];
    return sol;
  }

  SparseMatrix<double> A = Xw.transpose() * Xw;
  VectorXd b = Xw.transpose() * zw;
  for (int a = 0; a < k; ++a) {
    const int pos = active[a];
    A.coeffRef(a, a) += lambdas[pos];
    b[a] += lambdas[pos] * prior_mean[pos];
  }
  A.makeCompressed();

  VectorXd sol;
  if (natural_ordering) {
    Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::NaturalOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) return start;
    sol = solver.solve(b);
    if (solver.info() != Eigen::Success) return start;
  } else {
    Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::AMDOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) return start;
    sol = solver.solve(b);
    if (solver.info() != Eigen::Success) return start;
  }

  sol = constrain_vec(sol);
  for (int a = 0; a < k; ++a) sol[a] /= scaling[active[a]];
  return sol;
}

template <typename Design>
FitResult fit_channel_impl(const Design& X, const MatrixXd& Y,
                           const MatrixXd& prior_weights,
                           const MatrixXd& Offset,
                           const MatrixXd& prior_mean_full,
                           const MatrixXd* start_mat,
                           int channel,
                           const std::vector<int>& vars,
                           const std::vector<char>& always_include,
                           double lambda,
                           bool normalize,
                           int maxit,
                           double gamma,
                           double thresh,
                           double tol,
                           bool natural_ordering,
                           bool hard_threshold);

template <>
FitResult fit_channel_impl<Map<MatrixXd> >(const Map<MatrixXd>& X,
                                           const MatrixXd& Y,
                                           const MatrixXd& prior_weights,
                                           const MatrixXd& Offset,
                                           const MatrixXd& prior_mean_full,
                                           const MatrixXd* start_mat,
                                           int channel,
                                           const std::vector<int>& vars,
                                           const std::vector<char>& always_include,
                                           double lambda,
                                           bool normalize,
                                           int maxit,
                                           double gamma,
                                           double thresh,
                                           double tol,
                                           bool natural_ordering,
                                           bool hard_threshold) {
  (void) natural_ordering;
  const int n = X.rows();
  const int q = static_cast<int>(vars.size());
  const double eps = std::numeric_limits<double>::epsilon();

  VectorXd y = Y.col(channel);
  VectorXd prior_w = prior_weights.col(channel);
  VectorXd offset = Offset.col(channel);
  VectorXd prior_mean(q);
  for (int j = 0; j < q; ++j) {
    const int col = vars[j];
    prior_mean[j] = prior_mean_full(col, channel);
  }

  VectorXd scaling = normalize ? compute_scaling_dense(X, vars) :
    VectorXd::Ones(q);

  VectorXd beta(q);
  bool have_start = start_mat != NULL;
  if (have_start) {
    bool any_bad = false;
    bool any_nonzero = false;
    for (int j = 0; j < q; ++j) {
      beta[j] = (*start_mat)(j, channel);
      any_bad = any_bad || !std::isfinite(beta[j]);
      any_nonzero = any_nonzero || beta[j] != 0.0;
    }
    have_start = !any_bad && any_nonzero;
  }
  if (!have_start) {
    beta = init_start_dense(X, y, prior_w, offset, prior_mean, vars, lambda);
  }
  beta = constrain_vec(beta);

  VectorXd mu = y.array() + 0.1;
  VectorXd beta_old = beta;
  VectorXd weights = VectorXd::Zero(n);
  VectorXd z = y - offset;
  VectorXd lambdas = VectorXd::Constant(q, R_PosInf);

  for (int iter = 0; iter < maxit; ++iter) {
    for (int i = 0; i < n; ++i) {
      const double wi = prior_w[i] / mu[i];
      weights[i] = finite_pos(wi) ? wi : 0.0;
    }
    z = y - offset;

    for (int j = 0; j < q; ++j) {
      lambdas[j] = (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / std::pow(std::abs(beta[j]), gamma);
    }

    std::vector<int> active = active_positions(beta);
    if (active.empty()) break;
    bool bad_beta = false;
    for (int j = 0; j < q; ++j) bad_beta = bad_beta || !std::isfinite(beta[j]);
    if (bad_beta) break;

    VectorXd start_active(active.size());
    for (int a = 0; a < static_cast<int>(active.size()); ++a) {
      const int pos = active[a];
      start_active[a] = beta[pos] * scaling[pos];
    }

    VectorXd beta_sub = solve_ridge_dense(X, z, weights, vars, active, lambdas,
                                          prior_mean, scaling,
                                          start_active);
    for (int a = 0; a < static_cast<int>(active.size()); ++a) {
      beta[active[a]] = beta_sub[a];
    }

    const double relative_cutoff = hard_threshold ? beta.cwiseAbs().maxCoeff() * thresh : 0.0;
    for (int j = 0; j < q; ++j) {
      if (hard_threshold && (std::abs(beta[j]) < relative_cutoff || std::abs(beta[j]) < 1e-20)) {
        beta[j] = 0.0;
      }
      if (always_include[vars[j]]) prior_mean[j] = beta[j];
    }

    VectorXd eta = predict_dense(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    const double rel = (beta - beta_old).norm() /
      std::max(1e-12, beta.norm());
    if (rel < tol) break;
    beta_old = beta;
  }

  FitResult out;
  out.beta = beta;
  out.weights = weights;
  out.z = z;
  out.lambdas = lambdas;
  return out;
}

template <>
FitResult fit_channel_impl<Eigen::MappedSparseMatrix<double> >(
    const Eigen::MappedSparseMatrix<double>& X,
    const MatrixXd& Y,
    const MatrixXd& prior_weights,
    const MatrixXd& Offset,
    const MatrixXd& prior_mean_full,
    const MatrixXd* start_mat,
    int channel,
    const std::vector<int>& vars,
    const std::vector<char>& always_include,
    double lambda,
    bool normalize,
    int maxit,
    double gamma,
    double thresh,
    double tol,
    bool natural_ordering,
    bool hard_threshold) {
  const int n = X.rows();
  const int q = static_cast<int>(vars.size());
  const double eps = std::numeric_limits<double>::epsilon();

  VectorXd y = Y.col(channel);
  VectorXd prior_w = prior_weights.col(channel);
  VectorXd offset = Offset.col(channel);
  VectorXd prior_mean(q);
  for (int j = 0; j < q; ++j) {
    const int col = vars[j];
    prior_mean[j] = prior_mean_full(col, channel);
  }

  VectorXd scaling = normalize ? compute_scaling_sparse(X, vars) :
    VectorXd::Ones(q);

  VectorXd beta(q);
  bool have_start = start_mat != NULL;
  if (have_start) {
    bool any_bad = false;
    bool any_nonzero = false;
    for (int j = 0; j < q; ++j) {
      beta[j] = (*start_mat)(j, channel);
      any_bad = any_bad || !std::isfinite(beta[j]);
      any_nonzero = any_nonzero || beta[j] != 0.0;
    }
    have_start = !any_bad && any_nonzero;
  }
  if (!have_start) {
    beta = init_start_sparse(X, y, prior_w, offset, prior_mean, vars, lambda);
  }
  beta = constrain_vec(beta);

  VectorXd mu = y.array() + 0.1;
  VectorXd beta_old = beta;
  VectorXd weights = VectorXd::Zero(n);
  VectorXd z = y - offset;
  VectorXd lambdas = VectorXd::Constant(q, R_PosInf);

  for (int iter = 0; iter < maxit; ++iter) {
    for (int i = 0; i < n; ++i) {
      const double wi = prior_w[i] / mu[i];
      weights[i] = finite_pos(wi) ? wi : 0.0;
    }
    z = y - offset;

    for (int j = 0; j < q; ++j) {
      lambdas[j] = (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / std::pow(std::abs(beta[j]), gamma);
    }

    std::vector<int> active = active_positions(beta);
    if (active.empty()) break;
    bool bad_beta = false;
    for (int j = 0; j < q; ++j) bad_beta = bad_beta || !std::isfinite(beta[j]);
    if (bad_beta) break;

    VectorXd start_active(active.size());
    for (int a = 0; a < static_cast<int>(active.size()); ++a) {
      const int pos = active[a];
      start_active[a] = beta[pos] * scaling[pos];
    }

    VectorXd beta_sub = solve_ridge_sparse(X, z, weights, vars, active, lambdas,
                                           prior_mean, scaling,
                                           start_active, natural_ordering);
    for (int a = 0; a < static_cast<int>(active.size()); ++a) {
      beta[active[a]] = beta_sub[a];
    }

    const double relative_cutoff = hard_threshold ? beta.cwiseAbs().maxCoeff() * thresh : 0.0;
    for (int j = 0; j < q; ++j) {
      if (hard_threshold && (std::abs(beta[j]) < relative_cutoff || std::abs(beta[j]) < 1e-20)) {
        beta[j] = 0.0;
      }
      if (always_include[vars[j]]) prior_mean[j] = beta[j];
    }

    VectorXd eta = predict_sparse(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    const double rel = (beta - beta_old).norm() /
      std::max(1e-12, beta.norm());
    if (rel < tol) break;
    beta_old = beta;
  }

  FitResult out;
  out.beta = beta;
  out.weights = weights;
  out.z = z;
  out.lambdas = lambdas;
  return out;
}

template <typename Design>
Rcpp::List fit_all_impl(const Design& X,
                        const Map<MatrixXd>& Ymap,
                        const Map<MatrixXd>& Wmap,
                        const Map<MatrixXd>& Offsetmap,
                        const Map<MatrixXd>& Priormap,
                        Rcpp::Nullable<Rcpp::NumericMatrix> start_nullable,
                        Rcpp::IntegerVector always_include_r,
                        double lambda,
                        bool normalize,
                        bool grouped,
                        bool relaxed,
                        int group_maxit,
                        int maxit,
                        double gamma,
                        double thresh,
                        double tol,
                        int nthreads,
                        bool natural_ordering) {
  const int n = X.rows();
  const int p = X.cols();
  const int m = Ymap.cols();

  MatrixXd Y = Ymap;
  MatrixXd prior_weights = Wmap;
  MatrixXd Offset = Offsetmap;
  MatrixXd prior_mean_full = Priormap;

  MatrixXd start_full;
  bool has_start = start_nullable.isNotNull();
  if (has_start) {
    Rcpp::NumericMatrix sr(start_nullable);
    start_full = Rcpp::as<Map<MatrixXd> >(sr);
  }

  std::vector<char> always_include(p, 0);
  for (int i = 0; i < always_include_r.size(); ++i) {
    const int idx = always_include_r[i] - 1;
    if (idx >= 0 && idx < p) always_include[idx] = 1;
  }

  const int threads =
#ifdef _OPENMP
    std::max(1, nthreads);
#else
    1;
#endif

  const int effective_group_maxit = (grouped && m > 1) ?
    std::max(1, group_maxit) : 1;
  std::vector<int> vars(p);
  std::iota(vars.begin(), vars.end(), 0);

  MatrixXd coef = MatrixXd::Zero(p, m);
  MatrixXd lambdas = MatrixXd::Zero(p, m);
  MatrixXd weights_out = MatrixXd::Zero(n, m);
  MatrixXd z_out = MatrixXd::Zero(n, m);
  MatrixXd start_current;
  MatrixXd* start_ptr = has_start ? &start_full : NULL;

  auto fit_for_vars = [&](double lambda_fit, MatrixXd* start_for_fit,
                          bool hard_threshold) {
    const int q = static_cast<int>(vars.size());
    std::vector<FitResult> fits(m);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads) if(threads > 1 && m > 1)
#endif
    for (int col = 0; col < m; ++col) {
      fits[col] = fit_channel_impl<Design>(
        X, Y, prior_weights, Offset, prior_mean_full, start_for_fit, col, vars,
        always_include, lambda_fit, normalize, maxit,
        gamma, thresh, tol, natural_ordering, hard_threshold);
    }

    coef.setZero();
    lambdas.setZero();
    for (int col = 0; col < m; ++col) {
      for (int a = 0; a < q; ++a) {
        coef(vars[a], col) = fits[col].beta[a];
        lambdas(vars[a], col) = fits[col].lambdas[a];
      }
      weights_out.col(col) = fits[col].weights;
      z_out.col(col) = fits[col].z;
    }
  };

  for (int iter = 0; iter < effective_group_maxit; ++iter) {
    if (iter > 0) {
      group_support_start(coef, vars, start_current);
      if (vars.empty()) break;
      start_ptr = &start_current;
    }
    fit_for_vars(lambda, start_ptr, true);
  }

  if (relaxed) {
    group_support_start(coef, vars, start_current);
    if (!vars.empty()) {
      fit_for_vars(1e-30, &start_current, false);
    }
  }

  Rcpp::NumericVector dispersion(m, 1.0);
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coef,
    Rcpp::Named("weights") = weights_out,
    Rcpp::Named("dispersion") = dispersion,
    Rcpp::Named("lambdas") = lambdas,
    Rcpp::Named("z") = z_out,
    Rcpp::Named("openmp") =
#ifdef _OPENMP
      true,
#else
      false,
#endif
    Rcpp::Named("nthreads") = threads);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List L0nnpois_dense_cpp(const Eigen::Map<Eigen::MatrixXd>& X,
                              const Eigen::Map<Eigen::MatrixXd>& y,
                              const Eigen::Map<Eigen::MatrixXd>& weights,
                              const Eigen::Map<Eigen::MatrixXd>& offset,
                              const Eigen::Map<Eigen::MatrixXd>& prior_mean,
                              Rcpp::Nullable<Rcpp::NumericMatrix> start,
                              Rcpp::IntegerVector always_include,
                              double lambda,
                              bool normalize,
                              bool grouped,
                              bool relaxed,
                              int group_maxit,
                              int maxit,
                              double gamma,
                              double thresh,
                              double tol,
                              int nthreads,
                              bool natural_ordering) {
  return fit_all_impl<Map<MatrixXd> >(
    X, y, weights, offset, prior_mean, start, always_include,
    lambda, normalize, grouped, relaxed, group_maxit, maxit, gamma,
    thresh, tol, nthreads, natural_ordering);
}

// [[Rcpp::export]]
Rcpp::List L0nnpois_sparse_cpp(const Eigen::MappedSparseMatrix<double>& X,
                               const Eigen::Map<Eigen::MatrixXd>& y,
                               const Eigen::Map<Eigen::MatrixXd>& weights,
                               const Eigen::Map<Eigen::MatrixXd>& offset,
                               const Eigen::Map<Eigen::MatrixXd>& prior_mean,
                               Rcpp::Nullable<Rcpp::NumericMatrix> start,
                               Rcpp::IntegerVector always_include,
                               double lambda,
                               bool normalize,
                               bool grouped,
                               bool relaxed,
                               int group_maxit,
                               int maxit,
                               double gamma,
                               double thresh,
                               double tol,
                               int nthreads,
                               bool natural_ordering) {
  return fit_all_impl<Eigen::MappedSparseMatrix<double> >(
    X, y, weights, offset, prior_mean, start, always_include,
    lambda, normalize, grouped, relaxed, group_maxit, maxit, gamma,
    thresh, tol, nthreads, natural_ordering);
}
