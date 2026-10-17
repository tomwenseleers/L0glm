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
#include "ridge_profile.h"

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

inline std::vector<int> active_positions(const VectorXd& beta,
                                          const std::vector<int>& vars,
                                          const std::vector<char>& always_include) {
  std::vector<int> active;
  active.reserve(beta.size());
  for (int i = 0; i < beta.size(); ++i) {
    if (beta[i] != 0.0 || always_include[vars[i]]) active.push_back(i);
  }
  return active;
}

inline void group_support_start(const MatrixXd& coef,
                                std::vector<int>& vars,
                                MatrixXd& start_current,
                                const std::vector<char>& always_include) {
  const int p = coef.rows();
  const int m = coef.cols();
  vars.clear();
  std::vector<double> starts;
  starts.reserve(p);

  for (int j = 0; j < p; ++j) {
    bool active = always_include[j];
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

inline double coefficient_change(const VectorXd& beta, const VectorXd& previous,
                                  const VectorXd& scaling, const std::vector<int>& vars,
                                  const std::vector<char>& mandatory,
                                  const VectorXd& prediction, const VectorXd& old_prediction,
                                  const VectorXd& weights, bool normalize) {
  if (!normalize) return (beta - previous).norm() / std::max(1e-12, beta.norm());
  double change = 0.0, magnitude = 0.0;
  bool nuisance = false;
  for (int j = 0; j < beta.size(); ++j) {
    if (mandatory[vars[j]]) { nuisance = true; continue; }
    change += std::pow((beta[j] - previous[j]) * scaling[j], 2);
    magnitude += std::pow(beta[j] * scaling[j], 2);
  }
  double relative = std::sqrt(change) / std::max(1e-12, std::sqrt(magnitude));
  if (!nuisance) return relative;
  const double total = weights.sum();
  const double mean = total > 0.0 ? weights.dot(prediction) / total : 0.0;
  const double fitted_change = std::sqrt((weights.array() * (prediction - old_prediction).array().square()).sum());
  const double fitted_size = std::sqrt((weights.array() * (prediction.array() - mean).square()).sum());
  return std::max(relative, fitted_change / std::max(1e-12, std::max(std::sqrt(total), fitted_size)));
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
                                 const VectorXd& penalty) {
  const int n = X.rows();
  const int q = static_cast<int>(vars.size());
  VectorXd numer = VectorXd::Zero(q);
  VectorXd denom = penalty;
  for (int j = 0; j < q; ++j) numer[j] = penalty[j] * prior_mean[j];

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
                                  const VectorXd& penalty) {
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
    numer[j] = penalty[j] * prior_mean[j];
    denom[j] = penalty[j];
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
  VectorXd penalty(k), center(k), previous = start;
  for (int a = 0; a < k; ++a) {
    const double scale = scaling[active[a]];
    penalty[a] = lambdas[active[a]] / (scale * scale);
    center[a] = prior_mean[active[a]] * scale;
    previous[a] /= scale;
  }
  if (k > n && (penalty.array() >= 0.0).all() && ((penalty.array() == 0.0).all() || (penalty.array() == 0.0).count() <= n)) {
    VectorXd sol;
    if (!L0glm::ridge_dual(Xw, zw, penalty, center, sol, false)) return previous;
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
    A(a, a) += penalty[a];
    b[a] += penalty[a] * center[a];
  }

  Eigen::LLT<MatrixXd> solver(A);
  if (solver.info() != Eigen::Success) return previous;
  VectorXd sol = solver.solve(b);
  if (solver.info() != Eigen::Success) return previous;

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

  VectorXd penalty(k), center(k), previous = start;
  for (int a = 0; a < k; ++a) {
    const double scale = scaling[active[a]];
    penalty[a] = lambdas[active[a]] / (scale * scale);
    center[a] = prior_mean[active[a]] * scale;
    previous[a] /= scale;
  }
  if (k > n && (penalty.array() >= 0.0).all() && ((penalty.array() == 0.0).all() || (penalty.array() == 0.0).count() <= n)) {
    VectorXd sol;
    if (!L0glm::ridge_dual(Xw, zw, penalty, center, sol, natural_ordering)) return previous;
    sol = constrain_vec(sol);
    for (int a = 0; a < k; ++a) sol[a] /= scaling[active[a]];
    return sol;
  }
  SparseMatrix<double> A = Xw.transpose() * Xw;
  VectorXd b = Xw.transpose() * zw;
  for (int a = 0; a < k; ++a) {
    const int pos = active[a];
    A.coeffRef(a, a) += penalty[a];
    b[a] += penalty[a] * center[a];
  }
  A.makeCompressed();

  VectorXd sol;
  if (natural_ordering) {
    Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::NaturalOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) return previous;
    sol = solver.solve(b);
    if (solver.info() != Eigen::Success) return previous;
  } else {
    Eigen::SimplicialLLT<SparseMatrix<double>, Lower, Eigen::AMDOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) return previous;
    sol = solver.solve(b);
    if (solver.info() != Eigen::Success) return previous;
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
                           int normalize,
                           const VectorXd& scaling_full,
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
                                           int normalize,
                                           const VectorXd& scaling_full,
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

  VectorXd scaling(q), initial_penalty = VectorXd::Constant(q, lambda);
  for (int j = 0; j < q; ++j) {
    scaling[j] = scaling_full[vars[j]];
    if (normalize) initial_penalty[j] = always_include[vars[j]] ? 0.0 : lambda * scaling[j] * scaling[j];
  }

  VectorXd beta(q);
  bool have_start = start_mat != NULL;
  if (have_start) {
    bool any_bad = false;

    for (int j = 0; j < q; ++j) {
      beta[j] = (*start_mat)(j, channel);
      any_bad = any_bad || !std::isfinite(beta[j]);

    }
    have_start = !any_bad;
  }
  if (!have_start) {
    beta = init_start_dense(X, y, prior_w, offset, prior_mean, vars, initial_penalty);
  }
  beta = constrain_vec(beta);

  VectorXd mu = y.array() + 0.1;
  VectorXd beta_old = beta;
  VectorXd prediction_old = normalize ? predict_dense(X, VectorXd::Zero(n), vars, beta) : VectorXd();
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
      lambdas[j] = always_include[vars[j]] ? 0.0 : (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / std::pow(std::abs(beta[j]), gamma);
    }

    std::vector<int> active = active_positions(beta, vars, always_include);
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

    double maximum = 0.0;
    for (int j = 0; j < q; ++j) if (!always_include[vars[j]])
      maximum = std::max(maximum, std::abs(beta[j]) * scaling[j]);
    const double relative_cutoff = maximum * thresh;
    for (int j = 0; j < q; ++j) {
      const double magnitude = std::abs(beta[j]) * scaling[j];
      if (hard_threshold && !always_include[vars[j]] && (magnitude < relative_cutoff || magnitude < 1e-20)) {
        beta[j] = 0.0;
      }

    }

    VectorXd eta = predict_dense(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    const VectorXd prediction = normalize ? (eta - offset).eval() : VectorXd();
    const double rel = coefficient_change(beta, beta_old, scaling, vars, always_include,
                                          prediction, prediction_old, prior_w, normalize);
    if (rel < tol) break;
    beta_old = beta;
    if (normalize) prediction_old = prediction;
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
    int normalize,
    const VectorXd& scaling_full,
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

  VectorXd scaling(q), initial_penalty = VectorXd::Constant(q, lambda);
  for (int j = 0; j < q; ++j) {
    scaling[j] = scaling_full[vars[j]];
    if (normalize) initial_penalty[j] = always_include[vars[j]] ? 0.0 : lambda * scaling[j] * scaling[j];
  }

  VectorXd beta(q);
  bool have_start = start_mat != NULL;
  if (have_start) {
    bool any_bad = false;

    for (int j = 0; j < q; ++j) {
      beta[j] = (*start_mat)(j, channel);
      any_bad = any_bad || !std::isfinite(beta[j]);

    }
    have_start = !any_bad;
  }
  if (!have_start) {
    beta = init_start_sparse(X, y, prior_w, offset, prior_mean, vars, initial_penalty);
  }
  beta = constrain_vec(beta);

  VectorXd mu = y.array() + 0.1;
  VectorXd beta_old = beta;
  VectorXd prediction_old = normalize ? predict_sparse(X, VectorXd::Zero(n), vars, beta) : VectorXd();
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
      lambdas[j] = always_include[vars[j]] ? 0.0 : (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / std::pow(std::abs(beta[j]), gamma);
    }

    std::vector<int> active = active_positions(beta, vars, always_include);
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

    double maximum = 0.0;
    for (int j = 0; j < q; ++j) if (!always_include[vars[j]])
      maximum = std::max(maximum, std::abs(beta[j]) * scaling[j]);
    const double relative_cutoff = maximum * thresh;
    for (int j = 0; j < q; ++j) {
      const double magnitude = std::abs(beta[j]) * scaling[j];
      if (hard_threshold && !always_include[vars[j]] && (magnitude < relative_cutoff || magnitude < 1e-20)) {
        beta[j] = 0.0;
      }

    }

    VectorXd eta = predict_sparse(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    const VectorXd prediction = normalize ? (eta - offset).eval() : VectorXd();
    const double rel = coefficient_change(beta, beta_old, scaling, vars, always_include,
                                          prediction, prediction_old, prior_w, normalize);
    if (rel < tol) break;
    beta_old = beta;
    if (normalize) prediction_old = prediction;
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
                        int normalize,
                        const Eigen::Map<Eigen::VectorXd>& scaling_full,
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
  const VectorXd scaling = scaling_full;

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
  std::vector<int> vars;
  for (int j = 0; j < p; ++j) {
    if (!has_start || always_include[j] || (start_full.row(j).array() != 0.0).any())
      vars.push_back(j);
  }

  MatrixXd coef = MatrixXd::Zero(p, m);
  MatrixXd lambdas = MatrixXd::Zero(p, m);
  MatrixXd weights_out = MatrixXd::Zero(n, m);
  MatrixXd z_out = MatrixXd::Zero(n, m);
  MatrixXd start_current;
  if (has_start) {
    start_current.resize(vars.size(), m);
    for (int a = 0; a < static_cast<int>(vars.size()); ++a)
      start_current.row(a) = start_full.row(vars[a]);
  }
  MatrixXd* start_ptr = has_start ? &start_current : NULL;

  auto preserve_start_mask = [&]() {
    if (has_start) for (int a = 0; a < static_cast<int>(vars.size()); ++a)
      for (int col = 0; col < m; ++col)
        if (!always_include[vars[a]] && start_full(vars[a], col) == 0.0)
          start_current(a, col) = 0.0;
  };

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
        always_include, lambda_fit, normalize, scaling, maxit,
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
      group_support_start(coef, vars, start_current, always_include);
      preserve_start_mask();
      if (vars.empty()) break;
      start_ptr = &start_current;
    }
    fit_for_vars(lambda, start_ptr, true);
  }

  if (relaxed) {
    group_support_start(coef, vars, start_current, always_include);
    preserve_start_mask();
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
                              int normalize,
                              const Eigen::Map<Eigen::VectorXd>& scaling,
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
    lambda, normalize, scaling, grouped, relaxed, group_maxit, maxit, gamma,
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
                               int normalize,
                               const Eigen::Map<Eigen::VectorXd>& scaling,
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
    lambda, normalize, scaling, grouped, relaxed, group_maxit, maxit, gamma,
    thresh, tol, nthreads, natural_ordering);
}
