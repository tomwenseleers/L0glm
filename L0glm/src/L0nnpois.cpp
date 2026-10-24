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
#include "poisson_validation.h"
#include "boundary.h"

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
  int iterations = 0;
  bool converged = false;
  bool step_failed = false;
  PoissonValidationTrace validation;
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

inline double poisson_deviance(const VectorXd& y, const VectorXd& eta,
                                const VectorXd& weights) {
  const double eps = std::numeric_limits<double>::epsilon();
  double value = 0.0;
  for (int i = 0; i < y.size(); ++i) if (weights[i] > 0.0) {
    if (!std::isfinite(eta[i]) || eta[i] < 0.0) return R_PosInf;
    const double mu = std::max(eta[i], eps);
    value += 2.0 * weights[i] * (y[i] == 0.0 ? mu : y[i] * std::log(y[i] / mu) - (y[i] - mu));
  }
  return value;
}

inline double ridge_objective(const VectorXd& beta, const VectorXd& center,
                                const VectorXd& penalty) {
  double value = 0.0;
  for (int j = 0; j < beta.size(); ++j) if (std::isfinite(penalty[j]) && penalty[j] > 0.0)
    value += penalty[j] * (beta[j] - center[j]) * (beta[j] - center[j]);
  return value;
}

// Same support-aware rule as the general R loop. One solve proposes a step;
// damping only interpolates coefficients/predictors and checks the frozen target.
template <typename Predict>
bool poisson_step(VectorXd& beta, const VectorXd& previous, VectorXd& eta,
                    const VectorXd& previous_eta, const VectorXd& y, const VectorXd& weights,
                    const VectorXd& penalty, const VectorXd& center,
                    double& previous_deviance, bool check_decrease,
                    const Predict& predict) {
  const double deviance = poisson_deviance(y, eta, weights);
  const double proposed_value = deviance + ridge_objective(beta, center, penalty);
  if (!check_decrease && std::isfinite(proposed_value)) { previous_deviance = deviance; return true; }
  VectorXd anchor = previous;
  VectorXd old_eta = previous_eta.size() ? previous_eta : predict(previous, false);
  if (!std::isfinite(previous_deviance) && !std::isfinite(poisson_deviance(y, old_eta, weights))) {
    const VectorXd null_eta = predict(VectorXd::Zero(beta.size()), false);
    if (std::isfinite(poisson_deviance(y, null_eta, weights))) { anchor.setZero(); old_eta = null_eta; }
  }
  const double baseline = check_decrease ? previous_deviance +
    ridge_objective(previous, center, penalty) : R_PosInf;
  const double limit = std::isfinite(baseline) ? baseline + 1e-10 * std::max(1.0, std::abs(baseline)) : R_PosInf;
  if (std::isfinite(proposed_value) && proposed_value <= limit) { previous_deviance = deviance; return true; }
  const VectorXd proposal = beta, proposed_eta = eta;
  double fraction = 1.0;
  for (int halving = 1; halving <= 50; ++halving) {
    fraction *= 0.5;
    beta = anchor + fraction * (proposal - anchor);
    eta = old_eta + fraction * (proposed_eta - old_eta);
    const double trial_deviance = poisson_deviance(y, eta, weights);
    const double value = trial_deviance + ridge_objective(beta, center, penalty);
    if (std::isfinite(value) && value <= limit) { previous_deviance = trial_deviance; return true; }
  }
  beta = anchor; eta = old_eta;
  return false;
}

// Thresholding is a discrete selection step after the smooth update. Retain
// the smooth point if deleting tiny terms would leave the mean domain.
template <typename Predict>
void poisson_threshold(VectorXd& beta, VectorXd& eta, double& deviance,
                         const VectorXd& y, const VectorXd& weights, const VectorXd& scaling,
                         const std::vector<int>& vars, const std::vector<char>& mandatory,
                         double thresh, const Predict& predict) {
  double maximum = 0.0;
  for (int j = 0; j < beta.size(); ++j) if (!mandatory[vars[j]])
    maximum = std::max(maximum, std::abs(beta[j]) * scaling[j]);
  VectorXd removed = VectorXd::Zero(beta.size());
  bool any_removed = false;
  for (int j = 0; j < beta.size(); ++j) if (!mandatory[vars[j]]) {
    const double magnitude = std::abs(beta[j]) * scaling[j];
    if ((magnitude < maximum * thresh || magnitude < 1e-20) && beta[j] != 0.0) {
      removed[j] = beta[j]; beta[j] = 0.0; any_removed = true;
    }
  }
  if (!any_removed) return;
  const VectorXd trial_eta = eta - predict(removed, true);
  const double trial_deviance = poisson_deviance(y, trial_eta, weights);
  if (std::isfinite(trial_deviance)) { eta = trial_eta; deviance = trial_deviance; } else beta += removed;
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
                           double thresh,
                           double tol,
                           bool natural_ordering,
                           bool hard_threshold,
                           const PoissonValidation* validation);

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
                                           double thresh,
                                           double tol,
                                           bool natural_ordering,
                                           bool hard_threshold,
                                           const PoissonValidation* validation) {
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
  VectorXd eta_previous;
  double previous_deviance = R_PosInf;
  VectorXd prediction_old = normalize ? predict_dense(X, VectorXd::Zero(n), vars, beta) : VectorXd();
  VectorXd weights = VectorXd::Zero(n);
  VectorXd z = y - offset;
  VectorXd lambdas = VectorXd::Constant(q, R_PosInf);
  int iterations = 0;
  bool converged = false;
  bool step_failed = false;
  PoissonCheckpoint checkpoint;
  if (validation) checkpoint.initialize(*validation, vars, always_include, channel, lambda);

  for (int iter = 0; iter < maxit; ++iter) {
    iterations = iter + 1;
    for (int i = 0; i < n; ++i) {
      const double wi = prior_w[i] / mu[i];
      weights[i] = finite_pos(wi) ? wi : 0.0;
    }
    z = y - offset;

    for (int j = 0; j < q; ++j) {
      lambdas[j] = always_include[vars[j]] ? 0.0 : (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / (beta[j] * beta[j]);
    }

    std::vector<int> active = active_positions(beta, vars, always_include);
    if (active.empty()) { converged = true; break; }
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

    VectorXd eta = predict_dense(X, offset, vars, beta);
    bool support_changed = false;
    for (int j = 0; j < q; ++j) if (!always_include[vars[j]] && ((beta[j] != 0.0) != (beta_old[j] != 0.0)))
      { support_changed = true; break; }
    auto predict = [&](const VectorXd& b, bool zero_offset) {
      return predict_dense(X, zero_offset ? VectorXd::Zero(n).eval() : offset, vars, b);
    };
    if (!poisson_step(beta, beta_old, eta, eta_previous, y, prior_w, lambdas, prior_mean,
        previous_deviance, iter > 0 && !support_changed, predict)) {
      step_failed = true; break;
    }
    if (hard_threshold) poisson_threshold(beta, eta, previous_deviance, y, prior_w,
      scaling, vars, always_include, thresh, predict);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    if (validation && checkpoint.update(*validation, beta, lambdas, vars, channel, iter + 1)) break;
    const VectorXd prediction = normalize ? (eta - offset).eval() : VectorXd();
    const double rel = coefficient_change(beta, beta_old, scaling, vars, always_include,
                                          prediction, prediction_old, prior_w, normalize);
    if (rel < tol && iter > 0) { converged = true; break; }
    beta_old = beta;
    eta_previous = eta;
    if (normalize) prediction_old = prediction;
  }

  if (validation) {
    beta = checkpoint.beta; lambdas = checkpoint.lambdas;
    const VectorXd restored = predict_dense(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) weights[i] = prior_w[i] == 0.0 ? 0.0 : prior_w[i] / restored[i];
  }
  FitResult out;
  out.validation = checkpoint.trace;
  out.beta = beta;
  out.weights = weights;
  out.z = z;
  out.lambdas = lambdas;
  out.iterations = iterations;
  out.converged = converged;
  out.step_failed = step_failed;
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
    double thresh,
    double tol,
    bool natural_ordering,
    bool hard_threshold,
                                           const PoissonValidation* validation) {
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
  VectorXd eta_previous;
  double previous_deviance = R_PosInf;
  VectorXd prediction_old = normalize ? predict_sparse(X, VectorXd::Zero(n), vars, beta) : VectorXd();
  VectorXd weights = VectorXd::Zero(n);
  VectorXd z = y - offset;
  VectorXd lambdas = VectorXd::Constant(q, R_PosInf);
  int iterations = 0;
  bool converged = false;
  bool step_failed = false;
  PoissonCheckpoint checkpoint;
  if (validation) checkpoint.initialize(*validation, vars, always_include, channel, lambda);

  for (int iter = 0; iter < maxit; ++iter) {
    iterations = iter + 1;
    for (int i = 0; i < n; ++i) {
      const double wi = prior_w[i] / mu[i];
      weights[i] = finite_pos(wi) ? wi : 0.0;
    }
    z = y - offset;

    for (int j = 0; j < q; ++j) {
      lambdas[j] = always_include[vars[j]] ? 0.0 : (beta[j] == 0.0) ? R_PosInf :
        (lambda / 4.0) / (beta[j] * beta[j]);
    }

    std::vector<int> active = active_positions(beta, vars, always_include);
    if (active.empty()) { converged = true; break; }
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

    VectorXd eta = predict_sparse(X, offset, vars, beta);
    bool support_changed = false;
    for (int j = 0; j < q; ++j) if (!always_include[vars[j]] && ((beta[j] != 0.0) != (beta_old[j] != 0.0)))
      { support_changed = true; break; }
    auto predict = [&](const VectorXd& b, bool zero_offset) {
      return predict_sparse(X, zero_offset ? VectorXd::Zero(n).eval() : offset, vars, b);
    };
    if (!poisson_step(beta, beta_old, eta, eta_previous, y, prior_w, lambdas, prior_mean,
        previous_deviance, iter > 0 && !support_changed, predict)) {
      step_failed = true; break;
    }
    if (hard_threshold) poisson_threshold(beta, eta, previous_deviance, y, prior_w,
      scaling, vars, always_include, thresh, predict);
    for (int i = 0; i < n; ++i) mu[i] = std::max(eta[i], eps);

    if (validation && checkpoint.update(*validation, beta, lambdas, vars, channel, iter + 1)) break;
    const VectorXd prediction = normalize ? (eta - offset).eval() : VectorXd();
    const double rel = coefficient_change(beta, beta_old, scaling, vars, always_include,
                                          prediction, prediction_old, prior_w, normalize);
    if (rel < tol && iter > 0) { converged = true; break; }
    beta_old = beta;
    eta_previous = eta;
    if (normalize) prediction_old = prediction;
  }

  if (validation) {
    beta = checkpoint.beta; lambdas = checkpoint.lambdas;
    const VectorXd restored = predict_sparse(X, offset, vars, beta);
    for (int i = 0; i < n; ++i) weights[i] = prior_w[i] == 0.0 ? 0.0 : prior_w[i] / restored[i];
  }
  FitResult out;
  out.validation = checkpoint.trace;
  out.beta = beta;
  out.weights = weights;
  out.z = z;
  out.lambdas = lambdas;
  out.iterations = iterations;
  out.converged = converged;
  out.step_failed = step_failed;
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
                        double thresh,
                        double tol,
                        int nthreads,
                        bool natural_ordering,
                        Rcpp::Nullable<Rcpp::List> validation_r) {
  const int n = X.rows();
  const int p = X.cols();
  const int m = Ymap.cols();
  PoissonValidation validation_data(validation_r);
  const PoissonValidation* validation = validation_data.enabled ? &validation_data : NULL;
  PoissonValidationTrace group_trace;
  Rcpp::List validation_passes, relaxed_validation;
  MatrixXd best_coef, best_lambdas, best_weights, best_z;
  double pass_loss = 0.0;

  MatrixXd Y = Ymap;
  MatrixXd prior_weights = Wmap;
  MatrixXd Offset = Offsetmap;
  MatrixXd prior_mean_full = Priormap;
  const VectorXd scaling = scaling_full;

  MatrixXd start_full;
  bool has_start = start_nullable.isNotNull();
  if (has_start) {
    Rcpp::NumericMatrix sr(start_nullable.get());
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
  std::vector<std::string> stages;
  std::vector<int> passes, channels, iterations;
  std::vector<bool> converged;

  auto preserve_start_mask = [&]() {
    if (has_start) for (int a = 0; a < static_cast<int>(vars.size()); ++a)
      for (int col = 0; col < m; ++col)
        if (!always_include[vars[a]] && start_full(vars[a], col) == 0.0)
          start_current(a, col) = 0.0;
  };

  auto fit_for_vars = [&](double lambda_fit, MatrixXd* start_for_fit,
                          bool hard_threshold, int pass) {
    const int q = static_cast<int>(vars.size());
    std::vector<FitResult> fits(m);
    int failed_steps = 0;
    Rcpp::List channel_validation(m);
    pass_loss = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads) if(threads > 1 && m > 1)
#endif
    for (int col = 0; col < m; ++col) {
      fits[col] = fit_channel_impl<Design>(
        X, Y, prior_weights, Offset, prior_mean_full, start_for_fit, col, vars,
        always_include, lambda_fit, normalize, scaling, maxit,
        thresh, tol, natural_ordering, hard_threshold, validation);
    }

    coef.setZero();
    lambdas.setZero();
    for (int col = 0; col < m; ++col) {
      stages.push_back(hard_threshold ? "selection" : "relaxed");
      passes.push_back(pass);
      channels.push_back(col + 1);
      iterations.push_back(fits[col].iterations);
      converged.push_back(fits[col].converged);
      if (fits[col].step_failed) ++failed_steps;
      for (int a = 0; a < q; ++a) {
        coef(vars[a], col) = fits[col].beta[a];
        lambdas(vars[a], col) = fits[col].lambdas[a];
      }
      if (validation) {
        channel_validation[col] = fits[col].validation.info(fits[col].iterations);
        if (validation_data.mass[col] > 0.0)
          pass_loss += validation_data.mass[col] * fits[col].validation.best;
      }
      weights_out.col(col) = fits[col].weights;
      z_out.col(col) = fits[col].z;
    }
    if (validation) {
      if (hard_threshold) validation_passes.push_back(channel_validation);
      else relaxed_validation = channel_validation;
    }
    if (failed_steps) Rcpp::warning("No decreasing finite Poisson step in %d channel(s); retaining previous coefficients.", failed_steps);
  };

  for (int iter = 0; iter < effective_group_maxit; ++iter) {
    if (iter > 0) {
      group_support_start(coef, vars, start_current, always_include);
      preserve_start_mask();
      if (vars.empty()) break;
      start_ptr = &start_current;
    }
    fit_for_vars(lambda, start_ptr, true, iter + 1);
    if (validation) {
      if (group_trace.update(validation_data, pass_loss / validation_data.mass.sum(),
                             (coef.array() != 0.0).count(), iter + 1)) {
        best_coef = coef; best_lambdas = lambdas; best_weights = weights_out; best_z = z_out;
      }
      if (group_trace.stopped) break;
    }
  }

  if (validation) {
    coef = best_coef; lambdas = best_lambdas; weights_out = best_weights; z_out = best_z;
    validation_data.baseline = coef;
  }
  if (relaxed) {
    group_support_start(coef, vars, start_current, always_include);
    preserve_start_mask();
    if (!vars.empty()) {
      fit_for_vars(1e-30, &start_current, false, 1);
    }
  }

  if (relaxed) {
    const VectorXd lower = VectorXd::Zero(p);
    const bool positive_design = L0glm::boundary_nonnegative(X);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads) if(threads > 1 && m > 1)
#endif
    for (int col = 0; col < m; ++col) {
      VectorXd beta = coef.col(col), previous = beta, eta = Offset.col(col);
      for (int j : vars) if (beta[j] != 0.0)
        L0glm::boundary_column(X,j,[&](int i,double x) { eta[i] += x*beta[j]; });
      int removed = L0glm::boundary_identity(X,Y.col(col),prior_weights.col(col),
                                             beta,eta,lower,1,positive_design && (Offset.col(col).array()>=0.0).all());
      if (removed && validation) {
        VectorXd old_active(vars.size()), new_active(vars.size());
        for (int a = 0; a < static_cast<int>(vars.size()); ++a) {
          old_active[a] = previous[vars[a]]; new_active[a] = beta[vars[a]];
        }
        const double old_loss = validation_data.loss(old_active,vars,col);
        const double new_loss = validation_data.loss(new_active,vars,col);
        if (!std::isfinite(new_loss) || new_loss > old_loss) removed = 0;
      }
      if (removed) {
        coef.col(col) = beta;
        // Rebuild the final signal: endpoint subtractions can leave roundoff
        // at exact zero means, where inverse-variance weights are sensitive.
        eta = Offset.col(col);
        for (int j : vars) if (beta[j] != 0.0)
          L0glm::boundary_column(X,j,[&](int i,double x) { eta[i] += x*beta[j]; });
        for (int i = 0; i < n; ++i) weights_out(i,col) = prior_weights(i,col) == 0.0 ?
          0.0 : prior_weights(i,col) / std::max(eta[i],std::numeric_limits<double>::epsilon());
      }
    }
  }
  Rcpp::NumericVector dispersion(m, 1.0);
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("coefficients") = coef,
    Rcpp::Named("weights") = weights_out,
    Rcpp::Named("dispersion") = dispersion,
    Rcpp::Named("lambdas") = lambdas,
    Rcpp::Named("z") = z_out,
    Rcpp::Named("convergence") = Rcpp::DataFrame::create(
      Rcpp::Named("stage") = stages, Rcpp::Named("pass") = passes,
      Rcpp::Named("channel") = channels, Rcpp::Named("iterations") = iterations,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("dispersion_updates") = Rcpp::IntegerVector(stages.size(), 0),
      Rcpp::Named("dispersion_restarts") = Rcpp::IntegerVector(stages.size(), 0)),
    Rcpp::Named("openmp") =
#ifdef _OPENMP
      true,
#else
      false,
#endif
    Rcpp::Named("nthreads") = threads);
  if (validation) {
    Rcpp::List selection;
    if (m == 1) selection = Rcpp::as<Rcpp::List>(validation_passes[0])[0];
    else {
      selection = group_trace.info(validation_passes.size(), false);
      selection["passes"] = validation_passes;
    }
    Rcpp::List info = Rcpp::List::create(Rcpp::Named("selection") = selection);
    if (relaxed_validation.size()) info["relaxed"] = m == 1 ? relaxed_validation[0] : Rcpp::wrap(relaxed_validation);
    result["validation"] = info;
  }
  return result;
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
                              double thresh,
                              double tol,
                              int nthreads,
                              bool natural_ordering,
                              Rcpp::Nullable<Rcpp::List> validation = R_NilValue) {
  return fit_all_impl<Map<MatrixXd> >(
    X, y, weights, offset, prior_mean, start, always_include,
    lambda, normalize, scaling, grouped, relaxed, group_maxit, maxit,
    thresh, tol, nthreads, natural_ordering, validation);
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
                               double thresh,
                               double tol,
                               int nthreads,
                               bool natural_ordering,
                              Rcpp::Nullable<Rcpp::List> validation = R_NilValue) {
  return fit_all_impl<Eigen::MappedSparseMatrix<double> >(
    X, y, weights, offset, prior_mean, start, always_include,
    lambda, normalize, scaling, grouped, relaxed, group_maxit, maxit,
    thresh, tol, nthreads, natural_ordering, validation);
}
