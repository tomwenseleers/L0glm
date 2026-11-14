// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
#include <type_traits>
#include <vector>

namespace {
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// Feasible active-set polish of x'Gx/2 - h'x, lower <= x <= upper.
// Diagonal equilibration changes coordinates, not the quadratic or bounds.
template<class Matrix>
Rcpp::List bounded_polish(Matrix G, VectorXd h, VectorXd x,
                         VectorXd lower, VectorXd upper) {
  const int p = h.size();
  if (G.rows()!=p || G.cols()!=p || x.size()!=p || lower.size()!=p || upper.size()!=p)
    Rcpp::stop("Incompatible bounded ridge dimensions.");
  const VectorXd original = x;
  VectorXd scale(p);
  std::vector<int> state(p), free;
  bool valid = h.allFinite() && x.allFinite();
  for (int j=0; j<p; ++j) {
    const double diagonal = G.coeff(j,j);
    valid = valid && std::isfinite(diagonal) && diagonal>=0 &&
      !std::isnan(lower[j]) && !std::isnan(upper[j]) && lower[j]<=upper[j];
    scale[j] = diagonal>0 ? std::sqrt(diagonal) : 1;
    x[j] = std::max(lower[j],std::min(upper[j],x[j]));
    // A zero-information Gram column has no identifiable update.
    if (diagonal==0) {
      valid = valid && h[j]==0;
      lower[j] = upper[j] = x[j];
    }
  }
  int iterations = 0;
  double error = R_PosInf;
  bool success = valid && p==0;
  if (valid && p) {
    if constexpr (std::is_same_v<Matrix,MatrixXd>) {
      for (int j=0; j<p; ++j) for (int i=0; i<p; ++i)
        G(i,j) /= scale[i]*scale[j];
    } else {
      for (int j=0; j<G.outerSize(); ++j)
        for (typename Matrix::InnerIterator entry(G,j); entry; ++entry)
          entry.valueRef() /= scale[entry.row()]*scale[entry.col()];
    }
    h.array() /= scale.array();
    x.array() *= scale.array();
    lower.array() *= scale.array();
    upper.array() *= scale.array();
    const double denominator = std::max(1.,h.cwiseAbs().maxCoeff());
    const double tolerance = 1e-10*denominator;
    for (int j=0; j<p; ++j)
      state[j] = x[j]<=lower[j] ? -1 : x[j]>=upper[j] ? 1 : 0;

    for (; iterations<5*p+20 && x.allFinite(); ++iterations) {
      if (iterations%64==0) Rcpp::checkUserInterrupt();
      free.clear();
      for (int j=0; j<p; ++j) if (state[j]==0) free.push_back(j);
      VectorXd fixed = x;
      for (int j:free) fixed[j] = 0;
      const VectorXd rhs = h-G*fixed;
      const int q = free.size();
      VectorXd target(q);
      for (int k=0; k<q; ++k) target[k] = rhs[free[k]];
      if (q) {
        if constexpr (std::is_same_v<Matrix,MatrixXd>) {
          MatrixXd A(q,q);
          for (int j=0; j<q; ++j) for (int i=0; i<q; ++i)
            A(i,j) = G(free[i],free[j]);
          Eigen::LLT<MatrixXd> factor(A);
          if (factor.info()!=Eigen::Success) break;
          target = factor.solve(target).eval();
        } else {
          std::vector<int> index(p,-1);
          for (int k=0; k<q; ++k) index[free[k]] = k;
          std::vector<Eigen::Triplet<double>> entries;
          for (int k=0; k<q; ++k)
            for (typename Matrix::InnerIterator entry(G,free[k]); entry; ++entry)
              if (index[entry.row()]>=0)
                entries.emplace_back(index[entry.row()],k,entry.value());
          SparseMatrix<double> A(q,q);
          A.setFromTriplets(entries.begin(),entries.end());
          Eigen::SimplicialLLT<SparseMatrix<double>> factor(A);
          if (factor.info()!=Eigen::Success) break;
          target = factor.solve(target).eval();
        }
        if (!target.allFinite()) break;
        // Stop at the first bound, retain feasibility, then solve the new face.
        double step = 1;
        int hit = -1, side = 0;
        for (int k=0; k<q; ++k) {
          const int j = free[k];
          double fraction = 1;
          int bound = 0;
          if (target[k]<lower[j]) {
            fraction = (lower[j]-x[j])/(target[k]-x[j]); bound = -1;
          } else if (target[k]>upper[j]) {
            fraction = (upper[j]-x[j])/(target[k]-x[j]); bound = 1;
          }
          if (bound && fraction<=step) {
            step = std::max(0.,fraction); hit = j; side = bound;
          }
        }
        for (int k=0; k<q; ++k) x[free[k]] += step*(target[k]-x[free[k]]);
        if (hit>=0) {
          state[hit] = side;
          x[hit] = side<0 ? lower[hit] : upper[hit];
          continue;
        }
      }
      const VectorXd gradient = G*x-h;
      double violation = tolerance;
      int release = -1;
      for (int j=0; j<p; ++j) if (lower[j]<upper[j]) {
        const double value = state[j]<0 ? -gradient[j] : state[j]>0 ? gradient[j] : 0;
        if (value>violation) { violation = value; release = j; }
      }
      if (release<0) { success = true; break; }
      state[release] = 0;
    }
    // Fresh KKT residual, including stationarity on the free coefficients.
    for (int j=0; j<p; ++j) x[j] = std::max(lower[j],std::min(upper[j],x[j]));
    const VectorXd gradient = G*x-h;
    error = 0;
    for (int j=0; j<p; ++j) if (lower[j]<upper[j]) {
      double value = gradient[j];
      if ((x[j]<=lower[j] && value>0) || (x[j]>=upper[j] && value<0)) value = 0;
      error = std::max(error,std::abs(value));
    }
    error /= denominator;
    success = success && x.allFinite() && std::isfinite(error) && error<=1e-10;
    x.array() /= scale.array();
  }
  return Rcpp::List::create(Rcpp::_["coefficients"]=success ? x : original,
    Rcpp::_["success"]=success,Rcpp::_["iterations"]=iterations,
    Rcpp::_["error"]=p ? error : 0);
}
}

// [[Rcpp::export]]
Rcpp::List solve_bounded_gram(SEXP G, Eigen::VectorXd h, Eigen::VectorXd x,
                            Eigen::VectorXd lower, Eigen::VectorXd upper) {
  if (Rf_isS4(G))
    return bounded_polish(Rcpp::as<SparseMatrix<double>>(G),h,x,lower,upper);
  return bounded_polish(Rcpp::as<MatrixXd>(G),h,x,lower,upper);
}
