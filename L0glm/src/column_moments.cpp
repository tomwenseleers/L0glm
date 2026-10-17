// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// Column information without constructing X^2 or changing design storage.
// [[Rcpp::export]]
Rcpp::List L0glm_column_moments_cpp(SEXP X,
                                  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
  const bool s4 = Rf_isS4(X);
  Rcpp::IntegerVector dim = s4 ? Rcpp::as<Rcpp::IntegerVector>(Rcpp::S4(X).slot("Dim")) :
    Rcpp::as<Rcpp::IntegerVector>(Rf_getAttrib(X, R_DimSymbol));
  const int n = dim[0], p = dim[1];
  Rcpp::NumericVector values = s4 ? Rcpp::as<Rcpp::NumericVector>(Rcpp::S4(X).slot("x")) :
    Rcpp::as<Rcpp::NumericVector>(X);
  Rcpp::NumericVector w = weights.isNotNull() ? Rcpp::NumericVector(weights) : Rcpp::NumericVector(n, 1.0);
  if (w.size() != n) Rcpp::stop("weights must have one value per row.");
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(w[i]) || w[i] < 0.0) Rcpp::stop("weights must be finite and nonnegative.");
    total += w[i];
  }
  Rcpp::NumericVector mean(p), ss(p), centered(p);
  const bool sparse = s4 && Rcpp::S4(X).hasSlot("p");
  Rcpp::IntegerVector pointers, rows;
  if (sparse) { pointers = Rcpp::S4(X).slot("p"); rows = Rcpp::S4(X).slot("i"); }
  for (int j = 0; j < p; ++j) {
    double sum = 0.0, represented = 0.0;
    if (sparse) {
      for (int k = pointers[j]; k < pointers[j + 1]; ++k) {
        const double wi = w[rows[k]], x = values[k];
        sum += wi * x; ss[j] += wi * x * x; represented += wi;
      }
      mean[j] = total > 0.0 ? sum / total : 0.0;
      for (int k = pointers[j]; k < pointers[j + 1]; ++k) {
        const double difference = values[k] - mean[j];
        centered[j] += w[rows[k]] * difference * difference;
      }
      centered[j] += std::max(0.0, total - represented) * mean[j] * mean[j];
    } else {
      const R_xlen_t offset = static_cast<R_xlen_t>(j) * n;
      for (int i = 0; i < n; ++i) {
        const double x = values[offset + i];
        sum += w[i] * x; ss[j] += w[i] * x * x;
      }
      mean[j] = total > 0.0 ? sum / total : 0.0;
      for (int i = 0; i < n; ++i) {
        const double difference = values[offset + i] - mean[j];
        centered[j] += w[i] * difference * difference;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("mean") = mean, Rcpp::Named("ss") = ss,
                             Rcpp::Named("centered_ss") = centered,
                             Rcpp::Named("weight_sum") = total);
}
