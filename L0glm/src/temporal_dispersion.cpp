#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

double difference_median(std::vector<double>& x) {
  if (x.empty()) return NA_REAL;
  const size_t k = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + k, x.end());
  const double upper = x[k];
  return x.size() % 2 ? upper : .5 * (upper + *std::max_element(x.begin(), x.begin() + k));
}

// All R allocation occurs outside the parallel region. Channels are independent.
// [[Rcpp::export]]
Rcpp::List difference_scale_cpp(Rcpp::NumericMatrix residual,
                                Rcpp::NumericMatrix variance,
                                Rcpp::NumericMatrix level,
                                Rcpp::LogicalMatrix eligible,
                                Rcpp::NumericVector coefficient,
                                double minimum_mean, int nthreads) {
  const int n = residual.nrow(), m = residual.ncol(), r = coefficient.size() - 1;
  Rcpp::NumericVector energy(m), normalizer(m), cells(m), moment(m), robust(m);
  std::vector<std::vector<double>> ratios(m);
  const double median_chisq1 = 0.4549364231195727;
#ifdef _OPENMP
  const int threads = std::max(1, std::min(nthreads, m));
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
  for (int j = 0; j < m; ++j) {
    auto& values = ratios[j];values.reserve(n - r);
    double numerator = 0, denominator = 0;
    for (int i = 0; i < n - r; ++i) {
      double delta = 0, v = 0, local_level = 0, coefficient_mass = 0;
      bool valid = true;
      for (int k = 0; k <= r; ++k) {
        if (eligible(i + k, j) != TRUE) {valid = false;break;}
        const double a = coefficient[k], a2 = a * a;
        delta += a * residual(i + k, j);
        v += a2 * variance(i + k, j);
        local_level += a2 * level(i + k, j);coefficient_mass += a2;
      }
      if (!valid || v <= 0 || !std::isfinite(v) || !std::isfinite(delta) ||
          local_level < minimum_mean * coefficient_mass) continue;
      numerator += delta * delta;denominator += v;
      values.push_back(delta * delta / v);
    }
    energy[j] = numerator;normalizer[j] = denominator;cells[j] = values.size();
    moment[j] = denominator > 0 ? numerator / denominator : NA_REAL;
    robust[j] = difference_median(values) / median_chisq1;
  }
  std::vector<double> pooled;
  size_t count = 0;for (const auto& x : ratios) count += x.size();
  pooled.reserve(count);
  double total_energy = 0, total_normalizer = 0;
  for (int j = 0; j < m; ++j) {
    pooled.insert(pooled.end(), ratios[j].begin(), ratios[j].end());
    total_energy += energy[j];total_normalizer += normalizer[j];
  }
  return Rcpp::List::create(
    Rcpp::_ ["channel"] = Rcpp::DataFrame::create(
      Rcpp::_ ["channel"] = Rcpp::seq_len(m), Rcpp::_ ["cells"] = cells,
      Rcpp::_ ["energy"] = energy, Rcpp::_ ["normalizer"] = normalizer,
      Rcpp::_ ["moment"] = moment, Rcpp::_ ["robust"] = robust),
    Rcpp::_ ["pooled"] = Rcpp::NumericVector::create(
      Rcpp::_ ["cells"] = static_cast<double>(count),
      Rcpp::_ ["moment"] = total_normalizer > 0 ? total_energy / total_normalizer : NA_REAL,
      Rcpp::_ ["robust"] = difference_median(pooled) / median_chisq1));
}
