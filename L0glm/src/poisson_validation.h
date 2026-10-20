#ifndef L0GLM_POISSON_VALIDATION_H
#define L0GLM_POISSON_VALIDATION_H

// Constructed on the calling thread. Workers only read the shared design and
// response matrices; their traces and coefficient checkpoints are private.
struct PoissonValidation {
  bool enabled = false, sparse = false;
  Eigen::MatrixXd X, y, weights, offset, baseline;
  Eigen::RowVectorXd mass;
  Eigen::SparseMatrix<double> Xs;
  int patience = 3, minit = 3;
  double tol = 1e-8, floor = 1e-8;

  explicit PoissonValidation(Rcpp::Nullable<Rcpp::List> supplied) {
    if (supplied.isNull()) return;
    Rcpp::List v(supplied);
    enabled = true;
    SEXP design = v["X"];
    sparse = Rf_isS4(design) && Rcpp::S4(design).hasSlot("i");
    if (sparse) Xs = Rcpp::as<Eigen::SparseMatrix<double> >(design);
    else X = Rcpp::as<Eigen::MatrixXd>(design);
    y = Rcpp::as<Eigen::MatrixXd>(v["y"]);
    weights = Rcpp::as<Eigen::MatrixXd>(v["weights"]);
    mass = weights.colwise().sum();
    offset = Rcpp::as<Eigen::MatrixXd>(v["offset"]);
    baseline = Rcpp::as<Eigen::MatrixXd>(v["baseline"]);
    patience = Rcpp::as<int>(v["patience"]); minit = Rcpp::as<int>(v["minit"]);
    tol = Rcpp::as<double>(v["tol"]); floor = Rcpp::as<double>(v["mu_floor"]);
  }

  double loss(const Eigen::VectorXd& beta, const std::vector<int>& vars, int channel) const {
    Eigen::VectorXd eta = offset.col(channel);
    for (int j = 0; j < beta.size(); ++j) if (beta[j] != 0.0) {
      if (sparse) eta += Xs.col(vars[j]) * beta[j];
      else eta.noalias() += X.col(vars[j]) * beta[j];
    }
    double value = 0.0;
    for (int i = 0; i < y.rows(); ++i) {
      const double w = weights(i, channel);
      if (w == 0.0) continue;
      if (!std::isfinite(eta[i])) return std::numeric_limits<double>::infinity();
      const double mu = std::max(floor, eta[i]), yi = y(i, channel);
      const double term = yi == 0.0 ? mu : yi * std::log(yi / mu) - (yi - mu);
      if (!std::isfinite(term)) return std::numeric_limits<double>::infinity();
      value += 2.0 * w * std::max(0.0, term);
    }
    return mass[channel] > 0.0 ? value / mass[channel] : std::numeric_limits<double>::infinity();
  }
};

struct PoissonValidationTrace {
  std::vector<int> iteration, selected;
  std::vector<double> deviance;
  int best_iteration = 0, stalls = 0;
  bool stopped = false;
  double best = std::numeric_limits<double>::infinity();

  bool update(const PoissonValidation& v, double loss, int size, int iter) {
    const bool improved = iteration.empty() || (std::isfinite(loss) &&
      (!std::isfinite(best) || loss < best - v.tol * std::max(1.0, std::abs(best))));
    iteration.push_back(iter); selected.push_back(size); deviance.push_back(loss);
    if (improved) { best = loss; best_iteration = iter; stalls = 0; }
    else ++stalls;
    stopped = iter >= v.minit && stalls >= v.patience;
    return improved;
  }

  Rcpp::List info(int iterations, bool channel = true) const {
    return Rcpp::List::create(Rcpp::Named("best_iteration") = best_iteration,
      Rcpp::Named("iterations") = iterations, Rcpp::Named("stopped") = stopped,
      Rcpp::Named("deviance") = best,
      Rcpp::Named("history") = Rcpp::DataFrame::create(
        Rcpp::Named("iteration") = iteration, Rcpp::Named("deviance") = deviance,
        Rcpp::Named("selected") = selected,
        Rcpp::Named("dispersion") = Rcpp::NumericVector(iteration.size(), channel ? 1.0 : NA_REAL)));
  }
};

struct PoissonCheckpoint {
  PoissonValidationTrace trace;
  Eigen::VectorXd beta, lambdas;

  void initialize(const PoissonValidation& v, const std::vector<int>& vars,
                  const std::vector<char>& mandatory, int channel, double lambda, double gamma) {
    beta.resize(vars.size()); lambdas.resize(vars.size());
    for (int j = 0; j < beta.size(); ++j) {
      beta[j] = v.baseline(vars[j], channel);
      lambdas[j] = mandatory[vars[j]] ? 0.0 : beta[j] == 0.0 ?
        std::numeric_limits<double>::infinity() : lambda / 4.0 / std::pow(std::abs(beta[j]), gamma);
    }
    trace.update(v, v.loss(beta, vars, channel), (beta.array() != 0.0).count(), 0);
  }

  bool update(const PoissonValidation& v, const Eigen::VectorXd& b, const Eigen::VectorXd& d,
              const std::vector<int>& vars, int channel, int iter) {
    if (trace.update(v, v.loss(b, vars, channel), (b.array() != 0.0).count(), iter)) {
      beta = b; lambdas = d;
    }
    return trace.stopped;
  }
};
#endif
