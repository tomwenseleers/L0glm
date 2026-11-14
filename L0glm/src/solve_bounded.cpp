// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
#include "bounded.h"

namespace {
template<class Matrix>
Rcpp::List bounded_result(const Matrix& G,const Eigen::VectorXd& h,const Eigen::VectorXd& x,
  const Eigen::VectorXd& lower,const Eigen::VectorXd& upper,bool batched) {
  const int p=h.size();
  if(G.rows()!=p||G.cols()!=p||x.size()!=p||lower.size()!=p||upper.size()!=p)
    Rcpp::stop("Incompatible bounded ridge dimensions.");
  auto result=batched?L0glm::bounded_block(G,h,x,lower,upper,Rcpp::checkUserInterrupt):
    L0glm::bounded_polish(G,h,x,lower,upper,Rcpp::checkUserInterrupt);
  if(batched&&result.success) {
    const Eigen::VectorXd clipped=x.cwiseMax(lower).cwiseMin(upper);
    result.success=L0glm::bounded_improves(G,h,clipped,result.coefficients);
    if(!result.success)result.coefficients=x;
  }
  return Rcpp::List::create(Rcpp::_["coefficients"]=result.coefficients,
    Rcpp::_["success"]=result.success,Rcpp::_["iterations"]=result.iterations,Rcpp::_["error"]=result.error);
}
}

// [[Rcpp::export]]
Rcpp::List solve_bounded_gram(SEXP G, Eigen::VectorXd h, Eigen::VectorXd x,
                            Eigen::VectorXd lower, Eigen::VectorXd upper, bool batched = false) {
  if(Rf_isS4(G))return bounded_result(Rcpp::as<Eigen::SparseMatrix<double>>(G),h,x,lower,upper,batched);
  return bounded_result(Rcpp::as<Eigen::MatrixXd>(G),h,x,lower,upper,batched);
}
