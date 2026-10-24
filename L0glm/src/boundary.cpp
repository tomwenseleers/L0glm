#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include "boundary.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List boundary_dense_cpp(const Eigen::Map<Eigen::MatrixXd>& X,
  Eigen::VectorXd y,Eigen::VectorXd weights,Eigen::VectorXd beta,
  Eigen::VectorXd eta,Eigen::VectorXd lower,int family,bool nonnegative_components=false) {
  const int removed=L0glm::boundary_identity(X,y,weights,beta,eta,lower,family,nonnegative_components);
  return Rcpp::List::create(Rcpp::Named("beta")=beta,Rcpp::Named("eta")=eta,
    Rcpp::Named("removed")=removed);
}
// [[Rcpp::export]]
Rcpp::List boundary_sparse_cpp(const Eigen::MappedSparseMatrix<double>& X,
  Eigen::VectorXd y,Eigen::VectorXd weights,Eigen::VectorXd beta,
  Eigen::VectorXd eta,Eigen::VectorXd lower,int family,bool nonnegative_components=false) {
  const int removed=L0glm::boundary_identity(X,y,weights,beta,eta,lower,family,nonnegative_components);
  return Rcpp::List::create(Rcpp::Named("beta")=beta,Rcpp::Named("eta")=eta,
    Rcpp::Named("removed")=removed);
}
