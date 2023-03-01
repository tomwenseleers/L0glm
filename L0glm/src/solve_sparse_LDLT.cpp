#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Map;
using Eigen::SparseMatrix;
using Eigen::SimplicialLDLT;
using Eigen::Success;

// [[Rcpp::export]]
Rcpp::List solve_sparse_LDLT(Rcpp::List input) {

  const Map<SparseMatrix<double> > m1 = input[0]; // A
  const Map<VectorXd>              v1 = input[1]; // b
  SparseMatrix<double>             m2(m1.cols(), m1.cols());
  m2.selfadjointView<Lower>().rankUpdate(m1.adjoint());
  Eigen::SimplicialLDLT<SparseMatrix<double> > ff;
  ff.compute(m2);
  VectorXd                        res = ff.solve(m1.adjoint() * v1);
  return Rcpp::List::create(Rcpp::Named("res") = res);
}
