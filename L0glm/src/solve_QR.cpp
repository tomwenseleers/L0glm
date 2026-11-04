// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
#include <Eigen/SparseQR>
// Solver wrappers adapted from nk027/sanic.
namespace L0glm {
template<class Solver> SolverResult dense_qr(DenseRef A,DenseRef b,bool ranked) {
  SolverResult result;Solver solver;solver.compute(A);
  if(ranked && solver.rank()<A.cols())return result;
  result.coefficients=solver.solve(b);result.success=result.finite();return result;
}
SolverResult qr_solve(DenseRef A,DenseRef b,unsigned int pivot,unsigned int) {
  if(pivot==0) {
    SolverResult result;Eigen::HouseholderQR<Eigen::MatrixXd> solver(A);
    result.coefficients=solver.solve(b);result.success=result.finite();return result;
  }
  if(pivot==2)return dense_qr<Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>(A,b,true);
  return dense_qr<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>(A,b,true);
}
SolverResult qr_solve(SparseRef A,DenseRef b,unsigned int,unsigned int ordering) {
  if(ordering==0)return factor_solve<Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::AMDOrdering<int>>>(A,b);
  if(ordering==2)return factor_solve<Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::NaturalOrdering<int>>>(A,b);
  return factor_solve<Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>(A,b);
}
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_QR(const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int ordering=1) {
  if(ordering>2)Rcpp::warning("No valid ordering requested -- using default.");
  return L0glm::qr_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),1,ordering).list();
}
// [[Rcpp::export]]
Rcpp::List solve_dense_QR(const Eigen::Map<Eigen::MatrixXd> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int pivot=1) {
  if(pivot>2)Rcpp::warning("No valid pivoting scheme requested -- using default ColPivHouseholderQR.");
  return L0glm::qr_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),pivot).list(2);
}
