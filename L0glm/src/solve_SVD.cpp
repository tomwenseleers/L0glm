// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
namespace L0glm {
template<class Solver> SolverResult svd_factor(DenseRef A,DenseRef b,unsigned int flags,double threshold) {
  SolverResult result;Solver solver;if(threshold>0)solver.setThreshold(threshold);
  solver.compute(A,flags);result.coefficients=solver.solve(b);result.rank=solver.rank();
  result.success=result.finite();return result;
}
SolverResult svd_solve(DenseRef A,DenseRef b,unsigned int svd,unsigned int precond,double threshold) {
  const unsigned int thin=Eigen::ComputeThinU|Eigen::ComputeThinV;
  const unsigned int full=Eigen::ComputeFullU|Eigen::ComputeFullV;
  if(svd==0)return svd_factor<Eigen::BDCSVD<Eigen::MatrixXd>>(A,b,thin,threshold);
  if(precond==1)return svd_factor<Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::FullPivHouseholderQRPreconditioner>>(A,b,full,threshold);
  if(precond==2)return svd_factor<Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::HouseholderQRPreconditioner>>(A,b,thin,threshold);
  if(precond==3) {
    if(A.rows()!=A.cols())return SolverResult();
    return svd_factor<Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::NoQRPreconditioner>>(A,b,full,threshold);
  }
  return svd_factor<Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::ColPivHouseholderQRPreconditioner>>(A,b,thin,threshold);
}
SolverResult svd_solve(SparseRef A,DenseRef b,unsigned int svd,unsigned int precond,double threshold) {
  Eigen::MatrixXd dense(A);return svd_solve(DenseRef(dense),b,svd,precond,threshold);
}
}
// [[Rcpp::export]]
Rcpp::List solve_dense_SVD(const Eigen::Map<Eigen::MatrixXd> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int svd=0,unsigned int precond=0,double threshold=0) {
  if(svd>1)Rcpp::warning("No valid SVD requested -- using JacobiSVD.");
  if(precond>3)Rcpp::warning("No valid JacobiSVD preconditioner requested -- using ColPivHouseholderQRPreconditioner.");
  if(svd!=0 && precond==3 && A.rows()!=A.cols())return Rcpp::List::create(
    Rcpp::_["status"]=false,Rcpp::_["message"]="JacobiSVD NoQRPreconditioner requires a square matrix.");
  return L0glm::svd_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),svd,precond,threshold).list(3);
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_SVD(const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int svd=0,unsigned int precond=0,double threshold=0) {
  if(svd>1)Rcpp::warning("No valid SVD requested -- using JacobiSVD.");
  if(precond>3)Rcpp::warning("No valid JacobiSVD preconditioner requested -- using ColPivHouseholderQRPreconditioner.");
  if(svd!=0 && precond==3 && A.rows()!=A.cols())return Rcpp::List::create(
    Rcpp::_["status"]=false,Rcpp::_["message"]="JacobiSVD NoQRPreconditioner requires a square matrix.");
  return L0glm::svd_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),svd,precond,threshold).list(3);
}
