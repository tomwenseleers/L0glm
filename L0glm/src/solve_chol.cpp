// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
// Solver wrappers adapted from nk027/sanic.
namespace L0glm {
SolverResult chol_solve(DenseRef A,DenseRef b,unsigned int pivot,unsigned int) {
  if(pivot==0)return factor_solve<Eigen::LLT<Eigen::MatrixXd>>(A,b);
  return factor_solve<Eigen::LDLT<Eigen::MatrixXd>>(A,b);
}
SolverResult chol_solve(SparseRef A,DenseRef b,unsigned int pivot,unsigned int ordering) {
  using Eigen::SparseMatrix;
  if(pivot==0) {
    if(ordering==1)return factor_solve<Eigen::SimplicialLLT<SparseMatrix<double>,
      Eigen::Lower,Eigen::NaturalOrdering<int>>>(A,b);
    return factor_solve<Eigen::SimplicialLLT<SparseMatrix<double>,
      Eigen::Lower,Eigen::AMDOrdering<int>>>(A,b);
  }
  if(ordering==1)return factor_solve<Eigen::SimplicialLDLT<SparseMatrix<double>,
    Eigen::Lower,Eigen::NaturalOrdering<int>>>(A,b);
  return factor_solve<Eigen::SimplicialLDLT<SparseMatrix<double>,
    Eigen::Lower,Eigen::AMDOrdering<int>>>(A,b);
}
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_chol(const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int pivot=1,unsigned int ordering=0) {
  if(ordering>1){Rcpp::warning("No valid ordering requested -- using default.");ordering=0;}
  if(pivot>1)Rcpp::warning("No valid pivoting scheme requested -- using default.");
  return L0glm::chol_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),pivot,ordering).list();
}
// [[Rcpp::export]]
Rcpp::List solve_dense_chol(const Eigen::Map<Eigen::MatrixXd> &A,
    const Eigen::Map<Eigen::MatrixXd> b,unsigned int pivot=1) {
  if(pivot>1)Rcpp::warning("No valid pivoting scheme requested -- using default.");
  return L0glm::chol_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),pivot).list();
}
