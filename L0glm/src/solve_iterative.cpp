// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
// Solver wrappers adapted from nk027/sanic.
namespace L0glm {
SolverResult cg_solve(DenseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return iterative_solve<Eigen::ConjugateGradient<Eigen::MatrixXd,
    Eigen::Lower|Eigen::Upper,Eigen::IdentityPreconditioner>>(A,b,start,tol,maxit);
  return iterative_solve<Eigen::ConjugateGradient<Eigen::MatrixXd,
    Eigen::Lower|Eigen::Upper,Eigen::DiagonalPreconditioner<double>>>(A,b,start,tol,maxit);
}
SolverResult cg_solve(SparseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  using Eigen::SparseMatrix;
  if(precond==0)return iterative_solve<Eigen::ConjugateGradient<SparseMatrix<double>,
    Eigen::Lower|Eigen::Upper,Eigen::IdentityPreconditioner>>(A,b,start,tol,maxit);
  if(precond==2)return iterative_solve<Eigen::ConjugateGradient<SparseMatrix<double>,
    Eigen::Lower|Eigen::Upper,Eigen::IncompleteCholesky<double>>>(A,b,start,tol,maxit);
  return iterative_solve<Eigen::ConjugateGradient<SparseMatrix<double>,
    Eigen::Lower|Eigen::Upper,Eigen::DiagonalPreconditioner<double>>>(A,b,start,tol,maxit);
}
SolverResult lscg_solve(DenseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return iterative_solve<Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd,
    Eigen::IdentityPreconditioner>>(A,b,start,tol,maxit);
  return iterative_solve<Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd,
    Eigen::LeastSquareDiagonalPreconditioner<double>>>(A,b,start,tol,maxit);
}
SolverResult lscg_solve(SparseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return iterative_solve<Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
    Eigen::IdentityPreconditioner>>(A,b,start,tol,maxit);
  return iterative_solve<Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
    Eigen::LeastSquareDiagonalPreconditioner<double>>>(A,b,start,tol,maxit);
}
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_LSCG(const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,const Eigen::Map<Eigen::MatrixXd> x0,
    double tol=0,int maxit=0,unsigned int precond=1) {
  if(precond>1)Rcpp::warning("No valid preconditioner requested -- using default.");
  return L0glm::lscg_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_dense_LSCG(const Eigen::MatrixXd &A,const Eigen::MatrixXd b,const Eigen::MatrixXd x0,
    double tol=0,int maxit=0,unsigned int precond=1) {
  if(precond>1)Rcpp::warning("No valid preconditioner requested -- using default.");
  return L0glm::lscg_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_CG(Eigen::MappedSparseMatrix<double> &A,Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,double tol=0,int maxit=0,unsigned int precond=1) {
  if(precond>2)Rcpp::warning("No valid preconditioner requested -- using default.");
  return L0glm::cg_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_dense_CG(Eigen::Map<Eigen::MatrixXd> &A,Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,double tol=0,int maxit=0,unsigned int precond=1) {
  if(precond>2)Rcpp::warning("No valid preconditioner requested for dense CG -- using default.");
  return L0glm::cg_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
