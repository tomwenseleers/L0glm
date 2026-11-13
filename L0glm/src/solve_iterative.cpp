// [[Rcpp::depends(RcppEigen)]]
#include "solver_native.h"
#include <type_traits>
// Solver wrappers adapted from nk027/sanic.
namespace L0glm {
// Keep the dense CG operator; only the IC factorization needs sparse storage.
// Preserve every nonzero entry, including small adaptive-ridge diagonals.
class DenseIncompleteCholesky : public Eigen::IncompleteCholesky<double> {
public:
  template<class Matrix> DenseIncompleteCholesky& compute(const Matrix& A) {
    Eigen::SparseMatrix<double> sparse=A.sparseView(0.,0.);
    Eigen::IncompleteCholesky<double>::compute(sparse);
    return *this;
  }
};

// LSCG applies its preconditioner to the normal residual, so factor A' A
// rather than A. Matrix-vector products still use the original operator.
class LeastSquaresIncompleteCholesky : public Eigen::IncompleteCholesky<double> {
public:
  template<class Matrix> LeastSquaresIncompleteCholesky& compute(const Matrix& A) {
    Eigen::SparseMatrix<double> gram;
    if constexpr(std::is_same_v<typename Matrix::StorageKind,Eigen::Sparse>) {
      gram=A.transpose()*A;
      gram.prune(0.,0.);
    }else {
      Eigen::MatrixXd dense=A.transpose()*A;
      gram=dense.sparseView(0.,0.);
    }
    Eigen::IncompleteCholesky<double>::compute(gram);
    return *this;
  }
};

SolverResult cg_solve(DenseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return iterative_solve<Eigen::ConjugateGradient<Eigen::MatrixXd,
    Eigen::Lower|Eigen::Upper,Eigen::IdentityPreconditioner>>(A,b,start,tol,maxit);
  if(precond==2)return iterative_solve<Eigen::ConjugateGradient<Eigen::MatrixXd,
    Eigen::Lower|Eigen::Upper,DenseIncompleteCholesky>>(A,b,start,tol,maxit);
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
// Check every RHS using a fresh normal residual. If recursive convergence
// is inaccurate, IC-LSCG may restart once, reusing the factorization and the
// remaining per-RHS iteration budget. No extra factorization is needed.
template<class Preconditioner,class Matrix> SolverResult checked_lscg(
    const Matrix& A,DenseRef b,DenseRef start,double tol,int maxit,int restarts=0) {
  Eigen::LeastSquaresConjugateGradient<typename Matrix::PlainObject,Preconditioner> solver;
  if(tol>0)solver.setTolerance(tol);
  if(maxit>0)solver.setMaxIterations(maxit);
  solver.compute(A);SolverResult result;result.info=solver.info();
  if(solver.info()!=Eigen::Success)return result;
  const int budget=solver.maxIterations();const double target=solver.tolerance();
  result.coefficients=start;
  for(int j=0;j<b.cols();++j) {
    Eigen::VectorXd guess=start.col(j),x;
    const double denominator=(A.transpose()*b.col(j)).norm();
    double error=0;int used=0;
    for(int attempt=0;;++attempt) {
      const int remaining=budget-used;solver.setMaxIterations(remaining);
      x=solver.solveWithGuess(b.col(j),guess);
      result.iterations+=solver.iterations();
      // Eigen reports zero for convergence on the first update. Charge that
      // update too, conservatively including an already-converged guess.
      used+=std::min(remaining,int(solver.iterations())+1);
      const double numerator=(A.transpose()*(b.col(j)-A*x)).norm();
      error=denominator>0?numerator/denominator:
        (numerator==0?0:std::numeric_limits<double>::infinity());
      if(!std::isfinite(error))error=std::numeric_limits<double>::infinity();
      if(error<=target||solver.info()!=Eigen::Success||attempt>=restarts||
         used>=budget||!x.allFinite())break;
      guess=x;
    }
    result.coefficients.col(j)=x;
    result.error=std::max(result.error,error);
  }
  result.success=result.finite()&&result.error<=target;
  result.info=result.success?Eigen::Success:Eigen::NoConvergence;
  return result;
}
SolverResult lscg_solve(DenseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return checked_lscg<Eigen::IdentityPreconditioner>(A,b,start,tol,maxit);
  if(precond==2)return checked_lscg<LeastSquaresIncompleteCholesky>(A,b,start,tol,maxit,1);
  return checked_lscg<Eigen::LeastSquareDiagonalPreconditioner<double>>(A,b,start,tol,maxit);
}
SolverResult lscg_solve(SparseRef A,DenseRef b,DenseRef start,double tol,int maxit,unsigned int precond) {
  if(precond==0)return checked_lscg<Eigen::IdentityPreconditioner>(A,b,start,tol,maxit);
  if(precond==2)return checked_lscg<LeastSquaresIncompleteCholesky>(A,b,start,tol,maxit,1);
  return checked_lscg<Eigen::LeastSquareDiagonalPreconditioner<double>>(A,b,start,tol,maxit);
}
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_LSCG(const Eigen::MappedSparseMatrix<double> &A,
    const Eigen::Map<Eigen::MatrixXd> b,const Eigen::Map<Eigen::MatrixXd> x0,
    double tol=0,int maxit=0,unsigned int precond=2) {
  if(precond>2){Rcpp::warning("No valid preconditioner requested -- using incomplete Cholesky.");precond=2;}
  return L0glm::lscg_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_dense_LSCG(const Eigen::MatrixXd &A,const Eigen::MatrixXd b,const Eigen::MatrixXd x0,
    double tol=0,int maxit=0,unsigned int precond=2) {
  if(precond>2){Rcpp::warning("No valid preconditioner requested -- using incomplete Cholesky.");precond=2;}
  return L0glm::lscg_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_sparse_CG(Eigen::MappedSparseMatrix<double> &A,Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,double tol=0,int maxit=0,unsigned int precond=2) {
  if(precond>2){Rcpp::warning("No valid preconditioner requested -- using incomplete Cholesky.");precond=2;}
  return L0glm::cg_solve(L0glm::SparseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
// [[Rcpp::export]]
Rcpp::List solve_dense_CG(Eigen::Map<Eigen::MatrixXd> &A,Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,double tol=0,int maxit=0,unsigned int precond=2) {
  if(precond>2){Rcpp::warning("No valid preconditioner requested -- using incomplete Cholesky.");precond=2;}
  return L0glm::cg_solve(L0glm::DenseRef(A),L0glm::DenseRef(b),L0glm::DenseRef(x0),tol,maxit,precond).list(1);
}
