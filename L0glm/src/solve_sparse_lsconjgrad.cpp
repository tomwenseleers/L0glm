#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Map;
using Eigen::SparseMatrix;
using Eigen::IdentityPreconditioner;
using Eigen::IncompleteCholesky;
using Eigen::LeastSquaresConjugateGradient;
using Eigen::Success;

// [[Rcpp::export]]
Rcpp::List solve_sparse_lsconjgrad(Rcpp::List input) {

const Map<SparseMatrix<double> > m1 = input[0]; // A
const Map<VectorXd>              v1 = input[1]; // b
const Map<VectorXd>              x = input[2]; // initial guess x
SparseMatrix<double>             m2(m1.cols(), m1.cols());
// typedef ConjugateGradient<SparseMatrix<double>,Lower, IncompleteCholesky<double> > ICCG;
m2.selfadjointView<Lower>().rankUpdate(m1.adjoint());
// typedef Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>, IncompleteCholesky<double> > ICCG; // (m2); // Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner<double>
  Eigen::LeastSquaresConjugateGradient<SparseMatrix<double> > ff; //
  ff.setMaxIterations(100);
ff.setTolerance(1E-20);
ff.compute(m1);
VectorXd                        res = ff.solveWithGuess(v1, x); // or ff.solve(v1);
return Rcpp::List::create(Rcpp::Named("res") = res);
}
