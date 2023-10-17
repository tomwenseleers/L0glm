#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// Solve Ax = b using iterative methods (Least Squares Conjugate Gradient or
// Preconditioned Conjugate Gradient) for sparse or dense covariate matrix A
// adapted from https://github.com/nk027/sanic/blob/master/src/solve_iterative.cpp

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MappedSparseMatrix;
using Eigen::DiagonalPreconditioner;
using Eigen::LeastSquaresConjugateGradient;
using Eigen::Success;

// Solve Ax = b using sparse LeastSquaresConjugateGradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_sparse_LSCG(const Eigen::MappedSparseMatrix<double> &A,
                             const Eigen::Map<Eigen::MatrixXd> b,
                             const Eigen::Map<Eigen::MatrixXd> x0,
                             double tol = 0,
                             int maxit = 0,
                             unsigned int precond = 1
                             ) {

  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
                                       Eigen::LeastSquareDiagonalPreconditioner<double> >
                                         solver;
  if(precond == 0) {
    Eigen::LeastSquaresConjugateGradient < Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > solver;
  } else if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  if(tol) {
    solver.setTolerance(tol);
  }
  if(maxit) {
    solver.setMaxIterations(maxit);
  }

  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    // decomposition failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}

// Solve Ax = b using dense LeastSquaresConjugateGradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_dense_LSCG(const Eigen::MatrixXd &A,
                            const Eigen::MatrixXd b,
                            const Eigen::MatrixXd x0,
                            double tol = 0,
                            int maxit = 0,
                            unsigned int precond = 1) {

  Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd,
                                       Eigen::LeastSquareDiagonalPreconditioner<double> >
                                         solver;

  if(precond == 0) {
    Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd, Eigen::IdentityPreconditioner > solver;
  } else if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  if(tol) {
    solver.setTolerance(tol);
  }
  if(maxit) {
    solver.setMaxIterations(maxit);
  }

  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    // decomposition failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}


// Solve Ax = b using sparse Conjugate Gradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_sparse_CG(
    Eigen::MappedSparseMatrix<double> &A,
    Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,
    double tol = 0,
    int maxit = 0,
    unsigned int precond = 1) {

  // Set-up the Conjugate Gradient solver
  Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double> > solver;
  if(precond == 0) {
    Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner > solver;
  } else if(precond == 2) {
    Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > solver;
  } else if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  if(tol) {
    solver.setTolerance(tol);
  }
  if(maxit) {
    solver.setMaxIterations(maxit);
  }

  solver.compute(A);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}

// Solve Ax = b using dense Conjugate Gradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_dense_CG(
    Eigen::Map<Eigen::MatrixXd> &A,
    Eigen::Map<Eigen::MatrixXd> b,
    Eigen::Map<Eigen::MatrixXd> x0,
    double tol = 0,
    int maxit = 0,
    unsigned int precond = 1) {

  // Set-up the Conjugate Gradient solver
  Eigen::ConjugateGradient < Eigen::MatrixXd, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double> > solver;
  if(precond == 0) {
    Eigen::ConjugateGradient < Eigen::MatrixXd, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner > solver;
  } else if(precond == 2) {
    Eigen::ConjugateGradient < Eigen::MatrixXd, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > solver;
  } else if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  if(tol) {
    solver.setTolerance(tol);
  }
  if(maxit) {
    solver.setMaxIterations(maxit);
  }

  solver.compute(A);

  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  if(solver.info() != Eigen::Success) {
    // solving failed
    return Rcpp::List::create(Rcpp::Named("status") = false);
  }

  return Rcpp::List::create(Rcpp::Named("status") = true,
                            Rcpp::Named("coefficients") = x);
}
