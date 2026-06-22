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

template <typename SolverType, typename MatrixType, typename RhsType, typename GuessType>
Rcpp::List run_iterative_solver(const MatrixType &A,
                                const RhsType &b,
                                const GuessType &x0,
                                double tol,
                                int maxit) {
  SolverType solver;

  if(tol > 0) {
    solver.setTolerance(tol);
  }
  if(maxit > 0) {
    solver.setMaxIterations(maxit);
  }

  solver.compute(A);
  if(solver.info() != Eigen::Success) {
    return Rcpp::List::create(
      Rcpp::Named("status") = false,
      Rcpp::Named("info") = static_cast<int>(solver.info()),
      Rcpp::Named("iterations") = solver.iterations(),
      Rcpp::Named("error") = solver.error());
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  bool converged = solver.info() == Eigen::Success;

  if(!x.allFinite()) {
    return Rcpp::List::create(
      Rcpp::Named("status") = false,
      Rcpp::Named("info") = static_cast<int>(solver.info()),
      Rcpp::Named("iterations") = solver.iterations(),
      Rcpp::Named("error") = solver.error());
  }

  return Rcpp::List::create(
    Rcpp::Named("status") = converged,
    Rcpp::Named("info") = static_cast<int>(solver.info()),
    Rcpp::Named("iterations") = solver.iterations(),
    Rcpp::Named("error") = solver.error(),
    Rcpp::Named("coefficients") = x);
}

// Solve Ax = b using sparse LeastSquaresConjugateGradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_sparse_LSCG(const Eigen::MappedSparseMatrix<double> &A,
                             const Eigen::Map<Eigen::MatrixXd> b,
                             const Eigen::Map<Eigen::MatrixXd> x0,
                             double tol = 0,
                             int maxit = 0,
                             unsigned int precond = 1
                             ) {

  if(precond == 0) {
    typedef Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
                                                 Eigen::IdentityPreconditioner> Solver;
    return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
  }

  if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  typedef Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
                                               Eigen::LeastSquareDiagonalPreconditioner<double> > Solver;
  return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
}

// Solve Ax = b using dense LeastSquaresConjugateGradient with given preconditioner with initial guess x0
// [[Rcpp::export]]
Rcpp::List solve_dense_LSCG(const Eigen::MatrixXd &A,
                            const Eigen::MatrixXd b,
                            const Eigen::MatrixXd x0,
                            double tol = 0,
                            int maxit = 0,
                            unsigned int precond = 1) {

  if(precond == 0) {
    typedef Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd,
                                                 Eigen::IdentityPreconditioner> Solver;
    return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
  }

  if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  typedef Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXd,
                                               Eigen::LeastSquareDiagonalPreconditioner<double> > Solver;
  return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
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

  if(precond == 0) {
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                     Eigen::Lower|Eigen::Upper,
                                     Eigen::IdentityPreconditioner> Solver;
    return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
  }

  if(precond == 2) {
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                     Eigen::Lower|Eigen::Upper,
                                     Eigen::IncompleteCholesky<double> > Solver;
    return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
  }

  if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }

  typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                   Eigen::Lower|Eigen::Upper,
                                   Eigen::DiagonalPreconditioner<double> > Solver;
  return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
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

  if(precond == 0) {
    typedef Eigen::ConjugateGradient<Eigen::MatrixXd,
                                     Eigen::Lower|Eigen::Upper,
                                     Eigen::IdentityPreconditioner> Solver;
    return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
  }

  if(precond == 2) {
    Eigen::SparseMatrix<double> A_sparse = A.sparseView();
    typedef Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                     Eigen::Lower|Eigen::Upper,
                                     Eigen::IncompleteCholesky<double> > Solver;
    return run_iterative_solver<Solver>(A_sparse, b, x0, tol, maxit);
  }

  if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested for dense CG -- using default.");
  }

  typedef Eigen::ConjugateGradient<Eigen::MatrixXd,
                                   Eigen::Lower|Eigen::Upper,
                                   Eigen::DiagonalPreconditioner<double> > Solver;
  return run_iterative_solver<Solver>(A, b, x0, tol, maxit);
}
