controlSolverGen(solver="osqp")

# CG uses incomplete Cholesky for both dense and sparse systems.
controlSolverGen("eigen", method=3L)
# Diagonal preconditioning is a cheaper alternative for fully dense systems.
controlSolverGen("eigen", method=3L, precond=1L)
