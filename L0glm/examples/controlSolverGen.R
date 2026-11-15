controlSolverGen(solver="osqp")

# CG uses incomplete Cholesky for both dense and sparse systems.
controlSolverGen("eigen", method=3L)
# Diagonal preconditioning is a cheaper alternative for fully dense systems.
controlSolverGen("eigen", method=3L, precond=1L)

# LSCG also defaults to incomplete Cholesky, with checked residuals and an
# LLT fallback. Select diagonal preconditioning to avoid forming its normal
# matrix for a rectangular solve.
controlSolverGen("eigen", method=4L)
controlSolverGen("eigen", method=4L, precond=1L)

# Optional constrained polish of nonnegative working ridge steps. This reuses
# an existing predictor Gram matrix and is off by default.
controlSolverGen("eigen", method=1L, nonnegative.polish=TRUE)
controlSolverGen("eigen", method=2L, nonnegative.polish=TRUE)

# Compare external backends using the same clamping strategy as default Eigen.
controlSolverGen("osqp", nonnegative.clamp=TRUE)
controlSolverGen("bigGlm", nonnegative.clamp=TRUE)
