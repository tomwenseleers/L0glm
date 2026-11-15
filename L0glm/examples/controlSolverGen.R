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
controlSolverGen("eigen", method=1L, nonnegative.method="constrained")
controlSolverGen("eigen", method=2L, nonnegative.method="constrained")
# LSCG forms a memory-checked predictor Gram only when a bound correction
# is needed; its initial least-squares solve and LLT fallback are retained.
controlSolverGen("eigen", method=4L, nonnegative.method="constrained")

# Compare external backends using the same clamping strategy as default Eigen.
controlSolverGen("osqp", nonnegative.method="clamp")
controlSolverGen("bigGlm", nonnegative.method="clamp")
