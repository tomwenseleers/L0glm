# L0glm 0.2.1

## Development changes

* Removed `lower` and `upper` from `L0glm.fit()` and `L0nnpois()`. General box
  constraints are no longer part of either fitting interface.
* `L0glm.fit(nonnegative = TRUE)` constrains all coefficients, including any
  intercept and `always.include` terms. Binomial and quasibinomial intercepts
  remain unrestricted; columns of ones supplied in `X` are recognized too.
  The default `nonnegative = FALSE` allows coefficients of either sign.
* Removed `nonnegative` from `L0nnpois()`: its coefficients, including any
  intercept, are now always nonnegative. Remove this argument from existing
  calls, and use `L0glm.fit(nonnegative = FALSE)` for unrestricted fits.
* Applied the same constraint rules to initialization, grouped and relaxed
  fits, and serial and parallel execution. The NNLS ridge backend now handles
  unrestricted coefficients, including the binomial intercept.
* Updated help, examples and benchmarks to use the simplified interfaces.
