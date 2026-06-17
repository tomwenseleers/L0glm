if (requireNamespace("cobs", quietly = TRUE)) {
set.seed(1)
x      <- seq(-5, 5, length.out = 101)

## underlying "true" density on the grid
f_true <- dnorm(x)
f_true <- f_true / sum(f_true)          # ensure it integrates to 1

N      <- 100000                        # total expected counts
lambda <- N * f_true                    # Poisson means in each bin

y_cnt  <- rpois(length(x), lambda)      # Poisson-noisy counts
y      <- y_cnt / sum(y_cnt)            # rescale back to a density

g <- log_conc_smooth(y)

matplot(cbind(y, f_true, g), type = "l", lty = 1,
        col = c(2, 1, 4), ylab = "density",
        main = "Raw (Poisson) vs. log-concave smoothed curve + ground truth")
legend("topleft", legend = c("raw (Poisson)", "smoothed", "ground truth"),
       col = c(2, 4, 1), lty = 1, bty = "n")
}
