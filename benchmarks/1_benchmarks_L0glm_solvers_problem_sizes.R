# BENCHMARKS & TIMINGS OF L0GLM FITS BY PROBLEM SIZE AND SOLVER ####
# Uses automatic native fitting paths and repeated warm timings.
# Shared unit starts and disabled support search isolate the ridge solvers.
# "Dense" and "sparse" describe input storage: solvers may compress exact zeros.

library(L0glm)
library(ggplot2)
library(tidyr)
library(tidyselect)
library(microbenchmark)

repository_dir <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
output_dir <- file.path(repository_dir, "benchmarks", "1_benchmarks_solvers_problem_sizes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
                             
# Simulation sizes (with n = p)
sizes <- c(100, 500, 1000, 10000)

# Serial outcome fits give every solver the same thread budget.
ncores_fit <- 1L
blas_threads <- 1L
timing_repetitions <- 3L # median after one untimed warm fit per condition
# Generous benchmark-only ceiling: fits stop as soon as they converge.
# The package default remains 25 iterations.
fit_maxit <- 500L
fit_tol <- 1e-10
# Optional constrained polish for Eigen paths with an existing predictor Gram.
# Rectangular LSCG and Woodbury solves bypass this correction.
# Keep unit starts shared across all solvers when comparing this switch.
nonnegative_polish <- FALSE
# Solve OSQP/bigGlm ridge steps without bounds, then clamp as Eigen does.
# Keep FALSE to compare their exact constraints with default Eigen clamping.
nonnegative_clamp <- FALSE
RhpcBLASctl::blas_set_num_threads(blas_threads)

criterion <- "gic"
metric_columns <- c(
  "adj_R2_weighted",
  "GIC",
  "Time",
  "Precision",
  "Recall",
  "F1",
  "R2",
  "RMSE_coefficients",
  "deviance"
)

solver_grid <- data.frame(
  Solver = c(
    "OSQP + polish",
    "bigGlm + polish",
    "Eigen LLT",
    "Eigen LDLT",
    "Eigen CG-IC",
    "Eigen LSCG-IC"
  ),
  Solver_backend = c("osqp", "bigGlm", rep("eigen", 4)),
  Solver_method = c(1L, 1L, 1L, 2L, 3L, 4L),
  stringsAsFactors = FALSE
)

solver_colours <- c(
  "OSQP + polish" = "blue",
  "bigGlm + polish" = "red2",
  "Eigen LLT" = "green2",
  "Eigen LDLT" = "orange",
  "Eigen CG-IC" = "purple",
  "Eigen LSCG-IC" = "cyan2"
)
solver_shapes <- stats::setNames(c(1, 17, 0, 15, 18, 4), names(solver_colours))
solver_linetypes <- stats::setNames(
  c("dashed", "solid", "solid", "dotted", "longdash", "dotdash"),
  names(solver_colours))
solver_labels <- stats::setNames(
  c("OSQP + polish", "bigGlm + polish", "Eigen LLT", "Eigen LDLT", "Eigen CG-IC", "Eigen LSCG-IC"),
  solver_grid$Solver)
metric_labels <- c(adj_R2_weighted="Adjusted weighted prediction R-squared",
  GIC="GIC", Time="Time (s)", Precision="Precision", Recall="Recall",
  F1="F1", R2="Coefficient R-squared",
  RMSE_coefficients="Coefficient RMSE", deviance="True-signal deviance / n")

family_grid <- list(
  poisson = list(
    label = "Poisson",
    family = poisson(identity),
    simulation_method = 2L
  ),
  gaussian = list(
    label = "Gaussian",
    family = gaussian(identity),
    simulation_method = 1L
  )
)

matrix_grid <- list(
  sparse = list(label = "sparse", sparse = TRUE),
  dense = list(label = "dense", sparse = FALSE)
)

empty_metrics <- function() {
  as.list(stats::setNames(rep(NA_real_, length(metric_columns)), metric_columns))
}

fit_one_solver <- function(X, y, beta_true, family, lambda, solver_backend,
                           solver_method, solver_name, start_values) {
  fit <- NULL
  fit_time <- NA_real_
  fit_error <- NA_character_
  fit_warnings <- character()
  metric_values <- empty_metrics()

  run_fit <- function() L0glm::l0glmFit(
    X = X, y = y, start = start_values, family = family, lambda = lambda,
    nonnegative = TRUE, local.search = FALSE,
    maxit = fit_maxit, tol = fit_tol, thresh = 1e-5,ncores=ncores_fit, nthreads = blas_threads,
    solver = solver_backend, method = solver_method,
    control.solver = L0glm::controlSolverGen(solver_backend, solver_method,
      nonnegative.polish = solver_backend == "eigen" && nonnegative_polish,
      nonnegative.clamp = solver_backend != "eigen" && nonnegative_clamp))

  timed_fit <- withCallingHandlers(
    tryCatch(
      {
        fit <- run_fit()
        microbenchmark::microbenchmark({fit <- run_fit()}, times = timing_repetitions)
      },
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (!is.null(timed_fit)) {
    # microbenchmark uses a high-resolution timer; system.time() rounded some
    # successful fast fits to exactly zero, which cannot appear on a log axis.
    fit_time <- median(timed_fit$time) / 1e9
    metric_values$Time <- fit_time
  }

  if (!is.null(fit)) {
    metric_error <- NA_character_
    metrics <- tryCatch(
      L0glm::l0glmMetrics(fit, beta_true = as.vector(beta_true)),
      error = function(e) {
        metric_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(metrics)) {
      metric_values$adj_R2_weighted <- metrics$adj_R2_weighted
      metric_values$GIC <- metrics$gic
      metric_values$Precision <- metrics$precision
      metric_values$Recall <- metrics$recall
      metric_values$F1 <- metrics$F1_score
      metric_values$R2 <- metrics$R2_coefs
      metric_values$RMSE_coefficients <- metrics$RMSE_coefs
      metric_values$deviance <- metrics$deviance_true / nrow(X)
    } else {
      fit_error <- paste(
        na.omit(c(fit_error, paste("Metric error:", metric_error))),
        collapse = " | "
      )
    }
  }

  status <- if (is.null(fit)) {
    "fit_error"
  } else if (!is.na(fit_error) && nzchar(fit_error)) {
    "metric_error"
  } else {
    "ok"
  }

  c(
    list(
      Solver = solver_name,
      Solver_backend = solver_backend,
      Solver_method = solver_method,
      Native_driver = TRUE,
      Ncores = ncores_fit,
      BLAS_threads = blas_threads,
      Timing_repetitions = timing_repetitions,
      Warmup_fits = 1L,
      Timing_clock = "microbenchmark",
      Start = "1",
      Local_search = FALSE,
      Nonnegative_polish = solver_backend == "eigen" && nonnegative_polish,
      Nonnegative_clamp = solver_backend != "eigen" && nonnegative_clamp,
      Maxit = fit_maxit,
      Tolerance = fit_tol,
      Iterations = if (is.null(fit)) NA_integer_ else max(fit$convergence$iterations),
      Converged = if (is.null(fit)) NA else all(fit$convergence$converged),
      Status = status,
      Warnings = if (length(fit_warnings) > 0) {
        paste(unique(fit_warnings), collapse = " | ")
      } else {
        NA_character_
      },
      Error = fit_error
    ),
    metric_values
  )
}

solver_layers <- function(data, point_size = 1.5, point_stroke = 0.5) {
  llt <- data$Solver_backend == "eigen" & data$Solver_method == 1L
  osqp <- data$Solver_backend == "osqp"
  bigglm <- data$Solver_backend %in% c("bigGlm", "bigglm", "glmnet")
  # Solid bigGlm and dashed OSQP are drawn last, so coincident constrained
  # solutions remain visible. Routing identifiers are independent of labels.
  unlist(lapply(list(!llt & !osqp & !bigglm, llt, bigglm, osqp), function(rows) list(
    ggplot2::geom_line(data = data[rows, ], linewidth = 0.85, na.rm = TRUE),
    ggplot2::geom_point(data = data[rows, ], size = point_size,
                       stroke = point_stroke, na.rm = TRUE)
  )), recursive = FALSE)
}

solver_scales <- function() {
  list(
    ggplot2::scale_colour_manual(values = solver_colours, labels = solver_labels,
      limits = names(solver_colours), drop = FALSE),
    ggplot2::scale_shape_manual(values = solver_shapes, labels = solver_labels,
      limits = names(solver_shapes), drop = FALSE),
    ggplot2::scale_linetype_manual(values = solver_linetypes, labels = solver_labels,
      limits = names(solver_linetypes), drop = FALSE)
  )
}

save_plot <- function(plot, filename_stem, width = 7, height = 5, dpi = 600) {
  write_plot <- function(filename) {
    withCallingHandlers(
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi),
      warning = function(w) {
        # Some PNG devices signal a failed write as a warning, not an error.
        if (grepl("could not write|cannot open file", conditionMessage(w), ignore.case = TRUE))
          stop(conditionMessage(w), call. = FALSE)
      }
    )
  }
  for (extension in c("png", "pdf")) {
    filename <- file.path(output_dir, paste0(filename_stem, ".", extension))
    tryCatch(
      write_plot(filename),
      error = function(e) {
        # A PDF open in a Windows viewer may be locked against replacement.
        alternate <- file.path(output_dir, paste0(filename_stem, "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"), ".", extension))
        warning(sprintf("Could not save %s (%s); saving %s instead.",
                        basename(filename), conditionMessage(e), basename(alternate)),
                call. = FALSE)
        write_plot(alternate)
      }
    )
  }
}

result_rows <- list()

for (family_name in names(family_grid)) {
  family_config <- family_grid[[family_name]]
  family <- family_config$family

  for (matrix_name in names(matrix_grid)) {
    matrix_config <- matrix_grid[[matrix_name]]

    for (n in sizes) {
      message(
        sprintf(
          "Simulating %s noise, %s matrix, n = p = %s",
          family_config$label,
          matrix_config$label,
          n
        )
      )

      sim <- L0glm::simul(
        method = family_config$simulation_method,
        n = n,
        p = n,
        sparse = matrix_config$sparse,
        seed = 1L,
        options = list(sd_noise = 1)
      )

      X <- sim$X
      y <- sim$y
      beta_true <- sim$beta_true

      lambda_preset <- L0glm::presetLambda(
        IC = criterion,
        p = ncol(X),
        n = nrow(X)
      )
      start_values <- 1

      for (solver_idx in seq_len(nrow(solver_grid))) {
        solver_spec <- solver_grid[solver_idx, ]

        message(
          sprintf(
            "  Fitting %s, %s matrix, n = %s, solver = %s",
            family_name,
            matrix_name,
            n,
            solver_spec$Solver
          )
        )

        fit_result <- fit_one_solver(
          X = X,
          y = y,
          beta_true = beta_true,
          family = family,
          lambda = lambda_preset,
          solver_backend = solver_spec$Solver_backend,
          solver_method = solver_spec$Solver_method,
          solver_name = solver_spec$Solver,
          start_values = start_values
        )

        result_rows[[length(result_rows) + 1L]] <- data.frame(
          Noise_family = family_name,
          Link = family$link,
          Matrix_density = matrix_name,
          Sparse = matrix_config$sparse,
          Size = n,
          Lambda = lambda_preset,
          fit_result,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  }
}

results <- do.call(rbind, result_rows)
results$Noise_family <- factor(results$Noise_family, levels = names(family_grid))
results$Matrix_density <- factor(results$Matrix_density, levels = names(matrix_grid))
results$Solver <- factor(results$Solver, levels = solver_grid$Solver)
results <- results[order(
  results$Noise_family,
  results$Matrix_density,
  results$Size,
  results$Solver
), ]

metrics_file <- file.path(output_dir, "metrics_all_conditions.csv")
write.csv(results, metrics_file, row.names = FALSE)
saveRDS(results, file.path(output_dir, "metrics_all_conditions.rds"))

results_long <- tidyr::pivot_longer(
  results,
  cols = tidyselect::all_of(metric_columns),
  names_to = "metric",
  values_to = "value"
)
results_long$Solver <- factor(results_long$Solver, levels = solver_grid$Solver)
results_long$metric <- factor(results_long$metric, levels = metric_columns)

metrics_long_file <- file.path(output_dir, "metrics_all_conditions_long.csv")
write.csv(results_long, metrics_long_file, row.names = FALSE)

write.csv(
  data.frame(
    Solver = names(solver_colours),
    Colour = unname(solver_colours),
    Shape = unname(solver_shapes),
    Line_type = unname(solver_linetypes),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "solver_colours.csv"),
  row.names = FALSE
)

for (family_name in names(family_grid)) {
  family_label <- family_grid[[family_name]]$label

  for (matrix_name in names(matrix_grid)) {
    matrix_label <- matrix_grid[[matrix_name]]$label
    condition_long <- subset(
      results_long,
      Noise_family == family_name & Matrix_density == matrix_name
    )

    metrics_plot <- ggplot2::ggplot(
      condition_long,
      ggplot2::aes(x = Size, y = value, colour = Solver, shape = Solver,
                   linetype = Solver, group = Solver)
    ) +
      solver_layers(condition_long) +
      ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3,
                         labeller = ggplot2::as_labeller(metric_labels)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = sprintf("%s deconvolution (%s design)", family_label, matrix_label),
        subtitle = sprintf("start = 1; local.search = FALSE; GIC preset; maxit = %d; Eigen polish = %s",
                           fit_maxit, nonnegative_polish),
        caption = if(nonnegative_clamp) "OSQP/bigGlm: unconstrained ridge steps followed by clamping." else NULL,
        colour = "solver",
        shape = "solver", linetype = "solver",
        x = "Problem size (n = p)",
        y = NULL
      ) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3),
                     shape = ggplot2::guide_legend(ncol = 3),
                     linetype = ggplot2::guide_legend(ncol = 3)) +
      ggplot2::scale_x_continuous(trans = "log10", breaks = sizes,
                                labels = scales::label_number(big.mark = ",")) +
      solver_scales()

    save_plot(
      metrics_plot,
      sprintf("timings_metrics_%s_%s", matrix_name, family_name), width = 11, height = 8
    )

    time_data <- subset(condition_long, metric == "Time")
    time_data$value[time_data$value <= 0] <- NA_real_

    timing_plot <- ggplot2::ggplot(
      time_data,
      ggplot2::aes(x = Size, y = value, colour = Solver, shape = Solver,
                   linetype = Solver, group = Solver)
    ) +
      solver_layers(time_data) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = sprintf("%s deconvolution (%s design)", family_label, matrix_label),
        subtitle = sprintf("start = 1; local.search = FALSE; maxit = %d; median warm fitting time",
                           fit_maxit),
        caption = if(nonnegative_clamp) "OSQP/bigGlm: unconstrained ridge steps followed by clamping." else NULL,
        colour = "solver",
        shape = "solver", linetype = "solver",
        x = "Problem size (n = p)",
        y = "Time (s)"
      ) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3),
                     shape = ggplot2::guide_legend(ncol = 3),
                     linetype = ggplot2::guide_legend(ncol = 3)) +
      ggplot2::scale_x_continuous(trans = "log10", breaks = sizes,
                                labels = scales::label_number(big.mark = ",")) +
      ggplot2::scale_y_continuous(trans = "log10") +
      solver_scales()

    save_plot(
      timing_plot,
      sprintf("timings_%s_%s", matrix_name, family_name)
    )
  }
}

combined_timing_data <- subset(results_long, metric == "Time")
combined_timing_data$value[combined_timing_data$value <= 0] <- NA_real_
combined_timing_plot <- ggplot2::ggplot(
  combined_timing_data,
  ggplot2::aes(x = Size, y = value, colour = Solver, shape = Solver,
               linetype = Solver, group = Solver)
) +
  solver_layers(combined_timing_data, point_size = 2, point_stroke = 0.7) +
  ggplot2::facet_grid(Matrix_density ~ Noise_family) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "L0glm solver timings by problem size",
    subtitle = sprintf(paste0("start = 1; local.search = FALSE; maxit = %d; Eigen polish = %s\n",
                              "native drivers; ncores = %d; BLAS threads = %d; median of %d warm fits"),
                       fit_maxit, nonnegative_polish, ncores_fit, blas_threads, timing_repetitions),
    caption = if(nonnegative_clamp) "OSQP/bigGlm: unconstrained ridge steps followed by clamping." else NULL,
    colour = "solver",
    shape = "solver",
    linetype = "solver",
    x = "Problem size (n = p)",
    y = "Time (s)"
  ) +
  ggplot2::theme(legend.position = "bottom", panel.spacing = grid::unit(1, "lines")) +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2),
                  shape = ggplot2::guide_legend(ncol = 2),
                  linetype = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::scale_x_continuous(trans = "log10", breaks = sizes,
                            labels = scales::label_number(big.mark = ",")) +
  ggplot2::scale_y_continuous(trans = "log10") +
  solver_scales()

save_plot(combined_timing_plot, "timings_all_conditions", width = 10, height = 7)

run_summary <- c(
  paste("Output directory:", output_dir),
  paste("Wide metrics table:", metrics_file),
  paste("Long metrics table:", metrics_long_file),
  "Fitting driver: automatic native dispatch",
  paste("Outcome thread budget:", ncores_fit),
  paste("BLAS threads:", blas_threads),
  paste("High-resolution timing repetitions:", timing_repetitions),
  "Untimed warm fits per condition: 1",
  "Starting coefficients: 1 (working coefficient coordinates)",
  "Automatic local support search: FALSE",
  paste("Optional Eigen constrained polish:", nonnegative_polish),
  paste("Optional OSQP/bigGlm unconstrained solves followed by clamping:", nonnegative_clamp),
  "OSQP: eps_abs = eps_rel = 1e-9; diagonal equilibration; checked Cholesky polish",
  "Prediction R-squared: adjusted, final GLM working weights; Poisson mean floor 0.1",
  "GIC: fitted dispersion, coefficient-count degrees of freedom, prediction mean floor 1e-8",
  "Coincident constrained solutions: dashed OSQP overlays solid bigGlm",
  paste("Maximum IRLS iterations:", fit_maxit),
  paste("Converged fits:", sum(results$Converged, na.rm = TRUE), "of", nrow(results)),
  paste("Largest observed iteration count:", max(results$Iterations, na.rm = TRUE)),
  paste("Convergence tolerance:", fit_tol),
  paste("L0glm version:", packageVersion("L0glm")),
  paste("glmnet version:", if (isNamespaceLoaded("glmnet"))
    as.character(getNamespaceVersion("glmnet")) else "not loaded"),
  paste("R version:", R.version.string),
  paste("Run completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(run_summary, file.path(output_dir, "run_summary.txt"))
message(paste(run_summary, collapse = "\n"))
if (any(!results$Converged, na.rm = TRUE))
  warning("Some fits did not converge within the benchmark ceiling; inspect Converged and Iterations in the metrics table.",
          call. = FALSE)
