# BENCHMARKS & TIMINGS OF L0GLM FITS BY PROBLEM SIZE AND SOLVER ####
# Uses high-resolution timing and the general R fitting path for each solver.
# "Dense" and "sparse" describe input storage: solvers may compress exact zeros.

library(L0glm)
library(ggplot2)
library(tidyr)
library(tidyselect)
library(microbenchmark)

output_dir <- normalizePath("..\\benchmarks\\1_benchmarks_solvers_problem_sizes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
                             
# Simulation sizes (with n = p)
sizes <- c(100, 500, 1000, 10000)

# Compare the general R fitting path for every solver. Native/OpenMP LLT
# uses a different outer loop and can compress dense inputs automatically.
use_openmp <- FALSE
blas_threads <- 1L
timing_repetitions <- 1L # increase for repeated fits; Time is the median
RhpcBLASctl::blas_set_num_threads(blas_threads)

criterion <- "gic"
metric_columns <- c(
  "Time",
  "Precision",
  "Recall",
  "Specificity",
  "Accuracy",
  "F1",
  "R2",
  "Estimated_over_true_support_size",
  "deviance"
)

solver_grid <- data.frame(
  Solver = c(
    "osqp",
    "bigGlm",
    "eigen Cholesky LLT",
    "eigen Cholesky LDLT",
    "eigen PCG incomplete Cholesky precond",
    "eigen PLSCG"
  ),
  Solver_backend = c("osqp", "bigGlm", rep("eigen", 4)),
  Solver_method = c(1L, 1L, 1L, 2L, 3L, 4L),
  stringsAsFactors = FALSE
)

solver_colours <- c(
  "osqp" = "blue",
  "bigGlm" = "red2",
  "eigen Cholesky LLT" = "green2",
  "eigen Cholesky LDLT" = "orange",
  "eigen PCG incomplete Cholesky precond" = "purple",
  "eigen PLSCG" = "cyan2"
)
solver_shapes <- stats::setNames(c(16, 17, 0, 15, 18, 4), names(solver_colours))

family_grid <- list(
  poisson = list(
    label = "Poisson",
    family = poisson(identity),
    simulation_method = 2L,
    start_mu_offset = 0.1
  ),
  gaussian = list(
    label = "Gaussian",
    family = gaussian(identity),
    simulation_method = 1L,
    start_mu_offset = 0
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

  timed_fit <- withCallingHandlers(
    tryCatch(
      microbenchmark::microbenchmark({
        fit <- L0glm::L0glm.fit(
          X = X,
          y = y,
          start = start_values,
          family = family,
          lambda = lambda,
          nonnegative = TRUE,
          relaxed = FALSE,
          thresh = 1e-5,
          openmp = use_openmp,
          parallelize = FALSE,
          nthreads = blas_threads,
          solver = solver_backend,
          method = solver_method,
          control.solver = do.call(
            L0glm::control.solver.gen,
            list(solver = solver_backend, method = solver_method)
          )
        )
      }, times = timing_repetitions),
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
      L0glm:::get.metrics(
        beta = fit$coefficients,
        X = X,
        y = y,
        family = family,
        beta_true = as.vector(beta_true)
      ),
      error = function(e) {
        metric_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(metrics)) {
      metric_values$Precision <- metrics$precision
      metric_values$Recall <- metrics$recall
      metric_values$Specificity <- metrics$specificity
      metric_values$Accuracy <- metrics$accuracy
      metric_values$F1 <- metrics$F1_score
      metric_values$R2 <- metrics$R2_coefs
      metric_values$Estimated_over_true_support_size <-
        metrics$estimated_over_true_support_size
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
      OpenMP = use_openmp,
      BLAS_threads = blas_threads,
      Timing_repetitions = timing_repetitions,
      Timing_clock = "microbenchmark",
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
  llt <- data$Solver == "eigen Cholesky LLT"
  # Draw LLT after every other solver's lines AND points, irrespective of
  # factor levels or row order. Keep the legend order unchanged.
  list(
    ggplot2::geom_line(data = data[!llt, ], linewidth = 1, na.rm = TRUE),
    ggplot2::geom_point(data = data[!llt, ], size = point_size,
                       stroke = point_stroke, na.rm = TRUE),
    ggplot2::geom_line(data = data[llt, ], linewidth = 1, na.rm = TRUE),
    ggplot2::geom_point(data = data[llt, ], size = point_size,
                       stroke = point_stroke, na.rm = TRUE)
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

      lambda_preset <- L0glm:::preset.lambda(
        IC = criterion,
        p = ncol(X),
        n = nrow(X)
      )
      start_values <- L0glm:::init_coefs(
        X = X,
        y = y,
        mu = y + family_config$start_mu_offset,
        family = family,
        lambda = lambda_preset
      )

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

metrics_file <- normalizePath(file.path(output_dir, "metrics_all_conditions.csv"))
write.csv(results, metrics_file, row.names = FALSE)
saveRDS(results, normalizePath(file.path(output_dir, "metrics_all_conditions.rds")))

results_long <- tidyr::pivot_longer(
  results,
  cols = tidyselect::all_of(metric_columns),
  names_to = "metric",
  values_to = "value"
)
results_long$Solver <- factor(results_long$Solver, levels = solver_grid$Solver)

metrics_long_file <- normalizePath(file.path(output_dir, "metrics_all_conditions_long.csv"))
write.csv(results_long, metrics_long_file, row.names = FALSE)

write.csv(
  data.frame(
    Solver = names(solver_colours),
    Colour = unname(solver_colours),
    stringsAsFactors = FALSE
  ),
  normalizePath(file.path(output_dir, "solver_colours.csv")),
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
      ggplot2::aes(x = Size, y = value, colour = Solver, group = Solver)
    ) +
      solver_layers(condition_long) +
      ggplot2::facet_wrap(~metric, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0(
          "Benchmarks different L0glm solvers for varying problem sizes\n",
          "(deconvolution of blurred spike train with ",
          family_label,
          " noise, ",
          matrix_label,
          " covariate matrix)"
        ),
        colour = "solver",
        x = "Problem size (n = p)",
        y = NULL
      ) +
      ggplot2::theme(legend.position = "right") +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::scale_color_manual(
        values = solver_colours,
        limits = names(solver_colours),
        drop = FALSE
      )

    save_plot(
      metrics_plot,
      sprintf("timings_metrics_%s_%s", matrix_name, family_name)
    )

    time_data <- subset(condition_long, metric == "Time")
    time_data$value[time_data$value <= 0] <- NA_real_

    timing_plot <- ggplot2::ggplot(
      time_data,
      ggplot2::aes(x = Size, y = value, colour = Solver, group = Solver)
    ) +
      solver_layers(time_data) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0(
          "Timings different L0glm solvers for varying problem sizes\n",
          "(deconvolution of blurred spike train with ",
          family_label,
          " noise, ",
          matrix_label,
          " covariate matrix)"
        ),
        colour = "solver",
        x = "Problem size (n = p)",
        y = "Time (s)"
      ) +
      ggplot2::theme(legend.position = "right") +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::scale_y_continuous(trans = "log10") +
      ggplot2::scale_color_manual(
        values = solver_colours,
        limits = names(solver_colours),
        drop = FALSE
      )

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
  ggplot2::aes(x = Size, y = value, colour = Solver, shape = Solver, group = Solver)
) +
  solver_layers(combined_timing_data, point_size = 2, point_stroke = 0.7) +
  ggplot2::facet_grid(Matrix_density ~ Noise_family) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "L0glm solver timings by problem size",
    subtitle = sprintf("openmp = %s; BLAS threads = %d; timed fits per condition = %d",
                       use_openmp, blas_threads, timing_repetitions),
    colour = "solver",
    shape = "solver",
    x = "Problem size (n = p)",
    y = "Time (s)"
  ) +
  ggplot2::theme(legend.position = "bottom", panel.spacing = grid::unit(1, "lines")) +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2),
                  shape = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::scale_x_continuous(trans = "log10", breaks = sizes,
                            labels = scales::label_number(big.mark = ",")) +
  ggplot2::scale_y_continuous(trans = "log10") +
  ggplot2::scale_shape_manual(values = solver_shapes,
                             limits = names(solver_shapes), drop = FALSE) +
  ggplot2::scale_color_manual(
    values = solver_colours,
    limits = names(solver_colours),
    drop = FALSE
  )

save_plot(combined_timing_plot, "timings_all_conditions", width = 10, height = 7)

run_summary <- c(
  paste("Output directory:", output_dir),
  paste("Wide metrics table:", metrics_file),
  paste("Long metrics table:", metrics_long_file),
  paste("OpenMP/native dispatch:", use_openmp),
  paste("BLAS threads:", blas_threads),
  paste("High-resolution timing repetitions:", timing_repetitions),
  paste("Run completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(run_summary, normalizePath(file.path(output_dir, "run_summary.txt")))
message(paste(run_summary, collapse = "\n"))


# TODOS ####
#
# - For Gaussian noise, check whether the bigGlm solution quality can be
#   improved; in the original benchmark notes it was markedly worse than the
#   eigen and osqp solutions.
# - Check whether the PLSCG solution quality can be made to match the PCG
#   solver with incomplete Cholesky preconditioning, possibly by tuning solver
#   settings or adding an incomplete Cholesky preconditioner in the Rcpp solver.
# - Check whether all eigen solvers can be tuned to match osqp solution quality
#   more closely, especially in the Gaussian cases where osqp had lower
#   deviance but was much slower.
# - Check whether dense Gaussian eigen solver settings can be tuned to get
#   closer to the osqp solution while retaining the eigen solvers' speed.
