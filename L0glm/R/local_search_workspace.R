# A screened search changes coefficients only inside its retained pool. Work in
# those coordinates, preserving fit$p for every likelihood-complexity charge.
L0glm_local_workspace <- function(fit, args, control) {
  candidates <- fit$screening$candidates
  source <- inherits(args$X,'L0glm_data')
  if (!source && (control$scope == "all" || !length(candidates))) return(NULL)
  added <- as.integer(isTRUE(fit$added_intercept))
  beta <- coef(fit)
  rows <- sort(unique(c(fit$always.include, L0glm_local_pattern(fit)$rows,
                        candidates + added)))
  if(source && (control$scope=='all'||!length(candidates)))rows <- seq_len(nrow(beta))
  if (!source && length(rows) > nrow(beta) / 2) return(NULL)
  original <- fit
  fit$X <- fit$X[, rows, drop=FALSE]
  fit$coefficients <- beta[rows, , drop=FALSE]
  fit$always.include <- match(original$always.include, rows)
  fit$screening$candidates <- match(candidates + added, rows) - added
  for (name in c("center", "scale"))
    fit$preprocessing[[name]] <- fit$preprocessing[[name]][rows]
  fit$preprocessing$p <- length(rows) - added
  if (!is.null(fit$lambdas))
    fit$lambdas <- matrix(as.numeric(fit$lambdas), nrow(beta), ncol(beta))[rows, , drop=FALSE]
  if (!is.null(fit$selected)) fit$selected <- match(fit$selected, rows)
  if (!is.null(args$validation)) {
    context <- args$validation
    if (!inherits(context, "L0glm_validation"))
      context <- L0glm_validation_context(context, args$X, fit$y, fit$family, added == 1L)
    args$validation <- L0glm_subset_validation(context, columns=rows)
  }
  args$X <- L0glm_data_matrix(args$X[, (rows - added)[rows > added], drop=FALSE])
  if (!is.null(args$start) && length(args$start) > 1L) {
    start <- as.matrix(args$start)
    if (added && nrow(start) == nrow(beta) - 1L) start <- rbind(0, start)
    args$start <- start[rows, , drop=FALSE]
  }
  list(fit=fit, args=args, rows=rows, original=original)
}

L0glm_local_restore <- function(fit, workspace) {
  rows <- workspace$rows
  original <- workspace$original
  kind <- L0glm_ext_kind(original$family)
  quasi <- startsWith(original$family$family, "quasi")
  calibrated <- identical(original$dispersionCalibration$method,"cross-channel") &&
    isTRUE(original$dispersionCalibration$applied)
  d <- list(D=original$X, p=nrow(coef(original)), m=ncol(coef(original)),
    mandatory=original$always.include, grouped=isTRUE(original$grouped), kind=kind,
    estimated=original$dispersion.estimated,
    common.phi=if(is.null(kind) && (quasi||calibrated))
      pmax(abs(original$dispersion), 1e-8) else numeric())
  search <- fit$local.search
  fit <- if (nrow(search$history)>1L) L0glm_local_remap(fit, rows, d, original) else original
  fit$local.search <- search
  move <- function(value) vapply(strsplit(value, ":", fixed=TRUE), function(parts) {
    positions <- if (parts[1L] == "swap") 2:3 else
      if (parts[1L] %in% c("add", "drop", "add_entry", "drop_entry")) 2L else integer()
    parts[positions] <- as.character(rows[as.integer(parts[positions])])
    paste(parts, collapse=":")
  }, character(1), USE.NAMES=FALSE)
  for (name in c("history", "proposals")) if (nrow(fit$local.search[[name]]))
    fit$local.search[[name]]$move <- move(fit$local.search[[name]]$move)
  fit$local.search$pool <- rows[fit$local.search$pool]
  if (length(fit$local.refit$rows)) fit$local.refit$rows <- rows[fit$local.refit$rows]
  fit$local.search$workspace <- c(original.rows=nrow(coef(original)), working.rows=length(rows))
  if (inherits(fit, "L0glm_extended_fit")) {
    fit$metrics <- l0glmMetrics(fit, mu_floor=fit$local.search$scoring.floor)
    fit$objective <- fit$metrics$objective
  }
  fit
}
