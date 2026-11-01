#' Write a chromatogram matrix as an ANDI-MS CDF file
#'
#' Export scan-by-mass intensities for [GCMSdeco()] or external GC/MS software.
#' The NetCDF classic file includes scan indexing, total ion counts, mass ranges
#' and ANDI-MS administrative and instrument records.
#'
#' @param y Finite nonnegative numeric matrix, with scans in rows and nominal
#'   mass channels in columns. Zero channels are retained.
#' @param cdf_file Output `.cdf` path.
#' @param time Increasing scan times in minutes. Defaults to numeric row names.
#' @param mass Increasing nonnegative integer m/z values. Defaults to numeric
#'   column names.
#' @param sample_name Sample identity recorded in the header; defaults to the
#'   output filename without its extension.
#' @return Invisibly, the normalized output path.
#' @details Intensities are stored without rescaling or rounding. Integer data
#'   fitting in a signed 32-bit integer use ANDI `Long` storage; other data use
#'   `Double`. Scan times are converted to seconds in the file. Every channel
#'   is written in every scan, including exact zeros. The file's timestamp
#'   records export time; acquisition metadata are left unspecified.
#' @seealso [GCMSdeco()], [splogconNMF()]
#' @examples
#' \donttest{
#' if (requireNamespace("ncdf4", quietly = TRUE)) {
#'   data(GCMS_data1)
#'   cdf <- tempfile(fileext = ".cdf")
#'   write_GCMS_cdf(GCMS_data1, cdf)
#'   unlink(cdf)
#' }
#' }
#' @md
#' @export
write_GCMS_cdf <- function(y, cdf_file, time = as.numeric(rownames(y)),
                           mass = as.numeric(colnames(y)),
                           sample_name = tools::file_path_sans_ext(basename(cdf_file))) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("Install ncdf4 to write ANDI-MS files.", call. = FALSE)
  if (!is.matrix(y) || !is.numeric(y) || nrow(y) < 2L || !ncol(y) ||
      any(!is.finite(y) | y < 0)) stop("y must be a finite nonnegative scan-by-mass matrix with at least two scans.", call. = FALSE)
  n <- nrow(y);m <- ncol(y)
  if (!is.numeric(time) || length(time) != n || any(!is.finite(time) | time < 0) || any(diff(time) <= 0))
    stop("Supply increasing finite scan times in minutes.", call. = FALSE)
  if (!is.numeric(mass) || length(mass) != m || any(!is.finite(mass) | mass < 0 | mass != floor(mass)) || any(diff(mass) <= 0))
    stop("Supply increasing nonnegative integer m/z values.", call. = FALSE)
  if (!is.character(cdf_file) || length(cdf_file) != 1L || is.na(cdf_file) || !grepl("\\.cdf$", cdf_file, ignore.case = TRUE))
    stop("cdf_file must end in .cdf.", call. = FALSE)
  if (!is.character(sample_name) || length(sample_name) != 1L || is.na(sample_name) || !nzchar(sample_name))
    stop("sample_name must be a nonempty string.", call. = FALSE)
  points <- as.double(n) * m
  if (points > .Machine$integer.max) stop("The ANDI-MS scan index exceeds signed 32-bit storage.", call. = FALSE)
  seconds <- time * 60;dt <- c(diff(seconds), stats::median(diff(seconds)))
  integer_y <- all(y == round(y)) && max(y) <= .Machine$integer.max
  dim <- function(name, size, unlimited = FALSE)
    ncdf4::ncdim_def(name, "", seq_len(size), unlim = unlimited, create_dimvar = FALSE)
  scan <- dim("scan_number", n);point <- dim("point_number", points, TRUE)
  string32 <- dim("_32_byte_string", 32L);string64 <- dim("_64_byte_string", 64L)
  instrument <- dim("instrument_number", 1L);error <- dim("error_number", 1L)
  var <- function(name, dims, type = "double", units = "") ncdf4::ncvar_def(name, units, dims, prec = type)
  instruments <- c("instrument_name", "instrument_id", "instrument_mfr", "instrument_model",
    "instrument_sw_version", "instrument_os_version")
  scans_int <- c("scan_index", "point_count", "flag_count", "scan_type")
  scans_double <- c("a_d_sampling_rate", "scan_acquisition_time", "scan_duration", "mass_range_min",
    "mass_range_max", "resolution", "total_intensity")
  variables <- c(list(var("error_log", list(string64, error), "char")),
    lapply(instruments, var, dims = list(string32, instrument), type = "char"),
    lapply(scans_int, var, dims = list(scan), type = "integer"),
    lapply(scans_double, function(name) var(name, list(scan), units = switch(name,
      scan_acquisition_time = "seconds", total_intensity = "Total Counts", ""))),
    list(var("mass_values", list(point), units = "M/Z"),
      var("intensity_values", list(point), if (integer_y) "integer" else "double", "Arbitrary Intensity Units")))
  dir.create(dirname(cdf_file), recursive = TRUE, showWarnings = FALSE)
  nc <- ncdf4::nc_create(cdf_file, variables, force_v4 = FALSE)
  on.exit(ncdf4::nc_close(nc))
  stamp <- format(Sys.time(), "%Y%m%d%H%M%S%z")
  metadata <- list(dataset_completeness = "C1+C2", ms_template_revision = "1.0.1",
    netcdf_revision = "2.3.2", dataset_origin = "L0glm", dataset_owner = "",
    experiment_title = sample_name, experiment_type = "Centroided Mass Spectrum",
    netcdf_file_date_time_stamp = stamp, source_file_reference = sample_name,
    source_file_format = "ANDI-MS matrix export", operator_name = "", languages = "English",
    administrative_comments = "Intensity matrix exported by L0glm. Timestamp records export time; acquisition metadata unspecified.",
    sample_comments = sample_name, instrument_number = 1L, units = "Seconds",
    number_of_scans = as.integer(n), starting_scan_number = 0L,
    raw_data_mass_format = "Double", raw_data_time_format = "Double",
    raw_data_intensity_format = if (integer_y) "Long" else "Double",
    raw_data_uniform_sampling_flag = as.integer(max(abs(diff(seconds) - stats::median(diff(seconds)))) < 1e-6),
    actual_delay_time = min(seconds), actual_run_time = max(seconds),
    actual_run_time_length = max(seconds), global_mass_min = min(mass), global_mass_max = max(mass),
    mass_axis_label = "M/Z", intensity_axis_label = "Abundance")
  for (name in c("test_separation_type", "test_ms_inlet", "test_ionization_mode", "test_ionization_polarity",
      "test_detector_type", "test_scan_function", "test_scan_direction", "test_scan_law")) metadata[[name]] <- ""
  for (name in names(metadata)) ncdf4::ncatt_put(nc, 0, name, metadata[[name]])
  for (name in c("mass_values", "intensity_values")) ncdf4::ncatt_put(nc, name, "scale_factor", 1)
  values <- list(scan_index = as.integer((seq_len(n) - 1L) * m), point_count = rep.int(as.integer(m), n),
    flag_count = integer(n), scan_type = integer(n), a_d_sampling_rate = 1 / dt,
    scan_acquisition_time = seconds, scan_duration = dt, mass_range_min = rep(min(mass), n),
    mass_range_max = rep(max(mass), n), resolution = rep(1, n), total_intensity = rowSums(y),
    mass_values = rep(mass, n), intensity_values = as.vector(t(y)), error_log = "",
    instrument_name = "L0glm", instrument_id = "Matrix export", instrument_mfr = "",
    instrument_model = "", instrument_sw_version = "", instrument_os_version = "")
  for (name in names(values)) {
    if (name %in% c("mass_values", "intensity_values"))
      ncdf4::ncvar_put(nc, name, values[[name]], start = 1L, count = points) else
      ncdf4::ncvar_put(nc, name, values[[name]])
  }
  ncdf4::nc_sync(nc)
  invisible(normalizePath(cdf_file, winslash = "/"))
}
