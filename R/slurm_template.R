#' Generate a SLURM pipeline scaffold for cluster mode
#'
#' Emits a set of R + shell scripts into `scripts_dir` that run the
#' five cluster-mode functions in order with SLURM array job
#' dependencies. The user only needs to edit the SLURM resource
#' options at the top of each `.sh` and then run
#' `bash pipeline.sh`.
#'
#' @param scripts_dir Directory to write the scripts into (created if
#'   missing).
#' @param work_dir Working directory the scripts should use (becomes a
#'   variable in each script).
#' @param output_dir Final outputs directory used by step 04.
#' @param h5ad_file Path to the .h5ad file used by `geco_prepare()`.
#' @param n_jobs Number of SLURM array elements for the chunk step.
#'   Choose this so each chunk runs in a reasonable wall time (e.g.
#'   241 for ~6M pairs, mirroring the legacy script).
#' @param n_cell_types Number of cell types in the data. Used only to
#'   set the SLURM array size for the precision step
#'   (`1..(n_cell_types + 1)`, for the "shared" matrix plus each cell
#'   type). If `NULL`, the script leaves a placeholder array bound
#'   that the user must edit after running `geco_prepare()` (since the
#'   actual cell types are taken from the h5ad file).
#' @param prepare_args Named list of additional args injected into the
#'   `geco_prepare()` call (e.g. `subset_col = "SLE_status"`).
#' @param slurm Named list of SLURM options. Defaults are reasonable
#'   starting points; tune for your cluster.
#' @param r_executable Command used to run R, e.g. `"Rscript"` or
#'   `"micromamba run -n r Rscript"`.
#'
#' @return Invisibly returns `scripts_dir`.
#' @export
generate_slurm_pipeline <- function(scripts_dir,
                                    work_dir,
                                    output_dir,
                                    h5ad_file,
                                    n_jobs = 241,
                                    n_cell_types = NULL,
                                    prepare_args = list(),
                                    slurm = list(
                                      partition = "general",
                                      time      = "24:00:00",
                                      mem       = "32G",
                                      cpus      = 4
                                    ),
                                    r_executable = "Rscript") {

  dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)

  prep_args_str <- paste(
    sprintf('  %s = %s',
            names(prepare_args),
            vapply(prepare_args, function(x) deparse1(x), character(1))),
    collapse = ",\n"
  )
  prepare_args_block <- if (nzchar(prep_args_str)) paste0(",\n", prep_args_str) else ""

  ## ---------- R scripts ----------
  writeLines(c(
    "## 00_prepare.R - one-time setup",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "h5ad_file <- args[1]",
    "work_dir  <- args[2]",
    "GECO::geco_prepare(",
    "  h5ad_file = h5ad_file,",
    paste0("  work_dir  = work_dir", prepare_args_block),
    ")"
  ), file.path(scripts_dir, "00_prepare.R"))

  writeLines(c(
    "## 01_chunk.R - SLURM array element",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "work_dir <- args[1]",
    "n_jobs   <- as.integer(args[2])",
    "task_id  <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID',",
    "                                 unset = if (length(args) >= 3) args[3] else NA))",
    "if (is.na(task_id)) stop('No SLURM_ARRAY_TASK_ID set and no 3rd CLI arg.')",
    "GECO::geco_chunk(work_dir = work_dir,",
    "                 task_id  = task_id,",
    "                 n_jobs   = n_jobs)"
  ), file.path(scripts_dir, "01_chunk.R"))

  writeLines(c(
    "## 02_finalize_estimates.R",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "GECO::geco_finalize_estimates(work_dir = args[1])"
  ), file.path(scripts_dir, "02_finalize_estimates.R"))

  writeLines(c(
    "## 03_precision.R - SLURM array element (one per cell type + shared)",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "work_dir <- args[1]",
    "ct_idx <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID',",
    "                                unset = if (length(args) >= 2) args[2] else NA))",
    "if (is.na(ct_idx)) stop('No SLURM_ARRAY_TASK_ID set and no 2nd CLI arg.')",
    "GECO::geco_precision_chunk(work_dir = work_dir, ct_idx = ct_idx)"
  ), file.path(scripts_dir, "03_precision.R"))

  writeLines(c(
    "## 04_finalize_outputs.R",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "GECO::geco_finalize_outputs(work_dir = args[1], output_dir = args[2])"
  ), file.path(scripts_dir, "04_finalize_outputs.R"))

  ## ---------- Shell wrappers ----------
  hdr <- function(name, array = NULL) {
    out_pattern <- if (is.null(array)) {
      paste0(name, "-%j.out")
    } else {
      paste0(name, "-%A_%a.out")
    }
    lines <- c(
      "#!/bin/bash",
      paste0("#SBATCH --job-name=", name),
      paste0("#SBATCH --partition=", slurm$partition),
      paste0("#SBATCH --time=",      slurm$time),
      paste0("#SBATCH --mem=",       slurm$mem),
      paste0("#SBATCH --cpus-per-task=", slurm$cpus),
      paste0("#SBATCH --output=", out_pattern)
    )
    if (!is.null(array)) {
      lines <- c(lines, paste0("#SBATCH --array=", array))
    }
    lines
  }

  writeLines(c(
    hdr("geco_prepare"),
    "set -euo pipefail",
    paste0(r_executable, " 00_prepare.R \\"),
    paste0("  ", shQuote(h5ad_file), " \\"),
    paste0("  ", shQuote(work_dir))
  ), file.path(scripts_dir, "00_prepare.sh"))

  writeLines(c(
    hdr("geco_chunk", array = paste0("1-", n_jobs)),
    "set -euo pipefail",
    paste0(r_executable, " 01_chunk.R \\"),
    paste0("  ", shQuote(work_dir), " \\"),
    paste0("  ", n_jobs)
  ), file.path(scripts_dir, "01_chunk.sh"))

  writeLines(c(
    hdr("geco_fin_est"),
    "set -euo pipefail",
    paste0(r_executable, " 02_finalize_estimates.R \\"),
    paste0("  ", shQuote(work_dir))
  ), file.path(scripts_dir, "02_finalize_estimates.sh"))

  precision_array <- if (!is.null(n_cell_types)) {
    paste0("1-", as.integer(n_cell_types) + 1L)
  } else {
    "1-N   # edit to (n_cell_types + 1) after geco_prepare() runs"
  }
  writeLines(c(
    hdr("geco_precision", array = precision_array),
    "set -euo pipefail",
    paste0(r_executable, " 03_precision.R \\"),
    paste0("  ", shQuote(work_dir))
  ), file.path(scripts_dir, "03_precision.sh"))

  writeLines(c(
    hdr("geco_fin_out"),
    "set -euo pipefail",
    paste0(r_executable, " 04_finalize_outputs.R \\"),
    paste0("  ", shQuote(work_dir), " \\"),
    paste0("  ", shQuote(output_dir))
  ), file.path(scripts_dir, "04_finalize_outputs.sh"))

  ## ---------- Master orchestrator ----------
  writeLines(c(
    "#!/bin/bash",
    "# pipeline.sh - submits the entire GECO cluster pipeline via sbatch dependencies",
    "set -euo pipefail",
    "",
    "J0=$(sbatch --parsable 00_prepare.sh)",
    "echo \"prepare: $J0\"",
    "",
    "J1=$(sbatch --parsable --dependency=afterok:${J0} 01_chunk.sh)",
    "echo \"chunk:   $J1\"",
    "",
    "J2=$(sbatch --parsable --dependency=afterok:${J1} 02_finalize_estimates.sh)",
    "echo \"fin_est: $J2\"",
    "",
    "J3=$(sbatch --parsable --dependency=afterok:${J2} 03_precision.sh)",
    "echo \"precision: $J3\"",
    "",
    "J4=$(sbatch --parsable --dependency=afterok:${J3} 04_finalize_outputs.sh)",
    "echo \"fin_out: $J4\""
  ), file.path(scripts_dir, "pipeline.sh"))

  Sys.chmod(list.files(scripts_dir, "\\.sh$", full.names = TRUE), "0755")

  message("[generate_slurm_pipeline] Scripts written to ", scripts_dir)
  message("  Edit SLURM options at the top of each .sh, then: cd ", scripts_dir,
          " && bash pipeline.sh")
  invisible(scripts_dir)
}
