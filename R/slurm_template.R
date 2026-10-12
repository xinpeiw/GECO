## Internal: is `p` an absolute path? Handles Unix (starts with "/")
## and Windows (`C:/`, `C:\\`, `//server/share`) forms.
.is_abs_path <- function(p) {
  grepl("^(/|~|[A-Za-z]:[/\\\\]|//)", p)
}


#' Generate a SLURM pipeline scaffold for cluster mode
#'
#' Emits a set of R + shell scripts into `scripts_dir` that run the
#' five cluster-mode functions in order with SLURM array job
#' dependencies. The user only needs to edit the SLURM resource
#' options at the top of each `.sh` and then run
#' `bash pipeline.sh`.
#'
#' The number of cell types is taken straight from the .h5ad file —
#' the precision SLURM array is sized at runtime by step 02 (after
#' [geco_finalize_estimates()] has written `cell_type_order.rds`), so you never
#' have to tell the scaffolding how many cell types your data has.
#'
#' **A note on paths.** All five path arguments below (`scripts_dir`,
#' `work_dir`, `output_dir`, `h5ad_file`, `gene_list_file`) are
#' resolved to **absolute paths** at template-generation time via
#' `normalizePath(mustWork = FALSE)`. That means it doesn't matter
#' whether you pass them as absolute or relative to your R working
#' directory — the generated `.sh` files always contain the fully
#' resolved paths, so `sbatch` / `bash pipeline.sh` will work no
#' matter which directory you launch them from.
#'
#' @param scripts_dir Directory to write the scripts into (created if
#'   missing).
#' @param work_dir Working directory the scripts should use (becomes a
#'   variable in each script). Does not need to exist yet.
#' @param output_dir Final outputs directory used by step 04. Does
#'   not need to exist yet.
#' @param h5ad_file Path to the .h5ad file used by `geco_prepare()`.
#'   Must exist.
#' @param gene_list_file Path to a text file with one gene ID per line
#'   (loaded by [read_gene_list()]). Passed to `geco_prepare()` as the
#'   `genes` argument. Required. Must exist.
#' @param n_jobs Number of SLURM array elements for the chunk step.
#'   Choose this so each chunk runs in a reasonable wall time (e.g.
#'   300 is a reasonable default for ~10k genes / ~50M pairs).
#' @param prepare_args Named list of additional args injected into the
#'   `geco_prepare()` call (e.g. `subset_col = "SLE_status"`).
#' @param slurm Named list of SLURM options. Any subset of
#'   `partition`, `time`, `mem`, `cpus` may be given; entries you omit
#'   keep their defaults (`partition = "general"`, `time = "24:00:00"`,
#'   `mem = "32G"`, `cpus = 1`), so `slurm = list(mem = "64G")` changes
#'   only the memory. Unknown names are an error.
#' @param r_executable Command used to run R, e.g. `"Rscript"` or
#'   `"micromamba run -n r Rscript"`.
#'
#' @return Invisibly returns `scripts_dir`.
#' @export
generate_slurm_pipeline <- function(scripts_dir,
                                    work_dir,
                                    output_dir,
                                    h5ad_file,
                                    gene_list_file,
                                    n_jobs = 300,
                                    prepare_args = list(),
                                    slurm = list(),
                                    r_executable = "Rscript") {

  ## Merge user SLURM options over the defaults, so a partial list such
  ## as list(mem = "64G") does not leave the other #SBATCH lines empty.
  slurm_defaults <- list(partition = "general",
                         time      = "24:00:00",
                         mem       = "32G",
                         cpus      = 1)
  if (!is.list(slurm) ||
      (length(slurm) > 0L && (is.null(names(slurm)) || any(!nzchar(names(slurm)))))) {
    stop("`slurm` must be a named list, e.g. list(mem = \"64G\", time = \"12:00:00\").",
         call. = FALSE)
  }
  unknown <- setdiff(names(slurm), names(slurm_defaults))
  if (length(unknown) > 0L) {
    stop("Unknown `slurm` option(s): ", paste(unknown, collapse = ", "),
         ". Allowed: ", paste(names(slurm_defaults), collapse = ", "), ".",
         call. = FALSE)
  }
  slurm <- utils::modifyList(slurm_defaults, slurm)
  empty <- vapply(slurm, function(v) is.null(v) || length(v) == 0L ||
                                      is.na(v) || !nzchar(as.character(v)),
                  logical(1))
  if (any(empty)) {
    stop("`slurm` option(s) must not be empty/NA: ",
         paste(names(slurm)[empty], collapse = ", "), ".", call. = FALSE)
  }

  ## Resolve every path argument to an absolute path so the generated
  ## .sh files work no matter where the user launches sbatch from.
  ## normalizePath() with mustWork = FALSE may leave a non-existent
  ## relative path unchanged on some platforms, so we manually anchor
  ## relative paths at the R session's working directory.
  as_abs <- function(p) {
    p <- path.expand(p)
    if (!.is_abs_path(p)) p <- file.path(getwd(), p)
    normalizePath(p, mustWork = FALSE, winslash = "/")
  }
  scripts_dir    <- as_abs(scripts_dir)
  work_dir       <- as_abs(work_dir)
  output_dir     <- as_abs(output_dir)
  h5ad_file      <- as_abs(h5ad_file)
  gene_list_file <- as_abs(gene_list_file)

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
    "gene_file <- args[3]",
    "genes <- GECO::read_gene_list(gene_file)",
    "GECO::geco_prepare(",
    "  h5ad_file = h5ad_file,",
    "  work_dir  = work_dir,",
    paste0("  genes     = genes", prepare_args_block),
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
    "GECO::geco_finalize_estimates(work_dir = args[1],",
    "                              output_dir = args[2])"
  ), file.path(scripts_dir, "02_finalize_estimates.R"))

  writeLines(c(
    "## 03_precision.R - SLURM array element (one per cell type + shared)",
    "args <- commandArgs(trailingOnly = TRUE)",
    "library(GECO)",
    "work_dir   <- args[1]",
    "output_dir <- args[2]",
    "ct_idx <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID',",
    "                                unset = if (length(args) >= 3) args[3] else NA))",
    "if (is.na(ct_idx)) stop('No SLURM_ARRAY_TASK_ID set and no 3rd CLI arg.')",
    "GECO::geco_precision_chunk(work_dir   = work_dir,",
    "                           output_dir = output_dir,",
    "                           ct_idx     = ct_idx)"
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
    paste0("  ", shQuote(work_dir),  " \\"),
    paste0("  ", shQuote(gene_list_file))
  ), file.path(scripts_dir, "00_prepare.sh"))

  writeLines(c(
    hdr("geco_chunk", array = paste0("1-", n_jobs)),
    "set -euo pipefail",
    paste0(r_executable, " 01_chunk.R \\"),
    paste0("  ", shQuote(work_dir), " \\"),
    paste0("  ", n_jobs)
  ), file.path(scripts_dir, "01_chunk.sh"))

  ## Step 02 finalizes estimates AND then submits steps 03 + 04 with
  ## the correct precision array size (derived from cell_type_order.rds).
  writeLines(c(
    hdr("geco_fin_est"),
    "set -euo pipefail",
    "",
    paste0("WORK_DIR=",   shQuote(work_dir)),
    paste0("OUTPUT_DIR=", shQuote(output_dir)),
    "",
    paste0(r_executable,
           " 02_finalize_estimates.R \"${WORK_DIR}\" \"${OUTPUT_DIR}\""),
    "",
    "# How many precision jobs to launch? = (n_cell_types + 1).",
    "# Read it from cell_type_order.rds, which step 02 just wrote.",
    paste0("N_PREC=$(", r_executable,
           " -e \"cat(length(readRDS(file.path(commandArgs(TRUE)[1],",
           " 'cell_type_order.rds'))))\" --args \"${WORK_DIR}\")"),
    "echo \"n_precision_jobs = ${N_PREC}\"",
    "",
    "J3=$(sbatch --parsable --array=1-${N_PREC} 03_precision.sh)",
    "echo \"precision: ${J3}\"",
    "",
    "sbatch --parsable --dependency=afterok:${J3} 04_finalize_outputs.sh"
  ), file.path(scripts_dir, "02_finalize_estimates.sh"))

  ## 03_precision.sh has no #SBATCH --array directive — that is
  ## supplied at submit time by 02_finalize_estimates.sh.
  writeLines(c(
    hdr("geco_precision"),   # no array directive here
    "set -euo pipefail",
    paste0(r_executable, " 03_precision.R \\"),
    paste0("  ", shQuote(work_dir),   " \\"),
    paste0("  ", shQuote(output_dir))
  ), file.path(scripts_dir, "03_precision.sh"))

  writeLines(c(
    hdr("geco_fin_out"),
    "set -euo pipefail",
    paste0(r_executable, " 04_finalize_outputs.R \\"),
    paste0("  ", shQuote(work_dir),   " \\"),
    paste0("  ", shQuote(output_dir))
  ), file.path(scripts_dir, "04_finalize_outputs.sh"))

  ## ---------- Master orchestrator ----------
  writeLines(c(
    "#!/bin/bash",
    "# pipeline.sh - submits the GECO cluster pipeline via sbatch dependencies.",
    "# Steps 03 (precision array) and 04 (finalize outputs) are submitted",
    "# from inside step 02, after cell_type_order.rds tells us how many cell types we have.",
    "set -euo pipefail",
    "",
    "J0=$(sbatch --parsable 00_prepare.sh)",
    "echo \"prepare: $J0\"",
    "",
    "J1=$(sbatch --parsable --dependency=afterok:${J0} 01_chunk.sh)",
    "echo \"chunk:   $J1\"",
    "",
    "J2=$(sbatch --parsable --dependency=afterok:${J1} 02_finalize_estimates.sh)",
    "echo \"fin_est: $J2   (will submit 03 + 04 itself)\""
  ), file.path(scripts_dir, "pipeline.sh"))

  Sys.chmod(list.files(scripts_dir, "\\.sh$", full.names = TRUE), "0755")

  message("[generate_slurm_pipeline] Scripts written to ", scripts_dir)
  message("  Edit SLURM options at the top of each .sh, then: cd ", scripts_dir,
          " && bash pipeline.sh")
  invisible(scripts_dir)
}
