## ---------------------------------------------------------------------
##  Runtime dependency check.
##
##  The CRAN packages in DESCRIPTION's Imports are installed together
##  with GECO, so they never need checking. What CAN be missing is the
##  optional machinery the pipeline needs at run time:
##    * fst                               -- on-disk cell data (Suggests)
##    * an .h5ad reader: anndata (CRAN + the Python `anndata` module)
##                       or zellkonverter + SingleCellExperiment +
##                       SummarizedExperiment (Bioconductor)
##  check_dependencies() reports what is missing and, with the user's
##  consent, installs the R side. It is called by geco_prepare() and
##  generate_slurm_pipeline(); worker steps only verify, never install.
## ---------------------------------------------------------------------

.pkg_available <- function(p) requireNamespace(p, quietly = TRUE)

## What the pipeline needs for a given backend choice. Returns a list
## of "requirement groups"; each group is satisfied if all its R
## packages are installed. For backend = "auto", either h5ad group
## satisfies the h5ad requirement.
.geco_requirements <- function(backend = c("auto", "anndata", "zellkonverter")) {
  backend <- match.arg(backend)
  fst <- list(name = "fst", pkgs = "fst", source = "CRAN",
              why  = "column-selectable on-disk cell data (cell_data.fst)")
  ad  <- list(name = "anndata", pkgs = "anndata", source = "CRAN",
              why  = ".h5ad reader (uses the Python 'anndata' module via reticulate)")
  zk  <- list(name = "zellkonverter",
              pkgs = c("zellkonverter", "SingleCellExperiment", "SummarizedExperiment"),
              source = "Bioconductor",
              why  = ".h5ad reader (Bioconductor; bundles its own Python via basilisk)")
  list(fst = fst, h5ad = switch(backend,
                                auto          = list(ad, zk),
                                anndata       = list(ad),
                                zellkonverter = list(zk)))
}

.group_missing <- function(g) g$pkgs[!vapply(g$pkgs, .pkg_available, logical(1))]

## Install a set of packages from CRAN / Bioconductor into the first
## writable library. Stops with the manual command if it fails.
.install_pkgs <- function(pkgs, source) {
  repos <- getOption("repos")
  if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@") ||
      !nzchar(repos[1])) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }
  if (source == "CRAN") {
    utils::install.packages(pkgs, repos = repos)
  } else {
    if (!.pkg_available("BiocManager")) {
      utils::install.packages("BiocManager", repos = repos)
    }
    BiocManager::install(pkgs, ask = FALSE, update = FALSE)
  }
  still <- pkgs[!vapply(pkgs, .pkg_available, logical(1))]
  if (length(still) > 0L) {
    stop("Installation of ", paste(still, collapse = ", "), " did not ",
         "succeed. Install manually with:\n  ",
         .install_cmd(still, source), call. = FALSE)
  }
  invisible(TRUE)
}

.install_cmd <- function(pkgs, source) {
  q <- paste0('"', pkgs, '"', collapse = ", ")
  if (source == "CRAN") sprintf("install.packages(c(%s))", q)
  else sprintf('install.packages("BiocManager"); BiocManager::install(c(%s))', q)
}

#' Check (and optionally install) what the GECO pipeline needs at run time
#'
#' GECO's hard dependencies are installed together with the package.
#' Two optional pieces are needed to actually run the pipeline and are
#' checked here: the `fst` package, and an `.h5ad` reader (`anndata`
#' from CRAN, or `zellkonverter` + `SingleCellExperiment` +
#' `SummarizedExperiment` from Bioconductor).
#'
#' If something is missing, the behaviour depends on `install`:
#' - `"ask"` (default): in an interactive session, list the missing
#'   packages and ask for confirmation before installing them. In a
#'   non-interactive session (e.g. a SLURM job) nothing is installed;
#'   the function stops and prints the exact commands to run.
#' - `"always"`: install without asking (also non-interactively).
#' - `"never"`: only report; stop if anything is missing.
#'
#' The Python side of the `anndata` backend (`pip install anndata`)
#' is checked when `reticulate` is available but is never installed
#' automatically; the message tells you what to run.
#'
#' This is called automatically by [geco_prepare()] and
#' [generate_slurm_pipeline()]; call it yourself to prepare an
#' environment ahead of time.
#'
#' @param backend Which `.h5ad` reader the run will use (see
#'   [load_h5ad()]). With `"auto"`, having either reader is enough;
#'   if neither is present, `anndata` is the one offered for install.
#' @param install One of `"ask"`, `"always"`, `"never"` (see above).
#'
#' @return Invisibly `TRUE` if everything is available (possibly after
#'   installing). Stops otherwise.
#' @export
check_dependencies <- function(backend = c("auto", "anndata", "zellkonverter"),
                               install = c("ask", "always", "never")) {
  backend <- match.arg(backend)
  install <- match.arg(install)
  req <- .geco_requirements(backend)

  ## ---- what is missing --------------------------------------------
  todo <- list()                                 # groups to install
  if (length(.group_missing(req$fst)) > 0L) todo <- c(todo, list(req$fst))

  h5ad_ok <- any(vapply(req$h5ad, function(g) length(.group_missing(g)) == 0L,
                        logical(1)))
  if (!h5ad_ok) todo <- c(todo, list(req$h5ad[[1L]]))   # anndata first for "auto"

  if (length(todo) > 0L) {
    lines <- vapply(todo, function(g) {
      miss <- .group_missing(g)
      sprintf("  * %s (%s): %s", paste(miss, collapse = ", "), g$source, g$why)
    }, character(1))
    cmds <- vapply(todo, function(g) .install_cmd(.group_missing(g), g$source),
                   character(1))
    header <- paste0("GECO needs the following R package(s) that are not ",
                     "installed:\n", paste(lines, collapse = "\n"))
    manual <- paste0("Install them with:\n  ", paste(cmds, collapse = "\n  "),
                     "\nor run GECO::check_dependencies(install = \"always\").")

    do_install <- switch(install,
      never  = FALSE,
      always = TRUE,
      ask    = {
        if (!interactive()) {
          stop(header, "\nThis is a non-interactive session, so nothing was ",
               "installed. ", manual, call. = FALSE)
        }
        message(header)
        ans <- utils::askYesNo("Install them now into your R library?",
                               default = FALSE)
        isTRUE(ans)
      })
    if (!do_install) {
      stop(header, "\n", manual, call. = FALSE)
    }
    installed <- character(0)
    for (g in todo) {
      miss <- .group_missing(g)
      message("[GECO] Installing ", paste(miss, collapse = ", "),
              " from ", g$source, "...")
      .install_pkgs(miss, g$source)
      installed <- c(installed, miss)
    }
    message("[GECO] Installed: ", paste(installed, collapse = ", "), ".")
  }

  ## ---- Python side of the anndata backend --------------------------
  uses_anndata <- backend == "anndata" ||
    (backend == "auto" && .pkg_available("anndata"))
  if (uses_anndata && .pkg_available("reticulate")) {
    py_ok <- tryCatch(reticulate::py_module_available("anndata"),
                      error = function(e) NA)
    if (identical(py_ok, FALSE)) {
      msg <- paste0("The R package 'anndata' is installed but the Python ",
                    "module 'anndata' is not available to reticulate. In the ",
                    "Python environment R will use, run:\n  pip install anndata\n",
                    "and, if needed, point reticulate at it with\n  ",
                    "Sys.setenv(RETICULATE_PYTHON = \"/path/to/python\")")
      if (backend == "anndata") stop(msg, call. = FALSE)
      ## auto: fall back to zellkonverter if that is complete, else stop
      zk_ok <- length(.group_missing(.geco_requirements("zellkonverter")$h5ad[[1L]])) == 0L
      if (!zk_ok) stop(msg, call. = FALSE)
      message("[GECO] ", msg, "\n[GECO] Falling back to the zellkonverter backend.")
    }
  }
  invisible(TRUE)
}

## Internal: verify (never install) that a package is present. Used by
## the worker-side steps, which run non-interactively on compute nodes.
.require_pkg <- function(pkg, where) {
  if (!.pkg_available(pkg)) {
    stop(where, "(): package '", pkg, "' is not installed in this R ",
         "library (", .libPaths()[1L], "). Install it with ",
         .install_cmd(pkg, "CRAN"), " or run GECO::check_dependencies() ",
         "on this machine.", call. = FALSE)
  }
  invisible(TRUE)
}
