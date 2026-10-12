## ---------------------------------------------------------------------
##  Runtime dependency check.
##
##  The CRAN packages in DESCRIPTION's Imports are installed together
##  with GECO, so they never need checking. Two optional packages are
##  needed at run time and are checked here:
##    * fst    (CRAN)         -- on-disk cell data (cell_data.fst)
##    * rhdf5  (Bioconductor) -- reading the .h5ad file
##  check_dependencies() reports what is missing and, with the user's
##  consent, installs it. It is called by geco_prepare() and
##  generate_slurm_pipeline(); worker steps only verify, never install.
## ---------------------------------------------------------------------

.pkg_available <- function(p) requireNamespace(p, quietly = TRUE)

.geco_requirements <- function() {
  list(
    list(pkgs = "fst",   source = "CRAN",
         why  = "column-selectable on-disk cell data (cell_data.fst)"),
    list(pkgs = "rhdf5", source = "Bioconductor",
         why  = "reading the .h5ad file (HDF5)")
  )
}

.group_missing <- function(g) g$pkgs[!vapply(g$pkgs, .pkg_available, logical(1))]

.install_cmd <- function(pkgs, source) {
  q <- paste0('"', pkgs, '"', collapse = ", ")
  if (source == "CRAN") sprintf("install.packages(c(%s))", q)
  else sprintf('install.packages("BiocManager"); BiocManager::install(c(%s))', q)
}

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

#' Check (and optionally install) what the GECO pipeline needs at run time
#'
#' GECO's hard dependencies are installed together with the package.
#' Two optional packages are needed to actually run the pipeline and
#' are checked here: `fst` (CRAN) and `rhdf5` (Bioconductor).
#'
#' If something is missing, the behaviour depends on `install`:
#' - `"ask"` (default): in an interactive session, list the missing
#'   packages and ask for confirmation before installing them. In a
#'   non-interactive session (e.g. a SLURM job) nothing is installed;
#'   the function stops and prints the exact commands to run.
#' - `"always"`: install without asking (also non-interactively).
#' - `"never"`: only report; stop if anything is missing.
#'
#' This is called automatically by [geco_prepare()] and
#' [generate_slurm_pipeline()]; call it yourself to prepare an
#' environment ahead of time.
#'
#' @param install One of `"ask"`, `"always"`, `"never"` (see above).
#'
#' @return Invisibly `TRUE` if everything is available (possibly after
#'   installing). Stops otherwise.
#' @export
check_dependencies <- function(install = c("ask", "always", "never")) {
  install <- match.arg(install)
  todo <- Filter(function(g) length(.group_missing(g)) > 0L, .geco_requirements())
  if (length(todo) == 0L) return(invisible(TRUE))

  lines <- vapply(todo, function(g) {
    sprintf("  * %s (%s): %s", paste(.group_missing(g), collapse = ", "),
            g$source, g$why)
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
      isTRUE(utils::askYesNo("Install them now into your R library?",
                             default = FALSE))
    })
  if (!do_install) stop(header, "\n", manual, call. = FALSE)

  installed <- character(0)
  for (g in todo) {
    miss <- .group_missing(g)
    message("[GECO] Installing ", paste(miss, collapse = ", "),
            " from ", g$source, "...")
    .install_pkgs(miss, g$source)
    installed <- c(installed, miss)
  }
  message("[GECO] Installed: ", paste(installed, collapse = ", "), ".")
  invisible(TRUE)
}

## Internal: verify (never install) that a package is present. Used by
## the worker-side steps, which run non-interactively on compute nodes.
.require_pkg <- function(pkg, where, source = "CRAN") {
  if (!.pkg_available(pkg)) {
    stop(where, "(): package '", pkg, "' is not installed in this R ",
         "library (", .libPaths()[1L], "). Install it with ",
         .install_cmd(pkg, source), " or run GECO::check_dependencies() ",
         "on this machine.", call. = FALSE)
  }
  invisible(TRUE)
}
