## Declare this package as data.table-aware.
##
## `data.table`'s `[.data.table` method sets up special symbols such
## as `.SD`, `.SDcols`, `.N`, `.I`, and `.GRP` inside the `j`
## expression, but ONLY when it's called from an R session or package
## that has been declared "data.table-aware". Without this flag,
## calls like `dt[, lapply(.SD, mean), by = k, .SDcols = ...]` made
## from within a package function silently fall through to
## `[.data.frame`, which doesn't recognize `.SD` and errors with
## `object '.SD' not found`. See ?data.table::`datatable-package`
## and the "Importing data.table" vignette.
##
## Setting this once in any of the package's R files is enough.
.datatable.aware <- TRUE

## Silence R CMD check "no visible binding for global variable" for
## the special data.table symbols we reference inside [.data.table.
utils::globalVariables(c(".SD", ".N", ".I", ".GRP", ".SDcols"))
