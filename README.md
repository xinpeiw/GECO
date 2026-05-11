# GECO

**GECO** (**Ge**ne **Co**variance) is an R package that decomposes
gene-by-gene covariance single-cell expression data
(`.h5ad`) into cell-type-shared and cell-type-specific components, and
turns the result into shared and per-cell-type sparse precision matrices and node
degrees.

The entire legacy SLURM pipeline (data process → VS/Y generation →
estimation → reformatting → p-values → PSD transform → glasso →
eBIC selection → degree) is exposed in two ways:

- **In-memory mode** — one call to `run_geco()` runs everything on a
  single machine. Good for up to a few hundred genes.
- **Cluster mode** — a set of `geco_prepare()` / `geco_chunk()` /
  `geco_finalize_estimates()` / `geco_precision_chunk()` /
  `geco_finalize_outputs()` functions, plus a
  `generate_slurm_pipeline()` helper that emits ready-to-`sbatch`
  scripts. Use this for ~10k-gene runs that don't fit on one node.

## Inputs and outputs

**Input.** A `.h5ad` file. Cell-level metadata is read from
`obs`/`colData`; the expression matrix is read from the chosen assay.
Column names for the individual ID and cell-type are configurable.

**Outputs.** For every cell type plus the "shared" (across-individual)
component:

- `cov_<celltype>.rds`   — covariance matrix (`n_gene x n_gene`)
- `se_<celltype>.rds`    — standard error matrix
- `p_<celltype>.rds`     — two-sided p-value matrix
- `cov_psd_<celltype>.rds` — PSD-clipped covariance
- `precision_<celltype>.rds` — sparse precision matrix (glassoFast + eBIC)

Plus a single workbook:

- `degree.xlsx` — node degrees, with one column per cell type.

## Installation

```r
# Dependencies (CRAN)
install.packages(c("Matrix", "data.table", "glassoFast", "openxlsx"))

# For reading .h5ad: either zellkonverter (recommended, Bioconductor)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment",
                       "SummarizedExperiment"))
# ...or anndata (CRAN, uses reticulate)
# install.packages("anndata")

# Install GECO from local source
install.packages("path/to/GECO", repos = NULL, type = "source")
# or, if you've cloned the package and have devtools:
# devtools::install("path/to/GECO")
```

## Quick start (in-memory mode)

```r
library(GECO)

result <- run_geco(
  h5ad_file        = "subset_data.h5ad",
  ind_col          = "ind_cov",
  ct_col           = "cg_cov",
  # optional: restrict to a subset of cells
  # subset_col     = "SLE_status", subset_value = "Healthy",
  # pick which covariates to regress out (see available_covariates()):
  fixed_covariates = c("age_group", "sex", "PC1"),
  output_dir       = "./results"
)
# Cell-type names in the output are taken directly from
# unique(data[[ct_col]]) in the h5ad file.

# Everything also returned in-memory:
str(result$cov)        # named list of covariance matrices
str(result$precision)  # named list of sparse precision matrices
head(result$degree)    # gene x cell_type degree table
```

After `run_geco()` finishes, the `output_dir` will contain all the
`.rds` matrices and `degree.xlsx`.

### Choosing which fixed covariates to regress out

`run_geco()` and `geco_prepare()` both take a `fixed_covariates`
argument that decides which per-individual covariates are regressed
out before the GeCo decomposition runs.

The built-in catalog (see `GECO_COVARIATE_CATALOG` or call
`available_covariates()`) provides these names:

| Name | Source column | What it becomes |
|------|---------------|-----------------|
| `age_group` | `Age`        | binned to 1..15 |
| `sex`       | `Sex`        | 1 if Male, else 0 |
| `pop`       | `pop_cov`    | 1 if not European, else 0 |
| `cases`     | `SLE_status` | 1 if not Healthy, else 0 |
| `PC1`       | (derived)    | PC1 of per-individual mean expression |

Pick whichever subset you want — the most common cases:

```r
# No regression at all (default):
run_geco(..., fixed_covariates = NULL)

# Just age and sex:
run_geco(..., fixed_covariates = c("age_group", "sex"))

# The full legacy-script set:
run_geco(..., fixed_covariates = c("age_group", "sex", "pop", "cases", "PC1"))

# Custom covariate from a "Smoker" column, alongside age:
run_geco(..., fixed_covariates = list(
  age_group = GECO_COVARIATE_CATALOG$age_group,
  smoker    = list(col = "Smoker", type = "binary", positive = "Yes")
))

# Or compute your own per-individual matrix and pass it directly:
my_mat <- ...  # n_ind x p matrix, rownames = individual IDs
run_geco(..., fixed_covariates = my_mat)
```

## Step-by-step usage

Every step of the pipeline is also exposed in case you want to run only
part of it or inject custom logic:

```r
data <- load_h5ad("subset_data.h5ad")

genes <- setdiff(colnames(data),
                 c("ind_cov", "cg_cov", "SLE_status",
                   "Sex", "Age", "pop_cov"))

opms  <- compute_OPms(data, genes)
pc1   <- compute_PC1 (data, genes)
fcov  <- build_fixed_cov(data,
                         covariates = c("age_group", "sex", "PC1"),
                         pc1 = pc1)

est   <- estimate_geco(data, genes, opms, fixed_cov = fcov)
mats  <- reformat_to_matrices(est, genes)

p_mat   <- compute_p_matrix(mats$cov$T4, mats$se$T4)
cov_psd <- transform_psd   (mats$cov$T4)

fit <- fit_precision(cov_psd, p_mat, n_ind = 243)
deg <- compute_degree(fit$precision)
```

## Cluster (HPC) mode

For ~10k-gene runs that don't fit on a single node, use the
chunked pipeline. The cluster mode requires `fst`
(`install.packages("fst")`) for column-selectable on-disk cell data.

The simplest path is to let GECO generate ready-to-`sbatch` scripts
for you:

```r
library(GECO)

generate_slurm_pipeline(
  scripts_dir  = "/path/to/scripts",
  work_dir     = "/path/to/work",
  output_dir   = "/path/to/results",
  h5ad_file    = "/path/to/subset_data.h5ad",
  n_jobs       = 241,   # number of SLURM array chunks
  n_cell_types = 7,     # so the precision SLURM array is sized as 1..(n+1)
  prepare_args = list(fixed_covariates = c("age_group", "sex", "PC1")),
  slurm = list(partition = "general", time = "24:00:00",
               mem = "32G", cpus = 4)
)
```

This writes five `.R` + `.sh` files plus a `pipeline.sh` orchestrator
into `scripts_dir`. Adjust the SLURM headers if needed, then:

```bash
cd /path/to/scripts
bash pipeline.sh
```

`pipeline.sh` chains the five SLURM jobs with `--dependency=afterok`:

```
00_prepare.sh   →  one-time setup (h5ad → fst + OPms + fixed_cov + gene_pairs)
01_chunk.sh     →  SLURM array, runs geco_chunk per task_id
02_finalize_estimates.sh
                →  merge chunks, reformat into per-cell-type matrices,
                   compute p-values, PSD-transform
03_precision.sh →  SLURM array (one per "shared" + cell type), glasso+eBIC
04_finalize_outputs.sh
                →  build degree.xlsx and write all matrices to output_dir
```

You can also call the functions directly without the script generator:

```r
# Once
geco_prepare(h5ad_file = "...", work_dir = "work/",
             fixed_covariates = c("age_group", "sex", "PC1"))

# Per chunk (SLURM array, 1..n_jobs)
geco_chunk(work_dir = "work/", task_id = task_id, n_jobs = 241)

# After all chunks finish
geco_finalize_estimates(work_dir = "work/")

# Per cell type (SLURM array, 1..(n_ct + 1) for shared)
geco_precision_chunk(work_dir = "work/", ct_idx = i)

# Finally
geco_finalize_outputs(work_dir = "work/", output_dir = "results/")
```

## Notes and caveats

- In-memory mode keeps the cell-level expression matrix and all
  per-pair tensors in RAM. For ~10k genes this is too large; switch
  to cluster mode.
- Cluster mode stores the cell-level data once as `cell_data.fst`
  and each chunk only reads the two gene columns it needs, so each
  worker's memory footprint stays small.
- `fixed_covariates = c("age_group", "sex", "pop", "cases", "PC1")`
  reproduces the legacy script. See **Choosing which fixed covariates
  to regress out** above for the full menu (and how to add custom
  ones). `fixed_covariates = NULL` skips the regression entirely.
- Cell-type names in the output are taken directly from
  `unique(data[[ct_col]])` in the h5ad file — they are not provided
  by the user. The first matrix in every output list is always named
  `shared` (the across-individual component).
