# GECO

**GECO** (**Ge**ne **Co**expression) is an R package that unbiasedly decomposes
interindividual gene-by-gene coexpression into cell-type-shared and cell-type-specific components, and turns the result into cell-type-shared and per-cell-type-specific sparse precision matrices and node
degrees.

The pipeline (data process → Y/VS generation → covariance/se
estimation → p-values → glasso → degree) runs as a SLURM job chain:
`geco_prepare()` / `geco_chunk()` / `geco_finalize_estimates()` /
`geco_precision_chunk()` / `geco_finalize_outputs()`, plus a
`generate_slurm_pipeline()` helper that emits ready-to-`sbatch`
scripts. A genome-wide run (~10k genes, ~50M gene pairs) does not fit
on a single node, so this is the only supported way to run it.

## Contents

- [Inputs and outputs](#inputs-and-outputs)
- [Installation](#installation)
- [Running on a cluster](#running-on-a-cluster)
  - [Running it](#running-it)
  - [Choosing which covariates to adjust](#choosing-which-covariates-to-adjust)
  - [What ends up where on disk](#what-ends-up-where-on-disk)
  - [Wall-time breakdown](#wall-time-breakdown)
  - [Memory & disk budget](#memory--disk-budget)
  - [Resuming after a partial failure](#resuming-after-a-partial-failure)
  - [Running a single SLURM chunk manually](#running-a-single-slurm-chunk-manually)
  - [Common pitfalls](#common-pitfalls)
- [Input requirements and the "error instead of NA" policy](#input-requirements-and-the-error-instead-of-na-policy)
- [Notes and caveats](#notes-and-caveats)

## Inputs and outputs

**Input.** A `.h5ad` file. Cell-level metadata is read from
`obs`/`colData`; the expression matrix is read from the chosen assay.
Column names for the individual ID and cell-type are configurable.

**Outputs.** For every cell-type-specific plus the cell-type-shared component:

- `cov_<celltype>.rds`   — covariance matrix (`n_gene x n_gene`)
- `se_<celltype>.rds`    — standard error matrix
- `p_<celltype>.rds`     — two-sided p-value matrix
- `precision_<celltype>.rds` — sparse precision matrix

Plus a single workbook:

- `degree.xlsx` — node degrees, with one column per cell type.

## Installation

```r
# Dependencies (CRAN)
install.packages(c("Matrix", "data.table", "glassoFast", "openxlsx"))

# For reading .h5ad: one of the two backends below.
#
# (Recommended, especially on HPC) anndata — CRAN, uses whatever
# Python is already on your PATH via reticulate. Fast install, no
# private Python compile:
install.packages("anndata")
# And in your Python environment: pip install anndata

# Alternatively: zellkonverter (Bioconductor). Bundles its own private
# Python via basilisk. Works cleanly on workstations, but on HPC
# nodes the private Python build often fails (missing Tk / OpenSSL /
# build headers). Prefer the anndata backend there.
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("zellkonverter", "SingleCellExperiment",
#                        "SummarizedExperiment"))

# Install GECO from GitHub (recommended)
# install.packages("remotes")            # if you don't have it
remotes::install_github("xinpeiw/GECO")
# or with devtools:
# devtools::install_github("xinpeiw/GECO")

# Alternatively, from a local clone:
# install.packages("path/to/GECO", repos = NULL, type = "source")
# devtools::install("path/to/GECO")
```

Repository: <https://github.com/xinpeiw/GECO> (private — collaborators
will need a GitHub PAT with `repo` scope to `remotes::install_github()`,
e.g. via `Sys.setenv(GITHUB_PAT = "...")`).

## Running on a cluster

The pipeline requires `fst` (`install.packages("fst")`) for
column-selectable on-disk cell data.

### Running it

The simplest path is to let GECO generate ready-to-`sbatch` scripts
for you:

```r
library(GECO)

generate_slurm_pipeline(
  scripts_dir  = "/path/to/scripts",
  work_dir       = "/path/to/work",
  output_dir     = "/path/to/results",
  h5ad_file      = "/path/to/data_post_qc.h5ad",
  gene_list_file = "/path/to/genelist.txt",  # one gene ID per line
  n_jobs         = 300,   # number of SLURM array chunks (chunk step)
  fixed_covariates_file = "/path/to/covariates.rds",  # optional, see below
  slurm = list(partition = "general", time = "24:00:00",
               mem = "32G", cpus = 1)
)
```

This writes five `.R` + `.sh` files plus a `pipeline.sh` orchestrator
into `scripts_dir`. Adjust the SLURM headers if needed, then:

```bash
cd /path/to/scripts
bash pipeline.sh
```

`pipeline.sh` chains the first three jobs with `--dependency=afterok`.
Step 02 then submits steps 03 and 04 itself once it knows the number
of cell types (read from `cell_type_order.rds`):

```
00_prepare.sh   →  one-time setup (h5ad → fst + OPms + fixed_cov + gene_pairs)
01_chunk.sh     →  SLURM array, runs geco_chunk per task_id
02_finalize_estimates.sh
                →  merge chunks, reformat into per-cell-type matrices,
                   compute p-values, PSD-transform; then sbatch:
   03_precision.sh →  SLURM array (one per "shared" + cell type),
                      glasso+eBIC
   04_finalize_outputs.sh
                   →  build degree.xlsx and write all matrices to output_dir
```

You can also call the functions directly without the script generator:

```r
# Once
genes <- read_gene_list("genelist.txt")
geco_prepare(h5ad_file = "...", work_dir = "work/",
             genes = genes,
             fixed_covariates = read_covariates("covariates.rds"))  # or NULL

# Per chunk (SLURM array, 1..n_jobs)
geco_chunk(work_dir = "work/", task_id = task_id, n_jobs = 300)

# After all chunks finish — writes cov/se/p directly to output_dir
geco_finalize_estimates(work_dir = "work/", output_dir = "results/")

# Per cell type (SLURM array, 1..(n_ct + 1) for shared)
geco_precision_chunk(work_dir = "work/", output_dir = "results/", ct_idx = i)

# Finally — builds degree.xlsx in output_dir
geco_finalize_outputs(work_dir = "work/", output_dir = "results/")
```

### Choosing which covariates to adjust

You supply the per-individual covariates yourself, as a numeric
matrix (or `data.frame`) with one row per individual and one column
per covariate:

- row names = individual IDs exactly as they appear in `ind_col`;
  every individual in the data must have a row (extra rows are
  ignored, and rows are reordered to the data's individual order);
- all columns numeric: encode sex, batch, case/control etc. as 0/1
  (or dummy) columns yourself; bin or scale age however you prefer;
- no intercept column (one is fitted per gene and cell type) and no
  constant column;
- no NA.

Anything else is an error that names the offending individuals or
columns. `fixed_covariates = NULL` (the default) skips the
regression entirely.

With `generate_slurm_pipeline()`, save the matrix to a file and
point to it with `fixed_covariates_file`; the generated
`00_prepare.R` reads it with `read_covariates()` and hands it to
`geco_prepare()`. Accepted formats: `.rds` (a matrix or `data.frame`
with row names) or `.csv` / `.tsv` (first column = individual ID,
header row = covariate names).

```r
# e.g. build it from the h5ad obs table once, in R:
obs  <- as.data.frame(anndata::read_h5ad("data_post_qc.h5ad")$obs)
obs  <- obs[!duplicated(obs$ind_cov), ]
fcov <- data.frame(row.names = obs$ind_cov,
                   age  = as.numeric(obs$Age),
                   male = as.numeric(obs$Sex == "Male"))
saveRDS(fcov, "covariates.rds")

generate_slurm_pipeline(..., fixed_covariates_file = "covariates.rds")

# or, calling geco_prepare() directly:
geco_prepare(..., fixed_covariates = fcov)
```

If you want an expression principal component as a covariate, load
the data with `load_h5ad()`, call `compute_PC1(data, genes)` and
`cbind()` its `PC1` column (rows named by individual) to the matrix.

### What ends up where on disk

After a full run, your two directories look like this:

**`work_dir/`** — intermediate state, can be deleted once everything is done:

```
work_dir/
├── cell_data.fst         # raw cell-level data (fst, column-selectable)
├── prep.rds              # metadata: genes, individuals, cell_types,
│                         # ind_col, ct_col, opms, fixed_cov, n_pairs
├── gene_pairs.rds        # 2 x M matrix of all pairs to estimate
├── Y.rds                 # [n_ind, n_gene, n_ct] precomputed means
├── VS_diag.rds           # [n_ind, n_ct, n_gene] precomputed variances
├── chunks/               # transient — deleted automatically by step 02
│   ├── chunk_00001.rds   # (pass cleanup_chunks = FALSE to keep)
│   ├── chunk_00002.rds
│   └── ...               # one per SLURM array task
├── cell_type_order.rds   # ordered names of "shared" + each cell type
└── precision/
    ├── fit_shared.rds
    ├── fit_CD4NC.rds
    └── ...               # per-cell-type glasso lambda + eBIC paths
```

**`output_dir/`** — the deliverables:

```
output_dir/
├── cov_shared.rds        cov_CD4NC.rds       cov_CD4ET.rds       ...
├── se_shared.rds         se_CD4NC.rds        se_CD4ET.rds        ...
├── p_shared.rds          p_CD4NC.rds         p_CD4ET.rds         ...
├── precision_shared.rds  precision_CD4NC.rds precision_CD4ET.rds ...
└── degree.xlsx
```

### Wall-time breakdown

The chunk step splits ~`choose(n_gene, 2)` pairs across `n_jobs`
SLURM array tasks. Pick `n_jobs` so each task runs in a reasonable
wall time — typically 1–6 hours, well under your queue's time limit.
Numbers assume `n_ind ≈ 1000, n_ct ≈ 7, ~100 cells per (ind, ct)`;
scale linearly with total cell count.

| n_gene | One-time precompute (Y + VS_diag) | n_pairs |  per-pair: n_jobs = 100 | 300 | 500 | 1000 |
|------:|------:|--------:|------:|------:|------:|------:|
| 1,000   |  ~1 m    |    500,000 |   ~2 m / task  |  <1 m   |  <1 m  |  <1 m  |
| 5,000   |  ~5 m    | 12,500,000 |  ~45 m / task  | ~15 m   | ~10 m  |  ~5 m  |
| 10,000  | ~10 m    | 49,995,000 |   ~3 h / task  |  ~1 h   | ~35 m  | ~20 m  |
| 20,000  | ~20 m    |     ~200M  |  ~11 h / task  |  ~4 h   | ~2.5 h |  ~1 h  |

### Memory & disk budget

| Item | Size | Where |
|------|------|------|
| `cell_data.fst` | ~ `n_cells × n_genes × 8` bytes, compressed ~3-4× | disk (one-time) |
| `Y.rds` | `n_ind × n_gene × n_ct × 8` bytes | disk + memory in each chunk |
| `VS_diag.rds` | same as Y | disk + memory in each chunk |
| Per-pair working memory | `O(n_ind × n_ct)` | per chunk worker |
| `chunk_*.rds` per task | `n_pairs_in_chunk × (3 + 3·n_ct) × 8` bytes ×2 (cov+se) | disk — auto-deleted by step 02 unless `cleanup_chunks = FALSE` |
| `cov_*.rds` / `se_*.rds` / `p_*.rds` | `n_gene² × 8` per cell type | disk (output_dir) |
| `precision_*.rds` | `n_gene² × 8` per cell type | disk (output_dir) |

For `n_ind = 1000, n_gene = 10000, n_ct = 7`:

- `cell_data.fst`: depends on `n_cells`, e.g. 1M cells × 11k genes × 4 (compression) ≈ 22 GB
- `Y.rds` and `VS_diag.rds`: each ~620 MB (loaded once per chunk worker)
- `cov_*` / `se_*` / `p_*` / `precision_*` per cell type: ~1 GB each
  (`n_gene² × 8 bytes`, e.g. 11k² × 8 ≈ 970 MB), so ~30 GB total
  across 8 cell types × 4 matrix kinds

Per-chunk worker memory is therefore dominated by `Y.rds + VS_diag.rds + meta` ≈ **2 GB**. Set `#SBATCH --mem=4G` or higher to be safe.

### Resuming after a partial failure

Each stage is idempotent and reads everything it needs from
`work_dir/`. To resume:

- **A few `chunk_*.rds` files failed.** Look at the SLURM output to
  find which task IDs failed, then resubmit just those (note: step 02
  deletes chunks on success, so this only works if step 02 hasn't run
  or you passed `cleanup_chunks = FALSE`):

  ```bash
  sbatch --array=37,42,108 01_chunk.sh
  # then re-run step 02 by hand or via sbatch --dependency
  bash -c 'sbatch --parsable --dependency=afterok:<chunk_jobid> 02_finalize_estimates.sh'
  ```

- **Step 02 finished, one or more `precision_<ct>.rds` failed.**
  Resubmit only that cell type:

  ```bash
  sbatch --array=3 03_precision.sh                # ct_idx = 3
  sbatch --dependency=afterok:<above> 04_finalize_outputs.sh
  ```

- **Step 04 failed but everything upstream succeeded.** Just run it
  again: `bash 04_finalize_outputs.sh` (or
  `Rscript 04_finalize_outputs.R "$WORK_DIR" "$OUTPUT_DIR"`).

If you want to start over, `rm -rf work_dir/` and rerun
`bash pipeline.sh`.

### Running a single SLURM chunk manually

Useful when you want to time-profile one chunk before submitting the
full array, or to reproduce a failure interactively:

```bash
# Allocate an interactive node
salloc --time=2:00:00 --mem=8G --cpus-per-task=1

# Then run task_id = 1 of an n_jobs = 300 split:
SLURM_ARRAY_TASK_ID=1 Rscript 01_chunk.R /path/to/work_dir 300
```

The `SLURM_ARRAY_TASK_ID` env var is what `01_chunk.R` reads; you can
override it for any task ID you want.

### Common pitfalls

- **"`genes` is required" or "missing argument"**. You forgot to pass
  a gene list. Either pass a character vector or use
  `read_gene_list("genelist.txt")`.
- **"X requested genes are not in the h5ad"**. Your gene IDs in the
  list don't match `rowData(adata)$gene_ids` (or whatever
  `gene_id_col` resolves to). Check with `rownames(SummarizedExperiment::assay(zellkonverter::readH5AD("file.h5ad")))[1:10]`.
- **"individual(s) in the data have no row in `fixed_covariates`"**.
  The row names of your covariate matrix don't match the values of
  `ind_col` in the h5ad (check for type coercion, e.g. numeric IDs
  read as `"1"` vs `"1.0"`, or trailing whitespace).
- **Step 03 array bound off**. The auto-sized array uses the number
  of cell types in `work_dir/cell_type_order.rds`. If something looks
  wrong there, open it with `readRDS()` to see the names.
- **`fst` not installed on worker nodes**. The R install used by the
  worker (often a module load or conda env) must have `fst`. Set
  `r_executable = "micromamba run -n my_env Rscript"` in
  `generate_slurm_pipeline()` if you use conda/micromamba.
- **Out-of-memory in a chunk**. The chunk worker holds `Y` +
  `VS_diag` + the cell-data metadata (~2 GB for n_ind=1000,
  n_gene=11k, n_ct=7). Bump `slurm$mem` to 8 GB if you see OOM kills.
- **`zellkonverter` tries to build Python 3.x and fails on the HPC**
  (`BUILD FAILED (RedHatEnterprise ... _tkinter... No such file or
  directory`). `zellkonverter` uses `basilisk`, which downloads and
  compiles its own private Python via `pyenv` on first use — that
  build often fails on locked-down cluster nodes missing `tk-devel`,
  `openssl-devel`, or similar. Switch to the `anndata` backend, which
  uses whatever Python is already on `PATH` (via `reticulate`):

  ```r
  # R side
  install.packages("anndata")

  # Python side (in your loaded module / conda env)
  # pip install anndata

  # Then in R, force the backend:
  geco_prepare(h5ad_file, work_dir, genes,
               backend = "anndata", ...)
  # (with generate_slurm_pipeline(): prepare_args = list(backend = "anndata"))
  ```

  `backend = "auto"` (the default) already prefers `anndata` when both
  packages are installed, so simply having `anndata` present is
  usually enough. Point reticulate at your Python with
  `Sys.setenv(RETICULATE_PYTHON = "/path/to/python")` or
  `reticulate::use_python("/path/to/python", required = TRUE)` if it
  picks the wrong one.

## Input requirements and the "error instead of NA" policy

GECO's model needs a complete individual x cell-type grid, and
several downstream steps (the covariate regression in `GeCo()`, the
`glassoFast` C code) cannot consume NA at all. Rather than letting a
single NA propagate into an all-NA output or a crash deep inside
glasso, **every step stops as soon as it would produce an NA and
tells you why** (which individuals, cell types, genes or gene pairs
are affected, and what to do). There is no "skip and continue"
default anywhere in the pipeline.

Concretely, before running GECO your data must satisfy:

| Requirement | Checked by | What the error names |
|---|---|---|
| Every individual has **>= 2 cells in every cell type** | `compute_Y()`, `compute_VS_diag()`, `geco_prepare()` | the (individual, cell type) groups with 0 or 1 cells, and the individuals to drop |
| At least 2 cell types | `estimate_geco()`, `geco_prepare()`, `GeCo()` | - |
| `ind_col` / `ct_col` have no NA | all steps that group by them | number of cells with a missing label |
| No NA expression values | `compute_OPms()`, `compute_Y()`, `compute_VS_diag()` | the gene(s) |
| No gene is constant across individuals (`OPsd > 0`) | `compute_OPms()`, `compute_PC1()` | the gene(s) |
| `fixed_covariates` is numeric, complete (no NA), has a row for every individual (row names matched to individual IDs and reordered) and no constant column | `geco_prepare()` | the individuals / columns at fault |
| `Y.rds` / `VS_diag.rds` in `work_dir` match `prep.rds` and contain no NA | `geco_chunk()` | the dims / the NA positions |

And at the merge / finalize steps:

- `geco_finalize_estimates()` verifies that the `chunk_*.rds` files
  cover every gene pair exactly once and stops listing the uncovered
  pair ranges and the task_ids it found, so a SLURM task killed by OOM
  or a time limit cannot silently drop its pairs.
- `geco_finalize_outputs()` stops if any `precision_<ct>.rds` is
  missing, naming the cell type and the step-03 `ct_idx` to resubmit.
- `reformat_to_matrices()`, `compute_p_matrix()`, `transform_psd()`
  and `fit_precision()` all refuse NA input. `transform_psd()` keeps an
  explicit opt-in `na_action = "zero"` if you know why entries are
  missing.

The usual fix for the most common error (incomplete grid) is to
filter individuals before loading, e.g.

```r
data <- load_h5ad("data_post_qc.h5ad", genes = genes)
counts <- table(data$ind_cov, data$cg_cov)
keep   <- rownames(counts)[apply(counts >= 2, 1, all)]
data   <- data[data$ind_cov %in% keep, ]
```

## Notes and caveats

- Covariate effects are modelled per gene and shared across cell
  types (one coefficient per covariate per gene), with a separate
  intercept per gene and cell type. `fixed_covariates = NULL` skips
  the regression entirely.
- Cell-type names in the output are taken directly from
  `unique(data[[ct_col]])` in the h5ad file — they are not provided
  by the user. The first matrix in every output list is always named
  `shared` (the across-individual component).
