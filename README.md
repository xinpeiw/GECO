# GECO

**GECO** (**Ge**ne **Co**expression) is an R package that unbiasedly decomposes
interindividual gene-by-gene coexpression into cell-type-shared and cell-type-specific components, and turns the result into cell-type-shared and per-cell-type-specific sparse precision matrices and node
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

# For reading .h5ad: either zellkonverter (recommended, Bioconductor)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment",
                       "SummarizedExperiment"))
# ...or anndata (CRAN, uses reticulate)
# install.packages("anndata")

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

## Quick start (in-memory mode)

```r
library(GECO)

# Required: a vector of gene IDs to analyse. They must exist in the
# h5ad's expression matrix; any other columns in the .h5ad are ignored.
genes <- read_gene_list("genelist.txt")   # or just: c("ENSG...", "ENSG...", ...)

result <- run_geco(
  h5ad_file        = "subset_data.h5ad",
  genes            = genes,
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
out before the GECO decomposition runs.

The built-in catalog (see `GECO_COVARIATE_CATALOG` or call
`available_covariates()`) provides these names:

| Name | Source column | What it becomes |
|------|---------------|-----------------|
| `age_group` | `Age`        | binned to 1..15 |
| `sex`       | `Sex`        | 1 if Male, else 0 |
| `PC1`       | (derived)    | PC1 of per-individual mean expression |

Pick whichever subset you want — the most common cases:

```r
# No regression at all (default):
run_geco(..., fixed_covariates = NULL)

# Just age and sex:
run_geco(..., fixed_covariates = c("age_group", "sex"))

# All catalog covariates:
run_geco(..., fixed_covariates = c("age_group", "sex", "PC1"))

# Catalog covariates + a study-specific binary one ("Smoker"):
run_geco(..., fixed_covariates = list(
  age_group = GECO_COVARIATE_CATALOG$age_group,
  sex       = GECO_COVARIATE_CATALOG$sex,
  smoker    = list(col = "Smoker",     type = "binary", positive = "Yes"),
  cases     = list(col = "SLE_status", type = "binary", positive_not = "Healthy")
))

# Or compute your own per-individual matrix and pass it directly:
my_mat <- ...  # n_ind x p matrix, rownames = individual IDs
run_geco(..., fixed_covariates = my_mat)
```

## Step-by-step usage

Every step of the pipeline is also exposed in case you want to run only
part of it or inject custom logic:

```r
genes <- read_gene_list("genelist.txt")
data  <- load_h5ad("subset_data.h5ad", genes = genes)

opms  <- compute_OPms(data, genes)
pc1   <- compute_PC1 (data, genes)
fcov  <- build_fixed_cov(data,
                         covariates = c("age_group", "sex", "PC1"),
                         pc1 = pc1)

est   <- estimate_geco(data, genes, opms, fixed_cov = fcov)
mats  <- reformat_to_matrices(est, genes)

p_mat   <- compute_p_matrix(mats$cov$CD4NC, mats$se$CD4NC)
cov_psd <- transform_psd   (mats$cov$CD4NC)

n_ind <- length(unique(data$ind_cov))         # taken from the h5ad
fit   <- fit_precision(cov_psd, p_mat, n_ind = n_ind)
deg   <- compute_degree(fit$precision)
```

## Cluster (HPC) mode

For ~10k-gene runs that don't fit on a single node, use the chunked
pipeline. The cluster mode requires `fst`
(`install.packages("fst")`) for column-selectable on-disk cell data.

### How cluster mode differs from in-memory mode

In-memory mode (`run_geco()`) does everything in one R process: load,
estimate every gene pair, reformat, fit precision matrices, write
outputs. The full cell-level matrix and all per-pair tensors live in
RAM.

Cluster mode breaks the same pipeline into five disk-backed stages
because each SLURM job is a separate R process that doesn't share
memory with the others:

1. **Prepare** — one process loads the .h5ad, writes a compressed
   `cell_data.fst` (column-selectable), computes per-gene OPmean/OPsd,
   builds the fixed covariate matrix, and saves the gene-pair list.
   Everything that downstream workers might need lands in `work_dir/`.
2. **Chunk** — a SLURM array (default 300 tasks). Each worker reads
   `prep.rds` + just the two gene columns it needs for its slice of
   gene pairs from the fst file, runs `GeCo()` on every pair in the
   slice, and writes `work_dir/chunks/chunk_<id>.rds`.
3. **Merge & reformat** — one process concatenates all
   `chunk_*.rds` into a single long table, reformats it into
   per-cell-type matrices, computes p-values, and writes
   `cov_<ct>.rds` / `se_<ct>.rds` / `p_<ct>.rds` for each cell type
   directly to `output_dir`. A small `work_dir/cell_type_order.rds`
   records the cell-type ordering for step 4.
4. **Precision** — a second SLURM array (one task per cell type
   plus `shared`). Each worker reads `cov_<ct>.rds` and `p_<ct>.rds`
   from `output_dir`, PSD-clips the covariance on the fly, runs
   glassoFast + eBIC, and writes `precision_<ct>.rds` back into
   `output_dir`.
5. **Finalize** — one process reads the `precision_*.rds` files,
   computes degrees, and writes `degree.xlsx` to `output_dir`. All
   other deliverables are already in place from steps 3 and 4.

The extra disk hops (especially the chunk-merging step) are why
cluster mode has more files and more invocations than in-memory mode.
In exchange, no single worker ever has to hold the full expression
matrix or all gene-pair tensors at once.

### Running it

The simplest path is to let GECO generate ready-to-`sbatch` scripts
for you:

```r
library(GECO)

generate_slurm_pipeline(
  scripts_dir  = "/path/to/scripts",
  work_dir       = "/path/to/work",
  output_dir     = "/path/to/results",
  h5ad_file      = "/path/to/subset_data.h5ad",
  gene_list_file = "/path/to/genelist.txt",  # one gene ID per line
  n_jobs         = 300,   # number of SLURM array chunks (chunk step)
  prepare_args   = list(fixed_covariates = c("age_group", "sex", "PC1")),
  slurm = list(partition = "general", time = "24:00:00",
               mem = "32G", cpus = 1)
)
# The precision SLURM array is sized automatically by step 02 (after
# cell_type_order.rds is written), so no `n_cell_types` argument is needed.
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
             fixed_covariates = c("age_group", "sex", "PC1"))

# Per chunk (SLURM array, 1..n_jobs)
geco_chunk(work_dir = "work/", task_id = task_id, n_jobs = 300)

# After all chunks finish — writes cov/se/p directly to output_dir
geco_finalize_estimates(work_dir = "work/", output_dir = "results/")

# Per cell type (SLURM array, 1..(n_ct + 1) for shared)
geco_precision_chunk(work_dir = "work/", output_dir = "results/", ct_idx = i)

# Finally — builds degree.xlsx in output_dir
geco_finalize_outputs(work_dir = "work/", output_dir = "results/")
```

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

There are two cost components in the cluster pipeline:

1. **One-time precompute** in `geco_prepare()`: compute `Y` (means per
   `[ind, gene, ct]`) and `VS_diag` (variances per `[ind, ct, gene]`)
   over **all** genes, once. This is paid by a single process at the
   start of the run. Cost ≈ `n_gene × n_cells × constant` —
   roughly **70 ms/gene** at `n_ind = 1000, n_ct = 7, ~100 cells per
   (ind, ct)`.
2. **Per-pair work** inside each chunk worker: one column read from
   `cell_data.fst`, the off-diagonal `cov(g1, g2)` aggregation, and
   the `GeCo()` linear algebra. At the same shape this is
   **~0.16 s per pair**, split roughly as:

   | sub-step | time | what it does |
   |---------:|-----:|--------------|
   | fst read of 2 gene columns | ~5 ms | I/O |
   | Y slice from precomputed tensor | <1 ms | array indexing |
   | per-pair `cov(g1, g2)` per (ind, ct) | ~150 ms | data.table group-by — dominant |
   | `GeCo()` QR + decomposition | ~9 ms | core math |

Total cluster wall time = one-time precompute + per-pair × n_pairs / n_jobs.

#### Total time, all-in-one view

The chunk step splits ~`choose(n_gene, 2)` pairs across `n_jobs`
SLURM array tasks. Pick `n_jobs` so each task runs in a reasonable
wall time — typically 6–24 hours, well under your queue's time limit.
Numbers assume `n_ind ≈ 1000, n_ct ≈ 7, ~100 cells per (ind, ct)`;
scale linearly with total cell count.

| n_gene | One-time precompute (Y + VS_diag) | n_pairs |  per-pair: n_jobs = 100 | 300 | 500 | 1000 |
|------:|------:|--------:|------:|------:|------:|------:|
| 1,000   |  ~1.2 m  |    500,000 |   ~12 m / task |  ~4 m  |  ~3 m  |  ~1 m  |
| 5,000   |  ~6 m    | 12,500,000 |   ~5 h  / task |  ~1.7 h |  ~1 h  | ~30 m  |
| 10,000  | ~12 m    | 49,995,000 |  ~21 h  / task |  ~7 h   |  ~4 h  |  ~2 h  |
| 20,000  | ~24 m    |     ~200M  |  ~83 h  / task | ~28 h   | ~17 h  |  ~8 h  |

The precompute is small relative to the chunked work for any real
gene count, and it's only paid once. If a chunk runs out of time, just
resubmit only the failed task IDs (see "Resuming after a partial
failure" below).

### Supplying your own Y or VS_diag

`geco_prepare()` precomputes two large tensors once before the chunk
step and reuses them for every gene pair:

- `Y.rds` — `[n_ind, n_gene, n_ct]` per-(individual, gene, cell-type) means
- `VS_diag.rds` — `[n_ind, n_ct, n_gene]` per-(individual, cell-type) variances

If you already have these (e.g. precomputed on a different machine
or from a non-`.h5ad` source) you can hand them to `geco_prepare()`
via `Y_path` / `VS_diag_path` and skip the recompute:

```r
geco_prepare(
  h5ad_file        = "data.h5ad",
  work_dir         = "work/",
  genes            = read_gene_list("genelist.txt"),
  fixed_covariates = c("age_group", "sex", "PC1"),
  Y_path           = "/path/to/my_Y.rds",
  VS_diag_path     = "/path/to/my_VS_diag.rds"
)
```

The two `.rds` files are simply copied into `work_dir/Y.rds` and
`work_dir/VS_diag.rds`. Make sure their shapes/orderings match the
gene and individual order in `prep$genes` / `prep$individuals` — the
chunk step indexes them positionally.

For in-memory mode, `run_geco()` exposes the same hand-off via
`Y_all = ...` and `VS_diag = ...` arguments (in-memory arrays, not
file paths).

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

### Smoke-testing on your laptop first

Before submitting to SLURM, run the whole thing locally on a tiny
gene list (say 50 genes) to verify your inputs are wired correctly.
This takes a few minutes and shakes out the common configuration
issues (wrong column names, missing covariate metadata, etc.):

```r
library(GECO)
genes <- read_gene_list("genelist.txt")[1:50]

work_dir   <- tempfile("geco_test_")
output_dir <- tempfile("geco_test_out_")

# 1. Prepare
geco_prepare(h5ad_file        = "subset_data.h5ad",
             work_dir         = work_dir,
             genes            = genes,
             fixed_covariates = c("age_group", "sex", "PC1"))

# 2. Run all chunks serially (no SLURM)
n_jobs <- 4
for (tid in seq_len(n_jobs)) {
  geco_chunk(work_dir, task_id = tid, n_jobs = n_jobs)
}

# 3-5. Finalize
geco_finalize_estimates(work_dir, output_dir)
order <- readRDS(file.path(work_dir, "cell_type_order.rds"))
for (i in seq_along(order)) {
  geco_precision_chunk(work_dir, output_dir, ct_idx = i)
}
geco_finalize_outputs(work_dir, output_dir)

list.files(output_dir)
```

If that works, swap `subset_data.h5ad` and the full gene list into
`generate_slurm_pipeline()` and submit.

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
- **Covariate column missing**. `build_fixed_cov()` complains that
  e.g. `"Age"` isn't a column. Either rename your h5ad metadata
  columns or pass a custom spec (see "Choosing which fixed covariates
  to regress out").
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

## Notes and caveats

- In-memory mode keeps the cell-level expression matrix and all
  per-pair tensors in RAM. For ~10k genes this is too large; switch
  to cluster mode.
- Cluster mode stores the cell-level data once as `cell_data.fst`
  and each chunk only reads the two gene columns it needs, so each
  worker's memory footprint stays small.
- `fixed_covariates = c("age_group", "sex", "PC1")` is the
  default-recommended set. See **Choosing which fixed covariates to
  regress out** above for the full menu (and how to add study-specific
  binary covariates like population or case/control). `fixed_covariates
  = NULL` skips the regression entirely.
- Cell-type names in the output are taken directly from
  `unique(data[[ct_col]])` in the h5ad file — they are not provided
  by the user. The first matrix in every output list is always named
  `shared` (the across-individual component).
