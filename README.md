# GECO

**GECO** (**Ge**ne **Co**expression) is an R package that unbiasedly decomposes
interindividual gene-by-gene coexpression into cell-type-shared and cell-type-specific components, and turns the result into cell-type-shared and per-cell-type-specific sparse precision matrices and node
degrees.

## Contents

- [Inputs and outputs](#inputs-and-outputs)
- [Installation](#installation)
- [Running on a cluster](#running-on-a-cluster)
  - [Running it](#running-it)
  - [Choosing which covariates to adjust](#choosing-which-covariates-to-adjust)
  - [What ends up where on disk](#what-ends-up-where-on-disk)
  - [Wall-time breakdown](#wall-time-breakdown)
  - [Resuming after a partial failure](#resuming-after-a-partial-failure)
  - [Common pitfalls](#common-pitfalls)
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
# install.packages("remotes")            # if you don't have it
remotes::install_github("xinpeiw/GECO")
```

Repository: <https://github.com/xinpeiw/GECO> (private — collaborators
will need a GitHub PAT with `repo` scope to `remotes::install_github()`,
e.g. via `Sys.setenv(GITHUB_PAT = "...")`).

## Running on a cluster

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
- **"package 'fst' is not installed in this R library" in a chunk
  job**. The R the workers run (often a module load or conda env) is
  not the one where you ran `check_dependencies()`. Run
  `GECO::check_dependencies()` from that R, or set
  `r_executable = "micromamba run -n my_env Rscript"` in
  `generate_slurm_pipeline()` so every step uses the same R.
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
  # R side (asks before installing the 'anndata' R package):
  GECO::check_dependencies(backend = "anndata")

  # Python side (in your loaded module / conda env):
  # pip install anndata

  # Then force the backend:
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

## Notes and caveats

- Covariate effects are modelled per gene and shared across cell
  types (one coefficient per covariate per gene), with a separate
  intercept per gene and cell type. `fixed_covariates = NULL` skips
  the regression entirely.
- Cell-type names in the output are taken directly from
  `unique(data[[ct_col]])` in the h5ad file — they are not provided
  by the user. The first matrix in every output list is always named
  `shared` (the across-individual component).
