# Inputs: raw/networkMAGMA/magma/results/<trait>.genes.out from raw/networkMAGMA/01_MAGMA.sh; the networks built by
#   raw/networkMAGMA/networks.R from deposited matrices -- GECO and CLUES partial correlations from the precision matrices in
#   data/onek1k_precision/ and data/others/precision_clues/ (rho_ij = -omega_ij / sqrt(omega_ii omega_jj)); binary (non-zero = 1)
#   networks for the glasso-01 variant (data/others/glasso_precision/, GECO precision matrices) and for the perturb-seq variants
#   (data/others/pertseq_union_network_matrix.rds; overlap = edge non-zero in both the GECO and the perturb-seq network, on the
#   6,603 intersection genes); Entrez-to-Ensembl mapping from Bioconductor org.Hs.eg.db.
# Usage: Rscript raw/networkMAGMA/02_networkMAGMA.R <variant> <task_id>, variant in geco, clues, perturbseq, perturbseq_ct,
#   glasso01; one job per (network, trait) with task_id = (network_id - 1) * 78 + trait_id (run as a cluster array);
#   the seed is network_id * 1000 + trait_id with networks in the order of VARIANTS in networks.R and traits sorted by name.
# Output: raw/networkMAGMA/results_<variant>/<network>_<trait>.txt; raw/networkMAGMA/03_collect.R assembles them into data/.
args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]; task_id <- as.integer(args[2])
source("raw/networkMAGMA/networkMAGMA_fns.R")
source("raw/networkMAGMA/networks.R")

cts <- VARIANTS[[variant]]
traits <- sort(sub("[.]sumstats[.]gz$", "", list.files("raw/sLDSC/sumstats", pattern = "[.]sumstats[.]gz$")))
stopifnot(length(traits) == 78, task_id >= 1, task_id <= length(cts) * 78)

network_id <- (task_id - 1) %/% 78 + 1
trait_id   <- (task_id - 1) %% 78 + 1
trait <- traits[trait_id]
network <- build_network(variant, cts[network_id])
seed <- network_id * 1000 + trait_id

d <- as.data.frame(networkMAGMA(network, file.path("raw/networkMAGMA/magma/results", paste0(trait, ".genes.out")),
                                n_perm = 100000, seed = seed, null = "global", n_bins = 10, degree_type = "unweighted",
                                mhc_filter = "drop_internal_edges", mhc_chr = "6", mhc_start = 25000000, mhc_end = 35000000))
d$seed <- seed; d$trait <- trait; d$network <- cts[network_id]

outdir <- file.path("raw/networkMAGMA", paste0("results_", variant))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.table(d, file.path(outdir, paste0(cts[network_id], "_", trait, ".txt")),
            row.names = FALSE, col.names = TRUE, sep = " ", quote = FALSE)
