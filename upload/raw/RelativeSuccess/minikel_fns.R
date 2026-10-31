# Functions extracted verbatim (via parse/deparse) from Minikel et al. 2024,
# github.com/ericminikel/genetic_support, src/gensup_analysis.R (CC-BY-4.0):
# pipeline_best() and abs_or(). Nothing else from that script is run.

abs_or <-
function(odds_ratio) {
    abs_odds_ratio = odds_ratio
    flip_indices = odds_ratio < 1 & !is.na(odds_ratio)
    abs_odds_ratio[flip_indices] = 1/odds_ratio[flip_indices]
    return(abs_odds_ratio)
}

pipeline_best <-
function(merged_table, basis = "ti", phase = "combined", require_insight = TRUE, share_mode = "L2G", min_share = 0.5, max_share = 1, worst_rank = Inf, min_h4 = 0.9, include_missing = FALSE, associations = c("OMIM", "GWAS"), otg_subcat = c(""), genebass_subcat = NULL, mendelian_mechanism = "", min_year = 2005, max_year = 2022, firstyear = F, minusomim = F, lacking = NULL, andalso = NULL, minusothersubcat = F, mingenecount = 0, maxgenecount = Inf, mapping_basis = "all", min_beta = 0, max_beta = Inf, 
    min_or = 1, max_or = Inf, min_maf = 0, max_maf = 1, threshold = 0.8, network_list = NA, verbose = T) {
    start_time = Sys.time()
    mtable = merged_table
    if (verbose) {
        cat(file = stderr(), "Starting row count: ", nrow(mtable), "\n")
        flush.console()
    }
    if (basis %in% c("di_mesh", "drug-indication")) {
        mtable$di_uid = paste0(mtable$drugid, "-", mtable$indication_mesh_id)
        mtable$uid = mtable$di_uid
    }
    else if (basis %in% c("ti", "target-indication")) {
        mtable$ti_uid = paste0(mtable$gene, "-", mtable$indication_mesh_id)
        mtable$uid = mtable$ti_uid
    }
    else if (basis == "drug") {
        mtable$uid = mtable$drugid
    }
    if (phase == "active") {
        meta = meta_acat
        mtable$cat = mtable$acat
    }
    else if (phase == "historical") {
        meta = meta_hcat
        mtable$cat = mtable$hcat
    }
    else if (phase == "combined") {
        meta = meta_ccat
        mtable$cat = mtable$ccat
    }
    mtable$catnum = meta$num[match(mtable$cat, meta$cat)]
    mtable = mtable[mtable$cat != "Other", ]
    mtable$assoc_share = mtable$l2g_share
    mtable$assoc_rank = mtable$l2g_rank
    if (verbose) {
        cat(file = stderr(), "Selecting user-specified filters...")
        flush.console()
    }
    if (!include_missing) {
        mtable = mtable[mtable$gene != "" & mtable$indication_mesh_id != "" & !is.na(mtable$gene) & !is.na(mtable$indication_mesh_id), ]
    }
    if (require_insight) {
        mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id[indic$genetic_insight != "none"], ]
    }
    else {
        mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id, ]
    }
    if (minusomim) {
        omim_supported_ti <- mtable %>% filter(comb_norm >= threshold) %>% filter(assoc_source %in% c("OMIM")) %>% select(ti_uid)
        mtable <- mtable %>% filter(is.na(assoc_source) | !(ti_uid %in% omim_supported_ti$ti_uid))
    }
    if (!is.null(lacking)) {
        lackable_supported_ti <- mtable %>% filter(comb_norm >= threshold) %>% filter(assoc_source %in% lacking) %>% filter(!(assoc_source %in% "OTG" & l2g_share >= min_share)) %>% select(ti_uid)
        mtable <- mtable %>% filter(is.na(assoc_source) | !(ti_uid %in% lackable_supported_ti$ti_uid))
    }
    if (!is.null(andalso)) {
        andalso_supported_ti <- mtable %>% filter(comb_norm >= threshold) %>% filter(assoc_source %in% andalso) %>% filter(!(assoc_source %in% "OTG" & l2g_share < min_share)) %>% select(ti_uid)
        mtable <- mtable %>% filter(is.na(assoc_source) | (ti_uid %in% andalso_supported_ti$ti_uid))
    }
    source_map = tibble(source = c("OTG", "PICCOLO", "Genebass", "OMIM", "intOGen"), source_name = c("GWAS", "GWAS", "GWAS", "OMIM", "Somatic"))
    if (!(identical(associations, c("OMIM", "GWAS", "Somatic")))) {
        mtable <- mtable %>% left_join(source_map, by = c(assoc_source = "source")) %>% filter(is.na(source_name) | source_name %in% associations | assoc_source %in% associations)
    }
    if (otg_subcat != "") {
        mtable <- mtable %>% filter(is.na(assoc_source) | assoc_source %in% "OTG")
        mtable <- mtable %>% mutate(gwas_source = case_when(grepl("GCST", original_link) ~ "GWAS Catalog", grepl("FINNGEN", original_link) ~ "FinnGen", grepl("NEALE", original_link) ~ "Neale UKBB", TRUE ~ "Other")) %>% filter(is.na(assoc_source) | gwas_source %in% otg_subcat)
    }
    if (!is.null(genebass_subcat)) {
        grepstring = paste(genebass_subcat, collapse = "|")
        mtable <- mtable %>% filter(is.na(assoc_source) | !(assoc_source %in% "Genebass") | grepl(grepstring, assoc_info))
    }
    if (mendelian_mechanism != "") {
        mtable <- mtable %>% filter(!(mtable$assoc_source %in% "OMIM") | is.na(mtable$assoc_info) | grepl(mendelian_mechanism, mtable$assoc_info))
    }
    if (share_mode == "V2G") {
        mtable$assoc_share = mtable$v2g_share
        mtable$assoc_rank = mtable$v2g_rank
        assoc$assoc_share = assoc$v2g_share
    }
    else if (share_mode == "L2G") {
        mtable$assoc_share = mtable$l2g_share
        mtable$assoc_rank = mtable$l2g_rank
        assoc$assoc_share = assoc$l2g_share
    }
    if (worst_rank < Inf) {
        mtable <- mtable %>% filter(!(assoc_source %in% "OTG") | mtable$assoc_rank <= worst_rank)
    }
    mtable <- mtable %>% filter(!(assoc_source %in% "OTG") | (!is.na(assoc_share) & assoc_share >= min_share & assoc_share <= max_share))
    if (min_h4 > 0.9) {
        mtable <- mtable %>% filter(!(assoc_source %in% "PICCOLO") | (!is.na(mtable$pic_h4) & mtable$pic_h4 >= min_h4))
    }
    if (mingenecount > 0 | maxgenecount < Inf) {
        gene_counts <- assoc %>% filter(source %in% associations) %>% filter(source != "OTG" | (assoc_share >= min_share & assoc_share <= max_share)) %>% group_by(mesh_id) %>% summarize(.groups = "keep", n_genes = length(unique(gene)))
        mtable$gene_count = gene_counts$n_genes[match(mtable$assoc_mesh_id, gene_counts$mesh_id)]
        mtable <- mtable %>% filter(is.na(gene_count) | gene_count >= mingenecount & gene_count <= maxgenecount)
    }
    if (firstyear) {
        ti_first_sup <- mtable %>% filter(comb_norm >= threshold) %>% group_by(ti_uid) %>% summarize(.groups = "keep", min_assoc_year = min(assoc_year))
        mtable$min_assoc_year = ti_first_sup$min_assoc_year[match(mtable$ti_uid, ti_first_sup$ti_uid)]
        mtable <- mtable %>% filter(is.na(assoc_year) | assoc_year == min_assoc_year)
    }
    if (sum(!is.na(mtable$assoc_year)) == 0) {
        mtable_intrinsic_max_year = 2021
        mtable_intrinsic_min_year = 2000
    }
    else {
        mtable_intrinsic_max_year = max(mtable$assoc_year, na.rm = T)
        mtable_intrinsic_min_year = min(mtable$assoc_year, na.rm = T)
    }
    if (min_year > mtable_intrinsic_min_year | max_year < mtable_intrinsic_max_year) {
        mtable <- mtable %>% filter(is.na(assoc_source) | assoc_source != "OTG" | assoc_year >= min_year & assoc_year <= max_year)
    }
    mtable$abs_beta = abs(assoc$beta[match(mtable$arow, assoc$arow)])
    if (min_beta > 0 | max_beta < Inf) {
        stopifnot(associations == "OTG")
        mtable <- mtable %>% filter(is.na(assoc_source) | (!is.na(mtable$abs_beta) & mtable$abs_beta >= min_beta & mtable$abs_beta <= max_beta))
    }
    mtable$abs_or = abs_or(assoc$odds_ratio[match(mtable$arow, assoc$arow)])
    if (min_or > 1 | max_or < Inf) {
        stopifnot(associations == "OTG")
        mtable <- mtable %>% filter(is.na(assoc_source) | (!is.na(mtable$abs_or) & mtable$abs_or >= min_or & mtable$abs_or < max_or))
    }
    if (min_maf > 0 | max_maf < 1) {
        mtable$lead_maf = pmin(mtable$af_gnomad_nfe, 1 - mtable$af_gnomad_nfe)
        mtable$lead_maf[!is.na(mtable$lead_maf) & mtable$lead_maf < 0] = NA
        mtable <- mtable %>% filter(is.na(assoc_source) | !(assoc_source %in% c("OTG", "PICCOLO", "Genebass")) | (!is.na(lead_maf) & (lead_maf >= min_maf & lead_maf < max_maf)))
    }
    if (verbose) {
        cat(file = stderr(), "Selecting highest phase reached and best genetic similarity...")
        flush.console()
    }
    suppressWarnings(step1 <- mtable %>% group_by(uid) %>% summarize(maxsim = max(comb_norm, na.rm = T), maxcat = max(catnum, na.rm = T)))
    if (verbose) {
        cat(file = stderr(), nrow(step1), "rows remain.\n")
        flush.console()
    }
    if (verbose) {
        cat(file = stderr(), "Joining back in program details...")
        flush.console()
    }
    mtable <- mtable %>% filter(uid %in% step1$uid & comb_norm %in% unique(step1$maxsim) & catnum %in% unique(step1$maxcat))
    step2 <- step1 %>% left_join(mtable, by = c(uid = "uid", maxsim = "comb_norm", maxcat = "catnum")) %>% rename(similarity = maxsim, catnum = maxcat)
    if (verbose) {
        cat(file = stderr(), nrow(step2), "rows found after join.\n")
        flush.console()
    }
    if (verbose) {
        cat(file = stderr(), "Removing duplicates resulting from ties...")
        flush.console()
    }
    step2 <- step2 %>% mutate(highest_phase_with_known_outcome = case_when(!is.na(succ_3_a) ~ 3, !is.na(succ_2_3) ~ 2, !is.na(succ_1_2) ~ 1, !is.na(succ_p_1) ~ 0)) %>% arrange(uid, desc(highest_phase_with_known_outcome)) %>% group_by(uid) %>% slice(1) %>% ungroup()
    if (verbose) {
        cat(file = stderr(), nrow(step2), "rows remain.\n")
        flush.console()
    }
    step2$areas = indic$areas[match(step2$indication_mesh_id, indic$indication_mesh_id)]
    step2$genetic_insight = replace_na(indic$genetic_insight[match(step2$indication_mesh_id, indic$indication_mesh_id)], "none")
    step2$target_status = ""
    step2$target_status[step2$similarity >= threshold] = "genetically supported target"
    step2$target_status[step2$similarity < threshold] = "unsupported target"
    step2$target_status[require_insight & step2$genetic_insight == "none"] = "indication lacks genetic insight"
    step2$target_status[is.na(step2$gene) | step2$gene == ""] = "no target annotated"
    step2$target_status[is.na(step2$indication_mesh_id) | step2$indication_mesh_id == ""] = "no indication annotated"
    if (verbose) {
        cat(file = stderr(), paste0("Using sim threshold ", threshold, ", \"genetically supported target\" rows: ", sum(step2$target_status == "genetically supported target"), "....\n"))
        time_elapsed = (Sys.time() - start_time)
        cat(file = stderr(), "pipeline_best completed in", round(time_elapsed, 1), units(time_elapsed), ".\n")
        flush.console()
    }
    return(step2)
}

