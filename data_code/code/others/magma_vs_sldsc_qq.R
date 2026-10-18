rm(list=ls())
suppressMessages(library(data.table))

outdir <- "results/others"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

d <- fread(file.path(outdir, "magma_vs_sldsc_logp.tsv"))
stopifnot(nrow(d) == 624L)
qq <- melt(d[, .(trait, group, network, networkMAGMA = p_magma, `S-LDSC` = p_ldsc)],
           measure.vars = c("networkMAGMA","S-LDSC"),
           variable.name = "analysis", value.name = "p")
qq[, netclass := factor(fifelse(network == "Shared", "Shared network", "Cell type-specific networks"),
                        c("Shared network", "Cell type-specific networks"))]
qq[, group := factor(group, c("Immune","Non-immune"))]
setorder(qq, analysis, netclass, group, p)
qq[, obs := -log10(p)]
qq[, exp := -log10(ppoints(.N)), by = .(analysis, netclass, group)]
fwrite(qq[, .(analysis, netclass, group, network, trait, p, expected_neglog10p = exp, observed_neglog10p = obs)],
       file.path(outdir, "magma_vs_sldsc_qq_points.tsv"), sep = "\t")
cat("QQ points:", nrow(qq), "\n")
