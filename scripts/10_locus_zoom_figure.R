#!/usr/bin/env Rscript
# 10_locus_zoom_figure.R
# Replaces the P-vs-P scatter diagnostic (weak: doesn't show positional
# alignment) with a proper locus-zoom style plot: GWAS and eQTL -log10(P)
# stacked vertically, aligned by shared genomic position (hg19, from the
# already-harmonized rsID-matched table), lead GWAS SNP marked with a
# vertical guide line across both panels. Built for the two Tier-1 hits
# (RSPO3 @ locus_01, ZNF664 @ locus_03) referenced in the README.

library(data.table)
library(ggplot2)
library(patchwork)
library(ragg)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

harmonized <- fread("data/interim/04_harmonized_locus_eqtl.csv")
coloc <- fread("results/tables/04_coloc_results.csv")
loci <- fread("results/tables/02_lead_loci.csv")

symbol_map <- c(ENSG00000146374 = "RSPO3", ENSG00000179195 = "ZNF664")

make_locus_zoom <- function(locus, target_gene, gene_symbol, pp4, out_path) {
  gd <- harmonized[locus_id == locus & gene_id == target_gene]
  lead_pos <- loci[locus_id == locus, POS]
  lead_snp <- loci[locus_id == locus, SNP]

  gd[, is_lead := rsid == tstrsplit(lead_snp, ":", fixed = TRUE, keep = 1L)[[1]]]

  base_theme <- theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          panel.grid.minor = element_blank())

  p_top <- ggplot(gd, aes(x = POS, y = -log10(P_gwas))) +
    geom_point(color = "grey45", size = 1.2, alpha = 0.6) +
    geom_point(data = gd[is_lead == TRUE], color = "firebrick", size = 2.6) +
    geom_vline(xintercept = lead_pos, linetype = "dashed", color = "firebrick", alpha = 0.5) +
    labs(y = expression(-log[10](GWAS~italic(P))), x = NULL) +
    base_theme

  p_bottom <- ggplot(gd, aes(x = POS, y = -log10(P_eqtl))) +
    geom_point(color = "steelblue4", size = 1.2, alpha = 0.6) +
    geom_point(data = gd[is_lead == TRUE], color = "firebrick", size = 2.6) +
    geom_vline(xintercept = lead_pos, linetype = "dashed", color = "firebrick", alpha = 0.5) +
    labs(y = bquote(-log[10]("adipose eQTL"~italic(P))),
         x = paste0("Position, chr", loci[locus_id == locus, CHR], " (hg19)")) +
    base_theme

  combined <- p_top / p_bottom +
    plot_annotation(
      title = paste0(gene_symbol, " (", locus, "): GWAS and adipose eQTL signal by position"),
      subtitle = paste0("PP4 = ", pp4, "  |  red = lead GWAS variant  |  dashed line = lead variant position"),
      theme = theme(plot.title = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 10, color = "grey30"))
    )

  ragg::agg_png(out_path, width = 7, height = 6, units = "in", res = 300)
  print(combined)
  dev.off()
  cat("Saved:", out_path, "\n")
}

rspo3_pp4 <- coloc[gene_id == "ENSG00000146374" & locus_id == "locus_01", round(PP4, 3)]
znf664_pp4 <- coloc[gene_id == "ENSG00000179195" & locus_id == "locus_03", round(PP4, 3)]

make_locus_zoom("locus_01", "ENSG00000146374", "RSPO3", rspo3_pp4,
                "results/figures/02_coloc_locus01_rspo3.png")
make_locus_zoom("locus_03", "ENSG00000179195", "ZNF664", znf664_pp4,
                "results/figures/02_coloc_best_locus_diagnostic.png")
