# 08_final_prioritization.R
# Builds the final integrated prioritization table and summary figure.
#
# Transparent tier framework (thresholds stated explicitly, not a fabricated
# composite score):
#   Tier 1 (Strong):     PP4 > 0.5 AND cell-type fold-enrichment > 2x
#   Tier 2 (Suggestive):  PP4 0.2-0.5, OR PP4>0.5 with weak cell-type contrast
#   Tier 3 (Weak):        PP4 < 0.2
#   GWAS-only:            loci 4-8 -- eQTL/coloc/single-cell integration was
#                          not performed for these (config scoped coloc to the
#                          3 most significant loci; this is a stated project
#                          boundary, not a silent omission).
#
# Nearest-gene labels for loci 1-3 come from this pipeline's own coloc
# results. Nearest-gene labels for loci 4-8 are NOT included here, because
# this pipeline never ran a formal nearest-gene annotation step for them --
# only an informal literature cross-reference was done at Checkpoint 3, and
# blending that into the same column as pipeline-derived results would
# overstate what was actually computed.

library(data.table)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Final integrated prioritization ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))

loci <- fread("results/tables/02_lead_loci.csv")
coloc <- fread("results/tables/04_coloc_results.csv")
expr_summary <- fread("results/tables/07_candidate_gene_expression_summary.csv")
# --- 1. Best coloc gene per tested locus -------------------------------------
best_coloc <- coloc[, .SD[which.max(PP4)], by = locus_id]
log_msg("Best coloc gene per tested locus:")
log_msg(paste(capture.output(print(best_coloc[, .(locus_id, gene_id, nsnps, PP4)])), collapse = "\n"))

# Map Ensembl IDs to symbols for the genes we manually verified via Ensembl lookup.
symbol_map <- c(
  ENSG00000146374 = "RSPO3",
  ENSG00000179195 = "ZNF664",
  ENSG00000124587 = "PEX6"
)
best_coloc[, gene_symbol := symbol_map[gene_id]]
# --- 2. Cell-type fold-enrichment for each prioritized gene ------------------
compute_fold <- function(g) {
  sub <- expr_summary[gene == g]
  if (nrow(sub) == 0) return(list(top_cell_type = NA_character_, fold_enrichment = NA_real_))
  sub <- sub[order(-pct_expressing)]
  top <- sub[1]
  others_mean <- mean(sub[-1]$pct_expressing)
  fold <- if (others_mean > 0) top$pct_expressing / others_mean else NA_real_
  list(top_cell_type = top$cell_type, fold_enrichment = round(fold, 2))
}

fold_results <- rbindlist(lapply(best_coloc$gene_symbol, function(g) {
  r <- compute_fold(g)
  data.table(gene_symbol = g, top_cell_type = r$top_cell_type, fold_enrichment = r$fold_enrichment)
}))
best_coloc <- merge(best_coloc, fold_results, by = "gene_symbol")
# --- 3. Assign transparent tier ----------------------------------------------
best_coloc[, tier := fifelse(PP4 > 0.5 & fold_enrichment > 2, "Tier 1: Strong",
                       fifelse(PP4 >= 0.2, "Tier 2: Suggestive", "Tier 3: Weak"))]

log_msg("Tier assignment (PP4 > 0.5 & fold > 2x = Tier 1; PP4 0.2-0.5 = Tier 2; PP4 < 0.2 = Tier 3):")
log_msg(paste(capture.output(print(best_coloc[, .(locus_id, gene_symbol, PP4, top_cell_type, fold_enrichment, tier)])), collapse = "\n"))
# --- 4. Build the full 8-locus table -----------------------------------------
final_table <- merge(
  loci[, .(locus_id, SNP, CHR, POS, P, BETA)],
  best_coloc[, .(locus_id, gene_symbol, nsnps, PP4, top_cell_type, fold_enrichment, tier)],
  by = "locus_id", all.x = TRUE
)
final_table[is.na(tier), tier := "GWAS-only (not coloc-tested)"]
setorder(final_table, P)
fwrite(final_table, "results/tables/08_final_prioritization.csv")
log_msg("Saved final prioritization table (8 loci): results/tables/08_final_prioritization.csv")
log_msg(paste(capture.output(print(final_table)), collapse = "\n"))
# --- 5. Nearest-gene vs coloc-supported-gene contrast table ------------------
# The central "proximity != colocalization" comparison this project is built around.
contrast <- data.table(
  locus_id = c("locus_02", "locus_02", "locus_03", "locus_03"),
  gene = c("VEGFA", "PEX6", "NCOR2", "ZNF664"),
  role = c("Locus namesake / nearest gene", "Best coloc-supported gene at locus",
           "Locus namesake / nearest gene", "Best coloc-supported gene at locus"),
  PP4 = c(coloc[gene_id == "ENSG00000112715", PP4], coloc[gene_id == "ENSG00000124587", PP4],
          coloc[gene_id == "ENSG00000196498", PP4], coloc[gene_id == "ENSG00000179195", PP4]),
  PP3 = c(coloc[gene_id == "ENSG00000112715", PP3], coloc[gene_id == "ENSG00000124587", PP3],
          coloc[gene_id == "ENSG00000196498", PP3], coloc[gene_id == "ENSG00000179195", PP3])
)
fwrite(contrast, "results/tables/08_nearest_vs_coloc_contrast.csv")
log_msg("Saved nearest-gene vs coloc-gene contrast table: results/tables/08_nearest_vs_coloc_contrast.csv")
log_msg(paste(capture.output(print(contrast)), collapse = "\n"))
writeLines(log_lines, "results/logs/08_final_prioritization_summary.txt")
log_msg("Summary written: results/logs/08_final_prioritization_summary.txt")
