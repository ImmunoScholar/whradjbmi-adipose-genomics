# 07_candidate_gene_mapping.R
# Maps GWAS/eQTL/coloc-prioritized genes into the annotated adipose SVF
# single-cell data. Includes both coloc-supported genes (RSPO3, ZNF664) and
# their loci's nearest/eponymous genes that did NOT colocalize (VEGFA, NCOR2)
# plus the weak-coloc gene at locus_02 (PEX6) -- deliberately showing the
# contrast between physical proximity and colocalization evidence, which is
# the central distinction this project is built to demonstrate.

library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)
library(ragg)
seu <- readRDS("data/processed/seurat_annotated.rds")

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Candidate gene x cell-type mapping ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
gene_table <- data.table(
  gene    = c("RSPO3", "ZNF664", "PEX6", "VEGFA", "NCOR2"),
  locus_id = c("locus_01", "locus_03", "locus_02", "locus_02", "locus_03"),
  PP4     = c(0.768, 0.918, 0.252, 0.028, 0.047),
  role    = c("Coloc-supported (best PP4 at locus)", "Coloc-supported (best PP4 at locus)",
              "Weak/ambiguous coloc (best PP4 at locus, still low)",
              "Locus namesake gene -- distinct signal (PP3=0.764), NOT coloc-supported",
              "Locus namesake gene -- coloc evidence weak despite proximity")
)
log_msg("Candidate genes for mapping:")
log_msg(paste(capture.output(print(gene_table)), collapse = "\n"))
genes <- gene_table$gene
present <- genes %in% rownames(seu)
log_msg("Gene presence check in single-cell data:")
for (i in seq_along(genes)) log_msg("  ", genes[i], ": ", if (present[i]) "present" else "NOT FOUND -- check symbol/alias")

genes_present <- genes[present]
if (length(genes_present) < length(genes)) {
  log_msg("WARNING: ", sum(!present), " gene(s) missing from the single-cell matrix -- excluded from plots.")
}
# --- Minimum group-size check (Checkpoint 10: don't trust expression from tiny groups)
group_sizes <- table(seu$cell_type, seu$Tissue)
log_msg("Cell-type x Tissue group sizes (flagging any group < 30 cells as unreliable):")
log_msg(paste(capture.output(print(group_sizes)), collapse = "\n"))
small_groups <- which(group_sizes < 30, arr.ind = TRUE)
if (nrow(small_groups) > 0) {
  for (i in seq_len(nrow(small_groups))) {
    ct <- rownames(group_sizes)[small_groups[i, 1]]
    ts <- colnames(group_sizes)[small_groups[i, 2]]
    log_msg("  FLAG: ", ct, " x ", ts, " has only ", group_sizes[small_groups[i,1], small_groups[i,2]],
            " cells -- any expression pattern here should not be over-interpreted.")
  }
}
# --- Dot plot by broad cell type ---------------------------------------------
dp <- DotPlot(seu, features = genes_present, group.by = "cell_type") +
  RotatedAxis() +
  labs(title = "Candidate gene expression by adipose SVF cell type",
       subtitle = "Dot size = % cells expressing, color = scaled average expression") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ragg::agg_png("results/figures/05_candidate_gene_dotplot.png", width = 7, height = 5, units = "in", res = 300)
print(dp)
dev.off()
log_msg("Saved cell-type dot plot: results/figures/05_candidate_gene_dotplot.png")
# --- Dot plot split by depot (SAT/VAT) ---------------------------------------
seu$cell_type_tissue <- paste0(seu$cell_type, " (", seu$Tissue, ")")
dp_depot <- DotPlot(seu, features = genes_present, group.by = "cell_type_tissue") +
  RotatedAxis() +
  labs(title = "Candidate gene expression by cell type x depot",
       subtitle = "Dot size = % cells expressing, color = scaled average expression") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ragg::agg_png("results/figures/05b_candidate_gene_dotplot_by_depot.png", width = 8, height = 6, units = "in", res = 300)
print(dp_depot)
dev.off()
log_msg("Saved depot-split dot plot: results/figures/05b_candidate_gene_dotplot_by_depot.png")

# --- Numeric summary table: mean expression + pct expressing per gene x cell type
expr_data <- FetchData(seu, vars = c(genes_present, "cell_type", "Tissue"))
setDT(expr_data)  # FetchData returns a base data.frame; data.table syntax below needs setDT
summary_list <- list()
for (g in genes_present) {
  by_ct <- expr_data[, .(mean_expr = mean(get(g)), pct_expressing = 100 * mean(get(g) > 0), n_cells = .N),
                      by = cell_type]
  by_ct[, gene := g]
  summary_list[[g]] <- by_ct
}
gene_summary <- rbindlist(summary_list)
setcolorder(gene_summary, c("gene", "cell_type", "mean_expr", "pct_expressing", "n_cells"))
setorder(gene_summary, gene, -mean_expr)
fwrite(gene_summary, "results/tables/07_candidate_gene_expression_summary.csv")
log_msg("Saved expression summary table: results/tables/07_candidate_gene_expression_summary.csv")
log_msg(paste(capture.output(print(gene_summary)), collapse = "\n"))
# --- Dropout / single-cell-domination check ----------------------------------
# For each gene, report whether "expressed" cells are spread across many cells
# or concentrated in a handful (which would make a dot-plot signal misleading).
log_msg("Dropout / concentration check per gene (top cell type only):")
for (g in genes_present) {
  top_ct <- gene_summary[gene == g][1]
  n_expr <- round(top_ct$n_cells * top_ct$pct_expressing / 100)
  log_msg("  ", g, " in ", top_ct$cell_type, ": ", n_expr, " of ", top_ct$n_cells,
          " cells express it (", round(top_ct$pct_expressing, 1), "%) -- ",
          if (n_expr < 10) "WARNING: fewer than 10 expressing cells, signal may be noisy" else "adequate cell count")
}
writeLines(log_lines, "results/logs/07_candidate_gene_mapping_summary.txt")
log_msg("Summary written: results/logs/07_candidate_gene_mapping_summary.txt")
