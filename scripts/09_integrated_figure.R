# 09_integrated_figure.R
# Figure 7: integrated evidence summary for the 3 coloc-tested loci --
# Panel A shows PP4 across all genes tested at each locus (winner highlighted),
# Panel B shows the winning gene's cell-type expression across categories.

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ragg)
coloc <- fread("results/tables/04_coloc_results.csv")
final_table <- fread("results/tables/08_final_prioritization.csv")
expr_summary <- fread("results/tables/07_candidate_gene_expression_summary.csv")
symbol_map <- c(
  ENSG00000146374 = "RSPO3", ENSG00000179195 = "ZNF664", ENSG00000124587 = "PEX6"
)
winners <- final_table[!is.na(gene_symbol), gene_symbol]
# --- Panel A: PP4 distribution per locus, winner highlighted -----------------
coloc[, gene_symbol_lookup := symbol_map[gene_id]]

winner_lookup <- final_table[!is.na(gene_symbol), .(locus_id, gene_symbol)]
coloc[, is_winner := FALSE]
coloc[winner_lookup, is_winner := TRUE, on = .(locus_id, gene_symbol_lookup = gene_symbol)]

panel_a <- ggplot(coloc, aes(x = locus_id, y = PP4)) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.3, color = "grey50") +
  geom_point(data = coloc[is_winner == TRUE], aes(x = locus_id, y = PP4),
             color = "firebrick", size = 3) +
  geom_text_repel(data = coloc[is_winner == TRUE],
                   aes(x = locus_id, y = PP4, label = gene_symbol_lookup),
                   color = "firebrick", size = 3.2, nudge_y = 0.08) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  labs(title = "A. Colocalization evidence across all tested genes",
       subtitle = "Each point = one gene; red = best-PP4 gene per locus",
       x = "Locus", y = "PP4 (posterior probability of shared causal variant)") +
  ylim(0, 1) +
  theme_minimal(base_size = 11)
# --- Panel B: cell-type expression for winning genes -------------------------
winner_expr <- expr_summary[gene %in% winners]
winner_expr[, is_top := pct_expressing == max(pct_expressing), by = gene]

panel_b <- ggplot(winner_expr, aes(x = cell_type, y = pct_expressing, fill = is_top)) +
  geom_col() +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue"), guide = "none") +
  labs(title = "B. Cell-type expression of coloc-prioritized genes",
       subtitle = "GSE129363 adipose SVF; blue = highest-expressing cell type",
       x = NULL, y = "% cells expressing") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

combined <- panel_a / panel_b + plot_layout(heights = c(1, 1.1))

ragg::agg_png("results/figures/06_integrated_evidence.png", width = 8, height = 9, units = "in", res = 300)
print(combined)
dev.off()
cat("Saved integrated evidence figure: results/figures/06_integrated_evidence.png\n")
