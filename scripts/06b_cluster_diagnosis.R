# 06b_cluster_diagnosis.R
# Diagnostic follow-up on Checkpoint 9's finding: clusters 0 and 5 were
# assigned "Endothelial" by argmax over all-negative module scores, which is
# not real evidence. This script finds each cluster's actual top marker genes
# (data-driven) instead of relying on our small pre-chosen marker panel.
library(Seurat)
library(data.table)
seu <- readRDS("data/processed/seurat_annotated.rds")

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}
log_msg("=== Cluster identity diagnosis ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))

# Restore numeric cluster identities from the stored PCA/clustering (cell_type
# was derived from these but we want the original cluster numbers back).
Idents(seu) <- seu$seurat_clusters

markers_all <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
setDT(markers_all)

top10 <- markers_all[order(cluster, -avg_log2FC), .SD[1:10], by = cluster]
fwrite(top10, "results/tables/06b_cluster_top_markers.csv")
log_msg("Saved top-10 marker genes per cluster: results/tables/06b_cluster_top_markers.csv")

log_msg("Top 10 marker genes for cluster 0 (currently mislabeled Endothelial via weak argmax):")
log_msg(paste(top10[cluster == "0", gene], collapse = ", "))
log_msg("Top 10 marker genes for cluster 5 (currently mislabeled Endothelial via weak argmax):")
log_msg(paste(top10[cluster == "5", gene], collapse = ", "))

log_msg("Top 10 marker genes for cluster 8 (the confident Endothelial cluster, score=0.631, for comparison):")
log_msg(paste(top10[cluster == "8", gene], collapse = ", "))

# Also check for mesothelial markers specifically (MSLN, WT1, UPK3B), since
# the paper describes a mesothelial-derived VAT-specific progenitor subtype
# that our 3-gene Progenitor panel (PDGFRA/DCN/LUM) may have missed.
meso_genes <- intersect(c("MSLN", "WT1", "UPK3B", "KRT19"), rownames(seu))
log_msg("Checking mesothelial markers (", paste(meso_genes, collapse = ", "), ") by cluster:")
if (length(meso_genes) > 0) {
  meso_expr <- FetchData(seu, vars = c(meso_genes, "seurat_clusters"))
  meso_means <- aggregate(meso_expr[, meso_genes, drop = FALSE], by = list(cluster = meso_expr$seurat_clusters), mean)
  log_msg(paste(capture.output(print(meso_means)), collapse = "\n"))
}
writeLines(log_lines, "results/logs/06b_cluster_diagnosis_summary.txt")
log_msg("Summary written: results/logs/06b_cluster_diagnosis_summary.txt")
