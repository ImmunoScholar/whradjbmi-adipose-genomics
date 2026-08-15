# 06_singlecell_seurat.R
# Builds a Seurat object from GSE129363 (Vijay et al. 2020, Nat Metab), subsets
# to the authors' own QC-passed cells, clusters, and annotates broad cell
# types via canonical markers.
#
# Includes a Mesothelial category (MSLN/UPK3B/KRT19) alongside
# Progenitor/Immune/Endothelial -- added after Checkpoint 9 caught two large
# clusters (0, 5) being mislabeled "Endothelial" by argmax over all-negative
# module scores. Their real top marker genes (MSLN, ITLN1/omentin, KRT19,
# UPK3B) match the paper's own description of VAT-specific mesothelial-
# derived progenitor clusters (P1/P3) almost exactly. Cluster identity is now
# only assigned when the winning module score is genuinely positive; ties/
# all-negative clusters are labeled "Unclassified" rather than forced.
library(Seurat)
library(data.table)
library(yaml)
library(ggplot2)
library(ragg)

cfg <- yaml::read_yaml("config/config.yml")
sc_dir <- cfg$single_cell$matrix_dir
set.seed(cfg$seed)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Single-cell processing: Seurat pipeline (v2, post Checkpoint 9 fix) ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))

# --- 1. Build Seurat object --------------------------------------------------
mat <- ReadMtx(
  mtx    = file.path(sc_dir, "GSE129363_Discovery_Cohort_matrix.mtx.gz"),
  cells  = file.path(sc_dir, "GSE129363_Discovery_Cohort_barcodes.tsv.gz"),
  features = file.path(sc_dir, "GSE129363_Discovery_Cohort_genes.tsv.gz"),
  feature.column = 2
)
log_msg("Raw matrix loaded: ", nrow(mat), " genes x ", ncol(mat), " cells")

anno <- fread(file.path(sc_dir, "GSE129363_Discovery_Cohort_CellAnnotation.txt.gz"))
anno_df <- as.data.frame(anno)
rownames(anno_df) <- anno_df$CellID

# --- 2. Subset to annotated (QC-passed) cells only --------------------------
n_before <- ncol(mat)
mat <- mat[, colnames(mat) %in% anno_df$CellID]
n_after <- ncol(mat)
log_msg("Subset to annotated cells: ", n_after, " of ", n_before,
        " (confirmed via Vijay et al. 2020 Methods as the authors' own QC exclusions)")

seu <- CreateSeuratObject(counts = mat, project = "GSE129363", min.cells = 3, min.features = 200)
seu <- AddMetaData(seu, metadata = anno_df[colnames(seu), c("SampleName", "Condition", "Tissue")])
log_msg("Seurat object created: ", ncol(seu), " cells x ", nrow(seu), " genes")
log_msg("Tissue distribution: ", paste(capture.output(print(table(seu$Tissue))), collapse = " "))
log_msg("Condition distribution: ", paste(capture.output(print(table(seu$Condition))), collapse = " "))
# --- 3. QC metrics ------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

qc_before <- ncol(seu)
seu <- subset(seu, subset = nFeature_RNA >= cfg$single_cell$min_genes_per_cell &
                            percent.mt <= cfg$single_cell$max_mito_pct)
qc_after <- ncol(seu)
log_msg("QC filter: ", qc_after, " of ", qc_before, " cells retained")

qc_plot <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    group.by = "Tissue", pt.size = 0, ncol = 3)
ragg::agg_png("results/figures/03_singlecell_qc.png", width = 12, height = 4, units = "in", res = 300)
print(qc_plot)
dev.off()
log_msg("Saved QC violin plot: results/figures/03_singlecell_qc.png")
# --- 4. Normalize, cluster, UMAP ---------------------------------------------
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

log_msg("Clustering complete: ", length(unique(Idents(seu))), " clusters identified")
log_msg(paste(capture.output(print(table(Idents(seu)))), collapse = "\n"))
# --- 5. Broad-category annotation via canonical markers ---------------------
marker_sets <- list(
  Progenitor  = c("PDGFRA", "DCN", "LUM"),
  Immune      = c("PTPRC"),
  Endothelial = c("PECAM1", "VWF"),
  Mesothelial = c("MSLN", "UPK3B", "KRT19")
)
marker_sets <- lapply(marker_sets, function(g) intersect(g, rownames(seu)))
log_msg("Marker genes actually present in the dataset:")
for (ct in names(marker_sets)) log_msg("  ", ct, ": ", paste(marker_sets[[ct]], collapse = ", "))

for (ct in names(marker_sets)) {
  seu <- AddModuleScore(seu, features = list(marker_sets[[ct]]), name = paste0(ct, "_score"), verbose = FALSE)
}
score_cols <- paste0(names(marker_sets), "_score1")

cluster_scores <- sapply(score_cols, function(col) {
  tapply(seu[[col]][, 1], Idents(seu), mean)
})
colnames(cluster_scores) <- names(marker_sets)
log_msg("Mean module score per cluster:")
log_msg(paste(capture.output(print(round(cluster_scores, 3))), collapse = "\n"))
# Only assign a category when its score is the max AND genuinely positive
# (>0.1, a modest floor above zero) -- otherwise "Unclassified". This is the
# direct fix for Checkpoint 9's finding: argmax over all-negative scores was
# manufacturing confident-looking but meaningless labels.
score_floor <- 0.1
max_score <- apply(cluster_scores, 1, max)
max_cat   <- colnames(cluster_scores)[apply(cluster_scores, 1, which.max)]
cluster_call <- ifelse(max_score > score_floor, max_cat, "Unclassified")
names(cluster_call) <- rownames(cluster_scores)
log_msg("Cluster -> broad cell type assignment (score floor = ", score_floor, "):")
for (cl in names(cluster_call)) log_msg("  Cluster ", cl, " -> ", cluster_call[cl], " (max score = ", round(max_score[cl], 3), ")")
ct_vec <- cluster_call[as.character(Idents(seu))]
names(ct_vec) <- colnames(seu)
seu$cell_type <- ct_vec

log_msg("Final cell-type distribution:")
log_msg(paste(capture.output(print(table(seu$cell_type))), collapse = "\n"))
log_msg("Final cell-type distribution (%):")
log_msg(paste(capture.output(print(round(100 * table(seu$cell_type) / ncol(seu), 1))), collapse = "\n"))
log_msg("Cell-type x Tissue cross-tab (checking the paper's depot-specificity claim for mesothelial cells):")
log_msg(paste(capture.output(print(table(seu$cell_type, seu$Tissue))), collapse = "\n"))
# --- 6. UMAP figure -----------------------------------------------------------
umap_plot <- DimPlot(seu, group.by = "cell_type", label = TRUE, repel = TRUE) +
  ggtitle("GSE129363 adipose SVF: broad cell-type clusters") +
  theme_minimal(base_size = 11)

ragg::agg_png("results/figures/04_singlecell_umap.png", width = 7, height = 6, units = "in", res = 300)
print(umap_plot)
dev.off()
log_msg("Saved UMAP plot: results/figures/04_singlecell_umap.png")

umap_depot <- DimPlot(seu, group.by = "Tissue") +
  ggtitle("GSE129363 adipose SVF: by depot (SAT/VAT)") +
  theme_minimal(base_size = 11)
ragg::agg_png("results/figures/04b_singlecell_umap_by_depot.png", width = 7, height = 6, units = "in", res = 300)
print(umap_depot)
dev.off()
log_msg("Saved depot-colored UMAP: results/figures/04b_singlecell_umap_by_depot.png")
# --- 7. Save processed object -------------------------------------------------
saveRDS(seu, "data/processed/seurat_annotated.rds", compress = TRUE)
log_msg("Saved annotated Seurat object: data/processed/seurat_annotated.rds")

writeLines(log_lines, "results/logs/06_singlecell_seurat_summary.txt")
log_msg("Summary written: results/logs/06_singlecell_seurat_summary.txt")
