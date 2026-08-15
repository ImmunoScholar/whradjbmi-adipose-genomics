# 05_load_singlecell.R
# Loads the GSE129363 human adipose SVF single-cell matrix and inspects its
# real structure before assuming anything about column meaning. Deliberately
# stops short of building the final annotated Seurat object -- per Checkpoint
# 8, we don't assume a column named "cell_type" (or similar) is correct
# without seeing the actual values and cross-checking against the dataset's
# documentation.
library(Matrix)
library(data.table)
library(yaml)
cfg <- yaml::read_yaml("config/config.yml")
sc_dir <- cfg$single_cell$matrix_dir
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)
log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Single-cell data inspection: GSE129363 ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
files <- list.files(sc_dir, full.names = TRUE)
log_msg("Files found in ", sc_dir, ":")
for (f in files) log_msg("  ", basename(f), " (", format(file.info(f)$size / 1e6, digits = 3), " MB)")
# --- 1. Load the raw count matrix -------------------------------------------
mtx_path   <- file.path(sc_dir, "GSE129363_Discovery_Cohort_matrix.mtx.gz")
cells_path <- file.path(sc_dir, "GSE129363_Discovery_Cohort_barcodes.tsv.gz")
genes_path <- file.path(sc_dir, "GSE129363_Discovery_Cohort_genes.tsv.gz")
anno_path  <- file.path(sc_dir, "GSE129363_Discovery_Cohort_CellAnnotation.txt.gz")

for (p in c(mtx_path, cells_path, genes_path, anno_path)) {
  if (!file.exists(p)) stop("Expected file not found: ", p)
}
mat <- Matrix::readMM(mtx_path)
log_msg("Raw matrix dimensions (as stored): ", nrow(mat), " x ", ncol(mat))
barcodes <- fread(cells_path, header = FALSE)
genes    <- fread(genes_path, header = FALSE)
log_msg("Barcodes file: ", nrow(barcodes), " rows, ", ncol(barcodes), " columns")
log_msg("Genes file: ", nrow(genes), " rows, ", ncol(genes), " columns")
log_msg("Genes file first 3 rows:")
log_msg(paste(capture.output(print(head(genes, 3))), collapse = "\n"))
# Determine matrix orientation: genes-x-cells (10x convention) vs cells-x-genes.
if (nrow(mat) == nrow(genes) && ncol(mat) == nrow(barcodes)) {
  orientation <- "genes_x_cells (standard 10x convention)"
} else if (nrow(mat) == nrow(barcodes) && ncol(mat) == nrow(genes)) {
  orientation <- "cells_x_genes (transposed relative to 10x convention)"
} else {
  orientation <- "MISMATCH -- matrix dimensions do not match either barcodes or genes file length"
}
log_msg("Matrix orientation check: ", orientation)
# --- 2. Load and inspect the annotation file --------------------------------
# Do NOT assume column names here -- print everything and stop.
anno <- fread(anno_path)
log_msg("Annotation file dimensions: ", nrow(anno), " rows, ", ncol(anno), " columns")
log_msg("Annotation column names: ", paste(names(anno), collapse = ", "))
log_msg("Annotation file first 5 rows:")
log_msg(paste(capture.output(print(head(anno, 5))), collapse = "\n"))

# For every character/factor-like column with a manageable number of unique
# values, print the value counts -- this is how we'll actually identify which
# column is cell type, which is depot, which is sample/donor, etc.
for (col in names(anno)) {
  vals <- anno[[col]]
  if (is.character(vals) || is.factor(vals)) {
    n_unique <- length(unique(vals))
    if (n_unique <= 40) {
      log_msg("Column '", col, "' (", n_unique, " unique values):")
      log_msg(paste(capture.output(print(table(vals, useNA = "ifany"))), collapse = "\n"))
    } else {
      log_msg("Column '", col, "' has ", n_unique, " unique values (too many to print, likely an ID column) -- first 5: ",
              paste(head(unique(vals), 5), collapse = ", "))
    }
  } else {
    log_msg("Column '", col, "' is numeric, range: ", min(vals, na.rm = TRUE), " to ", max(vals, na.rm = TRUE))
  }
}
# --- 3. Barcode overlap check between annotation and matrix -----------------
# Try matching on whichever annotation column looks most like a barcode
# (highest overlap with the matrix's own barcode list) -- report the overlap,
# do not assume which column it is in advance.
matrix_barcodes <- barcodes$V1
best_overlap <- data.table(column = character(0), n_overlap = integer(0), pct_of_matrix = numeric(0))
for (col in names(anno)) {
  vals <- as.character(anno[[col]])
  ov <- length(intersect(vals, matrix_barcodes))
  if (ov > 0) {
    best_overlap <- rbind(best_overlap, data.table(
      column = col, n_overlap = ov,
      pct_of_matrix = round(100 * ov / length(matrix_barcodes), 1)
    ))
  }
}
log_msg("Annotation columns with nonzero overlap against matrix barcodes:")
log_msg(paste(capture.output(print(best_overlap[order(-n_overlap)])), collapse = "\n"))
writeLines(log_lines, "results/logs/05_singlecell_inspection.txt")
log_msg("Inspection log written: results/logs/05_singlecell_inspection.txt")
log_msg("STOPPING HERE deliberately -- next script builds the Seurat object")
log_msg("using whichever columns this inspection identifies as barcode/cell-type/depot,")
log_msg("not assumed column names.")
