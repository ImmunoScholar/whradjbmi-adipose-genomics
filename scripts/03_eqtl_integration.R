# 03_eqtl_integration.R
# For each lead locus (hg19), liftOver the cis-window to GRCh38 and retrieve
# adipose subcutaneous cis-eQTL associations (GTEx v8, eQTL Catalogue
# QTD000116) via remote tabix range queries (Rsamtools::scanTabix -- the
# system tabix binary on this machine lacks libcurl/https support, confirmed
# by a failed smoke test; Rsamtools bundles its own htslib build that does
# support remote HTTPS access).
#
# Scope: retrieval + liftOver only. GWAS/eQTL allele harmonisation and the
# coloc-ready merged table are built in 04_harmonize_coloc.R.
library(data.table)
library(yaml)
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
cfg <- yaml::read_yaml("config/config.yml")

dir.create("data/interim", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== eQTL retrieval (liftOver + Rsamtools tabix) ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
loci <- fread("results/tables/02_lead_loci.csv")
log_msg("Loci to process: ", nrow(loci))
# --- 1. liftOver: hg19 cis-windows -> GRCh38 --------------------------------
chain_path <- cfg$liftover$local_path
if (!file.exists(chain_path)) {
  stop("liftOver chain not found at ", chain_path, " -- expected the uncompressed .chain file.")
}
chain <- import.chain(chain_path)

window_bp <- cfg$eqtl$cis_window_kb * 1000
loci[, `:=`(
  win_start_hg19 = pmax(1, POS - window_bp),
  win_end_hg19   = POS + window_bp
)]

gr_hg19 <- GRanges(
  seqnames = paste0("chr", loci$CHR),
  ranges   = IRanges(start = loci$win_start_hg19, end = loci$win_end_hg19),
  locus_id = loci$locus_id
)

lifted <- liftOver(gr_hg19, chain)
liftover_diag <- data.table(
  locus_id       = loci$locus_id,
  chr_hg19       = loci$CHR,
  win_start_hg19 = loci$win_start_hg19,
  win_end_hg19   = loci$win_end_hg19,
  width_hg19     = loci$win_end_hg19 - loci$win_start_hg19 + 1,
  n_fragments    = elementNROWS(lifted)
)

get_span <- function(gr_list_elt) {
  if (length(gr_list_elt) == 0) {
    return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_, n_chr = 0L))
  }
  chr_tab <- table(as.character(seqnames(gr_list_elt)))
  modal_chr <- names(chr_tab)[which.max(chr_tab)]
  sub <- gr_list_elt[as.character(seqnames(gr_list_elt)) == modal_chr]
  list(chr = modal_chr, start = min(start(sub)), end = max(end(sub)), n_chr = length(chr_tab))
}
spans <- lapply(lifted, get_span)
liftover_diag[, `:=`(
  chr_hg38     = sapply(spans, `[[`, "chr"),
  start_hg38   = sapply(spans, `[[`, "start"),
  end_hg38     = sapply(spans, `[[`, "end"),
  n_target_chr = sapply(spans, `[[`, "n_chr")
)]
liftover_diag[, width_hg38 := end_hg38 - start_hg38 + 1]
liftover_diag[, pct_width_retained := round(100 * width_hg38 / width_hg19, 1)]

for (i in seq_len(nrow(liftover_diag))) {
  row <- liftover_diag[i]
  log_msg(row$locus_id, ": hg19 chr", row$chr_hg19, ":", row$win_start_hg19, "-", row$win_end_hg19,
          " (", row$n_fragments, " fragments) -> hg38 ", row$chr_hg38, ":", row$start_hg38,
          "-", row$end_hg38, " (", row$pct_width_retained, "% width retained",
          if (row$n_target_chr > 1) paste0(", WARNING: mapped to ", row$n_target_chr, " different chromosomes") else "",
          ")")
}

fwrite(liftover_diag, "results/tables/03_liftover_diagnostics.csv")
log_msg("Saved liftOver diagnostics: results/tables/03_liftover_diagnostics.csv")
poor_liftover <- liftover_diag[is.na(pct_width_retained) | pct_width_retained < 90 | n_target_chr > 1]
if (nrow(poor_liftover) > 0) {
  log_msg("WARNING: ", nrow(poor_liftover), " locus/loci have <90% width retention or split across chromosomes:")
log_msg(paste(capture.output(print(poor_liftover)), collapse = "\n"))
} else {
  log_msg("All loci lifted cleanly (>=90% width retained, single target chromosome).")
}
# --- 2. Validate the eQTL file's column header ------------------------------
eqtl_url <- paste0(cfg$eqtl$ftp_base, "/", cfg$eqtl$study_id, "/", cfg$eqtl$dataset_id,
                    "/", cfg$eqtl$dataset_id, ".all.tsv.gz")
expected_header <- c("molecular_trait_id", "chromosome", "position", "ref", "alt",
                      "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type",
                      "ac", "an", "r2", "molecular_trait_object_id", "gene_id",
                      "median_tpm", "rsid")

header_tmp <- tempfile(fileext = ".tsv.gz")
download.file(eqtl_url, destfile = header_tmp, mode = "wb", quiet = TRUE,
              headers = c(Range = "bytes=0-65535"))
header_line <- tryCatch({
  con <- gzfile(header_tmp, "rt")
  line1 <- readLines(con, n = 1)
  close(con)
  line1
}, error = function(e) NA_character_)
if (is.na(header_line)) {
  stop("Could not read header from remote eQTL file -- check connectivity/URL: ", eqtl_url)
}
actual_header <- strsplit(header_line, "\t")[[1]]

if (!identical(actual_header, expected_header)) {
  log_msg("EXPECTED header: ", paste(expected_header, collapse = ", "))
  log_msg("ACTUAL header:   ", paste(actual_header, collapse = ", "))
  stop("eQTL Catalogue column schema does not match what this script expects.")
}
log_msg("eQTL file header verified: matches expected 19-column schema.")

# --- 3. Remote tabix query per locus via Rsamtools --------------------------
tbx <- TabixFile(eqtl_url)
open(tbx)

query_locus <- function(chr, start, end) {
  chr_bare <- sub("^chr", "", chr)
  gr <- GRanges(chr_bare, IRanges(start, end))
  res <- tryCatch(scanTabix(tbx, param = gr),
                   error = function(e) list(list()))
  lines <- res[[1]]
  if (is.null(lines)) lines <- character(0)
  list(ok = TRUE, lines = lines)
}

eqtl_list <- vector("list", nrow(liftover_diag))
for (i in seq_len(nrow(liftover_diag))) {
  row <- liftover_diag[i]
  if (is.na(row$chr_hg38)) {
    log_msg(row$locus_id, ": SKIPPED -- liftOver produced no mapped region.")
    next
  }
  res <- query_locus(row$chr_hg38, row$start_hg38, row$end_hg38)
  if (length(res$lines) == 0) {
    log_msg(row$locus_id, ": query succeeded, 0 variants returned (no eQTL coverage in this region for this dataset).")
    next
  }
  dt <- fread(text = paste(res$lines, collapse = "\n"), header = FALSE, col.names = expected_header)
  dt[, locus_id := row$locus_id]
  eqtl_list[[i]] <- dt
  log_msg(row$locus_id, ": ", nrow(dt), " variant-gene association rows retrieved, ",
          length(unique(dt$gene_id)), " unique genes")
}
close(tbx)
eqtl_raw <- rbindlist(eqtl_list, use.names = TRUE, fill = TRUE)
log_msg("Total eQTL association rows across all loci: ", format(nrow(eqtl_raw), big.mark = ","))
log_msg("Loci with zero returned variants: ", sum(sapply(eqtl_list, is.null)), " of ", nrow(liftover_diag))

fwrite(eqtl_raw, "data/interim/03_eqtl_raw.csv")
log_msg("Saved raw eQTL associations: data/interim/03_eqtl_raw.csv")

writeLines(log_lines, "results/logs/03_eqtl_retrieval_summary.txt")
log_msg("Summary written: results/logs/03_eqtl_retrieval_summary.txt")
