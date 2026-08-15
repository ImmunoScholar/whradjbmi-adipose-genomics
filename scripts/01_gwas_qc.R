# 01_gwas_qc.R
# Level 1/2 quality control on the Pulit et al. 2019 WHRadjBMI GWAS summary
# statistics (GIANT+UKBB meta-analysis, hg19). Produces a cleaned processed
# table plus a QC summary log and diagnostic plots. Does NOT define loci --
# that is 02_locus_selection.R, kept separate so each gets its own checkpoint.
library(data.table)
library(yaml)
library(ggplot2)
library(patchwork)
library(ragg)
cfg <- yaml::read_yaml("config/config.yml")
gwas_path <- cfg$gwas$local_path
if (!file.exists(gwas_path)) {
  stop("GWAS file not found at ", gwas_path,
       " -- has the download finished? Run gzip -t on it first.")
}
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

log_msg("=== GWAS QC: ", gwas_path, " ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
# --- 1. Read ------------------------------------------------------------
# Explicit colClasses keeps memory use predictable across the ~27M rows.
expected_cols <- c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele",
                    "Freq_Tested_Allele", "BETA", "SE", "P", "N", "INFO")
gwas <- fread(gwas_path, showProgress = TRUE)

log_msg("Rows read: ", format(nrow(gwas), big.mark = ","))
log_msg("Columns present: ", paste(names(gwas), collapse = ", "))
missing_cols <- setdiff(expected_cols, names(gwas))
extra_cols   <- setdiff(names(gwas), expected_cols)
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}
if (length(extra_cols) > 0) {
  log_msg("NOTE: unexpected extra columns present: ", paste(extra_cols, collapse = ", "))
}
# --- 2. Missingness -------------------------------------------------------
na_counts <- sapply(gwas[, ..expected_cols], function(x) sum(is.na(x)))
log_msg("Missing values per column:")
for (nm in names(na_counts)) log_msg("  ", nm, ": ", na_counts[nm])
# --- 3. Duplicate variants -------------------------------------------------
dup_snp <- sum(duplicated(gwas$SNP))
log_msg("Duplicated SNP entries: ", dup_snp)

# --- 4. Chromosome representation -----------------------------------------
chr_tab <- table(gwas$CHR)
log_msg("Chromosomes present: ", paste(names(chr_tab), collapse = ", "))
unexpected_chr <- setdiff(names(chr_tab), as.character(1:22))
if (length(unexpected_chr) > 0) {
  log_msg("NOTE: non-autosomal or unexpected CHR values: ", paste(unexpected_chr, collapse = ", "))
}
# --- 5. Position sanity -----------------------------------------------------
bad_pos <- sum(gwas$POS <= 0 | !is.finite(gwas$POS))
log_msg("Rows with non-positive/non-finite POS: ", bad_pos)
# --- 6. Allele sanity: single-character ACGT for Tested/Other_Allele ------
valid_bases <- c("A", "C", "G", "T")
bad_tested <- sum(!(toupper(gwas$Tested_Allele) %in% valid_bases))
bad_other  <- sum(!(toupper(gwas$Other_Allele)  %in% valid_bases))
log_msg("Tested_Allele values outside {A,C,G,T} (indels expected here): ", bad_tested)
log_msg("Other_Allele values outside {A,C,G,T} (indels expected here): ", bad_other)
# Cross-check: alleles embedded in the SNP column (rsID:A1:A2) should match
# Tested_Allele/Other_Allele in some order. Mismatches would indicate a
# column-alignment problem upstream -- this is the single most important
# internal-consistency check for this file.
snp_parts <- tstrsplit(gwas$SNP, ":", fixed = TRUE)
snp_a1 <- toupper(snp_parts[[2]])
snp_a2 <- toupper(snp_parts[[3]])
allele_match <- (snp_a1 == toupper(gwas$Tested_Allele) & snp_a2 == toupper(gwas$Other_Allele)) |
                (snp_a2 == toupper(gwas$Tested_Allele) & snp_a1 == toupper(gwas$Other_Allele))
n_allele_mismatch <- sum(!allele_match, na.rm = TRUE)
log_msg("SNP-column alleles vs Tested/Other_Allele mismatches: ", n_allele_mismatch,
        " (", round(100 * n_allele_mismatch / nrow(gwas), 3), "% of rows)")
# --- 7. Value-range sanity -------------------------------------------------
bad_p <- sum(gwas$P < 0 | gwas$P > 1 | !is.finite(gwas$P))
log_msg("Rows with P outside [0,1] or non-finite: ", bad_p)
bad_freq <- sum(gwas$Freq_Tested_Allele < 0 | gwas$Freq_Tested_Allele > 1, na.rm = TRUE)
log_msg("Rows with Freq_Tested_Allele outside [0,1]: ", bad_freq)
bad_beta_se <- sum(!is.finite(gwas$BETA) | !is.finite(gwas$SE) | gwas$SE <= 0)
log_msg("Rows with non-finite BETA or non-positive/non-finite SE: ", bad_beta_se)
log_msg("N range: ", min(gwas$N, na.rm = TRUE), " to ", max(gwas$N, na.rm = TRUE),
        " (paper reports max N = 694,649)")
log_msg("INFO range: ", round(min(gwas$INFO, na.rm = TRUE), 3), " to ",
        round(max(gwas$INFO, na.rm = TRUE), 3))
# --- 8. Genome-wide significant count (sanity, not final locus selection) --
n_gws <- sum(gwas$P < cfg$gwas$pval_threshold, na.rm = TRUE)
log_msg("Variants at P < ", cfg$gwas$pval_threshold, ": ", format(n_gws, big.mark = ","))
# --- 9. Diagnostic plots ----------------------------------------------------
p_hist <- ggplot(gwas[sample(.N, min(.N, 2e6))], aes(P)) +
  geom_histogram(bins = 50, fill = "grey40") +
  labs(title = "P-value distribution (2M-row subsample)", x = "P", y = "Count") +
  theme_minimal(base_size = 11)

info_hist <- ggplot(gwas[sample(.N, min(.N, 2e6))], aes(INFO)) +
  geom_histogram(bins = 50, fill = "grey40") +
  labs(title = "Imputation INFO score", x = "INFO", y = "Count") +
  theme_minimal(base_size = 11)
freq_hist <- ggplot(gwas[sample(.N, min(.N, 2e6))], aes(Freq_Tested_Allele)) +
  geom_histogram(bins = 50, fill = "grey40") +
  labs(title = "Effect-allele frequency", x = "Freq_Tested_Allele", y = "Count") +
  theme_minimal(base_size = 11)
qc_panel <- p_hist + info_hist + freq_hist
ragg::agg_png("results/figures/00_gwas_qc_diagnostics.png",
              width = 12, height = 4, units = "in", res = 300)
print(qc_panel)
dev.off()
log_msg("Saved diagnostic plot: results/figures/00_gwas_qc_diagnostics.png")
# --- 10. Save cleaned/processed table --------------------------------------
# QC is diagnostic here, not filtering -- we flag problems rather than
# silently drop rows. Locus selection (02) applies its own explicit filters.
saveRDS(gwas, "data/processed/gwas_qc.rds", compress = TRUE)
log_msg("Saved processed table: data/processed/gwas_qc.rds")
writeLines(log_lines, "results/logs/01_gwas_qc_summary.txt")
log_msg("QC summary written: results/logs/01_gwas_qc_summary.txt")
log_msg("=== sessionInfo ===")
log_msg(paste(capture.output(sessionInfo()), collapse = "\n"))
writeLines(log_lines, "results/logs/01_gwas_qc_summary.txt")
