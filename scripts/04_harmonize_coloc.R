# 04_harmonize_coloc.R
# Harmonizes GWAS (hg19, Tested/Other_Allele) and eQTL Catalogue (GRCh38,
# ref/alt) summary statistics via shared rsID, explicitly handles ambiguous
# palindromic SNPs and allele-flip alignment, then runs coloc.abf for every
# gene present at the 3 most significant loci (config$coloc$n_loci).
#
# Convention: eQTL Catalogue 'alt' is the effect allele (beta direction).
# GWAS Tested_Allele is the effect allele. We align GWAS beta to be relative
# to the eQTL alt allele, flipping sign and frequency where needed.
library(data.table)
library(yaml)
library(coloc)
library(ggplot2)
library(ragg)
cfg <- yaml::read_yaml("config/config.yml")

dir.create("data/interim", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Allele harmonization + colocalization ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
gwas <- readRDS("data/processed/gwas_qc.rds")
liftover_diag <- fread("results/tables/03_liftover_diagnostics.csv")
eqtl <- fread("data/interim/03_eqtl_raw.csv")
loci <- fread("results/tables/02_lead_loci.csv")
gwas[, rsid := tstrsplit(SNP, ":", fixed = TRUE, keep = 1L)]
# --- 1. Subset GWAS to each locus's window, tag with locus_id --------------
gwas_locus_list <- vector("list", nrow(liftover_diag))
for (i in seq_len(nrow(liftover_diag))) {
  row <- liftover_diag[i]
  sub <- gwas[CHR == row$chr_hg19 & POS >= row$win_start_hg19 & POS <= row$win_end_hg19]
  sub[, locus_id := row$locus_id]
  gwas_locus_list[[i]] <- sub
}
gwas_windows <- rbindlist(gwas_locus_list)
log_msg("GWAS variants across all 8 locus windows: ", format(nrow(gwas_windows), big.mark = ","))
# --- 2. Merge on (locus_id, rsid) -------------------------------------------
merged <- merge(
  gwas_windows[, .(locus_id, rsid, CHR, POS, Tested_Allele, Other_Allele,
                    Freq_Tested_Allele, BETA_gwas = BETA, SE_gwas = SE, P_gwas = P, N_gwas = N)],
  eqtl[, .(locus_id, rsid, gene_id, ref, alt, maf, BETA_eqtl = beta, SE_eqtl = se,
           P_eqtl = pvalue, ma_samples)],
  by = c("locus_id", "rsid"),
  allow.cartesian = TRUE  # one GWAS row can match multiple genes at the same variant
)
log_msg("Rows after GWAS-eQTL merge on (locus_id, rsid): ", format(nrow(merged), big.mark = ","))
# --- 3. Duplicate check ------------------------------------------------------
dup_key <- merged[, .N, by = .(locus_id, gene_id, rsid)][N > 1]
n_dup_groups <- nrow(dup_key)
if (n_dup_groups > 0) {
  log_msg("WARNING: ", n_dup_groups, " (locus_id, gene_id, rsid) combinations have duplicate rows ",
          "(likely from overlapping liftOver fragments) -- keeping first occurrence only.")
}
setkey(merged, locus_id, gene_id, rsid)
merged <- unique(merged, by = c("locus_id", "gene_id", "rsid"))
log_msg("Rows after deduplication: ", format(nrow(merged), big.mark = ","))
# --- 4. Allele classification: palindromic / direct match / flipped / mismatch
is_palindromic <- function(a1, a2) {
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
  (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}
merged[, palindromic := is_palindromic(toupper(Tested_Allele), toupper(Other_Allele))]
merged[, direct_match  := toupper(Tested_Allele) == toupper(alt) & toupper(Other_Allele) == toupper(ref)]
merged[, flipped_match := toupper(Tested_Allele) == toupper(ref) & toupper(Other_Allele) == toupper(alt)]
merged[, allele_status := fifelse(palindromic, "palindromic_excluded",
                            fifelse(direct_match, "direct_match",
                            fifelse(flipped_match, "flipped_match", "mismatch_excluded")))]
status_tab <- merged[, .N, by = allele_status]
log_msg("Allele harmonization outcome across all merged rows:")
for (i in seq_len(nrow(status_tab))) {
  log_msg("  ", status_tab$allele_status[i], ": ", format(status_tab$N[i], big.mark = ","))
}
# Align GWAS beta/frequency to the eQTL alt allele for matched rows only.
harmonized <- merged[allele_status %in% c("direct_match", "flipped_match")]
harmonized[, `:=`(
  BETA_gwas_aligned = fifelse(allele_status == "direct_match", BETA_gwas, -BETA_gwas),
  Freq_aligned       = fifelse(allele_status == "direct_match", Freq_Tested_Allele, 1 - Freq_Tested_Allele)
)]

log_msg("Rows retained for coloc after excluding palindromic/mismatched alleles: ",
        format(nrow(harmonized), big.mark = ","),
        " (", round(100 * nrow(harmonized) / nrow(merged), 1), "% of merged rows)")
fwrite(harmonized, "data/interim/04_harmonized_locus_eqtl.csv")
log_msg("Saved harmonized table: data/interim/04_harmonized_locus_eqtl.csv")
# --- 5. Per-locus, per-gene coloc.abf ---------------------------------------
coloc_loci <- loci[order(P)][seq_len(cfg$coloc$n_loci), locus_id]
log_msg("Loci selected for colocalization (", cfg$coloc$n_loci, " most significant): ",
        paste(coloc_loci, collapse = ", "))

min_snps <- 30
results_list <- list()
for (lid in coloc_loci) {
  locus_data <- harmonized[locus_id == lid]
  genes <- unique(locus_data$gene_id)
  log_msg(lid, ": ", length(genes), " genes present, testing each (min ", min_snps, " SNPs required)")

for (g in genes) {
    gd <- locus_data[gene_id == g]
    gd <- gd[!is.na(BETA_gwas_aligned) & !is.na(SE_gwas) & SE_gwas > 0 &
             !is.na(BETA_eqtl) & !is.na(SE_eqtl) & SE_eqtl > 0 &
             !is.na(Freq_aligned) & Freq_aligned > 0 & Freq_aligned < 1 &
             !is.na(maf) & maf > 0 & maf < 1]

    if (nrow(gd) < min_snps) {
      log_msg("  ", g, ": SKIPPED -- only ", nrow(gd), " usable SNPs (< ", min_snps, ")")
      next
    }
dataset1 <- list(beta = gd$BETA_gwas_aligned, varbeta = gd$SE_gwas^2,
                      N = median(gd$N_gwas), MAF = gd$Freq_aligned,
                      type = "quant", snp = gd$rsid)
    dataset2 <- list(beta = gd$BETA_eqtl, varbeta = gd$SE_eqtl^2,
                      N = cfg$eqtl$sample_size, MAF = gd$maf,
                      type = "quant", snp = gd$rsid)

    res <- tryCatch(
      suppressWarnings(coloc.abf(dataset1 = dataset1, dataset2 = dataset2,
                                  p1 = cfg$coloc$p1, p2 = cfg$coloc$p2, p12 = cfg$coloc$p12)),
      error = function(e) NULL
    )

    if (is.null(res)) {
      log_msg("  ", g, ": coloc.abf FAILED on ", nrow(gd), " SNPs -- skipped")
      next
    }
s <- res$summary
    log_msg("  ", g, ": nsnps=", s["nsnps"], " PP0=", round(s["PP.H0.abf"], 3),
            " PP1=", round(s["PP.H1.abf"], 3), " PP2=", round(s["PP.H2.abf"], 3),
            " PP3=", round(s["PP.H3.abf"], 3), " PP4=", round(s["PP.H4.abf"], 3))

    results_list[[paste(lid, g)]] <- data.table(
      locus_id = lid, gene_id = g, nsnps = s["nsnps"],
      PP0 = s["PP.H0.abf"], PP1 = s["PP.H1.abf"], PP2 = s["PP.H2.abf"],
      PP3 = s["PP.H3.abf"], PP4 = s["PP.H4.abf"]
    )
  }
}

coloc_results <- rbindlist(results_list)
setorder(coloc_results, locus_id, -PP4)
fwrite(coloc_results, "results/tables/04_coloc_results.csv")
log_msg("Saved coloc results: results/tables/04_coloc_results.csv (",
        nrow(coloc_results), " locus-gene pairs tested successfully)")

log_msg("Top PP4 gene per locus:")
top_per_locus <- coloc_results[, .SD[which.max(PP4)], by = locus_id]
print(top_per_locus)
log_msg(paste(capture.output(print(top_per_locus)), collapse = "\n"))
# --- 6. Diagnostic plot for the single best locus-gene coloc hit -----------
if (nrow(coloc_results) > 0) {
  best <- coloc_results[which.max(PP4)]
  gd <- harmonized[locus_id == best$locus_id & gene_id == best$gene_id]

p_scatter <- ggplot(gd, aes(x = -log10(P_gwas), y = -log10(P_eqtl))) +
    geom_point(alpha = 0.5, size = 1.2, color = "grey30") +
    labs(title = paste0(best$locus_id, " x ", best$gene_id,
                         " (PP4 = ", round(best$PP4, 3), ", n = ", best$nsnps, " SNPs)"),
         subtitle = "Each point is one variant; visual overlap of strong signals should track high PP4",
         x = expression(-log[10]("GWAS P")), y = expression(-log[10]("eQTL P"))) +
    theme_minimal(base_size = 11)
ragg::agg_png("results/figures/02_coloc_best_locus_diagnostic.png",
                width = 6, height = 5, units = "in", res = 300)
  print(p_scatter)
  dev.off()
  log_msg("Saved diagnostic plot for best coloc hit: results/figures/02_coloc_best_locus_diagnostic.png")
}

writeLines(log_lines, "results/logs/04_harmonize_coloc_summary.txt")
log_msg("Summary written: results/logs/04_harmonize_coloc_summary.txt")
