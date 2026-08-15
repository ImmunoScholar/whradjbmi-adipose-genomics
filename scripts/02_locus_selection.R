# 02_locus_selection.R
# Distance-based greedy locus definition on QC'd WHRadjBMI GWAS.
# NOT true LD-clumping -- no LD reference panel is used. This is an explicit
# approximation: iteratively take the most significant remaining variant as
# a lead SNP, exclude everything within +/- locus_window_kb of it, repeat.
# coloc.abf (stage 4) does not require an LD panel, so this approximation is
# acceptable for locus definition but must be labelled as an approximation
# everywhere it appears (README, figure captions, methods write-up).
library(data.table)
library(yaml)
library(ggplot2)
library(ggrepel)
library(ragg)

cfg <- yaml::read_yaml("config/config.yml")
gwas <- readRDS("data/processed/gwas_qc.rds")
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

log_msg("=== Locus selection ===")
log_msg("Run timestamp: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
log_msg("Input rows: ", format(nrow(gwas), big.mark = ","))
# --- 1. Drop rows with missing CHR/POS/N ------------------------------------
# Confirmed at Checkpoint 2: 706 rows lack coordinates (different SNP-ID
# format, no embedded alleles). None reach P < 5e-8 (min P = 1.573e-05), so
# dropping them loses no genome-wide-significant signal. Logged explicitly.
n_before <- nrow(gwas)
gwas <- gwas[!is.na(CHR) & !is.na(POS) & !is.na(N)]
n_dropped <- n_before - nrow(gwas)
log_msg("Dropped rows missing CHR/POS/N: ", n_dropped,
        " (confirmed at QC: min P among these = 1.573e-05, none genome-wide significant)")
log_msg("Rows remaining: ", format(nrow(gwas), big.mark = ","))
# --- 2. Restrict to genome-wide significant variants ------------------------
sig <- gwas[P < cfg$gwas$pval_threshold]
setorder(sig, P)
log_msg("Variants at P < ", cfg$gwas$pval_threshold, ": ", format(nrow(sig), big.mark = ","))
# --- 3. Greedy distance-based pruning ---------------------------------------
window_bp <- cfg$gwas$locus_window_kb * 1000
loci <- data.table()
remaining <- copy(sig)
iter <- 0
while (nrow(remaining) > 0 && nrow(loci) < cfg$gwas$n_top_loci) {
  iter <- iter + 1
  lead <- remaining[1]
  loci <- rbind(loci, lead)
  remaining <- remaining[!(CHR == lead$CHR & abs(POS - lead$POS) <= window_bp)]
  log_msg("Locus ", iter, ": lead ", lead$SNP, " (chr", lead$CHR, ":", lead$POS,
          ", P=", format(lead$P, scientific = TRUE, digits = 3), ") -- ",
          nrow(remaining), " significant variants remain after exclusion")
}

log_msg("Independent loci selected: ", nrow(loci), " (target was ", cfg$gwas$n_top_loci, ")")
# --- 4. Save lead-locus table -------------------------------------------------
loci[, locus_id := paste0("locus_", sprintf("%02d", .I))]
fwrite(loci, "results/tables/02_lead_loci.csv")
log_msg("Saved lead-locus table: results/tables/02_lead_loci.csv")
print(loci[, .(locus_id, CHR, POS, SNP, Tested_Allele, Other_Allele, BETA, SE, P, N)])
# --- 5. Manhattan plot --------------------------------------------------------
set.seed(cfg$seed)
plot_data <- rbind(
  gwas[P >= cfg$gwas$pval_threshold][sample(.N, min(.N, 3e5))],
  gwas[P < cfg$gwas$pval_threshold]
)
plot_data[, CHR := as.integer(CHR)]
setorder(plot_data, CHR, POS)
chr_max <- plot_data[, max(POS), by = CHR]
setorder(chr_max, CHR)
chr_offset <- c(0, cumsum(as.numeric(chr_max$V1))[-nrow(chr_max)])
names(chr_offset) <- chr_max$CHR
plot_data[, pos_cum := POS + chr_offset[as.character(CHR)]]
axis_df <- plot_data[, .(center = (min(pos_cum) + max(pos_cum)) / 2), by = CHR]

lead_plot <- plot_data[SNP %in% loci$SNP]
manhattan <- ggplot(plot_data, aes(x = pos_cum, y = -log10(P), color = factor(CHR %% 2))) +
  geom_point(size = 0.4, alpha = 0.6) +
  geom_hline(yintercept = -log10(cfg$gwas$pval_threshold), linetype = "dashed", color = "firebrick") +
  geom_point(data = lead_plot, aes(x = pos_cum, y = -log10(P)),
             inherit.aes = FALSE, color = "black", size = 1.6) +
  ggrepel::geom_text_repel(data = lead_plot, aes(x = pos_cum, y = -log10(P), label = SNP),
                            inherit.aes = FALSE, size = 2.8, max.overlaps = 20) +
  scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$center) +
  scale_color_manual(values = c("grey60", "grey35"), guide = "none") +
  labs(title = "WHRadjBMI GWAS (Pulit et al. 2019)",
       subtitle = paste0(cfg$gwas$n_top_loci, " lead loci highlighted, distance-based pruning (",
                          cfg$gwas$locus_window_kb, " kb window)"),
       x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ragg::agg_png("results/figures/01_manhattan.png", width = 10, height = 5, units = "in", res = 300)
print(manhattan)
dev.off()
log_msg("Saved Manhattan plot: results/figures/01_manhattan.png")
writeLines(log_lines, "results/logs/02_locus_selection_summary.txt")
log_msg("Summary written: results/logs/02_locus_selection_summary.txt")
