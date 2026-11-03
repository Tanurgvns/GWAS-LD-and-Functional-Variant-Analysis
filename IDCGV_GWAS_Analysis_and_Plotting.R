# =============================================================================
# GWAS QC, association results and visualisation
# Reads PLINK 1.9 outputs produced by hwe_missingness_analysis.sh and
# gwas_ld_analysis.sh, summarises them and saves figures to figures/.
# Run from the repository root:  Rscript IDCGV_GWAS_Analysis_and_Plotting.R
# =============================================================================

library(qqman)

GENOME_WIDE <- 5e-8      # genome-wide significance threshold
SUGGESTIVE  <- 1e-5      # suggestive threshold
PHENOTYPES  <- c("dichotomous", "quantitative")
CHROMOSOMES <- c(1, 2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# PLINK pads its output with spaces, so read with whitespace separation
read_plink <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read.table(path, header = TRUE, stringsAsFactors = FALSE)
}

# Save a base-R plot straight to PNG
save_png <- function(filename, plot_fun) {
  png(file.path("figures", filename), width = 1600, height = 1000, res = 200)
  plot_fun()
  dev.off()
}

# Genomic inflation factor (lambda GC): ~1 means no systematic inflation
lambda_gc <- function(p) {
  p <- p[!is.na(p)]
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
}

# -----------------------------------------------------------------------------
# 1. Quality control: missingness and Hardy-Weinberg equilibrium
# -----------------------------------------------------------------------------
qc_plots <- function(pheno, chr) {
  tag   <- paste0(pheno, "_chr", chr)
  imiss <- read_plink(sprintf("%s_missing/missing_output_chr%s.imiss", pheno, chr))
  lmiss <- read_plink(sprintf("%s_missing/missing_output_chr%s.lmiss", pheno, chr))
  hwe   <- read_plink(sprintf("%s_hwe/hardy_output_chr%s.hwe", pheno, chr))

  save_png(paste0(tag, "_individual_missingness.png"), function()
    hist(imiss$F_MISS, breaks = 20, col = "skyblue",
         main = paste("Individual missingness -", tag),
         xlab = "Proportion of missing genotypes per individual"))

  save_png(paste0(tag, "_snp_missingness.png"), function()
    hist(lmiss$F_MISS, breaks = 20, col = "salmon",
         main = paste("SNP missingness -", tag),
         xlab = "Proportion of missing genotypes per SNP"))

  save_png(paste0(tag, "_hwe_pvalues.png"), function()
    hist(hwe$P, breaks = 50, col = "orchid",
         main = paste("HWE p-values -", tag),
         xlab = "HWE p-value"))

  data.frame(analysis = tag,
             individuals = nrow(imiss),
             snps = nrow(lmiss),
             mean_ind_missing = round(mean(imiss$F_MISS), 4),
             mean_snp_missing = round(mean(lmiss$F_MISS), 4),
             snps_hwe_p_below_1e6 = sum(hwe$P < 1e-6, na.rm = TRUE))
}

# -----------------------------------------------------------------------------
# 2. Association results: significant hits, QQ and Manhattan plots
# -----------------------------------------------------------------------------
gwas_results <- function(pheno, chr) {
  tag   <- paste0(pheno, "_chr", chr)
  model <- if (pheno == "dichotomous") "logistic" else "linear"
  gwas  <- read_plink(sprintf("%s_gwas/gwas_results_chr%s.assoc.%s", pheno, chr, model))

  # Keep the additive SNP effect only and drop failed tests
  if ("TEST" %in% names(gwas)) gwas <- gwas[gwas$TEST == "ADD", ]
  gwas <- gwas[!is.na(gwas$P), ]

  hits <- gwas[gwas$P < GENOME_WIDE, ]
  hits <- hits[order(hits$P), ]
  write.table(hits, file.path("results", paste0(tag, "_significant_hits.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  lam <- lambda_gc(gwas$P)

  save_png(paste0(tag, "_qq.png"), function()
    qq(gwas$P, main = sprintf("QQ plot - %s (lambda GC = %.3f)", tag, lam)))

  save_png(paste0(tag, "_manhattan.png"), function()
    manhattan(gwas, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
              main = paste("Manhattan plot -", tag),
              genomewideline = -log10(GENOME_WIDE),
              suggestiveline = -log10(SUGGESTIVE)))

  data.frame(analysis = tag,
             snps_tested = nrow(gwas),
             genome_wide_hits = nrow(hits),
             suggestive_hits = sum(gwas$P < SUGGESTIVE),
             lambda_gc = round(lam, 3),
             top_snp = if (nrow(hits) > 0) hits$SNP[1] else NA)
}

# -----------------------------------------------------------------------------
# 3. Population structure: PCA plot (if PCA was run)
# -----------------------------------------------------------------------------
if (file.exists("pca/pca.eigenvec")) {
  pcs <- read.table("pca/pca.eigenvec", header = FALSE)
  names(pcs) <- c("FID", "IID", paste0("PC", seq_len(ncol(pcs) - 2)))
  save_png("pca_pc1_pc2.png", function()
    plot(pcs$PC1, pcs$PC2, pch = 19, col = adjustcolor("steelblue", 0.6),
         xlab = "PC1", ylab = "PC2", main = "Population structure (PCA)"))
}

# -----------------------------------------------------------------------------
# Run everything and write summary tables
# -----------------------------------------------------------------------------
grid <- expand.grid(pheno = PHENOTYPES, chr = CHROMOSOMES, stringsAsFactors = FALSE)

qc_summary   <- do.call(rbind, Map(qc_plots,     grid$pheno, grid$chr))
gwas_summary <- do.call(rbind, Map(gwas_results, grid$pheno, grid$chr))

write.table(qc_summary,   "results/qc_summary.tsv",   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gwas_summary, "results/gwas_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

print(qc_summary)
print(gwas_summary)
cat("Figures saved to figures/, tables saved to results/\n")
