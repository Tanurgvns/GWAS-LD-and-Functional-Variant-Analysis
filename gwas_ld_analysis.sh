#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l cput=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M your.email@example.com

# =============================================================================
# Step 2: Population structure, GWAS and LD analysis with PLINK 1.9
#   - PCA on LD-pruned SNPs to capture population structure
#   - Logistic regression GWAS (binary trait) and linear regression GWAS
#     (quantitative trait), adjusted for the top principal components
#   - Multiple-testing correction (Bonferroni, FDR) via --adjust
#   - LD block estimation
# Run hwe_missingness_analysis.sh first.
# =============================================================================

set -e

# ---- Paths (edit these) -----------------------------------------------------
FILTERED="data/filtered_data"            # QC-filtered PLINK prefix from step 1
PHENO_BINARY="data/pheno_binary.txt"     # FID IID PHENO (1 = control, 2 = case)
PHENO_QUANT="data/pheno_quant.txt"       # FID IID PHENO (continuous)
LOG_FILE="gwas_log.txt"

N_PCS=5                                   # principal components used as covariates

echo "GWAS started at $(date)" > "$LOG_FILE"
mkdir -p pca dichotomous_gwas quantitative_gwas ld

# ---- 1. Population structure (PCA) ------------------------------------------
echo "LD pruning and PCA" >> "$LOG_FILE"
plink --bfile "$FILTERED" --indep-pairwise 50 5 0.2 --out pca/pruned
plink --bfile "$FILTERED" --extract pca/pruned.prune.in --pca 10 --out pca/pca

# ---- 2. Association testing, per chromosome ---------------------------------
for CHR in 1 2; do
    echo "Logistic GWAS (binary trait), chr$CHR" >> "$LOG_FILE"
    plink --bfile "$FILTERED" --chr "$CHR" \
          --pheno "$PHENO_BINARY" \
          --covar pca/pca.eigenvec --covar-number 1-"$N_PCS" \
          --logistic hide-covar --adjust \
          --out dichotomous_gwas/gwas_results_chr${CHR}

    echo "Linear GWAS (quantitative trait), chr$CHR" >> "$LOG_FILE"
    plink --bfile "$FILTERED" --chr "$CHR" \
          --pheno "$PHENO_QUANT" \
          --covar pca/pca.eigenvec --covar-number 1-"$N_PCS" \
          --linear hide-covar --adjust \
          --out quantitative_gwas/gwas_results_chr${CHR}
done

# ---- 3. Linkage disequilibrium ----------------------------------------------
echo "LD block estimation" >> "$LOG_FILE"
plink --bfile "$FILTERED" --blocks no-pheno-req --out ld/ld_results

echo "GWAS finished at $(date)" >> "$LOG_FILE"
