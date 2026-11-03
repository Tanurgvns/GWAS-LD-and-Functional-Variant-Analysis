#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l cput=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -m abe
#PBS -M your.email@example.com

# =============================================================================
# Step 1: Genotype quality control with PLINK 1.9
#   - Reports genotype missingness (per individual and per SNP)
#   - Reports Hardy-Weinberg equilibrium (HWE) p-values
#   - Writes a filtered dataset used by gwas_ld_analysis.sh
# =============================================================================

set -e

# ---- Paths (edit these) -----------------------------------------------------
GENO_FILE="data/genotype_data"        # PLINK binary prefix (.bed/.bim/.fam)
FILTERED="data/filtered_data"         # output prefix for QC-filtered data
LOG_FILE="qc_log.txt"

# ---- QC thresholds ----------------------------------------------------------
GENO=0.05     # remove SNPs missing in >5% of individuals
MIND=0.05     # remove individuals missing >5% of genotypes
MAF=0.01      # remove SNPs with minor allele frequency <1%
HWE=1e-6      # remove SNPs failing HWE at p < 1e-6

echo "QC started at $(date)" > "$LOG_FILE"

for PHENO in dichotomous quantitative; do
    mkdir -p "${PHENO}_missing" "${PHENO}_hwe"
    for CHR in 1 2; do
        echo "Missingness and HWE reports: $PHENO, chr$CHR" >> "$LOG_FILE"

        # Missingness report -> .imiss (individuals) and .lmiss (SNPs)
        plink --bfile "$GENO_FILE" --chr "$CHR" --missing \
              --out "${PHENO}_missing/missing_output_chr${CHR}"

        # HWE report -> .hwe
        plink --bfile "$GENO_FILE" --chr "$CHR" --hardy \
              --out "${PHENO}_hwe/hardy_output_chr${CHR}"
    done
done

# Apply all QC filters and write a clean dataset for the GWAS step
echo "Writing filtered dataset" >> "$LOG_FILE"
plink --bfile "$GENO_FILE" \
      --geno "$GENO" --mind "$MIND" --maf "$MAF" --hwe "$HWE" \
      --make-bed --out "$FILTERED"

echo "QC finished at $(date)" >> "$LOG_FILE"
