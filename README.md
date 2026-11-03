# GWAS-LD-and-Functional-Variant-Analysis

A reproducible workflow for genome-wide association studies (GWAS) and linkage disequilibrium (LD) analysis using **PLINK 1.9** and **R**, applied to simulated genotype data with both a binary (case/control) and a quantitative phenotype. It also includes a functional annotation of known disease variants with Ensembl VEP.

## Workflow

```
Raw genotypes (.bed/.bim/.fam)
        │
        ▼
1. QC ── missingness, HWE, MAF filtering           hwe_missingness_analysis.sh
        │
        ▼
2. PCA ── population structure on LD-pruned SNPs   gwas_ld_analysis.sh
        │
        ▼
3. GWAS ── logistic (binary) / linear (quantitative),
           adjusted for top 5 PCs; Bonferroni + FDR  gwas_ld_analysis.sh
        │
        ▼
4. LD ── haplotype block estimation                gwas_ld_analysis.sh
        │
        ▼
5. Summaries & plots ── QC histograms, QQ plots with
   genomic inflation (λGC), Manhattan plots, PCA    IDCGV_GWAS_Analysis_and_Plotting.R
```

## Scripts

| Script | Description |
|---|---|
| `hwe_missingness_analysis.sh` | Reports per-individual and per-SNP missingness and HWE p-values, then writes a QC-filtered dataset (`--geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-6`). |
| `gwas_ld_analysis.sh` | LD pruning and PCA, logistic and linear association tests with PC covariates, multiple-testing correction (`--adjust`) and LD block estimation. |
| `IDCGV_GWAS_Analysis_and_Plotting.R` | Reads PLINK outputs, extracts genome-wide significant hits (p < 5×10⁻⁸), computes λGC and saves QC, QQ, Manhattan and PCA plots plus summary tables. |

## Usage

```bash
# 1. Edit the paths at the top of each shell script, then:
qsub hwe_missingness_analysis.sh      # or: sh hwe_missingness_analysis.sh
qsub gwas_ld_analysis.sh

# 2. Summarise and plot
Rscript IDCGV_GWAS_Analysis_and_Plotting.R
```

Binary phenotypes follow PLINK coding (1 = control, 2 = case).

## Outputs

| Folder | Contents |
|---|---|
| `dichotomous_missing/`, `quantitative_missing/` | Individual (`.imiss`) and SNP (`.lmiss`) missingness |
| `dichotomous_hwe/`, `quantitative_hwe/` | Hardy-Weinberg equilibrium results (`.hwe`) |
| `dichotomous_gwas/`, `quantitative_gwas/` | Association results (`.assoc.logistic`, `.assoc.linear`) and corrected p-values (`.adjusted`) |
| `pca/` | LD-pruned SNP list and principal components |
| `ld/` | LD blocks |
| `figures/` | QC histograms, QQ, Manhattan and PCA plots |
| `results/` | Significant hits and QC/GWAS summary tables |

## Requirements

- PLINK 1.9
- R (≥ 4.0) with the `qqman` package: `install.packages("qqman")`

## Data

The project uses simulated genotype and phenotype data. Input datasets are not included in this repository.

## Functional variant analysis

As a separate component, two well-characterised **NOD2** variants associated with Crohn's disease susceptibility were annotated with the Ensembl Variant Effect Predictor (VEP) web interface.

| Variant | Consequence | Analysis |
|---|---|---|
| rs2066844 (R702W) | Missense | SIFT, PolyPhen-2 and CADD predictions of impact on protein function |
| rs2066847 (1007fs) | Frameshift | Predicted loss-of-function consequence |

Because the GWAS used simulated data, NOD2 serves as an illustrative candidate gene rather than a locus identified from these results.

## Acknowledgements

Developed as part of the MSc Bioinformatics programme at the University of Glasgow.
