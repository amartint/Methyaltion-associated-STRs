# Methylation associated STRs
Scripts used for performing methylation QTL mapping using STR genotypes.

## QTL mapping of STRs associated with DNA methylation levels (mSTRs)
This repository contains scripts for performing QTL mapping analysis using Short Tandem Repeats (STR) genotypes and genome-wide DNA methylation profiling generated from PCR-free Illumina sequencing using [HipSTR](https://github.com/HipSTR-ToolHipSTR) and the Illumina Infinium MethylationEPIC array.

### Genotyping of tandem repeats
Genotypes were generated for STRs included in [HipSTR-references](https://github.com/HipSTR-Tool/HipSTR-references). Resulting genotypes were then filtered using default parameters as indicated in [HipSTR-tool](https://hipstr-tool.github.io/HipSTR/#default-filtering), and selected according to the listed criteria

•	minimum genotyping rate of ≥50%

•	motif length of ≥2 bp

•	≥3 observed alleles across the tested cohort

•	non-modal allele frequency (NMAF) ≥0.1, _i.e._ the frequency of the major allele <0.9.

### Quality Control

Quality control (QC) of DNA methylation measurements (beta (β) values) and HipSTR genotypes were based on Principal Component Analysis (PCA) and density plots of β values and repeat length (average of allelic size), respectively. 

•	PCA was generated using prcomp function in R and outlier samples were removed based on PCs by manual inspection. 

•	Density plots were generated using density function and outliers from the distribution were removed.

### Additional QC steps

•	Removing low variance DNA methylation levels (standard deviation β values >0.02 across samples)

### QTL mapping

•	Methylation: Association test between HipSTR-derived genotypes and β values was performed using lm function in R.

•	Expression: Association test between HipSTR-derived genotypes and gene expression levels (RNA sequencing) was performed using lm function in R. 

Resulting associations were then filtered based on the number of observations per pairwise comparison, i.e, STR:CpG or STR:gene) and Bonferroni correction was applied to account for multiple testing using p.adjust function in R 

### Additional analysis

•	Fine-mapping of identified mSTRs. Identified mSTRs were fine-mapped against local SNVs using [CAVIAR](https://github.com/fhormoz/caviar).

•	Linkage disequilibrium (LD) between STRs and local SNVs (±250kb). Local LD between STR genotypes (average allelic size) and SNV dosages (alternate allele content, _i.e._ homozygous for reference allele 0, heterozygous 1 and homozygous for alternate allele 2) were computed using lm function in R. 

•	Population Stratification of STR genotypes with the Vst statistic ([Redon _et al._](https://www.nature.com/articles/nature05329)).
