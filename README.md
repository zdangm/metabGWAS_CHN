This repository contains codes from the following project:
Genetic mapping of serum metabolome to chronic diseases among Han Chinese

1.Heritability
The SNP-based heritability was estimated using GCTA-GREML.

2.GWAS
GWAS for the metabolites.
The GWAS was performed using the mixed linear model-based association analysis (MLMA) implemented in the GCTA software (version 1.94.0 beta)

3.Meta-analysis
For the three cohorts in the discovery set, we combined cohort-specific GWAS summary statistics-level data in a fixed-effect meta-analysis by inverse variance weighting (IVW) method using METAL (version 1.7)

4.Mendelian randomization
For metabolites with at least three IVs, inverse variance weighting (IVW) and MR-Egger were implemented as baseline methods. 
The intercept test in MR-Egger was implemented to detect potential horizontal pleiotropy. 
To ensure the robustness of associations, we employed another two MR methods: MR Robust Adjusted Profile Scoring (MR-RAPS) and Generalized Summary-data-based Mendelian Randomization (GSMR).
