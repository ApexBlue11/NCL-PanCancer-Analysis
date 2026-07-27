# Supplementary table index

Supplementary tables are referred to by S-number in the manuscript and stored under their analysis name in `results/tables/`.

| Supplementary | File | Contents |
|---|---|---|
| **S1** | `S1_cohorts_and_tissue_mapping.tsv` | TCGA study abbreviations, full names, sample counts, and TCGA-GTEx normal tissue mapping with match-quality annotation |
| **S2** | `T1_differential_expression.tsv` | NCL differential expression per cancer, both normal comparators |
| **S3** | `T2_stage_association.tsv` | Association between NCL expression and pathological stage |
| **S4** | `T3_survival.tsv` | Log-rank, univariate and multivariable Cox results with proportional-hazards diagnostics |
| **S5** | `T4_immune_infiltration.tsv` | All cancer x cell type x algorithm infiltration tests |
| **S6** | `T5_algorithm_concordance.tsv` | Cross-algorithm concordance per cancer and canonical cell type |
| **S7** | `T6_checkpoints.tsv` | Immune checkpoint correlations, unadjusted and adjusted for purity and proliferation |
| **S8** | `T8_genomic_scores.tsv` | Immune, stromal and microenvironment scores, TMB, MSI, aneuploidy, FGA |
| **S9** | `T7_cptac_validation.tsv` | CPTAC protein-level validation |
| **S10** | `T10_gsea_per_cancer.tsv.gz` | Per-cancer GSEA against Hallmark and Reactome |
| **S11** | `T9_ncl_gene_correlations.tsv.gz` | Within-cancer Spearman correlation of NCL against every gene (GSEA input) |
