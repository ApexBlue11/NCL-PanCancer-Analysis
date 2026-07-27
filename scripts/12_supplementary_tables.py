"""Emit Supplementary Table S1 and the S-number -> result-file mapping.

The manuscript refers to supplementary tables by S-number while the repository
stores them under analysis names (T1..T10). Without an explicit mapping a reader
cannot connect a claim in the text to the file that supports it, so the mapping
is generated here rather than left implicit.
"""
import os
import pandas as pd

import data_io as D
from cohorts import COHORTS, FULL_NAME

# S-number -> (repository file, description). S1 is generated below.
S_MAP = [
    ("S1", "S1_cohorts_and_tissue_mapping.tsv",
     "TCGA study abbreviations, full names, sample counts, and TCGA-GTEx normal "
     "tissue mapping with match-quality annotation"),
    ("S2", "T1_differential_expression.tsv",
     "NCL differential expression per cancer, both normal comparators"),
    ("S3", "T2_stage_association.tsv",
     "Association between NCL expression and pathological stage"),
    ("S4", "T3_survival.tsv",
     "Log-rank, univariate and multivariable Cox results with "
     "proportional-hazards diagnostics"),
    ("S5", "T4_immune_infiltration.tsv",
     "All cancer x cell type x algorithm infiltration tests"),
    ("S6", "T5_algorithm_concordance.tsv",
     "Cross-algorithm concordance per cancer and canonical cell type"),
    ("S7", "T6_checkpoints.tsv",
     "Immune checkpoint correlations, unadjusted and adjusted for purity and "
     "proliferation"),
    ("S8", "T8_genomic_scores.tsv",
     "Immune, stromal and microenvironment scores, TMB, MSI, aneuploidy, FGA"),
    ("S9", "T7_cptac_validation.tsv",
     "CPTAC protein-level validation"),
    ("S10", "T10_gsea_per_cancer.tsv.gz",
     "Per-cancer GSEA against Hallmark and Reactome"),
    ("S11", "T9_ncl_gene_correlations.tsv.gz",
     "Within-cancer Spearman correlation of NCL against every gene (GSEA input)"),
]


def build_s1():
    ann, gtex = D.sample_groups()
    tum = ann[ann.group == "tumour"].cohort.value_counts()
    nor = ann[ann.group == "tcga_normal"].cohort.value_counts()
    rows = []
    for code in sorted(COHORTS):
        detailed, site, quality = COHORTS[code]
        rows.append({
            "TCGA_abbreviation": code,
            "full_name": FULL_NAME[code],
            "n_tumour": int(tum.get(code, 0)),
            "n_TCGA_adjacent_normal": int(nor.get(code, 0)),
            "GTEx_tissue": site or "none available",
            "n_GTEx_normal": len(gtex.get(code, [])),
            "tissue_match_quality": quality,
        })
    return pd.DataFrame(rows)


def main():
    s1 = build_s1()
    dest = os.path.join(D.TABLES, "S1_cohorts_and_tissue_mapping.tsv")
    s1.to_csv(dest, sep="\t", index=False)
    print("wrote %s  (%d cohorts)" % (dest, len(s1)))
    print("  good=%d approximate=%d none=%d"
          % ((s1.tissue_match_quality == "good").sum(),
             (s1.tissue_match_quality == "approximate").sum(),
             (s1.tissue_match_quality == "none").sum()))

    lines = ["# Supplementary table index", "",
             "Supplementary tables are referred to by S-number in the manuscript "
             "and stored under their analysis name in `results/tables/`.", "",
             "| Supplementary | File | Contents |", "|---|---|---|"]
    missing = []
    for s, fn, desc in S_MAP:
        exists = os.path.exists(os.path.join(D.TABLES, fn))
        if not exists:
            missing.append(fn)
        lines.append("| **%s** | `%s`%s | %s |"
                     % (s, fn, "" if exists else " *(missing)*", desc))
    lines.append("")
    with open(os.path.join(D.TABLES, "SUPPLEMENTARY_INDEX.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("wrote %s" % os.path.join(D.TABLES, "SUPPLEMENTARY_INDEX.md"))
    print("missing files: %s" % (missing or "none"))


if __name__ == "__main__":
    main()
