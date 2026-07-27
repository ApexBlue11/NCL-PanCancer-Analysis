"""Build the supplementary tables as they are referred to in the manuscript.

The analysis writes result tables under working names (T1..T10). The manuscript
refers to them by S-number. A reader browsing the repository for "Supplementary
Table S5" should not have to consult an index to discover it is called
T4_immune_infiltration.tsv, so this script materialises the S-numbered set:

  results/supplementary/S1_*.tsv ... S11_*.tsv   one file per supplementary table
  results/supplementary/Supplementary_Tables_S1-S10.xlsx   multi-sheet workbook
  results/supplementary/README.md                          index with row counts

S11 (the gene-by-cancer correlation matrix, ~58,000 rows x 33 columns) is
distributed only as a compressed TSV: it is too large to sit comfortably in a
workbook alongside the others.

S1 is generated here from the cohort definitions; the rest are derived from the
analysis outputs, so the supplementary set cannot drift from the results.
"""
import os, shutil
import pandas as pd

import data_io as D
from cohorts import COHORTS, FULL_NAME

SUPP = os.path.join(D.RESULTS, "supplementary")
os.makedirs(SUPP, exist_ok=True)

# S-number -> (source file in results/tables, destination stem, description)
S_TABLES = [
    ("S1", None, "S1_cohorts_and_tissue_mapping",
     "TCGA study abbreviations, full names, tumour and normal sample counts, and "
     "TCGA-GTEx normal tissue mapping with match-quality annotation"),
    ("S2", "T1_differential_expression.tsv", "S2_differential_expression",
     "NCL differential expression per cancer against both normal comparators, with "
     "Cliff's delta, Hedges' g, confidence intervals and FDR-adjusted p-values"),
    ("S3", "T2_stage_association.tsv", "S3_stage_association",
     "Association between NCL expression and pathological stage: per-stage sample "
     "sizes, Kruskal-Wallis and Jonckheere-Terpstra statistics, late-versus-early "
     "effect sizes"),
    ("S4", "T3_survival.tsv", "S4_survival",
     "Log-rank, univariate Cox and multivariable Cox results for overall, "
     "disease-specific and progression-free survival, with covariates used and "
     "proportional-hazards diagnostics"),
    ("S5", "T4_immune_infiltration.tsv", "S5_immune_infiltration",
     "All cancer x cell type x algorithm immune infiltration tests, unadjusted and "
     "purity-adjusted, across seven deconvolution algorithms"),
    ("S6", "T5_algorithm_concordance.tsv", "S6_algorithm_concordance",
     "Cross-algorithm concordance per cancer and canonical cell type, including the "
     "number of algorithms resolving each cell type"),
    ("S7", "T6_checkpoints.tsv", "S7_immune_checkpoints",
     "NCL correlations with 16 immune checkpoint and immunomodulatory genes, "
     "unadjusted and adjusted for tumour purity and for proliferation"),
    ("S8", "T8_genomic_scores.tsv", "S8_genomic_and_immune_scores",
     "NCL correlations with immune, stromal and microenvironment scores, tumour "
     "mutational burden, MANTIS MSI score, aneuploidy score and fraction of genome "
     "altered"),
    ("S9", "T7_cptac_validation.tsv", "S9_cptac_validation",
     "Independent CPTAC protein-level validation across nine cohorts, rank-sum and "
     "paired tests"),
    ("S10", "T10_gsea_per_cancer.tsv.gz", "S10_gsea_per_cancer",
     "Per-cancer pre-ranked GSEA against Hallmark and Reactome collections"),
    ("S11", "T9_ncl_gene_correlations.tsv.gz", "S11_ncl_gene_correlations",
     "Within-cancer Spearman correlation of NCL against every gene; the GSEA "
     "ranking statistic"),
]

# S11 is excluded from the workbook on size grounds.
IN_WORKBOOK = {s for s, *_ in S_TABLES if s != "S11"}


def build_s1():
    ann, gtex = D.sample_groups()
    tum = ann[ann.group == "tumour"].cohort.value_counts()
    nor = ann[ann.group == "tcga_normal"].cohort.value_counts()
    rows = []
    for code in sorted(COHORTS):
        _detailed, site, quality = COHORTS[code]
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
    frames, written = {}, []

    for snum, src, stem, desc in S_TABLES:
        if snum == "S1":
            df = build_s1()
            dest = os.path.join(SUPP, stem + ".tsv")
            df.to_csv(dest, sep="\t", index=False)
            # keep the copy under results/tables for backwards compatibility
            df.to_csv(os.path.join(D.TABLES, "S1_cohorts_and_tissue_mapping.tsv"),
                      sep="\t", index=False)
        else:
            src_path = os.path.join(D.TABLES, src)
            if not os.path.exists(src_path):
                print("  MISSING SOURCE for %s: %s" % (snum, src))
                continue
            df = pd.read_csv(src_path, sep="\t", low_memory=False)
            gz = src.endswith(".gz")
            dest = os.path.join(SUPP, stem + (".tsv.gz" if gz else ".tsv"))
            shutil.copyfile(src_path, dest)

        frames[snum] = df
        written.append((snum, os.path.basename(dest), len(df), df.shape[1], desc))
        print("  %-4s %-38s %7d rows x %2d cols" % (snum, os.path.basename(dest),
                                                    len(df), df.shape[1]))

    # ---- multi-sheet workbook
    xlsx = os.path.join(SUPP, "Supplementary_Tables_S1-S10.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        idx = pd.DataFrame(
            [{"Table": s, "File": f, "Rows": r, "Columns": c, "Contents": d}
             for s, f, r, c, d in written])
        idx.to_excel(xw, sheet_name="Index", index=False)
        for snum, _f, _r, _c, _d in written:
            if snum not in IN_WORKBOOK:
                continue
            frames[snum].to_excel(xw, sheet_name=snum, index=False)
    print("\n  wrote %s (%.1f MB)" % (os.path.basename(xlsx),
                                      os.path.getsize(xlsx) / 1e6))

    # ---- index for the repository
    lines = ["# Supplementary tables", "",
             "Supplementary tables as referred to by S-number in the manuscript.",
             "",
             "`Supplementary_Tables_S1-S10.xlsx` contains S1 to S10 as separate "
             "sheets. S11 is distributed only as a compressed TSV because it is too "
             "large for a workbook.", "",
             "| Table | File | Rows | Cols | Contents |", "|---|---|---|---|---|"]
    for s, f, r, c, d in written:
        lines.append("| **%s** | `%s` | %s | %s | %s |" % (s, f, format(r, ","), c, d))
    lines += ["", "## Working names", "",
              "The analysis scripts write these same tables to `results/tables/` "
              "under working names (T1..T10). The mapping is:", "",
              "| Supplementary | Working name |", "|---|---|"]
    for snum, src, _stem, _desc in S_TABLES:
        lines.append("| %s | `%s` |" % (snum, src or "generated by 12_supplementary_tables.py"))
    lines.append("")
    with open(os.path.join(SUPP, "README.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("  wrote supplementary/README.md")

    # keep the older index in results/tables consistent
    with open(os.path.join(D.TABLES, "SUPPLEMENTARY_INDEX.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    missing = [s for s, *_ in S_TABLES if s not in frames]
    print("\n%d of %d supplementary tables written%s"
          % (len(written), len(S_TABLES),
             "" if not missing else "; MISSING: %s" % missing))


if __name__ == "__main__":
    main()
