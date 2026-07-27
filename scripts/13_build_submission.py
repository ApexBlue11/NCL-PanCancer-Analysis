"""Assemble a submission copy of the manuscript with figures and tables embedded.

The working manuscript carries only legends and table captions. Journals want
the actual figures and tables in the file, so this script renders publication-
formatted Tables 1-6 from the result tables, places each figure image directly
beneath its legend, and writes manuscript_submission.md for conversion to DOCX.

Figure files are matched to legends by number and the mapping is asserted, so a
renamed or missing figure fails loudly rather than silently producing a
manuscript with the wrong image under a legend.
"""
import os, re
import pandas as pd
import numpy as np

import data_io as D

MS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "manuscript")
SRC = os.path.join(MS_DIR, "manuscript_revised.md")
DEST = os.path.join(MS_DIR, "manuscript_submission.md")

# Figure number -> (file, keyword that must appear in that figure's legend)
FIGURES = {
    1: ("Figure1_expression.png", "expression"),
    2: ("Figure2_survival.png", "Survival"),
    3: ("Figure3_immune.png", "infiltration"),
    4: ("Figure4_checkpoints.png", "checkpoint"),
    5: ("Figure5_cptac.png", "CPTAC"),
    6: ("Figure6_gsea.png", "GSEA"),
}


def L(f):
    return pd.read_csv(os.path.join(D.TABLES, f), sep="\t", low_memory=False)


def fmt_q(q):
    if not np.isfinite(q):
        return "-"
    if q < 1e-4:
        return "%.1e" % q
    return "%.4f" % q


def md_table(df, align=None):
    cols = list(df.columns)
    align = align or (["---"] * len(cols))
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join(align) + " |"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def table1():
    de = L("T1_differential_expression.tsv")
    g = de[(de.comparator == "GTEx") & de.q.notna()].set_index("cohort")
    a = de[(de.comparator == "TCGA_adjacent") & de.q.notna()].set_index("cohort")
    rows = []
    for c in g.sort_values("cliffs_delta", ascending=False).index:
        r = g.loc[c]
        rows.append({
            "Cancer": c,
            "n tumour": int(r.n_tumour),
            "n GTEx": int(r.n_normal),
            "Cliff's delta (95% CI)": "%+.2f (%+.2f, %+.2f)" % (r.cliffs_delta, r.delta_lo, r.delta_hi),
            "q": fmt_q(r.q),
            "Tissue match": r.tissue_match,
            "Adjacent-normal delta": ("%+.2f" % a.loc[c].cliffs_delta) if c in a.index else "n/a",
            "Adjacent q": fmt_q(a.loc[c].q) if c in a.index else "n/a",
        })
    return pd.DataFrame(rows)


def table2():
    st = L("T2_stage_association.tsv").sort_values("trend_p")
    rows = []
    for _, r in st.iterrows():
        rows.append({
            "Cancer": r.cohort,
            "n": int(r.n_total),
            "I / II / III / IV": "%d / %d / %d / %d" % (r.n_I, r.n_II, r.n_III, r.n_IV),
            "Kruskal-Wallis q": fmt_q(r.kruskal_q),
            "Trend z": "%+.2f" % r.trend_z,
            "Trend q": fmt_q(r.trend_q),
            "Stage IV vs I delta (95% CI)": "%+.2f (%+.2f, %+.2f)" % (
                r.late_vs_early_delta, r.delta_lo, r.delta_hi),
        })
    return pd.DataFrame(rows)


def table3():
    sv = L("T3_survival.tsv")
    sv = sv[(sv.endpoint == "OS") & sv.adj_HR.notna()].sort_values("adj_HR", ascending=False)
    rows = []
    for _, r in sv.iterrows():
        ph = r.get("ph_p_NCL", np.nan)
        rows.append({
            "Cancer": r.cohort,
            "n": int(r.n_adj),
            "Events": int(r.events_adj),
            "Adjusted HR (95% CI)": "%.2f (%.2f, %.2f)" % (r.adj_HR, r.adj_HR_lo, r.adj_HR_hi),
            "Adjusted q": fmt_q(r.adj_q),
            "Univariate q": fmt_q(r.uni_q),
            "Covariates": r.adj_covariates.replace("male", "sex"),
            "PH violated": "yes" if (np.isfinite(ph) and ph < 0.05) else "no",
        })
    return pd.DataFrame(rows)


def table4():
    gs = L("T8_genomic_scores.tsv")
    rows = []
    for m in ["Immune score", "Stromal score", "Microenvironment score",
              "TMB", "MSI (MANTIS)", "Aneuploidy", "FGA"]:
        s = gs[gs.measure == m]
        sig = s[s.q < 0.05]
        rows.append({
            "Measure": m,
            "Cancers tested": len(s),
            "Significant (q<0.05)": len(sig),
            "Positive": int((sig.rho > 0).sum()),
            "Negative": int((sig.rho < 0).sum()),
            "Median rho": "%+.3f" % s.rho.median(),
            "Range of significant rho": ("%+.2f to %+.2f" % (sig.rho.min(), sig.rho.max())) if len(sig) else "-",
        })
    return pd.DataFrame(rows)


def table5():
    cp = L("T6_checkpoints.tsv")
    rows = []
    for alias, g in cp.groupby("alias"):
        rob = g[(g.q_purity_adj < 0.05) & (g.q_prolif_adj < 0.05)]
        pos = int((rob.rho_prolif_adj > 0).sum())
        neg = int((rob.rho_prolif_adj < 0).sum())
        rows.append({
            "Molecule": alias,
            "Gene": g.gene.iloc[0],
            "Median rho": "%+.3f" % g.rho.median(),
            "Significant unadjusted": int((g.q < 0.05).sum()),
            "Robust after both adjustments": len(rob),
            "Positive / negative": "%d / %d" % (pos, neg),
            "Robust rho range": ("%+.2f to %+.2f" % (rob.rho_prolif_adj.min(),
                                                     rob.rho_prolif_adj.max())) if len(rob) else "-",
        })
    df = pd.DataFrame(rows)
    return df.sort_values("Robust after both adjustments", ascending=False)


def table6():
    cv = L("T7_cptac_validation.tsv")
    cv = cv[cv.cliffs_delta.notna()].sort_values("cliffs_delta", ascending=False)
    rows = []
    for _, r in cv.iterrows():
        rows.append({
            "CPTAC cohort": r.cohort,
            "n tumour": int(r.n_tumour),
            "n normal": int(r.n_normal),
            "n paired": int(r.n_paired) if np.isfinite(r.n_paired) else 0,
            "Cliff's delta (95% CI)": "%+.2f (%+.2f, %+.2f)" % (r.cliffs_delta, r.delta_lo, r.delta_hi),
            "q (rank-sum)": fmt_q(r.q_ranksum),
            "q (paired)": fmt_q(r.q_paired) if np.isfinite(r.get("q_paired", np.nan)) else "-",
        })
    return pd.DataFrame(rows)


BUILDERS = {1: table1, 2: table2, 3: table3, 4: table4, 5: table5, 6: table6}


def main():
    t = open(SRC, encoding="utf-8").read()

    # ---- 1. place each figure image immediately after its legend
    legends = t[t.index("## Figure legends"):t.index("## Tables")]
    for num, (fname, keyword) in FIGURES.items():
        path = os.path.join(D.FIGURES, fname)
        assert os.path.exists(path), "missing figure file: %s" % path
        m = re.search(r"(\*\*Figure %d\..*?)(?=\n\n\*\*Figure |\Z)" % num, legends, re.S)
        assert m, "no legend found for Figure %d" % num
        legend_text = m.group(1)
        assert keyword.lower() in legend_text.lower(), (
            "Figure %d legend does not mention '%s'; figure/legend mapping may be wrong"
            % (num, keyword))
        img = "\n\n![](../results/figures/%s){width=100%%}\n" % fname
        legends = legends.replace(legend_text, legend_text.rstrip() + img)
    t = t[:t.index("## Figure legends")] + legends + t[t.index("## Tables"):]

    # ---- 2. replace each table caption with caption + rendered table
    tabsec = t[t.index("## Tables"):]
    for num in sorted(BUILDERS):
        df = BUILDERS[num]()
        m = re.search(r"(\*\*Table %d\.\*\*[^\n]*)" % num, tabsec)
        assert m, "no caption found for Table %d" % num
        cap = m.group(1)
        tabsec = tabsec.replace(cap, cap + "\n\n" + md_table(df) + "\n")
        print("  Table %d: %d rows x %d cols" % (num, len(df), df.shape[1]))
    t = t[:t.index("## Tables")] + tabsec

    # ---- 3. rename the figure-legend heading now that figures are inline
    t = t.replace("## Figure legends", "## Figures")

    with open(DEST, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("\nwrote %s" % DEST)
    print("figures embedded: %d" % len(FIGURES))


if __name__ == "__main__":
    main()
