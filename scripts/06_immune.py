"""Immune microenvironment analysis: infiltration, checkpoints, TMB/MSI and scores.

This replaces submitted Figure 3 and supplies the analysis the submitted title
and abstract asserted but never performed.

Three things are fixed here.

1. Algorithm consistency (Reviewer 1 point 3; Reviewer 2 point 5). The submitted
   Figure 3 used a different deconvolution method in each panel -- MCPcounter in
   3b, CIBERSORT-ABS in 3c, xCell in 3d, quanTIseq in 3e. Every algorithm that
   estimates a given cell type is applied here to every cancer, and concordance
   across algorithms is reported rather than one method being chosen per panel.

2. Immune checkpoints. The submitted abstract claimed NCL "correlates with
   immune checkpoint molecules such as PD-L1, CTLA-4 and TIM-3" and with IL-10
   and TGF-beta. No such analysis existed. It is run here.

3. Confounding (Reviewer 2 point 6). NCL is a ribosome-biogenesis and
   proliferation gene, and proliferation is itself associated with immune
   phenotype. Every checkpoint correlation is therefore reported raw, adjusted
   for tumour purity, and adjusted for a proliferation score, so that any
   immune-specific signal is separated from a generic proliferation signal.
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats

import data_io as D
import statsutil as S
from cohorts import COHORTS, FULL_NAME, CHECKPOINTS, PROLIFERATION

warnings.filterwarnings("ignore")
GENE = "NCL"
MIN_N = 30

# Cell types resolvable by more than one algorithm, used for the concordance test.
CANONICAL = ["T cell CD8+", "Macrophage", "Neutrophil", "B cell", "NK cell",
             "T cell regulatory (Tregs)", "Cancer associated fibroblast",
             "Myeloid dendritic cell", "Monocyte", "Endothelial cell"]


def load_panel():
    """NCL, purity proxy, proliferation score and TIMER2 estimates on shared samples."""
    ncl = D.gene_vector(GENE)
    ix = D.index()
    pos = {s: i for i, s in enumerate(ix["samples"])}
    ann, _ = D.sample_groups()
    tum = ann.index[ann.group == "tumour"]

    tim = D.timer2()
    shared = [s for s in tum if s in tim.index]
    panel = pd.DataFrame(index=shared)
    panel["cohort"] = ann.loc[shared, "cohort"].to_numpy()
    panel["NCL"] = [ncl[pos[s]] for s in shared]

    # Proliferation score: mean of within-cohort z-scores of canonical markers.
    prol, _found = D.gene_block(PROLIFERATION)
    prol = prol[shared].T
    panel["proliferation"] = np.nan
    for c, idx in panel.groupby("cohort").groups.items():
        sub = prol.loc[idx]
        z = (sub - sub.mean()) / sub.std(ddof=0)
        panel.loc[idx, "proliferation"] = z.mean(axis=1).to_numpy()

    # EPIC's "uncharacterized cell" fraction is its estimate of the non-immune,
    # non-stromal compartment, used here as the tumour-purity covariate.
    if "uncharacterized cell_EPIC" in tim.columns:
        panel["purity"] = tim.loc[shared, "uncharacterized cell_EPIC"].to_numpy()
    else:
        panel["purity"] = np.nan
    return panel, tim.loc[shared]


def infiltration(panel, tim):
    rows = []
    cols = [c for c in tim.columns if "_" in c]
    for code in sorted(COHORTS):
        idx = panel.index[panel.cohort == code]
        if len(idx) < MIN_N:
            continue
        x = panel.loc[idx, "NCL"].to_numpy()
        pur = panel.loc[idx, "purity"].to_numpy()
        for col in cols:
            cell, alg = col.rsplit("_", 1)
            y = pd.to_numeric(tim.loc[idx, col], errors="coerce").to_numpy()
            if np.isfinite(y).sum() < MIN_N or np.nanstd(y) == 0:
                continue
            rho, p = stats.spearmanr(x, y, nan_policy="omit")
            pr, pp, n_adj = S.partial_spearman(x, y, pur)
            lo, hi = S.spearman_ci(rho, int(np.isfinite(y).sum()))
            rows.append(dict(cohort=code, cell_type=cell, algorithm=alg,
                             n=int(np.isfinite(y).sum()), rho=rho, rho_lo=lo,
                             rho_hi=hi, p=p, rho_purity_adj=pr, p_purity_adj=pp))
    df = pd.DataFrame(rows)
    df["q"] = S.bh_fdr(df["p"].to_numpy())
    df["q_purity_adj"] = S.bh_fdr(df["p_purity_adj"].to_numpy())
    return df


def concordance(infil):
    """Per cancer and cell type, how many algorithms agree in sign and significance."""
    rows = []
    for (code, cell), g in infil[infil.cell_type.isin(CANONICAL)].groupby(
            ["cohort", "cell_type"]):
        n_alg = len(g)
        if n_alg < 2:
            continue
        sig = g[g.q_purity_adj < 0.05]
        pos_ = (sig.rho_purity_adj > 0).sum()
        neg_ = (sig.rho_purity_adj < 0).sum()
        rows.append(dict(cohort=code, cell_type=cell, n_algorithms=n_alg,
                         n_significant=len(sig), n_positive=pos_, n_negative=neg_,
                         median_rho=g.rho_purity_adj.median(),
                         concordant=(len(sig) >= 2 and (pos_ == 0 or neg_ == 0)),
                         direction=("positive" if pos_ > neg_ else
                                    "negative" if neg_ > pos_ else "mixed")))
    return pd.DataFrame(rows)


def checkpoints(panel):
    block, found = D.gene_block(list(CHECKPOINTS))
    block = block[panel.index].T
    rows = []
    for code in sorted(COHORTS):
        idx = panel.index[panel.cohort == code]
        if len(idx) < MIN_N:
            continue
        x = panel.loc[idx, "NCL"].to_numpy()
        pur = panel.loc[idx, "purity"].to_numpy()
        pro = panel.loc[idx, "proliferation"].to_numpy()
        for g in found:
            y = block.loc[idx, g].to_numpy()
            if np.nanstd(y) == 0:
                continue
            rho, p = stats.spearmanr(x, y, nan_policy="omit")
            lo, hi = S.spearman_ci(rho, len(idx))
            r_pur, p_pur, _ = S.partial_spearman(x, y, pur)
            r_pro, p_pro, _ = S.partial_spearman(x, y, pro)
            rows.append(dict(cohort=code, gene=g, alias=CHECKPOINTS[g], n=len(idx),
                             rho=rho, rho_lo=lo, rho_hi=hi, p=p,
                             rho_purity_adj=r_pur, p_purity_adj=p_pur,
                             rho_prolif_adj=r_pro, p_prolif_adj=p_pro))
    df = pd.DataFrame(rows)
    for c, q in (("p", "q"), ("p_purity_adj", "q_purity_adj"),
                 ("p_prolif_adj", "q_prolif_adj")):
        df[q] = S.bh_fdr(df[c].to_numpy())
    return df


def genomic_scores(panel):
    """NCL vs TMB, MSI, aneuploidy and the xCell immune/stromal scores."""
    clin = D.clinical()
    tim = D.timer2()
    rows = []
    targets = {"TMB_NONSYNONYMOUS": "TMB", "MSI_SCORE_MANTIS": "MSI (MANTIS)",
               "ANEUPLOIDY_SCORE": "Aneuploidy", "FRACTION_GENOME_ALTERED": "FGA"}
    xcell = {"immune score_XCELL": "Immune score", "stroma score_XCELL": "Stromal score",
             "microenvironment score_XCELL": "Microenvironment score"}
    for code in sorted(COHORTS):
        idx = panel.index[panel.cohort == code]
        if len(idx) < MIN_N:
            continue
        x = panel.loc[idx, "NCL"].to_numpy()
        for col, label in targets.items():
            sub = clin.reindex(idx)[col] if col in clin.columns else None
            if sub is None:
                continue
            y = pd.to_numeric(sub, errors="coerce").to_numpy()
            if np.isfinite(y).sum() < MIN_N:
                continue
            rho, p = stats.spearmanr(x, y, nan_policy="omit")
            lo, hi = S.spearman_ci(rho, int(np.isfinite(y).sum()))
            rows.append(dict(cohort=code, measure=label,
                             n=int(np.isfinite(y).sum()), rho=rho,
                             rho_lo=lo, rho_hi=hi, p=p))
        for col, label in xcell.items():
            if col not in tim.columns:
                continue
            y = pd.to_numeric(tim.reindex(idx)[col], errors="coerce").to_numpy()
            if np.isfinite(y).sum() < MIN_N:
                continue
            rho, p = stats.spearmanr(x, y, nan_policy="omit")
            lo, hi = S.spearman_ci(rho, int(np.isfinite(y).sum()))
            rows.append(dict(cohort=code, measure=label,
                             n=int(np.isfinite(y).sum()), rho=rho,
                             rho_lo=lo, rho_hi=hi, p=p))
    df = pd.DataFrame(rows)
    df["q"] = S.bh_fdr(df["p"].to_numpy())
    return df


def main():
    panel, tim = load_panel()
    print(f"tumours with expression + TIMER2: {len(panel)} across "
          f"{panel.cohort.nunique()} cohorts\n")

    infil = infiltration(panel, tim)
    infil.to_csv(os.path.join(D.TABLES, "T4_immune_infiltration.tsv"),
                 sep="\t", index=False)
    sig = (infil.q_purity_adj < 0.05).sum()
    print(f"=== Infiltration: {len(infil)} cancer x cell-type x algorithm tests, "
          f"{sig} significant at FDR<0.05 (purity-adjusted) ===")

    conc = concordance(infil)
    conc.to_csv(os.path.join(D.TABLES, "T5_algorithm_concordance.tsv"),
                sep="\t", index=False)
    print("\n=== Cross-algorithm concordance (canonical cell types) ===")
    tot = len(conc)
    print(f"  {conc.concordant.sum()}/{tot} cancer x cell-type combinations are "
          f"concordant (>=2 algorithms significant, no sign conflict)")
    print(f"  {(conc.direction=='mixed').sum()}/{tot} give conflicting directions "
          f"between algorithms")
    top = conc[conc.concordant].reindex(
        conc[conc.concordant].median_rho.abs().sort_values(ascending=False).index)
    print("\n  strongest concordant associations:")
    for _, r in top.head(12).iterrows():
        print(f"    {r.cohort:5s} {r.cell_type:28s} rho={r.median_rho:+.3f} "
              f"({r.n_significant}/{r.n_algorithms} algorithms, {r.direction})")

    cp = checkpoints(panel)
    cp.to_csv(os.path.join(D.TABLES, "T6_checkpoints.tsv"), sep="\t", index=False)
    print(f"\n=== Immune checkpoints: {len(cp)} cancer x gene tests ===")
    print(f"  raw:                 {(cp.q<0.05).sum()} significant at FDR<0.05")
    print(f"  purity-adjusted:     {(cp.q_purity_adj<0.05).sum()}")
    print(f"  proliferation-adj:   {(cp.q_prolif_adj<0.05).sum()}")
    print("\n  by gene (n cancers significant / n tested, purity-adjusted):")
    for g, grp in cp.groupby("alias"):
        s = (grp.q_purity_adj < 0.05)
        pos_ = ((grp.rho_purity_adj > 0) & s).sum()
        neg_ = ((grp.rho_purity_adj < 0) & s).sum()
        print(f"    {g:12s} {s.sum():2d}/{len(grp):2d}   "
              f"(+{pos_} / -{neg_})  median rho={grp.rho_purity_adj.median():+.3f}")

    gs = genomic_scores(panel)
    gs.to_csv(os.path.join(D.TABLES, "T8_genomic_scores.tsv"), sep="\t", index=False)
    print("\n=== NCL vs TMB / MSI / aneuploidy / immune scores ===")
    for m, grp in gs.groupby("measure"):
        s = (grp.q < 0.05)
        print(f"  {m:24s} {s.sum():2d}/{len(grp):2d} cancers significant, "
              f"median rho={grp.rho.median():+.3f} "
              f"(+{((grp.rho>0)&s).sum()} / -{((grp.rho<0)&s).sum()})")


if __name__ == "__main__":
    main()
