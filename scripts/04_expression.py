"""NCL expression: tumour vs normal, and association with pathological stage.

Replaces submitted Figures 1a/1b and corrects two errors:

  1. Figure 1a was TIMER2's TCGA-only differential-expression module while
     Methods described GEPIA2 with GTEx normals. Cancers such as PAAD therefore
     rested on 4 adjacent normals. Here GTEx normals are used as the primary
     comparator (both arms Toil-reprocessed, so comparable), with TCGA adjacent
     normals as a sensitivity analysis.

  2. Section 3.1 claimed NCL "progressively increased as the cancer advanced
     through stages I-IV". The UALCAN asterisks that claim rested on compare
     each stage against normal tissue, not against the preceding stage. A
     stage-ordered trend is tested here explicitly (Jonckheere-Terpstra).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

import data_io as D
import statsutil as S
from cohorts import COHORTS, FULL_NAME

GENE = "NCL"


def differential_expression():
    ncl = D.gene_vector(GENE)
    ann, gtex = D.sample_groups()
    ix = D.index()
    pos = {s: i for i, s in enumerate(ix["samples"])}

    rows = []
    for code in sorted(COHORTS):
        _detailed, gtex_site, quality = COHORTS[code]
        t_idx = np.array([pos[s] for s in ann.index[(ann.cohort == code) &
                                                    (ann.group == "tumour")]])
        n_tcga = np.array([pos[s] for s in ann.index[(ann.cohort == code) &
                                                     (ann.group == "tcga_normal")]])
        n_gtex = np.array([pos[s] for s in gtex.get(code, [])])

        if len(t_idx) < 10:
            continue
        tum = ncl[t_idx]

        for label, nidx in (("GTEx", n_gtex), ("TCGA_adjacent", n_tcga)):
            if len(nidx) < 3:
                rows.append(dict(cohort=code, comparator=label, n_tumour=len(t_idx),
                                 n_normal=len(nidx), note="insufficient normals"))
                continue
            nor = ncl[nidx]
            u = stats.mannwhitneyu(tum, nor, alternative="two-sided")
            d, dlo, dhi = S.cliffs_delta(tum, nor)
            g, glo, ghi = S.hedges_g(tum, nor)
            rows.append(dict(
                cohort=code, comparator=label,
                n_tumour=len(t_idx), n_normal=len(nidx),
                median_tumour=float(np.median(tum)), median_normal=float(np.median(nor)),
                log2FC=float(np.median(tum) - np.median(nor)),   # already log2 scale
                p=float(u.pvalue),
                cliffs_delta=d, delta_lo=dlo, delta_hi=dhi,
                hedges_g=g, g_lo=glo, g_hi=ghi,
                tissue_match=quality, gtex_tissue=gtex_site, note=""))

    df = pd.DataFrame(rows)
    # FDR is applied within each comparator family across cancers.
    for comp in df["comparator"].dropna().unique():
        m = (df["comparator"] == comp) & df["p"].notna()
        df.loc[m, "q"] = S.bh_fdr(df.loc[m, "p"].to_numpy())
    df["full_name"] = df["cohort"].map(FULL_NAME)
    return df


def stage_association():
    ncl = D.gene_vector(GENE)
    ix = D.index()
    pos = {s: i for i, s in enumerate(ix["samples"])}
    clin = D.clinical()
    ann, _ = D.sample_groups()

    def norm_stage(v):
        if not isinstance(v, str):
            return np.nan
        v = v.upper().replace("STAGE", "").strip()
        for k, lab in (("IV", "IV"), ("III", "III"), ("II", "II"), ("I", "I")):
            if v.startswith(k):
                return lab
        return np.nan

    order = ["I", "II", "III", "IV"]
    rows = []
    for code in sorted(COHORTS):
        tum = ann.index[(ann.cohort == code) & (ann.group == "tumour")]
        sub = clin.reindex([s for s in tum if s in clin.index])
        if sub.empty:
            continue
        stage = sub["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(norm_stage)
        vals = pd.Series({s: ncl[pos[s]] for s in sub.index})
        groups, ns = [], []
        for st in order:
            v = vals[stage[vals.index] == st].dropna().to_numpy()
            groups.append(v); ns.append(len(v))
        present = [g for g in groups if len(g) >= 5]
        if len(present) < 3 or sum(ns) < 40:
            continue
        kw = stats.kruskal(*present)
        _jt, z, p_trend = S.jonckheere_terpstra(present)
        first, last = present[0], present[-1]
        d, dlo, dhi = S.cliffs_delta(last, first)
        rows.append(dict(
            cohort=code, full_name=FULL_NAME[code],
            n_I=ns[0], n_II=ns[1], n_III=ns[2], n_IV=ns[3], n_total=sum(ns),
            median_I=float(np.median(groups[0])) if ns[0] else np.nan,
            median_IV=float(np.median(groups[3])) if ns[3] else np.nan,
            kruskal_p=float(kw.pvalue),
            trend_z=z, trend_p=p_trend,
            late_vs_early_delta=d, delta_lo=dlo, delta_hi=dhi))

    df = pd.DataFrame(rows)
    if not df.empty:
        df["kruskal_q"] = S.bh_fdr(df["kruskal_p"].to_numpy())
        df["trend_q"] = S.bh_fdr(df["trend_p"].to_numpy())
    return df


def main():
    de = differential_expression()
    de.to_csv(os.path.join(D.TABLES, "T1_differential_expression.tsv"),
              sep="\t", index=False)
    g = de[(de.comparator == "GTEx") & de.q.notna()].sort_values("cliffs_delta")
    print(f"=== Tumour vs GTEx normal ({len(g)} cancers) ===")
    print(f"{'':2}{'cohort':6} {'nT':>5} {'nN':>5} {'log2FC':>7} {'delta':>7} "
          f"{'95% CI':>16} {'q':>10}  match")
    for _, r in g.iterrows():
        sig = "*" if r.q < 0.05 else " "
        print(f"{sig} {r.cohort:6} {r.n_tumour:5d} {r.n_normal:5d} {r.log2FC:7.2f} "
              f"{r.cliffs_delta:7.3f} [{r.delta_lo:6.3f},{r.delta_hi:6.3f}] "
              f"{r.q:10.2e}  {r.tissue_match}")
    up = ((g.cliffs_delta > 0) & (g.q < 0.05)).sum()
    dn = ((g.cliffs_delta < 0) & (g.q < 0.05)).sum()
    print(f"\nsignificant at FDR<0.05: {up} up, {dn} down, "
          f"{len(g)-up-dn} not significant")

    st = stage_association()
    st.to_csv(os.path.join(D.TABLES, "T2_stage_association.tsv"), sep="\t", index=False)
    print(f"\n=== Stage association ({len(st)} cancers) ===")
    print(f"{'cohort':6} {'n':>5} {'KW q':>9} {'trend z':>8} {'trend q':>9} "
          f"{'IV vs I delta':>14}")
    for _, r in st.sort_values("trend_p").iterrows():
        print(f"{r.cohort:6} {r.n_total:5d} {r.kruskal_q:9.3f} {r.trend_z:8.2f} "
              f"{r.trend_q:9.3f} {r.late_vs_early_delta:8.3f} "
              f"[{r.delta_lo:.2f},{r.delta_hi:.2f}]")
    inc = ((st.trend_q < 0.05) & (st.trend_z > 0)).sum()
    dec = ((st.trend_q < 0.05) & (st.trend_z < 0)).sum()
    print(f"\nmonotonic stage trend at FDR<0.05: {inc} increasing, {dec} decreasing, "
          f"{len(st)-inc-dec} none")


if __name__ == "__main__":
    main()
