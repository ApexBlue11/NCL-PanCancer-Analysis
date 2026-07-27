"""Emit every number quoted in the manuscript, traced to its source table.

Written so that no value in the text is typed from memory: results/
MANUSCRIPT_NUMBERS.txt is regenerated from the result tables and the manuscript
quotes it. Any change to an analysis changes this file, and the discrepancy is
then visible rather than silently stale.
"""
import os
import numpy as np
import pandas as pd

import data_io as D

T = D.TABLES
OUT = os.path.join(D.RESULTS, "MANUSCRIPT_NUMBERS.txt")
lines = []


def w(s=""):
    lines.append(str(s))


def load(name):
    p = os.path.join(T, name)
    if not os.path.exists(p):
        return None
    return pd.read_csv(p, sep="\t", low_memory=False)


def sec(title):
    w(); w("=" * 78); w(title); w("=" * 78)


# ---------------------------------------------------------------- cohort sizes
sec("1. COHORT COMPOSITION")
ann, gtex = D.sample_groups()
tum = ann[ann.group == "tumour"].cohort.value_counts().sort_index()
nor = ann[ann.group == "tcga_normal"].cohort.value_counts()
w(f"TCGA tumours analysed: {int(tum.sum())} across {len(tum)} cohorts")
w(f"TCGA adjacent normals: {int(nor.sum())}")
w(f"GTEx normals mapped:   {sum(len(v) for v in gtex.values())}")
w()
w(f"{'cohort':8}{'tumour':>8}{'TCGA_N':>8}{'GTEx_N':>8}")
for c in tum.index:
    w(f"{c:8}{tum[c]:>8}{int(nor.get(c,0)):>8}{len(gtex.get(c,[])):>8}")

# ---------------------------------------------------------- differential expr
de = load("T1_differential_expression.tsv")
if de is not None:
    sec("2. DIFFERENTIAL EXPRESSION (Table 1 / Figure 1)")
    for comp in ["GTEx", "TCGA_adjacent"]:
        s = de[(de.comparator == comp) & de.q.notna()]
        up = ((s.cliffs_delta > 0) & (s.q < 0.05)).sum()
        dn = ((s.cliffs_delta < 0) & (s.q < 0.05)).sum()
        w(f"{comp}: {len(s)} cancers | up {up} | down {dn} | ns {len(s)-up-dn}")
    g = de[(de.comparator == "GTEx") & de.q.notna()].sort_values("cliffs_delta",
                                                                ascending=False)
    w()
    w("Largest effects vs GTEx (Cliff's delta [95% CI], q):")
    for _, r in g.head(8).iterrows():
        w(f"  {r.cohort:6} d={r.cliffs_delta:+.3f} [{r.delta_lo:+.3f},{r.delta_hi:+.3f}] "
          f"q={r.q:.2e}  nT={r.n_tumour:.0f} nN={r.n_normal:.0f}")
    w()
    w("Significantly DOWN vs GTEx:")
    for _, r in g[(g.cliffs_delta < 0) & (g.q < 0.05)].iterrows():
        w(f"  {r.cohort:6} d={r.cliffs_delta:+.3f} [{r.delta_lo:+.3f},{r.delta_hi:+.3f}] "
          f"q={r.q:.2e}")
    # comparator agreement
    a = de[de.comparator == "GTEx"].set_index("cohort")[["cliffs_delta", "q"]]
    b = de[de.comparator == "TCGA_adjacent"].set_index("cohort")[["cliffs_delta", "q"]]
    m = a.join(b, lsuffix="_g", rsuffix="_t", how="inner").dropna()
    both = m[(m.q_g < 0.05) & (m.q_t < 0.05)]
    opp = both[np.sign(both.cliffs_delta_g) != np.sign(both.cliffs_delta_t)]
    w()
    w(f"Significant in BOTH comparators: {len(both)}; concordant direction "
      f"{len(both)-len(opp)}; OPPOSITE {len(opp)} -> {opp.index.tolist()}")
    w("Comparator-discordant detail (GTEx delta vs adjacent-normal delta):")
    for c in ["KIRP", "KICH", "KIRC", "THCA"]:
        if c in m.index:
            r = m.loc[c]
            w(f"  {c:6} GTEx {r.cliffs_delta_g:+.3f} (q={r.q_g:.1e})   "
              f"adjacent {r.cliffs_delta_t:+.3f} (q={r.q_t:.1e})")
    w()
    w("Claims in the submitted manuscript that these results contradict:")
    for c in ["KIRP", "READ", "SKCM"]:
        if c in g.index.to_list() or c in set(g.cohort):
            r = g[g.cohort == c].iloc[0]
            w(f"  submitted: {c} 'lower than normal (P<0.01)'  -> here d={r.cliffs_delta:+.3f}, "
              f"q={r.q:.2e} ({'UP' if r.cliffs_delta>0 else 'DOWN'})")

# ------------------------------------------------------------------ stage
st = load("T2_stage_association.tsv")
if st is not None:
    sec("3. PATHOLOGICAL STAGE (Figure 1c)")
    inc = st[(st.trend_q < 0.05) & (st.trend_z > 0)]
    dec = st[(st.trend_q < 0.05) & (st.trend_z < 0)]
    w(f"{len(st)} cancers tested (Jonckheere-Terpstra, BH-FDR)")
    w(f"  monotonic INCREASE: {len(inc)} -> {inc.cohort.tolist()}")
    w(f"  monotonic DECREASE: {len(dec)} -> {dec.cohort.tolist()}")
    w(f"  no trend:           {len(st)-len(inc)-len(dec)}")
    w()
    for _, r in st.sort_values("trend_p").head(6).iterrows():
        w(f"  {r.cohort:6} z={r.trend_z:+.2f} q={r.trend_q:.3f} n={r.n_total:.0f} "
          f"(I={r.n_I:.0f} II={r.n_II:.0f} III={r.n_III:.0f} IV={r.n_IV:.0f})")
    w()
    named = ["BRCA", "COAD", "ESCA", "HNSC", "KICH", "LIHC", "STAD", "READ"]
    have = st[st.cohort.isin(named)]
    ok = have[(have.trend_q < 0.05) & (have.trend_z > 0)]
    w(f"Submitted manuscript claimed a stage I-IV increase in: {', '.join(named)}")
    w(f"  supported here in {len(ok)}/{len(have)} of those testable: {ok.cohort.tolist()}")

# --------------------------------------------------------------- survival
sv = load("T3_survival.tsv")
if sv is not None:
    sec("4. SURVIVAL (Table 2 / Figure 2)")
    for ep in ["OS", "DSS", "PFS"]:
        r = sv[(sv.endpoint == ep) & sv.adj_HR.notna()]
        u = sv[(sv.endpoint == ep) & sv.uni_q.notna()]
        sig = r[r.adj_q < 0.05]
        w(f"{ep}: {len(r)} cancers modelled | univariate sig {int((u.uni_q<0.05).sum())} "
          f"| ADJUSTED sig {len(sig)} -> {sig.cohort.tolist()}")
        for _, x in sig.iterrows():
            w(f"    {x.cohort} HR={x.adj_HR:.2f} [{x.adj_HR_lo:.2f},{x.adj_HR_hi:.2f}] "
              f"q={x.adj_q:.4f} n={x.n_adj:.0f} events={x.events_adj:.0f} "
              f"covars={x.adj_covariates}")
        if "ph_p_NCL" in r:
            ph = r[r.ph_p_NCL.notna()]
            bad = ph[ph.ph_p_NCL < 0.05]
            w(f"    PH tested {len(ph)}, violated {len(bad)} -> {bad.cohort.tolist()}")
    w()
    w("Submitted manuscript's headline OS claims, re-tested with adjustment:")
    for c, p in [("HNSC", "0.0027"), ("KICH", "0.017"), ("KIRP", "0.034"), ("LIHC", "1.4e-6")]:
        r = sv[(sv.endpoint == "OS") & (sv.cohort == c)]
        if len(r) and pd.notna(r.iloc[0].get("adj_HR")):
            x = r.iloc[0]
            w(f"  {c}: submitted P={p} -> adjusted HR={x.adj_HR:.2f} "
              f"[{x.adj_HR_lo:.2f},{x.adj_HR_hi:.2f}] q={x.adj_q:.3f} "
              f"({'RETAINED' if x.adj_q<0.05 else 'NOT SIGNIFICANT'}); "
              f"univariate q={x.uni_q:.4f}")
        else:
            w(f"  {c}: not modellable (insufficient events)")

# ---------------------------------------------------------------- immune
inf = load("T4_immune_infiltration.tsv")
conc = load("T5_algorithm_concordance.tsv")
if inf is not None and conc is not None:
    sec("5. IMMUNE INFILTRATION AND ALGORITHM CONCORDANCE (Figure 3)")
    w(f"{len(inf)} cancer x cell-type x algorithm tests; "
      f"{int((inf.q_purity_adj<0.05).sum())} significant at FDR<0.05 (purity-adjusted)")
    w(f"algorithms: {sorted(inf.algorithm.unique())}")
    w(f"concordance set: {len(conc)} cancer x cell-type combinations")
    w(f"  concordant (>=2 algorithms, no sign conflict): {int(conc.concordant.sum())}")
    w(f"  DIRECTION CONFLICT between algorithms:         {int((conc.direction=='mixed').sum())}")
    w()
    top = conc[conc.concordant].reindex(
        conc[conc.concordant].median_rho.abs().sort_values(ascending=False).index)
    w("Strongest concordant associations:")
    for _, r in top.head(10).iterrows():
        w(f"  {r.cohort:6} {r.cell_type:28} rho={r.median_rho:+.3f} "
          f"({r.n_significant}/{r.n_algorithms} algorithms, {r.direction})")
    thca = conc[(conc.cohort == "THCA") & (conc.cell_type == "Neutrophil")]
    if len(thca):
        r = thca.iloc[0]
        w()
        w(f"THCA neutrophil (submitted Fig 3b reported rho=0.539 by MCPcounter alone): "
          f"median rho={r.median_rho:+.3f} across {r.n_significant}/{r.n_algorithms} "
          f"algorithms after purity adjustment")

# ------------------------------------------------------------- checkpoints
cp = load("T6_checkpoints.tsv")
if cp is not None:
    sec("6. IMMUNE CHECKPOINTS (Figure 4) -- the analysis the submitted title claimed")
    w(f"{len(cp)} cancer x gene tests across {cp.cohort.nunique()} cancers, "
      f"{cp.gene.nunique()} genes")
    w(f"significant at FDR<0.05: raw {int((cp.q<0.05).sum())}, "
      f"purity-adjusted {int((cp.q_purity_adj<0.05).sum())}, "
      f"proliferation-adjusted {int((cp.q_prolif_adj<0.05).sum())}")
    w()
    w(f"{'molecule':12}{'median rho':>11}{'sig raw':>9}{'sig both adj':>14}{'dir(+/-)':>12}")
    rowsum = []
    for a, gp in cp.groupby("alias"):
        both = ((gp.q_purity_adj < 0.05) & (gp.q_prolif_adj < 0.05))
        pos = int(((gp.rho_prolif_adj > 0) & both).sum())
        neg = int(((gp.rho_prolif_adj < 0) & both).sum())
        rowsum.append((a, gp.rho.median(), int((gp.q < 0.05).sum()), int(both.sum()),
                       pos, neg, len(gp)))
    for a, med, sr, sb, pos, neg, n in sorted(rowsum, key=lambda x: -x[3]):
        w(f"{a:12}{med:>+11.3f}{sr:>6}/{n:<3}{sb:>11}/{n:<3}   +{pos}/-{neg}")
    w()
    w("Molecules the submitted abstract specifically claimed:")
    for a in ["PD-L1", "CTLA-4", "TIM-3", "IL-10", "TGF-beta1"]:
        gp = cp[cp.alias == a]
        both = ((gp.q_purity_adj < 0.05) & (gp.q_prolif_adj < 0.05)).sum()
        w(f"  {a:10} median rho={gp.rho.median():+.3f}, robust in {int(both)}/{len(gp)} cancers")
    w()
    b = cp[cp.alias == "B7-H3"].sort_values("rho", ascending=False)
    w(f"B7-H3 (CD276): robust in "
      f"{int(((b.q_purity_adj<0.05)&(b.q_prolif_adj<0.05)&(b.rho_prolif_adj>0)).sum())}"
      f"/{len(b)} cancers. Top:")
    for _, r in b.head(8).iterrows():
        w(f"  {r.cohort:6} rho={r.rho:+.3f} purity-adj={r.rho_purity_adj:+.3f} "
          f"prolif-adj={r.rho_prolif_adj:+.3f} q={r.q_prolif_adj:.1e} n={r.n:.0f}")

# ------------------------------------------------------------ genomic scores
gs = load("T8_genomic_scores.tsv")
if gs is not None:
    sec("7. IMMUNE / STROMAL SCORES, TMB, MSI")
    for m, gp in gs.groupby("measure"):
        s = gp.q < 0.05
        w(f"{m:26} sig {int(s.sum()):2}/{len(gp):2}  median rho={gp.rho.median():+.3f}  "
          f"(+{int(((gp.rho>0)&s).sum())}/-{int(((gp.rho<0)&s).sum())})")
    w()
    im = gs[(gs.measure == "Immune score") & (gs.q < 0.05)].sort_values("rho")
    w("Immune score, strongest negative:")
    for _, r in im.head(6).iterrows():
        w(f"  {r.cohort:6} rho={r.rho:+.3f} [{r.rho_lo:+.3f},{r.rho_hi:+.3f}] "
          f"n={r.n:.0f} q={r.q:.1e}")
    w(f"  positive: {im[im.rho>0].cohort.tolist()}")

# ------------------------------------------------------------------ CPTAC
cv = load("T7_cptac_validation.tsv")
if cv is not None and "p_ranksum" in cv:
    sec("8. CPTAC PROTEIN-LEVEL VALIDATION (Figure 5)")
    ok = cv[cv.p_ranksum.notna()]
    up = ((ok.cliffs_delta > 0) & (ok.q_ranksum < 0.05)).sum()
    w(f"{len(ok)} CPTAC cohorts; NCL protein significantly higher in tumour in {up}")
    for _, r in ok.sort_values("cliffs_delta", ascending=False).iterrows():
        tag = "UP" if r.cliffs_delta > 0 and r.q_ranksum < 0.05 else (
            "DOWN" if r.cliffs_delta < 0 and r.q_ranksum < 0.05 else "ns")
        w(f"  {r.cohort:7} d={r.cliffs_delta:+.3f} [{r.delta_lo:+.3f},{r.delta_hi:+.3f}] "
          f"nT={r.n_tumour:.0f} nN={r.n_normal:.0f} paired={r.n_paired:.0f} "
          f"q={r.q_ranksum:.2e} {tag}")

# ------------------------------------------------------------------- GSEA
gsea = load("T10_gsea_per_cancer.tsv.gz")
if gsea is not None:
    sec("9. GSEA (Figure 6)")
    qcol = next((c for c in ["FDR_q-val", "FDR_qval"] if c in gsea.columns), None)
    sig = gsea[(gsea[qcol] < 0.05) & gsea.NES.notna()]
    w(f"{len(gsea)} tests, {len(sig)} significant at FDR<0.05, "
      f"{gsea.cohort.nunique()} cancers")
    for coll in gsea.collection.unique():
        s = sig[sig.collection == coll]
        n_c = gsea[gsea.collection == coll].cohort.nunique()
        agg = (s.groupby("Term").agg(n_sig=("cohort", "nunique"),
                                     n_up=("NES", lambda v: (v > 0).sum()),
                                     n_dn=("NES", lambda v: (v < 0).sum()),
                                     med=("NES", "median"))
                 .sort_values("n_sig", ascending=False))
        cons = agg[(agg.n_up == 0) | (agg.n_dn == 0)]
        w(); w(f"--- {coll}: direction-consistent across cancers (of {n_c}) ---")
        for t, r in cons.head(20).iterrows():
            w(f"  {r.n_sig:2.0f}/{n_c}  {'UP' if r.n_up>0 else 'DOWN':4} "
              f"NES={r.med:+.2f}  {t[:70]}")
else:
    sec("9. GSEA -- NOT YET AVAILABLE (run 07_gsea.py)")

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\n\nwrote {OUT}")
