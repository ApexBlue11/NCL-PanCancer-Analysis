"""Independent protein-level validation of NCL in CPTAC tumour/normal cohorts.

Reviewers 1 and 2 both asked for validation beyond the discovery databases.
No new wet-lab data was generated for this revision; CPTAC is the strongest
available substitute because it is an independent patient series measured on a
different platform (mass spectrometry) at the protein rather than transcript
level, with adjacent normal tissue from the same patients.

Paired tumour/normal samples are analysed with a Wilcoxon signed-rank test where
pairing exists, and rank-sum otherwise.
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats

import data_io as D
import statsutil as S

warnings.filterwarnings("ignore")
GENE = "NCL"

# CPTAC cohorts with harmonised proteomics and adjacent normals.
COHORT_CLASS = {
    "BRCA": "Brca", "CCRCC": "Ccrcc", "COAD": "Coad", "GBM": "Gbm",
    "HNSCC": "Hnscc", "LSCC": "Lscc", "LUAD": "Luad", "OV": "Ov",
    "PDAC": "Pdac", "UCEC": "Ucec",
}
SOURCES = ["umich", "bcm", "broad"]


def get_proteomics(ca):
    for src in SOURCES:
        try:
            p = ca.get_proteomics(source=src)
            if p is not None and len(p):
                return p, src
        except Exception:
            continue
    return None, None


def main():
    import cptac
    rows = []
    for code, cls in COHORT_CLASS.items():
        try:
            ca = getattr(cptac, cls)()
            p, src = get_proteomics(ca)
            if p is None:
                rows.append(dict(cohort=code, note="no proteomics")); continue
            if isinstance(p.columns, pd.MultiIndex):
                p.columns = p.columns.get_level_values(0)
            if GENE not in p.columns:
                rows.append(dict(cohort=code, source=src, note="NCL not quantified"))
                continue
            v = p[GENE]
            if isinstance(v, pd.DataFrame):
                v = v.iloc[:, 0]
            v = pd.to_numeric(v, errors="coerce")

            sid = pd.Series(p.index.astype(str), index=p.index)
            is_norm = sid.str.endswith(".N")
            if is_norm.sum() < 3:
                rows.append(dict(cohort=code, source=src, n_tumour=int((~is_norm).sum()),
                                 n_normal=int(is_norm.sum()), note="too few normals"))
                continue

            tum = v[~is_norm].dropna()
            nor = v[is_norm].dropna()

            # Pair on patient id where the same patient contributes both.
            base_n = {s[:-2]: s for s in sid[is_norm]}
            pairs = [(v[sid[sid == t].index[0]], v[base_n[t]])
                     for t in sid[~is_norm] if t in base_n]
            pairs = [(a, b) for a, b in pairs if np.isfinite(a) and np.isfinite(b)]

            rec = dict(cohort=code, source=src, n_tumour=len(tum), n_normal=len(nor),
                       n_paired=len(pairs),
                       median_tumour=float(tum.median()), median_normal=float(nor.median()))
            u = stats.mannwhitneyu(tum, nor, alternative="two-sided")
            rec["p_ranksum"] = float(u.pvalue)
            d, lo, hi = S.cliffs_delta(tum.to_numpy(), nor.to_numpy())
            rec.update(cliffs_delta=d, delta_lo=lo, delta_hi=hi)
            g, glo, ghi = S.hedges_g(tum.to_numpy(), nor.to_numpy())
            rec.update(hedges_g=g, g_lo=glo, g_hi=ghi)
            if len(pairs) >= 8:
                a = np.array([x for x, _ in pairs]); b = np.array([y for _, y in pairs])
                rec["p_paired"] = float(stats.wilcoxon(a, b).pvalue)
                rec["median_paired_diff"] = float(np.median(a - b))
            rows.append(rec)
            print(f"  {code:6s} {src:6s} T={len(tum):4d} N={len(nor):4d} "
                  f"paired={len(pairs):4d} delta={d:+.3f} p={u.pvalue:.2e}", flush=True)
        except Exception as e:
            rows.append(dict(cohort=code, note=f"{type(e).__name__}: {e}"))
            print(f"  {code:6s} ERROR {type(e).__name__}: {e}", flush=True)

    df = pd.DataFrame(rows)
    for col, q in (("p_ranksum", "q_ranksum"), ("p_paired", "q_paired")):
        if col in df:
            df[q] = S.bh_fdr(df[col].to_numpy())
    df.to_csv(os.path.join(D.TABLES, "T7_cptac_validation.tsv"), sep="\t", index=False)

    ok = df[df.get("p_ranksum").notna()] if "p_ranksum" in df else pd.DataFrame()
    print(f"\n=== CPTAC protein-level validation ({len(ok)} cohorts) ===")
    if len(ok):
        print(f"{'':2}{'cohort':7} {'nT':>4} {'nN':>4} {'delta':>7} {'95% CI':>16} "
              f"{'q(rank)':>9} {'q(paired)':>10}")
        for _, r in ok.sort_values("cliffs_delta", ascending=False).iterrows():
            sig = "*" if r.q_ranksum < 0.05 else " "
            qp = r.get("q_paired", np.nan)
            print(f"{sig} {r.cohort:7} {r.n_tumour:4.0f} {r.n_normal:4.0f} "
                  f"{r.cliffs_delta:7.3f} [{r.delta_lo:6.3f},{r.delta_hi:6.3f}] "
                  f"{r.q_ranksum:9.2e} {qp:10.2e}")
        up = ((ok.cliffs_delta > 0) & (ok.q_ranksum < 0.05)).sum()
        print(f"\nNCL protein significantly higher in tumour: {up}/{len(ok)} cohorts")


if __name__ == "__main__":
    main()
