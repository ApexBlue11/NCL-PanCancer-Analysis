"""Within-cancer Spearman correlation of NCL against every gene.

Produces data/proc/ncl_correlations.tsv.gz, the ranking statistic consumed by
07_gsea.py. Split out from the GSEA step because it is the expensive part
(a full pass over the 4.6 GB expression matrix) and only needs running once.

Correlations are computed *within* each cancer. This matters: the submitted
manuscript's Figure 4c reported pan-cancer pooled correlations, which are driven
by tissue-of-origin differences in composition rather than by any within-tumour
relationship.

Implementation note: correlations are accumulated in a single sequential pass
over the memmap. Fancy-indexing the sample columns per cohort would re-read the
whole file once per cancer -- ~150 GB of I/O instead of 4.6 GB.
"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats

import data_io as D
from cohorts import COHORTS

warnings.filterwarnings("ignore")
GENE = "NCL"
BLOCK = 2000
MIN_N = 40
OUT = os.path.join(D.PROC, "ncl_correlations.tsv.gz")


def rank_correlations():
    ix = D.index()
    n_g, _ = ix["shape"]
    X = D.expr_matrix()
    ann, _ = D.sample_groups()
    pos = {s: i for i, s in enumerate(ix["samples"])}

    cohort_idx, ncl_ranks = {}, {}
    ncl = D.gene_vector(GENE)
    for code in sorted(COHORTS):
        sel = [pos[s] for s in ann.index[(ann.cohort == code) & (ann.group == "tumour")]]
        if len(sel) < MIN_N:
            continue
        cohort_idx[code] = np.array(sel)
        r = stats.rankdata(ncl[sel])
        # Pre-centre and scale so the correlation is a single dot product later.
        ncl_ranks[code] = (r - r.mean()) / (r.std() * len(r) ** 0.5)

    codes = list(cohort_idx)
    print(f"{len(codes)} cohorts, {n_g} genes", flush=True)
    out = np.full((n_g, len(codes)), np.nan, dtype=np.float32)
    for start in range(0, n_g, BLOCK):
        stop = min(start + BLOCK, n_g)
        blk = np.asarray(X[start:stop], dtype=np.float32)
        for j, code in enumerate(codes):
            sub = blk[:, cohort_idx[code]]
            rk = stats.rankdata(sub, axis=1)
            rk -= rk.mean(axis=1, keepdims=True)
            sd = rk.std(axis=1, keepdims=True)
            sd[sd == 0] = np.nan            # constant genes -> undefined correlation
            rk /= sd * sub.shape[1] ** 0.5
            out[start:stop, j] = rk @ ncl_ranks[code]
        if (start // BLOCK) % 5 == 0:
            print(f"  {stop}/{n_g} genes", flush=True)
    return pd.DataFrame(out, index=ix["genes"], columns=codes)


def collapse_to_symbols(corr):
    """Collapse gencode IDs to gene symbols, keeping the most variable row per symbol."""
    m = D.id2symbol()
    sym = pd.Series([m.get(g) for g in corr.index], index=corr.index)
    corr = corr[sym.notna().to_numpy()]
    corr = corr.assign(_sym=sym.dropna().to_numpy())
    corr = corr.assign(_spread=corr.drop(columns="_sym").std(axis=1).to_numpy())
    corr = corr.sort_values("_spread", ascending=False)
    corr = corr[~corr["_sym"].duplicated()]
    return corr.set_index("_sym").drop(columns="_spread")


def main():
    if os.path.exists(OUT):
        print(f"{OUT} exists; nothing to do")
        return
    corr = collapse_to_symbols(rank_correlations())
    corr.to_csv(OUT, sep="\t")
    corr.to_csv(os.path.join(D.TABLES, "T9_ncl_gene_correlations.tsv.gz"), sep="\t")
    print(f"\nwrote {OUT}  {corr.shape[0]} symbols x {corr.shape[1]} cancers")
    # Sanity check: NCL must correlate perfectly with itself in every cohort.
    self_corr = corr.loc[GENE]
    assert np.allclose(self_corr.dropna(), 1.0, atol=1e-6), "NCL self-correlation != 1"
    print("sanity check passed: NCL self-correlation = 1.0 in all cohorts")


if __name__ == "__main__":
    main()
