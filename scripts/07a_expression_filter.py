"""Per-cohort expressed-gene filter for the GSEA ranked lists.

The first GSEA run warned that ~10% of genes had tied ranking statistics. Those
ties come from genes that are essentially unexpressed in a given cohort: with
few distinct values their correlation with NCL is degenerate, and GSEA then
orders them arbitrarily. Ranking genes that are not expressed is not meaningful,
so each cohort's ranked list is restricted to genes detected in that cohort.

A gene is called expressed in a cohort when TPM > 1 -- i.e. log2(TPM + 0.001) > 0
in these units -- in at least 25% of that cohort's tumours.

Small blocks are used deliberately: this machine has little free RAM and large
blocks cause page thrashing against the 4.6 GB matrix.
"""
import os
import numpy as np
import pandas as pd

import data_io as D
from cohorts import COHORTS

BLOCK = 1000
MIN_FRACTION = 0.25
EXPRESSED_LOG2TPM = 0.0          # log2(TPM + 0.001) > 0  <=>  TPM > ~1
MIN_N = 40
OUT = os.path.join(D.PROC, "expressed_fraction.tsv.gz")


def main():
    if os.path.exists(OUT):
        print(f"{OUT} exists; nothing to do")
        return

    ix = D.index()
    n_g, _ = ix["shape"]
    X = D.expr_matrix()
    ann, _ = D.sample_groups()
    pos = {s: i for i, s in enumerate(ix["samples"])}

    cohort_idx = {}
    for code in sorted(COHORTS):
        sel = [pos[s] for s in ann.index[(ann.cohort == code) & (ann.group == "tumour")]]
        if len(sel) >= MIN_N:
            cohort_idx[code] = np.array(sel)
    codes = list(cohort_idx)
    print(f"{len(codes)} cohorts, {n_g} genes", flush=True)

    frac = np.zeros((n_g, len(codes)), dtype=np.float32)
    for start in range(0, n_g, BLOCK):
        stop = min(start + BLOCK, n_g)
        blk = np.asarray(X[start:stop], dtype=np.float32)
        for j, code in enumerate(codes):
            sub = blk[:, cohort_idx[code]]
            frac[start:stop, j] = (sub > EXPRESSED_LOG2TPM).mean(axis=1)
        if (start // BLOCK) % 10 == 0:
            print(f"  {stop}/{n_g}", flush=True)

    df = pd.DataFrame(frac, index=ix["genes"], columns=codes)
    # Collapse gencode IDs to symbols the same way the correlation matrix did,
    # so the two line up row-for-row.
    m = D.id2symbol()
    sym = pd.Series([m.get(g) for g in df.index], index=df.index)
    df = df[sym.notna().to_numpy()]
    df.index = sym.dropna().to_numpy()
    df = df.groupby(level=0).max()          # symbol counts as expressed if any ID is
    df.to_csv(OUT, sep="\t")

    keep = (df >= MIN_FRACTION)
    print(f"\nwrote {OUT}  shape={df.shape}")
    print(f"genes retained per cohort (>= {MIN_FRACTION:.0%} of tumours with TPM>1):")
    for c in df.columns:
        print(f"  {c:6s} {int(keep[c].sum()):6d} / {len(df)}")


if __name__ == "__main__":
    main()
