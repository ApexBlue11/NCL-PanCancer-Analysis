"""Per-cancer GSEA on genes ranked by their within-cancer correlation with NCL.

Replaces submitted Figure 4 and addresses Reviewer 2 comment 7 (GSEA with
Hallmark and Reactome rather than a single DAVID KEGG bar chart).

The submitted enrichment analysis had several problems this replaces outright:
Methods named miRDB while Results and the legend named miRNet; Figure 4b was a
miRNA-target network in which NCL was one node among many, not a set of
"NCL-interacting miRNAs"; the named miRNAs were assigned roles in EMT,
metastasis and apoptosis that were never tested; and Figure 4c reported
pan-cancer pooled correlations, which track tissue-of-origin composition rather
than regulation, with r = -0.19 described as a "moderately negative correlation".

Correlations here are computed *within* each cancer, so tissue composition
cannot drive them, and enrichment is run per cancer and then aggregated.

Run characteristics
-------------------
* The expensive step -- Spearman correlation of NCL against every gene in every
  cohort -- is cached by 07_gsea_correlations (data/proc/ncl_correlations.tsv.gz),
  so this script never reads the 4.6 GB expression matrix.
* Each (collection, cohort) unit is written to its own file under
  results/gsea_parts/ and skipped if already present, so an interrupted run
  resumes rather than restarting.
* Units are spread over a small process pool. The pool is deliberately small:
  this machine runs other workloads and has little free RAM.
"""
import os, sys, glob, time, argparse, warnings, logging
import numpy as np
import pandas as pd

import data_io as D
import statsutil as S

warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

PARTS = os.path.join(D.RESULTS, "gsea_parts")
os.makedirs(PARTS, exist_ok=True)
CORR = os.path.join(D.PROC, "ncl_correlations.tsv.gz")
EXPRESSED = os.path.join(D.PROC, "expressed_fraction.tsv.gz")
MIN_EXPRESSED_FRACTION = 0.25

GMTS = {"Hallmark": "h.all.v2024.1.Hs.symbols.gmt",
        "Reactome": "c2.cp.reactome.v2024.1.Hs.symbols.gmt"}


def ranked_list(corr, expressed, cohort):
    """Correlation-ranked genes for one cohort, restricted to expressed genes."""
    rnk = corr[cohort].dropna()
    if expressed is not None and cohort in expressed.columns:
        keep = expressed[cohort].reindex(rnk.index) >= MIN_EXPRESSED_FRACTION
        rnk = rnk[keep.fillna(False).to_numpy()]
    return rnk.sort_values(ascending=False)


def run_unit(args):
    """One (collection, cohort) GSEA. Returns a short status string."""
    coll, cohort, perms = args
    dest = os.path.join(PARTS, f"{coll}__{cohort}.tsv")
    if os.path.exists(dest):
        return f"skip {coll}/{cohort}"
    import gseapy
    corr = pd.read_csv(CORR, sep="\t", index_col=0)
    expressed = pd.read_csv(EXPRESSED, sep="\t", index_col=0) \
        if os.path.exists(EXPRESSED) else None
    rnk = ranked_list(corr, expressed, cohort)
    t0 = time.time()
    try:
        pre = gseapy.prerank(rnk=rnk, gene_sets=os.path.join(D.RAW, GMTS[coll]),
                             permutation_num=perms, min_size=15, max_size=500,
                             threads=1, seed=0, no_plot=True, verbose=False,
                             outdir=None)
        res = pre.res2d.copy()
        res["cohort"] = cohort
        res["collection"] = coll
        res["n_genes_ranked"] = len(rnk)
        # Write atomically so an interrupted write is never mistaken for a result.
        tmp = dest + ".tmp"
        res.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, dest)
        return f"done {coll}/{cohort} {len(res)} sets {time.time()-t0:.0f}s n={len(rnk)}"
    except Exception as e:
        return f"FAIL {coll}/{cohort} {type(e).__name__}: {e}"


def aggregate():
    files = sorted(glob.glob(os.path.join(PARTS, "*.tsv")))
    if not files:
        print("no parts to aggregate")
        return None
    res = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    res.columns = [c.replace(" ", "_") for c in res.columns]
    qcol = next((c for c in ["FDR_q-val", "FDR_qval", "fdr"] if c in res.columns), None)
    for c in ["NES", "ES", qcol, "NOM_p-val"]:
        if c and c in res.columns:
            res[c] = pd.to_numeric(res[c], errors="coerce")
    dest = os.path.join(D.TABLES, "T10_gsea_per_cancer.tsv.gz")
    res.to_csv(dest, sep="\t", index=False)
    print(f"\nwrote {dest}  ({len(res)} rows, {res.cohort.nunique()} cohorts)")

    sig = res[(res[qcol] < 0.05) & res["NES"].notna()]
    print(f"{len(sig)} of {len(res)} tests significant at FDR<0.05")
    for coll in res.collection.unique():
        s = sig[sig.collection == coll]
        n_cancers = res[res.collection == coll].cohort.nunique()
        agg = (s.groupby("Term")
                 .agg(n_sig=("cohort", "nunique"),
                      n_up=("NES", lambda v: (v > 0).sum()),
                      n_dn=("NES", lambda v: (v < 0).sum()),
                      median_NES=("NES", "median"))
                 .sort_values("n_sig", ascending=False))
        cons = agg[(agg.n_up == 0) | (agg.n_dn == 0)]
        print(f"\n--- {coll}: direction-consistent pathways (of {n_cancers} cancers) ---")
        for t, r in cons.head(18).iterrows():
            d = "UP" if r.n_up > 0 else "DOWN"
            print(f"  {int(r.n_sig):2d}/{n_cancers}  {d:4s} "
                  f"NES={r.median_NES:+.2f}  {t[:76]}")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--aggregate-only", action="store_true")
    a = ap.parse_args()

    if a.aggregate_only:
        aggregate(); return

    corr = pd.read_csv(CORR, sep="\t", index_col=0, nrows=2)
    cohorts = list(corr.columns)
    units = [(coll, c, a.perms) for coll in GMTS for c in cohorts
             if not os.path.exists(os.path.join(PARTS, f"{coll}__{c}.tsv"))]
    total = len(GMTS) * len(cohorts)
    print(f"{total - len(units)}/{total} units already complete; "
          f"{len(units)} to run on {a.workers} workers", flush=True)
    if not units:
        aggregate(); return

    # Cheap units first so partial progress is maximally useful.
    units.sort(key=lambda u: (u[0] != "Hallmark", u[1]))
    t0 = time.time()
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, msg in enumerate(ex.map(run_unit, units), 1):
            print(f"[{i}/{len(units)}] {msg}  (elapsed {time.time()-t0:.0f}s)",
                  flush=True)
    aggregate()


if __name__ == "__main__":
    main()
