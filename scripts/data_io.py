"""Shared loaders for the expression memmap, sample annotation and clinical table."""
import json, os, functools
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "raw")
PROC = os.path.join(HERE, "..", "data", "proc")
RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")
for d in (TABLES, FIGURES):
    os.makedirs(d, exist_ok=True)


@functools.lru_cache(maxsize=1)
def index():
    with open(os.path.join(PROC, "expr_index.json")) as fh:
        return json.load(fh)


@functools.lru_cache(maxsize=1)
def id2symbol():
    pm = pd.read_csv(os.path.join(RAW, "gencode.v23.annotation.gene.probemap"), sep="\t")
    return dict(zip(pm["id"], pm["gene"]))


@functools.lru_cache(maxsize=1)
def symbol_rows():
    """symbol -> row indices in the memmap (a symbol can map to several IDs)."""
    ix = index()
    m = id2symbol()
    out = {}
    for i, gid in enumerate(ix["genes"]):
        s = m.get(gid)
        if s:
            out.setdefault(s, []).append(i)
    return out


def expr_matrix(mode="r"):
    ix = index()
    return np.memmap(os.path.join(PROC, ix["memmap"]), dtype=ix["dtype"],
                     mode=mode, shape=tuple(ix["shape"]))


def gene_vector(symbol):
    """log2(TPM+0.001) across all samples for one gene symbol.

    Where a symbol maps to multiple gencode IDs the highest-expressing row is
    used, which is the standard collapse and avoids averaging a real transcript
    with a lowly-expressed readthrough.
    """
    rows = symbol_rows().get(symbol)
    if not rows:
        raise KeyError(f"gene symbol not found: {symbol}")
    X = expr_matrix()
    if len(rows) == 1:
        return np.asarray(X[rows[0]], dtype=float)
    cand = np.vstack([np.asarray(X[r], dtype=float) for r in rows])
    return cand[np.nanmean(cand, axis=1).argmax()]


def gene_block(symbols):
    """DataFrame of genes x samples for a list of symbols (missing ones skipped)."""
    ix = index()
    data, found = {}, []
    for s in symbols:
        try:
            data[s] = gene_vector(s)
            found.append(s)
        except KeyError:
            pass
    return pd.DataFrame(data, index=ix["samples"]).T, found


@functools.lru_cache(maxsize=1)
def phenotype():
    p = pd.read_csv(os.path.join(RAW, "TcgaTargetGTEX_phenotype.txt.gz"),
                    sep="\t", encoding="latin-1")
    return p.set_index("sample")


@functools.lru_cache(maxsize=1)
def clinical():
    c = pd.read_csv(os.path.join(PROC, "clinical_pancan.tsv"), sep="\t",
                    low_memory=False)
    return c.set_index("sample_barcode")


@functools.lru_cache(maxsize=1)
def timer2():
    t = pd.read_csv(os.path.join(RAW, "infiltration_estimation_for_tcga.csv.gz"),
                    index_col=0)
    t.index = t.index.astype(str)
    return t


def sample_groups():
    """Map each expression column to (cohort, group) with group in tumour/tcga_normal/gtex_normal."""
    from cohorts import COHORTS, TUMOUR_TYPES, TCGA_NORMAL, GTEX_NORMAL
    ph = phenotype()
    ix = index()
    samples = pd.Index(ix["samples"])
    ph = ph.reindex(samples)

    cohort = pd.Series(index=samples, dtype=object)
    group = pd.Series(index=samples, dtype=object)

    for code, (detailed, gtex_site, _q) in COHORTS.items():
        is_tcga = (ph["_study"] == "TCGA") & (ph["detailed_category"] == detailed)
        t = is_tcga & ph["_sample_type"].isin(TUMOUR_TYPES)
        n = is_tcga & (ph["_sample_type"] == TCGA_NORMAL)
        cohort[t] = code; group[t] = "tumour"
        cohort[n] = code; group[n] = "tcga_normal"

    # GTEx normals are assigned per cohort separately (one tissue can serve several).
    gtex = {}
    for code, (_d, site, _q) in COHORTS.items():
        if site is None:
            continue
        sel = ((ph["_study"] == "GTEX") & (ph["_sample_type"] == GTEX_NORMAL)
               & (ph["_primary_site"] == site))
        gtex[code] = samples[sel.fillna(False).to_numpy()]

    return pd.DataFrame({"cohort": cohort, "group": group}), gtex
