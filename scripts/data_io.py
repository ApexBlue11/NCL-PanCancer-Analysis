"""Shared loaders for the expression memmap, sample annotation and clinical table.

Also resolves where the raw data lives, and — more importantly — makes it
obvious when it does not.

The result tables in `results/tables/` are committed so that the figures and
every quoted number can be inspected without a 6 GB download. That convenience
has a cost: a script reading those tables produces output that *looks* like a
reproduction of the analysis when it is really a re-plot of numbers computed
elsewhere. `table()` therefore prints a provenance banner whenever it is
serving distributed tables rather than tables this checkout computed, and
`require_raw()` stops the scripts that genuinely need source data with a
message naming what is missing instead of a FileNotFoundError from inside gzip.
"""
import json, os, sys, functools
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")
for d in (TABLES, FIGURES):
    os.makedirs(d, exist_ok=True)

# The six files 01_download.py fetches. Anything downstream of 03 needs them.
RAW_INPUTS = (
    "TcgaTargetGtex_rsem_gene_tpm.gz",
    "TcgaTargetGTEX_phenotype.txt.gz",
    "gencode.v23.annotation.gene.probemap",
    "infiltration_estimation_for_tcga.csv.gz",
    "h.all.v2024.1.Hs.symbols.gmt",
    "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
)


def _candidate_roots():
    """Where a `data/` directory might reasonably be, best guess first.

    A clone that already holds the download should be used in place rather than
    downloaded again, and the download is large enough that people put it on a
    different disk. NCL_DATA wins outright; otherwise the repository's own
    data/ is preferred, then the working directory, then alongside the scripts.
    """
    env = os.environ.get("NCL_DATA")
    if env:
        yield os.path.abspath(env), "NCL_DATA"
    yield os.path.join(HERE, "..", "data"), "repository default"
    yield os.path.join(os.getcwd(), "data"), "working directory"
    yield os.path.join(HERE, "data"), "beside the scripts"
    yield os.getcwd(), "working directory itself"


def _resolve_data_root():
    """First candidate that actually holds raw inputs; else the default."""
    default = os.path.join(HERE, "..", "data")
    for root, how in _candidate_roots():
        raw = os.path.join(root, "raw")
        if os.path.isdir(raw) and any(
                os.path.exists(os.path.join(raw, f)) for f in RAW_INPUTS):
            return os.path.abspath(root), how
    return os.path.abspath(default), "repository default (not present)"


DATA_ROOT, DATA_ROOT_ORIGIN = _resolve_data_root()
RAW = os.path.join(DATA_ROOT, "raw")
PROC = os.path.join(DATA_ROOT, "proc")


def missing_raw():
    """Which of RAW_INPUTS are not on disk."""
    return [f for f in RAW_INPUTS if not os.path.exists(os.path.join(RAW, f))]


def have_raw():
    return not missing_raw()


def require_raw(what):
    """Stop with an actionable message when source data is needed but absent."""
    missing = missing_raw()
    if not missing:
        return
    sys.exit(
        "\n%s needs the raw source data, and it is not present.\n"
        "  looked in : %s  (%s)\n"
        "  missing   : %s\n"
        "\nRun `python scripts/01_download.py` first (about 6 GB), or point\n"
        "NCL_DATA at an existing copy:  NCL_DATA=/path/to/data python ...\n"
        % (what, RAW, DATA_ROOT_ORIGIN, ", ".join(missing)))


_ANNOUNCED = False


def table(name, **kw):
    """Load a result table from results/tables, stating where the numbers came from.

    Called by the scripts downstream of the analysis (figures, quoted numbers).
    Those can run from the committed tables alone, which is useful but easy to
    mistake for a reproduction, so say plainly which it is — once per process.
    """
    global _ANNOUNCED
    if not _ANNOUNCED:
        _ANNOUNCED = True
        if have_raw():
            print("  [data] source data found in %s (%s); tables in results/tables\n"
                  "         were computed by steps 01-08 of this checkout."
                  % (RAW, DATA_ROOT_ORIGIN))
        else:
            print(
                "  [data] NOTE: no raw source data found (looked in %s).\n"
                "         Falling back to the result tables distributed with the\n"
                "         repository. Output below is RE-PLOTTED from committed\n"
                "         numbers, not recomputed from TCGA/GTEx/CPTAC.\n"
                "         To reproduce from source, run scripts 01-08 first." % RAW)
    p = os.path.join(TABLES, name)
    if not os.path.exists(p):
        sys.exit("result table not found: %s\n"
                 "It is produced by an earlier pipeline step; see README." % p)
    return pd.read_csv(p, sep="\t", low_memory=True, **kw)


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
