"""Record exactly what produced the results: input checksums, versions, shapes.

Public data resources are not immutable. UCSC Xena, TIMER2.0 and MSigDB
periodically reissue files under the same URL, cBioPortal serves a live API, and
the `cptac` package downloads data at run time. Without a record of what was
actually retrieved, "reproducible" means only "the code runs", not "the code
reproduces these numbers".

This writes results/MANIFEST.json and results/MANIFEST.md containing, for every
raw input: SHA-256, byte size, modification time and source URL; for every
result table: shape and SHA-256; and the resolved version of every package that
affects a numeric result.

Re-run after the pipeline. If a future run disagrees with the published numbers,
diffing two manifests localises the cause to an input file, a package version or
the code.
"""
import os, json, hashlib, platform, sys, datetime
import importlib.metadata as md

import data_io as D

RAW_SOURCES = {
    "TcgaTargetGtex_rsem_gene_tpm.gz":
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "TcgaTargetGTEX_phenotype.txt.gz":
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz",
    "gencode.v23.annotation.gene.probemap":
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
    "infiltration_estimation_for_tcga.csv.gz":
        "http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz",
    "h.all.v2024.1.Hs.symbols.gmt":
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
    "c2.cp.reactome.v2024.1.Hs.symbols.gmt":
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
}

# Packages whose version can change a reported number.
NUMERIC_PACKAGES = ["numpy", "scipy", "pandas", "statsmodels", "lifelines",
                    "gseapy", "cptac", "matplotlib", "requests"]

# Seeds fixed in the code, restated here so the manifest is self-contained.
SEEDS = {
    "cliffs_delta_bootstrap": {"seed": 0, "n_boot": 2000,
                               "where": "statsutil.cliffs_delta"},
    "gsea_permutation": {"seed": 0, "permutation_num": 1000,
                         "where": "07_gsea.py -> gseapy.prerank"},
    "jonckheere_permutation": {"seed": 0, "n_perm": 0,
                               "where": "statsutil.jonckheere_terpstra "
                                        "(asymptotic p reported; permutation optional)"},
}


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def describe(path):
    st = os.stat(path)
    return {
        "bytes": st.st_size,
        "modified_utc": datetime.datetime.fromtimestamp(
            st.st_mtime, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sha256": sha256(path),
    }


def main():
    # A manifest whose whole job is recording what produced the results is worse
    # than useless if it records nothing: run this without the inputs on disk and
    # every checksum is replaced by "absent", silently destroying the provenance
    # of a completed run. Refuse, unless the caller says that is what they meant.
    missing = D.missing_raw()
    if missing and "--allow-missing" not in sys.argv:
        sys.exit(
            "\nrefusing to write the manifest: %d of %d raw inputs are missing.\n"
            "  looked in : %s  (%s)\n"
            "  missing   : %s\n"
            "\nWriting now would overwrite the existing checksums with 'absent'\n"
            "and lose the record of the run that produced results/.\n"
            "Fetch the inputs (`python scripts/01_download.py`), point NCL_DATA at\n"
            "an existing copy, or pass --allow-missing if a partial manifest is\n"
            "genuinely what you want.\n"
            % (len(missing), len(D.RAW_INPUTS), D.RAW, D.DATA_ROOT_ORIGIN,
               ", ".join(missing)))

    man = {
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "packages": {},
        "seeds": SEEDS,
        "raw_inputs": {},
        "live_apis": {
            "cBioPortal": {
                "endpoint": "https://www.cbioportal.org/api",
                "studies": "*_tcga_pan_can_atlas_2018",
                "note": "Live API with no version pin. Survival endpoints are the "
                        "TCGA-CDR definitions; clinical values are stable in "
                        "practice but the service is not archival. The retrieved "
                        "table is checksummed below as clinical_pancan.tsv.",
            },
            "CPTAC": {
                "accessed_via": "cptac Python package",
                "note": "Datasets are downloaded at run time by the package; the "
                        "package version below determines the release retrieved.",
            },
        },
        "derived": {},
        "results": {},
    }

    for pkg in NUMERIC_PACKAGES:
        try:
            man["packages"][pkg] = md.version(pkg)
        except Exception:
            man["packages"][pkg] = None

    for name, url in RAW_SOURCES.items():
        path = os.path.join(D.RAW, name)
        if os.path.exists(path):
            man["raw_inputs"][name] = dict(describe(path), url=url)
        else:
            man["raw_inputs"][name] = {"url": url, "status": "absent"}

    for name in ["clinical_pancan.tsv", "expr_index.json",
                 "ncl_correlations.tsv.gz", "expressed_fraction.tsv.gz"]:
        path = os.path.join(D.PROC, name)
        if os.path.exists(path):
            man["derived"][name] = describe(path)

    import pandas as pd
    for fn in sorted(os.listdir(D.TABLES)):
        path = os.path.join(D.TABLES, fn)
        if not os.path.isfile(path):
            continue
        entry = describe(path)
        try:
            df = pd.read_csv(path, sep="\t", low_memory=False)
            entry["rows"], entry["cols"] = int(df.shape[0]), int(df.shape[1])
        except Exception:
            pass
        man["results"][fn] = entry

    dest = os.path.join(D.RESULTS, "MANIFEST.json")
    with open(dest, "w") as fh:
        json.dump(man, fh, indent=2, sort_keys=True)

    lines = ["# Run manifest", "",
             "Generated %s on %s, Python %s."
             % (man["generated_utc"], man["platform"], man["python"]), "",
             "Regenerate with `python scripts/11_manifest.py`. Diff two manifests "
             "to localise any disagreement in results to an input file, a package "
             "version or the code.", "",
             "## Package versions", "", "| package | version |", "|---|---|"]
    for k, v in man["packages"].items():
        lines.append("| %s | %s |" % (k, v or "not installed"))
    lines += ["", "## Random seeds", "", "| procedure | seed | draws | location |",
              "|---|---|---|---|"]
    for k, v in SEEDS.items():
        lines.append("| %s | %s | %s | `%s` |"
                     % (k, v["seed"], v.get("n_boot") or v.get("permutation_num")
                        or v.get("n_perm"), v["where"]))
    lines += ["", "## Raw inputs", "", "| file | bytes | sha256 (first 16) |",
              "|---|---|---|"]
    for k, v in man["raw_inputs"].items():
        lines.append("| `%s` | %s | `%s` |"
                     % (k, format(v["bytes"], ",") if "bytes" in v else "absent",
                        v.get("sha256", "-")[:16]))
    lines += ["", "## Result tables", "", "| table | rows | cols | sha256 (first 16) |",
              "|---|---|---|---|"]
    for k, v in man["results"].items():
        lines.append("| `%s` | %s | %s | `%s` |"
                     % (k, v.get("rows", "-"), v.get("cols", "-"),
                        v.get("sha256", "-")[:16]))
    lines += ["", "## Live services (not checksummable at source)", ""]
    for k, v in man["live_apis"].items():
        lines.append("- **%s** — %s" % (k, v["note"]))
    lines.append("")

    with open(os.path.join(D.RESULTS, "MANIFEST.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    print("wrote %s" % dest)
    print("wrote %s" % os.path.join(D.RESULTS, "MANIFEST.md"))
    print("\n%d raw inputs, %d derived, %d result tables"
          % (len(man["raw_inputs"]), len(man["derived"]), len(man["results"])))


if __name__ == "__main__":
    main()
