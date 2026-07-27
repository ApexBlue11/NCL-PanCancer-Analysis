"""Probe every planned data source with a ranged GET before committing to downloads."""
import urllib.request as U

SOURCES = {
    "xena_tpm":       "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "xena_pheno":     "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz",
    "xena_probemap":  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
    "tcga_cdr":       "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
    "timer2_infil":   "http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz",
    "thorsson_immune":"https://api.gdc.cancer.gov/data/1c6174d9-8ffb-466e-b5ee-07b204c15cf8",
    "hallmark_gmt":   "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
    "reactome_gmt":   "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
}

for name, url in SOURCES.items():
    try:
        req = U.Request(url, headers={"User-Agent": "Mozilla/5.0", "Range": "bytes=0-255"})
        with U.urlopen(req, timeout=60) as r:
            size = r.headers.get("Content-Range") or r.headers.get("Content-Length")
            head = r.read(64)
            print(f"{name:16s} OK   {r.status:3d}  size={size}  magic={head[:8]!r}")
    except Exception as e:
        print(f"{name:16s} FAIL      {type(e).__name__}: {e}")
