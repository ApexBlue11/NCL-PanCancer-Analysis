"""Download raw data with resume + retry. Safe to re-run; skips complete files."""
import os, sys, time, requests

RAW = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
os.makedirs(RAW, exist_ok=True)

FILES = [
    ("TcgaTargetGtex_rsem_gene_tpm.gz",
     "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz"),
    ("TcgaTargetGTEX_phenotype.txt.gz",
     "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz"),
    ("gencode.v23.annotation.gene.probemap",
     "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap"),
    ("infiltration_estimation_for_tcga.csv.gz",
     "http://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz"),
    ("h.all.v2024.1.Hs.symbols.gmt",
     "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"),
    ("c2.cp.reactome.v2024.1.Hs.symbols.gmt",
     "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt"),
]

UA = {"User-Agent": "Mozilla/5.0"}


def total_size(url):
    """Advertised byte length, or None when it cannot be trusted.

    Servers that apply transfer compression report the *encoded* length, which
    will never match the bytes we write to disk -- so treat that as unknown
    rather than as a size mismatch.
    """
    try:
        r = requests.head(url, headers=UA, timeout=60, allow_redirects=True)
        if r.headers.get("Content-Encoding"):
            return None
        n = r.headers.get("Content-Length")
        return int(n) if n else None
    except Exception:
        return None


def fetch(name, url, tries=6):
    dest = os.path.join(RAW, name)
    want = total_size(url)
    for attempt in range(1, tries + 1):
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if want and have == want:
            print(f"[ok]   {name} ({have/1e6:.1f} MB) already complete", flush=True)
            return True
        headers = dict(UA)
        mode = "wb"
        # Only resume when the server advertised a length we can compare against.
        if have and want:
            headers["Range"] = f"bytes={have}-"
            mode = "ab"
        try:
            with requests.get(url, headers=headers, stream=True, timeout=(30, 300)) as r:
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"HTTP {r.status_code}")
                if r.status_code == 200:
                    mode, have = "wb", 0  # server ignored Range; restart clean
                with open(dest, mode) as fh:
                    last = time.time()
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                        have += len(chunk)
                        if time.time() - last > 20:
                            pct = f"{100*have/want:.1f}%" if want else "?"
                            print(f"       {name}: {have/1e6:.0f} MB ({pct})", flush=True)
                            last = time.time()
            got = os.path.getsize(dest)
            if want and got != want:
                raise RuntimeError(f"size mismatch {got} != {want}")
            print(f"[done] {name} ({got/1e6:.1f} MB)", flush=True)
            return True
        except Exception as e:
            print(f"[retry {attempt}/{tries}] {name}: {type(e).__name__}: {e}", flush=True)
            time.sleep(min(60, 5 * attempt))
    print(f"[FAIL] {name}", flush=True)
    return False


if __name__ == "__main__":
    # Attempt every file; one failure must not skip the rest.
    results = [(n, fetch(n, u)) for n, u in FILES]
    bad = [n for n, ok in results if not ok]
    print("ALL DOWNLOADS COMPLETE" if not bad else f"FAILED: {', '.join(bad)}")
    sys.exit(0 if not bad else 1)
