"""Stream the Xena Toil TPM matrix into a flat float32 array plus an index file.

The source is log2(TPM+0.001), genes x samples, uniformly reprocessed across
TCGA/TARGET/GTEx (Vivian et al. 2017), which is what makes tumour-vs-normal
comparisons against GTEx defensible -- unlike mixing separately processed
TCGA and GTEx pipelines.

Written with sequential file I/O rather than a writable memmap: on a machine
with little free RAM a 4.6 GB writable mapping accumulates dirty pages faster
than they can be flushed and the build becomes I/O-bound. Reading it back as a
read-only memmap later is fine, since those pages are clean and evictable.
"""
import os, gzip, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "raw")
PROC = os.path.join(HERE, "..", "data", "proc")
os.makedirs(PROC, exist_ok=True)

SRC = os.path.join(RAW, "TcgaTargetGtex_rsem_gene_tpm.gz")
OUT = os.path.join(PROC, "expr.f32.memmap")


def main():
    genes = []
    t0 = time.time()
    with gzip.open(SRC, "rt") as fh, open(OUT, "wb", buffering=1 << 22) as out:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        n_s = len(samples)
        print(f"samples: {n_s}", flush=True)
        for i, line in enumerate(fh):
            tab = line.index("\t")
            genes.append(line[:tab])
            row = np.array(line[tab + 1:].split("\t"), dtype=np.float32)
            if row.size != n_s:
                raise ValueError(f"row {i} ({genes[-1]}) has {row.size} values, "
                                 f"expected {n_s}")
            out.write(row.tobytes())
            if i and i % 5000 == 0:
                print(f"  {i} genes  ({time.time()-t0:.0f}s)", flush=True)

    n_g = len(genes)
    size = os.path.getsize(OUT)
    assert size == n_g * n_s * 4, f"size {size} != {n_g}*{n_s}*4"
    with open(os.path.join(PROC, "expr_index.json"), "w") as fh:
        json.dump({"genes": genes, "samples": samples,
                   "shape": [n_g, n_s], "dtype": "float32",
                   "units": "log2(TPM+0.001)", "memmap": os.path.basename(OUT)}, fh)
    print(f"\nwrote {OUT} ({size/1e9:.2f} GB)")
    print(f"{n_g} genes x {n_s} samples in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
