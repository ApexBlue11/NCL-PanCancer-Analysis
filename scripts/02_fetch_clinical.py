"""Pull clinical, survival, TMB and MSI for all TCGA PanCancer Atlas studies from cBioPortal.

Produces data/proc/clinical_pancan.tsv, one row per sample, wide format.
Survival endpoints here are the curated TCGA-CDR values (Liu et al., Cell 2018),
which is what the revision needs for multivariable Cox.
"""
import os, time, sys
import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
PROC = os.path.join(HERE, "..", "data", "proc")
os.makedirs(PROC, exist_ok=True)

API = "https://www.cbioportal.org/api"
S = requests.Session()
S.headers.update({"User-Agent": "Mozilla/5.0"})

PATIENT_ATTRS = [
    "OS_MONTHS", "OS_STATUS", "DSS_MONTHS", "DSS_STATUS",
    "PFS_MONTHS", "PFS_STATUS", "DFS_MONTHS", "DFS_STATUS",
    "AGE", "SEX", "RACE", "AJCC_PATHOLOGIC_TUMOR_STAGE",
    "SUBTYPE", "HISTORY_NEOADJUVANT_TRTYN", "RADIATION_THERAPY",
]
# GRADE, TMB and MSI are sample-level attributes in cBioPortal, not patient-level.
SAMPLE_ATTRS = [
    "SAMPLE_TYPE", "GRADE", "TMB_NONSYNONYMOUS",
    "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE", "MUTATION_COUNT",
    "FRACTION_GENOME_ALTERED", "ANEUPLOIDY_SCORE",
    "PATH_T_STAGE", "PATH_N_STAGE", "PATH_M_STAGE",
]


def get(url, params, tries=5):
    for a in range(1, tries + 1):
        try:
            r = S.get(url, params=params, timeout=(30, 300))
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f"    retry {a}/{tries}: {type(e).__name__}", flush=True)
            time.sleep(5 * a)
    raise RuntimeError(f"failed: {url}")


def pull(study, kind):
    rows = get(f"{API}/studies/{study}/clinical-data",
               {"clinicalDataType": kind, "projection": "SUMMARY", "pageSize": 500000})
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    keep = PATIENT_ATTRS if kind == "PATIENT" else SAMPLE_ATTRS
    df = df[df["clinicalAttributeId"].isin(keep)]
    idx = "patientId" if kind == "PATIENT" else "sampleId"
    wide = df.pivot_table(index=idx, columns="clinicalAttributeId",
                          values="value", aggfunc="first")
    if kind == "SAMPLE":
        wide["patientId"] = df.groupby("sampleId")["patientId"].first()
    return wide.reset_index()


def main():
    studies = [s["studyId"] for s in get(f"{API}/studies", {})
               if s["studyId"].endswith("_tcga_pan_can_atlas_2018")]
    print(f"{len(studies)} PanCancer Atlas studies", flush=True)

    out = []
    for i, st in enumerate(sorted(studies), 1):
        pat = pull(st, "PATIENT")
        smp = pull(st, "SAMPLE")
        if pat.empty or smp.empty:
            print(f"  [{i:2d}/{len(studies)}] {st:45s} SKIPPED (empty)", flush=True)
            continue
        m = smp.merge(pat, on="patientId", how="left")
        m["studyId"] = st
        out.append(m)
        print(f"  [{i:2d}/{len(studies)}] {st:45s} {len(m):5d} samples", flush=True)

    all_df = pd.concat(out, ignore_index=True)
    # Xena/TIMER sample barcodes are 15-char (TCGA-XX-XXXX-01); cBioPortal matches already.
    all_df["sample_barcode"] = all_df["sampleId"].str.slice(0, 15)
    # cBioPortal has no CANCER_TYPE_ACRONYM here; derive it from the study id.
    # coadread and gbmlgg are merged studies and get split later using the Xena phenotype.
    all_df["cohort"] = (all_df["studyId"]
                        .str.replace("_tcga_pan_can_atlas_2018", "", regex=False)
                        .str.upper())
    for c in ["OS_MONTHS", "DSS_MONTHS", "PFS_MONTHS", "DFS_MONTHS", "AGE",
              "TMB_NONSYNONYMOUS", "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE",
              "MUTATION_COUNT", "FRACTION_GENOME_ALTERED", "ANEUPLOIDY_SCORE"]:
        if c in all_df:
            all_df[c] = pd.to_numeric(all_df[c], errors="coerce")

    dest = os.path.join(PROC, "clinical_pancan.tsv")
    all_df.to_csv(dest, sep="\t", index=False)
    print(f"\nwrote {dest}  shape={all_df.shape}", flush=True)
    print("cohorts:", all_df["cohort"].nunique())
    print("\ncoverage of key fields (non-null):")
    for c in ["OS_MONTHS", "OS_STATUS", "DSS_MONTHS", "PFS_MONTHS", "AGE", "SEX",
              "AJCC_PATHOLOGIC_TUMOR_STAGE", "GRADE", "TMB_NONSYNONYMOUS",
              "MSI_SCORE_MANTIS", "SUBTYPE"]:
        n = all_df[c].notna().sum() if c in all_df else 0
        print(f"   {c:32s} {n:6d} / {len(all_df)}")
    print("\nGRADE availability by cohort (top 12):")
    g = (all_df[all_df["GRADE"].notna()].groupby("cohort").size()
         .sort_values(ascending=False))
    print(g.head(12).to_string())


if __name__ == "__main__":
    main()
