"""Survival analysis: univariate log-rank plus multivariable Cox per cancer.

Addresses Reviewer 2 comment 4 ("only univariate Kaplan-Meier analyses are
presented; multivariate Cox proportional hazards analyses incorporating age,
stage, grade and other clinical covariates are necessary") and comment 3
(multiple-testing correction), neither of which the submitted manuscript did.

It also replaces submitted Figure 2, in which three of the eight Kaplan-Meier
panels carry cohort sizes incompatible with the cancer they are labelled with
(the "KICH" panel shows 877 patients; TCGA KICH has 65 primary tumours).

Endpoints are the curated TCGA-CDR definitions (Liu et al., Cell 2018):
overall survival, disease-specific survival and progression-free survival.
"""
import os, warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test

import data_io as D
import statsutil as S
from cohorts import COHORTS, FULL_NAME

warnings.filterwarnings("ignore")
GENE = "NCL"
ENDPOINTS = {"OS": ("OS_MONTHS", "OS_STATUS"),
             "DSS": ("DSS_MONTHS", "DSS_STATUS"),
             "PFS": ("PFS_MONTHS", "PFS_STATUS")}


def build_frame():
    """One row per TCGA primary tumour with NCL expression and clinical covariates."""
    ncl = D.gene_vector(GENE)
    ix = D.index()
    pos = {s: i for i, s in enumerate(ix["samples"])}
    ann, _ = D.sample_groups()
    clin = D.clinical()

    tum = ann.index[ann.group == "tumour"]
    keep = [s for s in tum if s in clin.index]
    df = clin.reindex(keep).copy()
    df["cohort"] = ann.loc[keep, "cohort"].to_numpy()
    df["NCL"] = [ncl[pos[s]] for s in keep]

    def ev(v):
        if not isinstance(v, str):
            return np.nan
        return 1.0 if v.startswith("1") else (0.0 if v.startswith("0") else np.nan)

    for ep, (tcol, scol) in ENDPOINTS.items():
        df[f"{ep}_time"] = pd.to_numeric(df[tcol], errors="coerce")
        df[f"{ep}_event"] = df[scol].map(ev)

    def stage_num(v):
        if not isinstance(v, str):
            return np.nan
        v = v.upper().replace("STAGE", "").strip()
        for k, n in (("IV", 4), ("III", 3), ("II", 2), ("I", 1)):
            if v.startswith(k):
                return n
        return np.nan

    def grade_num(v):
        if not isinstance(v, str):
            return np.nan
        v = v.upper().strip()
        for k, n in (("G4", 4), ("G3", 3), ("G2", 2), ("G1", 1)):
            if v.startswith(k):
                return n
        return np.nan

    df["stage_num"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(stage_num)
    df["grade_num"] = df["GRADE"].map(grade_num)
    df["age"] = pd.to_numeric(df["AGE"], errors="coerce")
    df["male"] = df["SEX"].map(lambda v: 1.0 if str(v).upper().startswith("M")
                               else (0.0 if str(v).upper().startswith("F") else np.nan))
    return df


def analyse(df):
    rows = []
    for code in sorted(COHORTS):
        sub = df[df.cohort == code]
        if len(sub) < 30:
            continue
        # z-score NCL within cancer so HRs are per SD and comparable across cohorts.
        z = (sub["NCL"] - sub["NCL"].mean()) / sub["NCL"].std(ddof=0)
        for ep in ENDPOINTS:
            d = pd.DataFrame({
                "T": sub[f"{ep}_time"], "E": sub[f"{ep}_event"], "NCL_z": z,
                "age": sub["age"], "male": sub["male"],
                "stage": sub["stage_num"], "grade": sub["grade_num"]})
            d = d[(d["T"] > 0) & d["T"].notna() & d["E"].notna()]
            if len(d) < 30 or d["E"].sum() < 10:
                continue

            rec = dict(cohort=code, full_name=FULL_NAME[code], endpoint=ep,
                       n=len(d), events=int(d["E"].sum()))

            # --- univariate: median split log-rank (what the paper reported) ---
            hi = d["NCL_z"] > d["NCL_z"].median()
            lr = logrank_test(d["T"][hi], d["T"][~hi], d["E"][hi], d["E"][~hi])
            rec["logrank_p"] = float(lr.p_value)
            km = KaplanMeierFitter()
            try:
                km.fit(d["T"][hi], d["E"][hi]); m_hi = km.median_survival_time_
                km.fit(d["T"][~hi], d["E"][~hi]); m_lo = km.median_survival_time_
                rec["median_surv_high"], rec["median_surv_low"] = m_hi, m_lo
            except Exception:
                pass

            # --- univariate Cox on continuous NCL ---
            try:
                cph = CoxPHFitter().fit(d[["T", "E", "NCL_z"]], "T", "E")
                rec["uni_HR"] = float(np.exp(cph.params_["NCL_z"]))
                ci = cph.confidence_intervals_.loc["NCL_z"]
                rec["uni_HR_lo"], rec["uni_HR_hi"] = float(np.exp(ci.iloc[0])), float(np.exp(ci.iloc[1]))
                rec["uni_p"] = float(cph.summary.loc["NCL_z", "p"])
            except Exception as e:
                rec["uni_note"] = type(e).__name__

            # --- multivariable Cox: adjust for what is available in this cohort ---
            covars = ["NCL_z"]
            for c in ["age", "male", "stage", "grade"]:
                v = d[c]
                # Require good coverage and actual variation to avoid degenerate fits.
                if v.notna().mean() > 0.8 and v.nunique(dropna=True) > 1:
                    covars.append(c)
            dm = d[["T", "E"] + covars].dropna()
            rec["adj_covariates"] = "+".join(c for c in covars if c != "NCL_z") or "none"
            rec["n_adj"] = len(dm)
            rec["events_adj"] = int(dm["E"].sum()) if len(dm) else 0
            if len(dm) >= 30 and dm["E"].sum() >= 10 and len(covars) > 1:
                try:
                    cph = CoxPHFitter().fit(dm, "T", "E")
                    rec["adj_HR"] = float(np.exp(cph.params_["NCL_z"]))
                    ci = cph.confidence_intervals_.loc["NCL_z"]
                    rec["adj_HR_lo"] = float(np.exp(ci.iloc[0]))
                    rec["adj_HR_hi"] = float(np.exp(ci.iloc[1]))
                    rec["adj_p"] = float(cph.summary.loc["NCL_z", "p"])
                    # Proportional-hazards check (Schoenfeld residuals) on the
                    # NCL term and on the model as a whole.
                    try:
                        zph = proportional_hazard_test(cph, dm, time_transform="rank")
                        s = zph.summary
                        rec["ph_p_NCL"] = float(s.loc["NCL_z", "p"])
                        rec["ph_p_min_any"] = float(s["p"].min())
                    except Exception as e:
                        rec["ph_note"] = type(e).__name__
                except Exception as e:
                    rec["adj_note"] = type(e).__name__
            rows.append(rec)

    res = pd.DataFrame(rows)
    for ep in ENDPOINTS:
        m = res["endpoint"] == ep
        for col, q in (("logrank_p", "logrank_q"), ("uni_p", "uni_q"), ("adj_p", "adj_q")):
            if col in res:
                res.loc[m, q] = S.bh_fdr(res.loc[m, col].to_numpy())
    return res


def main():
    df = build_frame()
    print(f"samples with expression + clinical: {len(df)}")
    res = analyse(df)
    res.to_csv(os.path.join(D.TABLES, "T3_survival.tsv"), sep="\t", index=False)

    for ep in ["OS", "DSS", "PFS"]:
        r = res[(res.endpoint == ep) & res.adj_HR.notna()].sort_values("adj_HR")
        print(f"\n=== {ep}: multivariable Cox, HR per SD of NCL ===")
        print(f"{'':2}{'cohort':6} {'n':>5} {'ev':>4} {'adjHR':>6} {'95% CI':>15} "
              f"{'adj q':>9} {'uni q':>9}  covariates")
        for _, x in r.iterrows():
            sig = "*" if x.adj_q < 0.05 else " "
            uq = x.uni_q if pd.notna(x.uni_q) else np.nan
            print(f"{sig} {x.cohort:6} {x.n_adj:5.0f} {x.events_adj:4.0f} "
                  f"{x.adj_HR:6.2f} [{x.adj_HR_lo:5.2f},{x.adj_HR_hi:5.2f}] "
                  f"{x.adj_q:9.4f} {uq:9.4f}  {x.adj_covariates}")
        s = (r.adj_q < 0.05)
        print(f"  significant after adjustment: {s.sum()}/{len(r)} "
              f"({(r.adj_HR[s] > 1).sum()} adverse, {(r.adj_HR[s] < 1).sum()} protective)")
        if "ph_p_NCL" in r:
            ph = r[r.ph_p_NCL.notna()]
            print(f"  proportional hazards on NCL term tested in {len(ph)} models; "
                  f"violated (p<0.05) in {(ph.ph_p_NCL < 0.05).sum()}"
                  + (f" -> {ph[ph.ph_p_NCL<0.05].cohort.tolist()}"
                     if (ph.ph_p_NCL < 0.05).any() else ""))


if __name__ == "__main__":
    main()
