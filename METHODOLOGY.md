# Methodology and Computational Workflow

Pan-cancer analysis of Nucleolin (NCL): complete specification of data sources, analytical decisions, script order and outputs.

This document is written so that a reader who has never seen the project can reproduce every number in the manuscript from public data. Each analysis step names the script that performs it, the inputs it consumes, the outputs it writes, and the manuscript element it supports. Where an analytical choice could reasonably have been made differently, the reasoning is stated rather than left implicit.

---

## 1. Scope and purpose

This is a revision of a previously submitted manuscript. The revision was prompted by peer review, and it does two things:

1. **Adds the analyses the reviewers required** — multiple-testing correction, multivariable survival modelling, multi-algorithm immune deconvolution, immune-checkpoint testing, GSEA, effect sizes with confidence intervals, and independent external validation.
2. **Corrects errors in the previous version.** Several of the previous conclusions did not survive reanalysis. Section 8 lists every correction and the evidence for it.

The workflow is deliberately linear and file-based. Every step writes its output to disk before the next step reads it, so any step can be re-run in isolation, and an interrupted run resumes without repeating completed work.

---

## 2. Environment

| Component | Version |
|---|---|
| OS | Windows 11 (10.0.26200) |
| Python | 3.14.3 |
| pandas | 2.3.3 |
| numpy | 2.4.3 |
| scipy | 1.17.1 |
| statsmodels | 0.14.6 |
| lifelines | 0.30.3 |
| gseapy | 1.3.1 |
| cptac | 1.5.14 |
| matplotlib | 3.10.8 |
| seaborn | 0.13.2 |

Install with:

```bash
pip install -r requirements.txt
```

**Determinism.** All stochastic procedures are seeded: bootstrap resampling for Cliff's delta confidence intervals (`seed=0`, 2,000 resamples) and GSEA permutation (`seed=0`, 1,000 permutations). Re-running the pipeline on the same inputs reproduces the same numbers.

**Resource requirements.** Peak disk use is approximately 6 GB, of which 4.6 GB is the expression matrix. Peak memory is under 1 GB — the pipeline is written to stream rather than load, because it was developed on a machine with limited free RAM (see §4.1). Total runtime is roughly 90 minutes, dominated by GSEA.

---

## 3. Data sources

| Dataset | Source | Retrieved by | Purpose |
|---|---|---|---|
| `TcgaTargetGtex_rsem_gene_tpm` | UCSC Xena Toil recompute | `01_download.py` | Expression, 60,498 genes × 19,131 samples, log₂(TPM+0.001) |
| `TcgaTargetGTEX_phenotype` | UCSC Xena | `01_download.py` | Sample → study, tissue, sample type |
| `gencode.v23.annotation.gene.probemap` | UCSC Xena | `01_download.py` | Ensembl ID → gene symbol |
| TCGA PanCancer Atlas clinical | cBioPortal REST API | `02_fetch_clinical.py` | Survival endpoints, age, sex, stage, grade, TMB, MSI |
| `infiltration_estimation_for_tcga` | TIMER2.0 | `01_download.py` | Immune infiltration, 7 algorithms |
| MSigDB Hallmark v2024.1.Hs | Broad Institute | `01_download.py` | GSEA gene sets |
| MSigDB Reactome v2024.1.Hs | Broad Institute | `01_download.py` | GSEA gene sets |
| CPTAC proteomics | `cptac` package | `08_cptac_validation.py` | Independent protein-level validation |

### 3.1 Why the Xena Toil recompute

TCGA and GTEx were originally processed with different pipelines. Comparing them directly introduces batch effects that can be larger than the biological difference being measured. The Toil recompute applies one alignment and quantification pipeline to TCGA, TARGET and GTEx together, which is what makes a tumour-versus-GTEx comparison defensible.

This matters concretely for the previous version of this manuscript. Its Methods described GEPIA2 with GTEx normals, but the figure presented was produced by TIMER2's differential-expression module, which uses TCGA adjacent normals only. Several cancers therefore rested on very few normals — pancreatic adenocarcinoma on four. The present analysis uses 167 GTEx pancreas samples for that comparison.

### 3.2 Why cBioPortal for clinical data

Survival endpoints are the curated definitions of the TCGA Clinical Data Resource (Liu et al., *Cell* 2018), which harmonises endpoint definitions across studies. Raw TCGA follow-up fields are inconsistent between projects and are not recommended for cross-cancer comparison. cBioPortal also supplies tumour mutational burden and MANTIS microsatellite-instability scores, both of which reviewers asked for, from the same harmonised release.

---

## 4. Script order

Run in this order. Each script is idempotent: it skips work whose output already exists.

```
scripts/
  00_check_sources.py        # probe all data URLs before downloading
  01_download.py             # download raw data (resumable)
  02_fetch_clinical.py       # cBioPortal clinical/survival/TMB/MSI
  03_build_matrix.py         # stream expression into a float32 matrix
  04_expression.py           # tumour vs normal; stage association
  05_survival.py             # log-rank, univariate and multivariable Cox
  06_immune.py               # infiltration, checkpoints, scores, TMB/MSI
  07a_expression_filter.py   # per-cohort expressed-gene filter
  07_gsea_correlations.py    # NCL vs all genes, within each cancer
  07_gsea.py                 # pre-ranked GSEA (resumable, parallel)
  08_cptac_validation.py     # CPTAC protein-level validation
  09_manuscript_numbers.py   # emit every number quoted in the manuscript
  10_figures.py              # publication figures at 600 dpi

  cohorts.py                 # cohort definitions, TCGA↔GTEx map, gene panels
  data_io.py                 # shared loaders
  statsutil.py               # effect sizes, CIs, FDR, trend and partial tests
  test_statsutil.py          # verification suite for statsutil
```

Steps 04, 05, 06 and 08 are mutually independent and may be run in any order once 03 completes. Step 07 requires 07a and 07_gsea_correlations.

### 4.1 Two implementation decisions worth knowing

**`03_build_matrix.py` writes sequentially, not through a writable memmap.** The first implementation opened a 4.6 GB writable `np.memmap` and assigned into it. On a machine with little free RAM this accumulates dirty pages faster than they can be flushed, and the build became I/O-bound — it processed roughly a third of the file in the time the current version takes to process all of it. Writing sequentially with an ordinary file handle and later *reading* through a read-only memmap avoids this: read-only pages are clean and evictable. Runtime fell from an extrapolated ~50 minutes to 601 seconds.

**`07_gsea_correlations.py` makes one pass over the matrix, not one per cancer.** Correlating NCL against every gene within each of 32 cohorts invites the pattern `X[:, cohort_columns]`, which for a memmap re-reads every row — that is, the entire 4.6 GB file — once per cohort, about 150 GB of I/O. Instead the file is read once in row blocks, and all cohorts are updated from each block.

---

## 5. Analysis specification

### 5.1 Cohort and comparator definition — `cohorts.py`, `data_io.sample_groups()`

Tumour samples are `Primary Tumor`, plus `Primary Blood Derived Cancer – Peripheral Blood` for LAML. Two normal comparators are defined per cancer:

- **GTEx normal**, matched by tissue of origin.
- **TCGA adjacent normal** (`Solid Tissue Normal`).

Every TCGA↔GTEx pairing carries a match-quality flag:

- `good` — the GTEx tissue is the tumour's tissue of origin (e.g. LIHC↔Liver).
- `approximate` — the closest available tissue, but not an exact counterpart: CHOL↔Liver (bile duct absent), HNSC↔Salivary Gland (squamous mucosa absent), READ↔Colon, SARC↔Adipose, DLBC↔Spleen, LAML↔Bone Marrow.
- `none` — no acceptable counterpart: MESO, THYM, UVM. These are analysed against adjacent normals only.

Flags are carried through to the results tables and reported in the manuscript, so that a reader can discount approximate comparisons.

**Why both comparators are reported.** Neither is unambiguously correct. Adjacent normal tissue shares the patient's genotype and exposures but is subject to field effects, inflammation and surgical ischaemia. GTEx tissue is unaffected by the tumour but is post-mortem and from different individuals. For ribosome-biogenesis genes such as NCL, which respond to metabolic and ischaemic state, this choice changes the answer — in kidney and thyroid it reverses the direction (§8.3). Reporting one comparator silently would have concealed that.

### 5.2 Differential expression — `04_expression.py` → `T1_differential_expression.tsv`

Two-sided Wilcoxon rank-sum test per cancer per comparator. Effect sizes:

- **Cliff's delta** with a 95% percentile bootstrap CI (2,000 stratified resamples). Non-parametric, because expression distributions are skewed and some cohorts have very few normals. A bootstrap CI is used rather than Cliff's asymptotic variance, which is unreliable at these sample sizes.
- **Hedges' *g*** with an analytic CI, for readers who prefer a standardised mean difference.

FDR is applied across cancers *within* each comparator family.

### 5.3 Stage association — `04_expression.py` → `T2_stage_association.tsv`

Two tests per cancer:

- **Kruskal–Wallis** — do stages differ at all?
- **Jonckheere–Terpstra** — is there a *monotonic trend* across ordered stages I–IV? Tie-corrected normal approximation.

**This distinction is the point of the analysis.** The previous version claimed NCL "progressively increased as the cancer advanced through stages I–IV". That claim came from web-tool plots whose significance annotations compare each stage against *normal tissue*, not against the preceding stage. Such annotations can all be significant while expression is flat or falling across stages — which is what several of the cited cancers actually show. Jonckheere–Terpstra is the test the claim required.

Inclusion: ≥3 stage groups with ≥5 patients each, and ≥40 staged patients.

### 5.4 Survival — `05_survival.py` → `T3_survival.tsv`

Endpoints: overall survival (OS), disease-specific survival (DSS), progression-free survival (PFS).

NCL is z-scored *within* each cancer, so hazard ratios are per standard deviation and comparable across cohorts.

Three models per cancer per endpoint:

1. **Log-rank** on a median split — the analysis reported in most prior work, retained for comparability.
2. **Univariate Cox** on continuous NCL.
3. **Multivariable Cox** adjusting for age, sex, stage and grade.

Covariates are included only where available for >80% of the cohort and showing more than one distinct value; the covariates actually used are recorded per model, because they differ (grade is unavailable for many cancers, stage for GBM/LGG/LAML).

Inclusion: ≥30 patients and ≥10 events.

**Proportional hazards** are tested for every multivariable model via Schoenfeld residuals with rank-transformed time (`lifelines.statistics.proportional_hazard_test`). Violations are reported per model. This matters: a hazard ratio from a model violating PH summarises a time-varying effect and must not be read as a constant risk multiplier. The assumption was violated for the NCL term in 10 of 74 models — but *not* in either cancer that retained significance, so the headline results are interpretable as stated.

> An earlier implementation attempted this check through `CoxPHFitter.check_assumptions()` and silently parsed nothing, so zero models were actually tested while the code appeared to succeed. It was caught by inspecting the output column rather than trusting the absence of an error. The lesson is recorded here because a manuscript claiming PH had been verified would otherwise have been wrong.

FDR is applied across cancers within each endpoint and model type.

### 5.5 Immune microenvironment — `06_immune.py` → `T4`, `T5`, `T6`, `T8`

**Infiltration.** For every cancer × cell type × algorithm, Spearman correlation with NCL and partial Spearman adjusted for tumour purity. Purity is represented by the EPIC "uncharacterised cell" fraction — EPIC's explicit estimate of the compartment that is neither immune nor stromal. Partial correlation is computed by regressing the ranks of both variables on the ranks of purity and correlating the residuals, with the t-test taking n−3 degrees of freedom.

**Cross-algorithm concordance** (`T5`). Ten cell types resolvable by more than one algorithm are designated canonical. For each cancer × canonical cell type we record how many of the seven algorithms give a significant purity-adjusted association and whether they agree in sign:

- *concordant* — ≥2 algorithms significant, no sign disagreement;
- *conflicting* — significant estimates of opposite sign both present.

This exists because the previous version used a different algorithm in each panel of its immune figure (MCP-counter, CIBERSORT-ABS, xCell, quanTIseq), which is indistinguishable from selecting the algorithm that gave the desired result. Quantifying concordance is the honest alternative. It turns out to matter: a third of combinations conflict in direction (§7).

**Checkpoints** (`T6`). Sixteen genes, each reported three ways: unadjusted, purity-adjusted, and **proliferation-adjusted**.

The proliferation adjustment is the analytical core of this section. NCL is a ribosome-biogenesis gene whose expression tracks proliferative rate, and proliferative tumours differ systematically in immune composition. A raw NCL–checkpoint correlation may therefore reflect nothing more specific than shared covariation with proliferation. The proliferation score is the mean within-cohort z-score of ten canonical markers (*MKI67*, *PCNA*, *TOP2A*, *CCNB1*, *BUB1*, *AURKA*, *CDK1*, *TYMS*, *RRM2*, *TK1*). Associations described as **robust** are those surviving *both* purity and proliferation adjustment.

**Genomic and composite scores** (`T8`). NCL versus xCell immune/stromal/microenvironment scores, TMB, MANTIS MSI, aneuploidy score, and fraction of genome altered.

FDR is applied across all tests within each family (all infiltration tests together; all checkpoint tests together; all score tests together).

### 5.6 Gene set enrichment — `07a`, `07_gsea_correlations.py`, `07_gsea.py` → `T9`, `T10`

**Ranking statistic.** Spearman correlation of NCL against every gene, computed **within each cancer**. Pooled pan-cancer correlations — as used in the previous version's Figure 4c — are dominated by tissue-of-origin differences in composition and describe no within-tumour relationship. Sanity check: NCL's self-correlation must equal 1.0 in every cohort, asserted in the script.

**Expressed-gene restriction** (`07a`). Ranked lists are restricted to genes with TPM > 1 in ≥25% of that cohort's tumours. Without this, ~10% of genes carried tied ranking statistics: genes that are essentially undetected produce degenerate correlations that GSEA then orders arbitrarily. The restriction reduced ties to below 0.2%. It also cut GSEA runtime roughly five-fold, since the ranked list shrinks from ~58,000 to ~15,000–18,000 genes.

**Enrichment.** Pre-ranked GSEA against Hallmark and Reactome (v2024.1.Hs), 1,000 permutations, gene sets restricted to 15–500 genes.

**Resumability.** Each (collection, cohort) unit writes its own file under `results/gsea_parts/` and is skipped if present. Writes are atomic (write to `.tmp`, then `os.replace`) so an interrupted write is never mistaken for a completed result. This was added after a machine shutdown lost an entire in-memory GSEA run; the correlation cache meant only the enrichment step had to be repeated.

### 5.7 Independent validation — `08_cptac_validation.py` → `T7`

CPTAC proteomics for ten cohorts, tumour versus adjacent normal.

CPTAC is the strongest validation available without new experimental work because it is simultaneously: an **independent patient series**, an **orthogonal platform** (mass spectrometry rather than sequencing), a **different analyte** (protein rather than transcript), and **paired** (adjacent normal tissue from the same patients).

Normal samples are identified by the `.N` sample-ID suffix. Tumour and normal are compared by Wilcoxon rank-sum with Cliff's delta, and by Wilcoxon signed-rank where the same patient contributes both.

This does **not** discharge the reviewers' request for functional experiments. See §9.

### 5.8 Number provenance — `09_manuscript_numbers.py`

Regenerates `results/MANUSCRIPT_NUMBERS.txt`, containing every value quoted in the manuscript, traced to its source table. The manuscript quotes this file rather than values typed from memory, so that changing an analysis surfaces any stale number rather than leaving it silently wrong.

---

## 6. Multiple testing

Benjamini–Hochberg FDR is applied within each family of tests:

| Family | Tests |
|---|---|
| Differential expression, per comparator | 29 / 21 |
| Stage trend | 17 |
| Survival, per endpoint × model | 21–29 each |
| Immune infiltration | 3,910 |
| Immune checkpoints | 495 |
| Genomic/composite scores | 230 |
| GSEA | per collection |

Adjusted values are reported as *q*; significance is q<0.05 throughout. Every association is reported with a point estimate, a 95% confidence interval, and its sample size.

---

## 7. Verification

`test_statsutil.py` verifies every statistical utility before use. Run:

```bash
python scripts/test_statsutil.py
```

Checks performed:

- **Cliff's delta** against brute-force pairwise dominance across four sample-size/effect combinations, and on tied data.
- **Hedges' *g*** against the closed-form expression; CI brackets the estimate.
- **BH-FDR** against `statsmodels.stats.multitest.multipletests` (agreement to 1×10⁻¹⁶); NaN inputs pass through.
- **Jonckheere–Terpstra**: strong increasing and decreasing trends recovered with correct sign; **null calibration** over 300 simulated datasets gives a 5.7% type-I error rate at α=0.05; a **non-monotonic (up-then-down) pattern is not reported as a trend**; permutation and asymptotic p-values agree.
- **Partial Spearman**: a correlation induced entirely by a confounder is removed (ρ 0.82 → 0.08); a genuine association independent of the confounder is retained.
- **Spearman CI**: brackets the estimate and narrows with n.

All checks pass. The null-calibration and non-monotonic checks are the important ones: they establish that the trend test does not manufacture the very finding the previous version over-claimed.

---

## 8. Corrections to the previous version

Each correction below is reproducible from the tables named.

### 8.1 Claims made without analysis

The previous title and abstract asserted correlations between NCL and PD-L1, CTLA-4, TIM-3, IL-10 and TGF-β. **No checkpoint analysis existed** in that work — no method, no result, no figure. The title also cited CancerSEA, and the text cited cBioPortal and "genetic alterations"; neither analysis existed either.

Resolution: the checkpoint analysis has now been performed (`T6`). Those specific associations are weak — CTLA-4 robust in 4 of 33 cancers, PD-1 in 3, TIM-3 in 9 — while **B7-H3 (CD276) is robust in 21 of 33** and is the finding that replaces them. CancerSEA references are removed. cBioPortal is now genuinely used, as the clinical-data source.

### 8.2 Direction of dysregulation

Previously reported as *reduced* in KIRP, READ and cutaneous melanoma. Re-analysis (`T1`): READ δ=+0.90 (q=2.2×10⁻³⁸) and SKCM δ=+0.65 (q=4.3×10⁻²⁵) are **elevated** against every available comparator. Only KIRP reproduces a reduction, and only against adjacent normals.

### 8.3 Comparator dependence (new finding)

In KIRP, KICH, KIRC and THCA the direction of NCL dysregulation **depends on the normal used**:

| Cancer | vs GTEx | vs adjacent normal |
|---|---|---|
| KIRP | +0.63 (q=5.1×10⁻⁸) | −0.29 (q=1.1×10⁻²) |
| KICH | +0.19 (ns) | −0.80 (q=1.6×10⁻⁸) |
| KIRC | +0.78 (q=7.3×10⁻¹²) | −0.07 (ns) |
| THCA | +0.37 (q=4.8×10⁻¹⁷) | −0.27 (q=1.4×10⁻³) |

CPTAC ccRCC protein (δ=+0.65) supports the GTEx direction in kidney.

### 8.4 Stage progression

Previously claimed as a progressive I–IV increase across eight cancers. A formal trend test (`T2`) supports **one** of them (LIHC). Two cancers show an increasing trend overall (LIHC, LUAD); one *decreases* (KIRC); 14 show none.

### 8.5 Survival

Previously reported adverse OS associations in HNSC, KICH, KIRP and LIHC. After covariate adjustment and FDR correction (`T3`):

| Cancer | Previously | Adjusted result |
|---|---|---|
| KIRP | P=0.034 | **HR 2.12 (1.37–3.29), q=0.019 — retained** |
| LIHC | P=1.4×10⁻⁶ | HR 1.32 (1.07–1.62), q=0.079 — not significant |
| HNSC | P=0.0027 | HR 1.06 (0.91–1.25), q=0.72 — null |
| KICH | P=0.017 | not modellable (too few events) |

### 8.6 Mislabelled survival panels

Three of the eight Kaplan–Meier panels in the previous version carried cohort sizes incompatible with the cancer named. Most clearly, the "KICH" panel showed 877 patients; TCGA KICH contains **65** primary tumours. The sizes correspond to pooled Human Protein Atlas cohorts (renal 877, lung 994) rather than individual TCGA projects.

### 8.7 Cancer nomenclature

The previous Table 1, presented as the TCGA cancer list, contained 39 entries for 33 projects, included entities that are not TCGA projects (AST, CML, HGG, NSCLC, ODG, MEL), used non-TCGA codes (HCC, RCC, CM, CCC, UEC, AML), and **defined ACC as adenoid cystic carcinoma when TCGA ACC is adrenocortical carcinoma**. It appears to have been a CancerSEA list. Corrected in `cohorts.py::FULL_NAME`.

### 8.8 Tool misattribution

Methods described GEPIA2 with GTEx normals; the figure was TIMER2's TCGA-only module. Consequence: PAAD rested on 4 adjacent normals. Now 167 GTEx pancreas samples.

### 8.9 Enrichment analysis

The previous miRNA analysis named miRDB in Methods and miRNet in Results and the legend; its figure was a miRNA–target network in which NCL was one node among many, presented as "NCL-interacting miRNAs"; and it assigned roles in EMT, metastasis and apoptosis to individual miRNAs that were never tested. It is **not replaced** — it supported no conclusion, and no substitute claim is made. GSEA replaces the pathway analysis.

---

## 9. What this workflow does not do

**It generates no new experimental data.** Both reviewers asked for validation in patient cohorts, tissue samples, or functional in vitro experiments. CPTAC addresses the *cohort* element — independent patients, orthogonal platform, protein-level, paired normals — and substantiates the overexpression claim. It does not address the *functional* element at all.

Specifically, the association between NCL and B7-H3 is the most interesting finding here, and nothing in this workflow establishes that NCL regulates B7-H3, or that either is causal for the immune-excluded phenotype. That requires perturbation: NCL knockdown or overexpression with measurement of B7-H3 and immune composition. It has not been done, and the manuscript says so.

Other constraints:

- Deconvolution infers composition from bulk expression; it is not a measurement of infiltration. The poor cross-algorithm concordance is itself evidence of this limit.
- Six TCGA↔GTEx pairings are approximate; three cancers have none.
- Small cohorts (CHOL n=36, UCS n=57, KICH n=65) are underpowered; absence of significance there is not evidence of absence.
- Covariate sets differ between survival models because stage and grade are not universally available.
- TCGA is predominantly treatment-naive primary resections from limited populations.

---

## 10. Output inventory

| File | Contents | Manuscript element |
|---|---|---|
| `T1_differential_expression.tsv` | Tumour vs normal, both comparators | Table 1, Figure 1a–b |
| `T2_stage_association.tsv` | Kruskal–Wallis and trend tests | Table 2, Figure 1c |
| `T3_survival.tsv` | Log-rank, univariate and multivariable Cox, PH tests | Table 3, Figure 2 |
| `T4_immune_infiltration.tsv` | 3,910 cancer × cell type × algorithm tests | Figure 3a |
| `T5_algorithm_concordance.tsv` | Cross-algorithm agreement | Figure 3b |
| `T6_checkpoints.tsv` | 495 checkpoint tests, three adjustments | Table 5, Figure 4 |
| `T7_cptac_validation.tsv` | CPTAC protein-level validation | Table 6, Figure 5 |
| `T8_genomic_scores.tsv` | Immune/stromal scores, TMB, MSI, aneuploidy | Table 4 |
| `T9_ncl_gene_correlations.tsv.gz` | NCL vs all genes, per cancer | GSEA input |
| `T10_gsea_per_cancer.tsv.gz` | Hallmark and Reactome enrichment | Figure 6 |
| `MANUSCRIPT_NUMBERS.txt` | Every quoted number, traced to source | whole manuscript |

---

## 11. Reviewer comments addressed

| Comment | Where addressed |
|---|---|
| R1① Correlation vs causation | Title changed to "association"; causal language removed throughout; §4.2 states explicitly that no regulatory relationship is demonstrated |
| R1② Lacks experimental validation | **Partially.** CPTAC protein validation (§5.7); functional work explicitly not done (§9) |
| R1③ Inconsistent analysis, cancer differences | Uniform 7-algorithm immune analysis with concordance reporting (§5.5); heterogeneity is now the paper's thesis |
| R1④ Extent of AI use unclear | Explicit AI-use declaration in the manuscript |
| R2① Novelty not justified | Reframed around findings that are new *and* corrective: B7-H3 association, comparator dependence, absence of independent prognostic value |
| R2② No independent validation | CPTAC, 9 cohorts (§5.7) |
| R2③ No multiple-testing correction | BH-FDR across every family (§6) |
| R2④ Only univariate survival | Multivariable Cox with PH diagnostics (§5.4) |
| R2⑤ Single deconvolution; no immune scores/TMB/MSI | 7 algorithms + concordance; immune/stromal/microenvironment scores, TMB, MSI, aneuploidy (§5.5) |
| R2⑥ Mechanism speculative | Proliferation-adjusted analysis separates immune-specific from generic proliferation signal (§5.5) |
| R2⑦ Enrichment lacks depth | Per-cancer GSEA on Hallmark and Reactome (§5.6) |
| R2⑧ Heterogeneity insufficiently addressed | Central to the revised conclusions; universal-biomarker framing withdrawn |
| R2⑨ No external datasets | CPTAC (§5.7) |
| R2⑩ Figures and reporting | 600 dpi figures; effect sizes with CIs and sample sizes throughout |
