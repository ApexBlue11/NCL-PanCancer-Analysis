# Methodology and Computational Workflow

Pan-cancer analysis of Nucleolin (NCL): complete specification of data sources, analytical decisions, script order and outputs.

This document is written so that a reader who has never seen the project can reproduce every number in the manuscript from public data. Each analysis step names the script that performs it, the inputs it consumes, the outputs it writes, and the manuscript element it supports. Where an analytical choice could reasonably have been made differently, the reasoning is stated rather than left implicit.

---

## 1. Scope and purpose

This document specifies the computational workflow behind the pan-cancer NCL
analysis. It is written so that a reader who has never seen the project can
reproduce every number in the manuscript from public data.

Each analysis step names the script that performs it, the inputs it consumes,
the outputs it writes, and the manuscript element it supports. Where an
analytical choice could reasonably have been made differently, which normal
tissue to compare against, which deconvolution algorithm to trust, whether to
adjust for proliferation, the reasoning is stated rather than left implicit,
and Section 8 sets out the failure modes those choices are designed to avoid.

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

**Inputs are not archival, so they are checksummed.** UCSC Xena, TIMER2.0 and MSigDB reissue files at stable URLs; cBioPortal serves a live API with no version pin; the `cptac` package downloads data at run time. `11_manifest.py` therefore records the SHA-256, byte size and retrieval date of every input, the shape and checksum of every result table, and the resolved version of every package affecting a numeric result, into `results/MANIFEST.json` and `MANIFEST.md`. If a future run disagrees with the published numbers, diffing two manifests localises the cause to an input file, a package version or the code. This is the difference between "the code runs" and "the code reproduces these numbers".

**Resource requirements.** Peak disk use is approximately 6 GB, of which 4.6 GB is the expression matrix. Peak memory is under 1 GB, the pipeline is written to stream rather than load, because it was developed on a machine with limited free RAM (see §4.1). Total runtime is roughly 90 minutes, dominated by GSEA.

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

This matters concretely for the previous version of this manuscript. Its Methods described GEPIA2 with GTEx normals, but the figure presented was produced by TIMER2's differential-expression module, which uses TCGA adjacent normals only. Several cancers therefore rested on very few normals, pancreatic adenocarcinoma on four. The present analysis uses 167 GTEx pancreas samples for that comparison.

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
  12_supplementary_tables.py # Supplementary Table S1 and S-number index
  13_build_submission.py     # manuscript with figures and tables embedded
  11_manifest.py             # input checksums, package versions, seeds

  cohorts.py                 # cohort definitions, TCGA↔GTEx map, gene panels
  data_io.py                 # shared loaders
  statsutil.py               # effect sizes, CIs, FDR, trend and partial tests
  test_statsutil.py          # verification suite for statsutil
```

Steps 04, 05, 06 and 08 are mutually independent and may be run in any order once 03 completes. Step 07 requires 07a and 07_gsea_correlations.

### 4.1 Two implementation decisions worth knowing

**`03_build_matrix.py` writes sequentially, not through a writable memmap.** The first implementation opened a 4.6 GB writable `np.memmap` and assigned into it. On a machine with little free RAM this accumulates dirty pages faster than they can be flushed, and the build became I/O-bound, it processed roughly a third of the file in the time the current version takes to process all of it. Writing sequentially with an ordinary file handle and later *reading* through a read-only memmap avoids this: read-only pages are clean and evictable. Runtime fell from an extrapolated ~50 minutes to 601 seconds.

**`07_gsea_correlations.py` makes one pass over the matrix, not one per cancer.** Correlating NCL against every gene within each of 32 cohorts invites the pattern `X[:, cohort_columns]`, which for a memmap re-reads every row (that is, the entire 4.6 GB file) once per cohort, about 150 GB of I/O. Instead the file is read once in row blocks, and all cohorts are updated from each block.

---

## 5. Analysis specification

### 5.1 Cohort and comparator definition: `cohorts.py`, `data_io.sample_groups()`

Tumour samples are `Primary Tumor`, plus `Primary Blood Derived Cancer – Peripheral Blood` for LAML. Two normal comparators are defined per cancer:

- **GTEx normal**, matched by tissue of origin.
- **TCGA adjacent normal** (`Solid Tissue Normal`).

Every TCGA↔GTEx pairing carries a match-quality flag:

- `good`, the GTEx tissue is the tumour's tissue of origin (e.g. LIHC↔Liver).
- `approximate`, the closest available tissue, but not an exact counterpart: CHOL↔Liver (bile duct absent), HNSC↔Salivary Gland (squamous mucosa absent), READ↔Colon, SARC↔Adipose, DLBC↔Spleen, LAML↔Bone Marrow.
- `none`, no acceptable counterpart: MESO, THYM, UVM. These are analysed against adjacent normals only.

Flags are carried through to the results tables and reported in the manuscript, so that a reader can discount approximate comparisons.

**Why both comparators are reported.** Neither is unambiguously correct. Adjacent normal tissue shares the patient's genotype and exposures but is subject to field effects, inflammation and surgical ischaemia. GTEx tissue is unaffected by the tumour but is post-mortem and from different individuals. For ribosome-biogenesis genes such as NCL, which respond to metabolic and ischaemic state, this choice changes the answer, in kidney and thyroid it reverses the direction (§8.3). Reporting one comparator silently would have concealed that.

### 5.2 Differential expression: `04_expression.py` → `T1_differential_expression.tsv`

Two-sided Wilcoxon rank-sum test per cancer per comparator. Effect sizes:

- **Cliff's delta** with a 95% percentile bootstrap CI (2,000 stratified resamples). Non-parametric, because expression distributions are skewed and some cohorts have very few normals. A bootstrap CI is used rather than Cliff's asymptotic variance, which is unreliable at these sample sizes.
- **Hedges' *g*** with an analytic CI, for readers who prefer a standardised mean difference.

FDR is applied across cancers *within* each comparator family.

### 5.3 Stage association: `04_expression.py` → `T2_stage_association.tsv`

Two tests per cancer:

- **Kruskal–Wallis**, do stages differ at all?
- **Jonckheere–Terpstra**, is there a *monotonic trend* across ordered stages I–IV? Tie-corrected normal approximation.

**This distinction is the point of the analysis.** The previous version claimed NCL "progressively increased as the cancer advanced through stages I–IV". That claim came from web-tool plots whose significance annotations compare each stage against *normal tissue*, not against the preceding stage. Such annotations can all be significant while expression is flat or falling across stages, which is what several of the cited cancers actually show. Jonckheere–Terpstra is the test the claim required.

Inclusion: ≥3 stage groups with ≥5 patients each, and ≥40 staged patients.

### 5.4 Survival: `05_survival.py` → `T3_survival.tsv`

Endpoints: overall survival (OS), disease-specific survival (DSS), progression-free survival (PFS).

NCL is z-scored *within* each cancer, so hazard ratios are per standard deviation and comparable across cohorts.

Three models per cancer per endpoint:

1. **Log-rank** on a median split, the analysis reported in most prior work, retained for comparability.
2. **Univariate Cox** on continuous NCL.
3. **Multivariable Cox** adjusting for age, sex, stage and grade.

Covariates are included only where available for >80% of the cohort and showing more than one distinct value; the covariates actually used are recorded per model, because they differ (grade is unavailable for many cancers, stage for GBM/LGG/LAML).

Inclusion: ≥30 patients and ≥10 events.

**Proportional hazards** are tested for every multivariable model via Schoenfeld residuals with rank-transformed time (`lifelines.statistics.proportional_hazard_test`). Violations are reported per model. This matters: a hazard ratio from a model violating PH summarises a time-varying effect and must not be read as a constant risk multiplier. The assumption was violated for the NCL term in 10 of 74 models, but *not* in either cancer that retained significance, so the headline results are interpretable as stated.

> An earlier implementation attempted this check through `CoxPHFitter.check_assumptions()` and silently parsed nothing, so zero models were actually tested while the code appeared to succeed. It was caught by inspecting the output column rather than trusting the absence of an error. The lesson is recorded here because a manuscript claiming PH had been verified would otherwise have been wrong.

FDR is applied across cancers within each endpoint and model type.

### 5.5 Immune microenvironment: `06_immune.py` → `T4`, `T5`, `T6`, `T8`

**Infiltration.** For every cancer × cell type × algorithm, Spearman correlation with NCL and partial Spearman adjusted for tumour purity. Purity is represented by the EPIC "uncharacterised cell" fraction, EPIC's explicit estimate of the compartment that is neither immune nor stromal. Partial correlation is computed by regressing the ranks of both variables on the ranks of purity and correlating the residuals, with the t-test taking n−3 degrees of freedom.

**Cross-algorithm concordance** (`T5`). Ten cell types resolvable by more than one algorithm are designated canonical. For each cancer × canonical cell type we record how many of the seven algorithms give a significant purity-adjusted association and whether they agree in sign:

- *concordant*, ≥2 algorithms significant, no sign disagreement;
- *conflicting*, significant estimates of opposite sign both present.

This exists because the previous version used a different algorithm in each panel of its immune figure (MCP-counter, CIBERSORT-ABS, xCell, quanTIseq), which is indistinguishable from selecting the algorithm that gave the desired result. Quantifying concordance is the honest alternative. It turns out to matter: a third of combinations conflict in direction (§7).

**Checkpoints** (`T6`). Sixteen genes, each reported three ways: unadjusted, purity-adjusted, and **proliferation-adjusted**.

The proliferation adjustment is the analytical core of this section. NCL is a ribosome-biogenesis gene whose expression tracks proliferative rate, and proliferative tumours differ systematically in immune composition. A raw NCL–checkpoint correlation may therefore reflect nothing more specific than shared covariation with proliferation. The proliferation score is the mean within-cohort z-score of ten canonical markers (*MKI67*, *PCNA*, *TOP2A*, *CCNB1*, *BUB1*, *AURKA*, *CDK1*, *TYMS*, *RRM2*, *TK1*). Associations described as **robust** are those surviving *both* purity and proliferation adjustment.

**Genomic and composite scores** (`T8`). NCL versus xCell immune/stromal/microenvironment scores, TMB, MANTIS MSI, aneuploidy score, and fraction of genome altered.

FDR is applied across all tests within each family (all infiltration tests together; all checkpoint tests together; all score tests together).

### 5.6 Gene set enrichment: `07a`, `07_gsea_correlations.py`, `07_gsea.py` → `T9`, `T10`

**Ranking statistic.** Spearman correlation of NCL against every gene, computed **within each cancer**. Pooled pan-cancer correlations (as used in the previous version's Figure 4c) are dominated by tissue-of-origin differences in composition and describe no within-tumour relationship. Sanity check: NCL's self-correlation must equal 1.0 in every cohort, asserted in the script.

**Expressed-gene restriction** (`07a`). Ranked lists are restricted to genes with TPM > 1 in ≥25% of that cohort's tumours. Without this, ~10% of genes carried tied ranking statistics: genes that are essentially undetected produce degenerate correlations that GSEA then orders arbitrarily. The restriction reduced ties to below 0.2%. It also cut GSEA runtime roughly five-fold, since the ranked list shrinks from ~58,000 to ~15,000–18,000 genes.

**Enrichment.** Pre-ranked GSEA against Hallmark and Reactome (v2024.1.Hs), 1,000 permutations, gene sets restricted to 15–500 genes.

**Resumability.** Each (collection, cohort) unit writes its own file under `results/gsea_parts/` and is skipped if present. Writes are atomic (write to `.tmp`, then `os.replace`) so an interrupted write is never mistaken for a completed result. This was added after a machine shutdown lost an entire in-memory GSEA run; the correlation cache meant only the enrichment step had to be repeated.

### 5.7 Independent validation: `08_cptac_validation.py` → `T7`

CPTAC proteomics for ten cohorts, tumour versus adjacent normal.

CPTAC is the strongest validation available without new experimental work because it is simultaneously: an **independent patient series**, an **orthogonal platform** (mass spectrometry rather than sequencing), a **different analyte** (protein rather than transcript), and **paired** (adjacent normal tissue from the same patients).

Normal samples are identified by the `.N` sample-ID suffix. Tumour and normal are compared by Wilcoxon rank-sum with Cliff's delta, and by Wilcoxon signed-rank where the same patient contributes both.

This does **not** discharge the reviewers' request for functional experiments. See §9.

### 5.8 Number provenance: `09_manuscript_numbers.py`

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

## 8. Analytical pitfalls this workflow is built to avoid

Each of the following is a way a pan-cancer expression analysis can produce a
confident but wrong answer. They are listed because the design decisions above
only make sense once the failure mode is visible.

**Reading stage-versus-normal annotations as a stage trend.** Several widely used
web tools plot expression by stage and annotate each stage against *normal*
tissue. Every stage can be significantly different from normal while expression
is flat or falling across stages. Testing a trend requires a test on the ordered
groups; see §5.3. Applied to NCL, a formal trend test supports a stage-ordered
increase in 2 of 17 cancers, where reading the annotations would suggest most.

**Reporting univariate survival only.** Proliferation is embedded in grade and
stage, so any proliferation-coupled gene will appear prognostic until those are
modelled. For NCL the count of significant cancers falls from 6/24 to 1/24 for
overall survival on adjustment (§5.4). Univariate screening is a reasonable first
pass; presenting it as evidence of independent prognostic value is not.

**Choosing a deconvolution algorithm per result.** Algorithms disagree in sign for
a third of cancer x cell-type combinations (§5.5). Selecting one algorithm per
comparison (particularly a different one for each) is indistinguishable from
selecting whichever gave the expected answer. Reporting concordance across all
available algorithms costs nothing and is far more interpretable.

**Ignoring the proliferation confounder in immune correlations.** A gene that
tracks proliferation will correlate with immune composition simply because
proliferative tumours differ immunologically. Any immune association claimed for
such a gene should survive explicit adjustment for a proliferation score
(§5.5); for NCL this removes most classical checkpoint associations while
leaving B7-H3, CD73 and CD39 intact.

**Treating the normal comparator as neutral.** Adjacent normal tissue is subject
to field effects and surgical ischaemia; GTEx tissue is post-mortem and from
different donors. For genes sensitive to metabolic state, which includes most
of the ribosome-biogenesis machinery, the choice can reverse the direction of
the result, as it does for NCL in kidney and thyroid (§5.1).

**Correlating across pooled cancers.** Pan-cancer pooled correlations are
dominated by tissue-of-origin differences in composition and describe no
within-tumour relationship. All correlations here are computed within cancer
(§5.6).

**Ranking genes that are not expressed.** Genes undetected in a cohort produce
degenerate correlations that GSEA orders arbitrarily; before filtering, ~10% of
ranked genes carried tied statistics (§5.6).

**Assuming an assumption was checked.** The proportional-hazards diagnostic in
this pipeline initially parsed nothing and silently reported zero tested models
while appearing to succeed. It was caught by inspecting the output column, not
by an error. Absence of an exception is not evidence that a check ran; §7
therefore verifies the statistical utilities against reference implementations
before they are used.


## 9. What this workflow does not do

**It generates no new experimental data.** Both reviewers asked for validation in patient cohorts, tissue samples, or functional in vitro experiments. CPTAC addresses the *cohort* element (independent patients, orthogonal platform, protein-level, paired normals) and substantiates the overexpression claim. It does not address the *functional* element at all.

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
| `S1_cohorts_and_tissue_mapping.tsv` | Cohort names, sample counts, TCGA-GTEx map | Supplementary Table S1 |
| `SUPPLEMENTARY_INDEX.md` | Manuscript S-number to filename mapping | all supplementary |

---
