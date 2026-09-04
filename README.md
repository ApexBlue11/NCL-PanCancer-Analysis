# Nucleolin (NCL) across 33 human cancers

Reproducible pan-cancer analysis of Nucleolin: expression, prognosis, immune microenvironment and pathway associations across 9,358 TCGA tumours, with GTEx normal references and independent CPTAC protein-level validation.

Code and result tables for the manuscript *"Nucleolin overexpression correlates with poor prognosis and immune checkpoint regulation across various cancer types: Insights from The Cancer Genome Atlas and GTEx analyses"*

---

## What this study found

**NCL co-expresses with B7-H3 (CD276), not with T-cell checkpoints.**
Across 33 cancers, NCL correlates positively with *CD276* in **21 of 33** cancers after adjusting for both tumour purity and proliferation, with no cancer showing a significant inverse association (rho up to +0.58). The adenosine-pathway ectoenzymes *NT5E* (CD73) and *ENTPD1* (CD39) follow the same pattern (17 of 33 each). By contrast, the T-cell checkpoint receptors are largely uncorrelated: PD-1 in 3 of 33 cancers, CTLA-4 in 4, TIGIT in 5, LAG-3 in 6.

This distinction is biologically coherent. B7-H3, CD73 and CD39 are expressed by tumour cells; PD-1, CTLA-4, LAG-3 and TIGIT are expressed by infiltrating lymphocytes. A tumour-cell-expressed gene such as NCL tracks the former, not the latter.

**NCL marks an immune-excluded phenotype.**
NCL correlates negatively with the xCell microenvironment score in **all 20** cancers where the association is significant, with immune score in 18 of 20, and with stromal score in 14 of 15. NCL-high tumours are immune- and stroma-poor.

**Total NCL transcript carries little independent prognostic information.**
Across 74 multivariable Cox models spanning overall, disease-specific and progression-free survival, NCL survived covariate adjustment and FDR correction in two: KIRP (OS HR 2.12, 95% CI 1.37–3.29; DSS HR 2.90, 1.60–5.27) and ACC (PFS HR 2.41, 1.40–4.16).

This is consistent with, rather than contrary to, the immunohistochemical literature. The largest prior meta-analysis (12 studies, 1,221 patients, 8 cancers) found *cytoplasmic* NCL adverse (HR 4.32) but *nuclear* NCL protective (HR 0.42). Bulk transcriptome measurement sums those compartments, so a near-null total-transcript effect is what the localisation hypothesis predicts. The practical implication: prognostic use of NCL requires an assay that resolves subcellular localisation.

**NCL is broadly overexpressed, and confirmed at protein level.**
Elevated in 24 of 29 evaluable cancers (Cliff's delta up to +0.93), reduced in ovarian carcinoma. Confirmed in **7 of 9** independent CPTAC proteomic cohorts with paired adjacent normals (LUAD delta +0.99, LSCC +0.97, GBM +0.99, OV +0.91, COAD +0.89, CCRCC +0.65, HNSCC +0.52). PDAC and UCEC did not confirm and are reported as such.

**NCL is proliferation-coupled in every cancer.**
GSEA on genes ranked by within-cancer correlation with NCL: the Hallmark G2M checkpoint signature is enriched in **32 of 32** cancers, MYC targets in 31, E2F targets in 30; Reactome results are dominated by mRNA processing, snRNP assembly and chromatin modification (31 of 32). No immune signature is consistently enriched, which is why every checkpoint association here is reported after explicit proliferation adjustment.

### Two methodological results

**Deconvolution algorithms disagree more than is generally reported.** Of 330 cancer × cell-type combinations, only 25% were concordant across algorithms, and **33% produced significant associations of opposite sign** depending on which algorithm was used. A pan-cancer immune association supported by one algorithm has roughly a one-in-three chance of being contradicted by another applied to the same data.

**The normal comparator can determine the direction of the result.** In KIRP, KICH, KIRC and THCA, whether NCL appears up- or down-regulated depends on whether GTEx or adjacent normal tissue is the reference. Reports of NCL being reduced in renal or thyroid tumours are reconcilable with reports of elevation; both can be correct given different references.

### What is new here

To our knowledge this is the first pan-cancer analysis of NCL (`nucleolin AND (pan-cancer OR TCGA) AND immune infiltration` returns no prior study), and the NCL–B7-H3 association has not previously been reported (`nucleolin AND (B7-H3 OR CD276)` returns one unrelated review). The immune-excluded phenotype, the comparator dependence, and the cross-algorithm discordance quantification are likewise not previously described for this gene.

### What this does not show

Nothing here demonstrates that NCL *regulates* B7-H3, or that either is causal for the immune phenotype. These are correlational analyses of bulk tumour data. Establishing a regulatory relationship requires perturbation (NCL knockdown or overexpression with measurement of B7-H3 and immune composition) which this study does not perform. See [METHODOLOGY.md §9](METHODOLOGY.md).

---

## Quick start

```bash
pip install -r requirements.txt
python scripts/test_statsutil.py       # verify statistical utilities (prints ALL PASS)
python scripts/00_check_sources.py     # confirm all data URLs are reachable
```

Then run the pipeline in order:

```bash
python scripts/01_download.py
python scripts/02_fetch_clinical.py
python scripts/03_build_matrix.py
python scripts/04_expression.py
python scripts/05_survival.py
python scripts/06_immune.py
python scripts/07a_expression_filter.py
python scripts/07_gsea_correlations.py
python scripts/07_gsea.py --workers 4
python scripts/08_cptac_validation.py
python scripts/09_manuscript_numbers.py
python scripts/10_figures.py
python scripts/12_supplementary_tables.py
python scripts/11_manifest.py
```

Every script is idempotent (it skips work whose output already exists) so an interrupted run resumes rather than restarting.

**Requirements:** Python 3.14, ~6 GB free disk, <1 GB RAM. Runtime ≈ 90 minutes, dominated by GSEA. All inputs are public; no credentials required.

---

## Repository layout

```
├── README.md                  you are here
├── METHODOLOGY.md             every analytical decision, script and output
├── requirements.txt           pinned environment
├── scripts/                   the pipeline
└── results/
    ├── tables/                all per-test results under working names
    ├── supplementary/         S1-S11 as cited in the article, plus a
    │                          multi-sheet workbook (Supplementary_Tables_S1-S11.xlsx)
    ├── figures/               600 dpi PNG, TIFF (LZW) and vector PDF
    ├── MANIFEST.json/.md      input checksums, package versions, seeds
    └── MANUSCRIPT_NUMBERS.txt every number quoted in the manuscript
```

`data/` is generated by the pipeline and not tracked (~6 GB, of which 4.6 GB is the expression matrix); it is fully reproducible from `scripts/01_download.py`.

### Where the numbers come from

The result tables in `results/tables/` are committed, so `10_figures.py` will
redraw every figure on a fresh clone without downloading anything. That is a
convenience for inspection, **not** a reproduction: it re-plots numbers that
were computed elsewhere. The scripts say which of the two you are getting —

```
  [data] NOTE: no raw source data found (looked in .../data/raw).
         Falling back to the result tables distributed with the
         repository. Output below is RE-PLOTTED from committed
         numbers, not recomputed from TCGA/GTEx/CPTAC.
         To reproduce from source, run scripts 01-08 first.
```

and the steps that cannot run on committed tables at all
(`09_manuscript_numbers.py`, `12_supplementary_tables.py`, and everything from
`03` to `08`) stop with a message naming the missing files rather than failing
somewhere inside a reader.

To go from the databases to the tables, figures and conclusions, run the
pipeline in order from `01_download.py`. Steps 01–08 do the analysis; 09–12
only report it.

**Locating the data.** Raw inputs are looked for in the first of these that
holds them: `$NCL_DATA`, `<repo>/data`, `./data`, `scripts/data`, then the
working directory itself. So an existing 6 GB download can be reused in place:

```bash
NCL_DATA=/mnt/big-disk/ncl-data python scripts/04_expression.py
```

| Script | Purpose |
|---|---|
| `00_check_sources.py` | Probe every data URL before downloading |
| `01_download.py` | Download raw data (resumable) |
| `02_fetch_clinical.py` | cBioPortal clinical, survival, TMB, MSI |
| `03_build_matrix.py` | Stream expression into a float32 matrix |
| `04_expression.py` | Tumour vs normal; stage-ordered trend |
| `05_survival.py` | Log-rank, univariate and multivariable Cox, PH diagnostics |
| `06_immune.py` | Infiltration (7 algorithms), checkpoints, scores, TMB/MSI |
| `07a_expression_filter.py` | Per-cohort expressed-gene filter for GSEA ranking |
| `07_gsea_correlations.py` | NCL vs all genes, within each cancer |
| `07_gsea.py` | Pre-ranked GSEA, resumable, parallel |
| `08_cptac_validation.py` | CPTAC protein-level validation |
| `09_manuscript_numbers.py` | Emit every manuscript number, traced to source |
| `10_figures.py` | Publication figures at 600 dpi |
| `12_supplementary_tables.py` | Supplementary Tables S1-S11 and the workbook |
| `11_manifest.py` | Input checksums, package versions, seeds |
| `cohorts.py` / `data_io.py` / `statsutil.py` | Cohort definitions, loaders, statistics |
| `figstyle.py` | Shared colour encoding and the legend-collision audit |
| `test_statsutil.py` / `test_figstyle.py` | Verification suites for `statsutil` and the legend audit |

---

## Data sources

All public; none require registration.

| Dataset | Source |
|---|---|
| Expression (TCGA/TARGET/GTEx, Toil recompute) | [UCSC Xena](https://xenabrowser.net/datapages/) |
| Clinical, survival, TMB, MSI | [cBioPortal](https://www.cbioportal.org) PanCancer Atlas |
| Immune infiltration (7 algorithms) | [TIMER2.0](http://timer.cistrome.org) |
| Hallmark and Reactome gene sets | [MSigDB](https://www.gsea-msigdb.org) v2024.1.Hs |
| Proteomics | [CPTAC](https://proteomics.cancer.gov) via the `cptac` package |

Survival endpoints are the curated TCGA Clinical Data Resource definitions (Liu et al., *Cell* 2018), not raw TCGA follow-up fields.

---

## Statistical approach

- **Effect sizes throughout.** Cliff's delta with bootstrap 95% CIs; Hedges' *g* alongside; hazard ratios per standard deviation with 95% CIs; Spearman rho with Fisher-*z* CIs.
- **Multiple testing.** Benjamini–Hochberg FDR within each family of tests; significance at q<0.05.
- **Confounder control.** Immune correlations are purity-adjusted; checkpoint correlations additionally adjusted for a ten-gene proliferation score, because NCL is proliferation-coupled and an unadjusted correlation may simply restate that.
- **Assumption checking.** Proportional hazards tested via Schoenfeld residuals for every Cox model; violations reported and flagged on the figures (10 of 74 models; neither significant result affected).
- **Trend testing.** Stage association uses Jonckheere–Terpstra, which tests monotonic trend across ordered stages rather than comparing each stage against normal tissue.

`test_statsutil.py` verifies each utility against brute-force or published reference implementations, including a null-calibration check (5.7% type-I rate at α=0.05) and confirmation that a non-monotonic pattern is not reported as a trend.

---

## Reproducibility

- All stochastic steps seeded (`seed=0`): bootstrap CIs (2,000 resamples), GSEA permutation (1,000 permutations).
- Each script writes to disk before the next reads, so any step can be re-run in isolation.
- `09_manuscript_numbers.py` regenerates every number quoted in the manuscript from the result tables, so a changed analysis surfaces stale text rather than leaving it silently wrong.
- Environment pinned in `requirements.txt`.
- `11_manifest.py` writes `results/MANIFEST.json` and `MANIFEST.md` recording the SHA-256, size and retrieval date of every input file, the shape and checksum of every result table, and the resolved version of every package that affects a number. This matters because the sources are not archival: UCSC Xena, TIMER2.0 and MSigDB reissue files at stable URLs, cBioPortal serves a live API, and `cptac` downloads at run time. Diffing two manifests localises any disagreement to a specific input, package version or code change.

## Citation

Prakash K, Balaji J, Babu S, et al. Nucleolin overexpression correlates with poor
prognosis and immune checkpoint regulation across various cancer types: Insights
from The Cancer Genome Atlas and GTEx analyses. *Eurasian J Med Oncol*. Published
online September 4, 2026. doi:10.36922/EJMO026220241

<https://doi.org/10.36922/EJMO026220241>

```bibtex
@article{Prakash2026Nucleolin,
  author  = {Kruthika Prakash and Janani Balaji and Surya Babu and
             Ramya Lakshmi Rajendran and Prakash Gangadaran and
             ArulJothi Kandasamy Nagarajan and Byeong-Cheol Ahn},
  title   = {Nucleolin overexpression correlates with poor prognosis and immune
             checkpoint regulation across various cancer types: Insights from
             The Cancer Genome Atlas and {GTEx} analyses},
  journal = {Eurasian Journal of Medicine and Oncology},
  year    = {2026},
  note    = {Published online 4 September 2026},
  doi     = {10.36922/EJMO026220241},
  url     = {https://accscience.com/journal/EJMO/articles/online_first/9151},
}
```



## License

Code released under the MIT License. Underlying data remain subject to the terms of TCGA, GTEx, cBioPortal, TIMER2.0, MSigDB and CPTAC.
