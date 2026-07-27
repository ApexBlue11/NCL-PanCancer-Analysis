# Pan-Cancer Analysis of Nucleolin (NCL): Broad Overexpression, Limited Independent Prognostic Value, and Association with an Immune-Excluded Phenotype and Tumour-Intrinsic Immunosuppressive Ligands

Kruthika Prakash^1^, Janani Balaji^1^, Surya Babu^1^, Ramya Lakshmi Rajendran^2,3,4^, Prakash Gangadaran^2,3,4,\*^, ArulJothi Kandasamy Nagarajan^1,\*^ and Byeong-Cheol Ahn^2,3,4,5,\*^

^1^Department of Genetic Engineering, College of Engineering and Technology, SRM Institute of Science and Technology, Kattankulathur, Chengalpattu, Tamil Nadu, India – 603203
^2^BK21 FOUR KNU Convergence Educational Program of Biomedical Sciences for Creative Future Talents, Department of Biomedical Sciences, School of Medicine, Kyungpook National University, Daegu 41944, Republic of Korea
^3^Department of Nuclear Medicine, School of Medicine, Kyungpook National University, Daegu 41944, Republic of Korea
^4^Cardiovascular Research Institute, Kyungpook National University, Daegu 41944, Republic of Korea
^5^Department of Nuclear Medicine, Kyungpook National University Hospital, Daegu 41944, Republic of Korea

**\*Correspondence:**
Prof. ArulJothi Kandasamy Nagarajan, aruljotn@srmist.edu.in; Tel: +91-784-520-6014
Prof. Prakash Gangadaran, prakashg@knu.ac.kr; Tel: +82-53-420-4914
Prof. Byeong-Cheol Ahn, abc2000@knu.ac.kr; Tel: +82-53-420-5583

---

## Significance statement

Nucleolin (NCL) is widely proposed as a pan-cancer biomarker and therapeutic target. In a systematic reanalysis of 9,358 tumours we confirm that NCL is broadly overexpressed at both transcript and protein level, but show that it is only rarely an independent prognostic factor once age, sex, stage and grade are accounted for. NCL instead marks a proliferative, immune-excluded tumour phenotype and co-expresses with tumour-intrinsic immunosuppressive ligands, most consistently B7-H3 (CD276), rather than with T-cell checkpoint receptors. These findings refine, and in places correct, the prevailing view of NCL as a universal prognostic biomarker.

## Abstract

**Objectives.** Nucleolin (NCL) is widely described as a pan-cancer prognostic biomarker and immunomodulatory target, but this has not been tested with adjustment for clinical covariates or multiple comparisons. We asked whether NCL carries independent prognostic information and how it relates to the tumour immune microenvironment.

**Methods.** We analysed 9,358 TCGA tumours across 33 cancer types against 7,262 GTEx and 727 adjacent normal samples using uniformly reprocessed transcriptomes. Differential expression used Wilcoxon tests with Cliff's delta and stage trends the Jonckheere–Terpstra test. Three survival endpoints were modelled by multivariable Cox regression adjusting for age, sex, stage and grade, with Schoenfeld-residual checks. Immune infiltration was estimated by seven deconvolution algorithms as purity-adjusted partial correlations; checkpoint correlations were additionally proliferation-adjusted. All test families were Benjamini–Hochberg corrected, and findings validated in nine CPTAC proteomic cohorts.

**Results.** NCL was overexpressed in 24 of 29 evaluable cancers (Cliff's delta up to +0.93, q<0.05) and reduced in ovarian carcinoma (−0.27, q=1.1×10⁻⁴). A stage-ordered increase occurred in only 2 of 17 cancers. After adjustment and FDR correction, NCL remained independently prognostic in one cancer per endpoint: KIRP (overall survival HR 2.12, 95% CI 1.37–3.29; disease-specific HR 2.90, 1.60–5.27) and ACC (progression-free HR 2.41, 1.40–4.16). Seven of nine CPTAC cohorts confirmed elevated NCL protein. NCL correlated negatively with microenvironment score (20/20 cancers) and positively with B7-H3 in 21 of 33 cancers, whereas CTLA-4, PD-1 and TIM-3 were largely uncorrelated.

**Conclusion.** NCL is broadly overexpressed but seldom an independent prognostic biomarker, associating with tumour-intrinsic immunosuppressive ligands and an immune-excluded phenotype rather than T-cell checkpoints.

**Keywords:** nucleolin; NCL; pan-cancer analysis; B7-H3; tumour immune microenvironment; prognostic biomarker; The Cancer Genome Atlas

---

## 1. Introduction

Cancer remains among the leading causes of morbidity and mortality worldwide, and the identification of molecules that are both mechanistically informative and clinically actionable continues to be a central objective of translational oncology.^1,2^ Large public consortia — most notably The Cancer Genome Atlas (TCGA), the Genotype-Tissue Expression project (GTEx) and the Clinical Proteomic Tumor Analysis Consortium (CPTAC) — have made it feasible to evaluate a candidate molecule across the full spectrum of human malignancy rather than one tumour type at a time.^3–5^ This capability has produced a large literature of "pan-cancer biomarker" reports. It has also produced a recognised methodological problem: such analyses frequently rely on univariate association testing, apply no correction for the many thousands of comparisons performed, and report statistical significance without effect sizes, so that findings which are real but negligible are presented as clinically meaningful.^6,7^

Nucleolin (NCL) is a highly conserved, multifunctional phosphoprotein that is predominantly nucleolar but also shuttles to the nucleoplasm, cytoplasm and cell surface.^8,9^ It participates in ribosome biogenesis, ribosomal RNA processing, chromatin remodelling, mRNA stability and the cellular stress response, and it is required for the high rates of ribosome production that proliferating cells demand.^10,11^ Cell-surface NCL has attracted particular attention because it is accessible to targeting agents: aptamers such as AS1411, the pseudopeptide N6L and the F3 peptide all exploit surface NCL, and several have entered early-phase evaluation.^12–14^ NCL has been reported as overexpressed and prognostically adverse in individual malignancies including endometrial carcinoma, lung cancer, hepatocellular carcinoma and pancreatic ductal adenocarcinoma, and a meta-analysis has argued that its subcellular localisation determines its prognostic direction.^15–18^

Two questions nevertheless remain open. First, it is unclear whether NCL carries prognostic information that is *independent* of the clinical variables oncologists already use. NCL expression is tightly coupled to proliferation, and proliferation is itself associated with tumour grade and stage; an association between NCL and survival may therefore be a restatement of established prognostic factors rather than an addition to them. Distinguishing these possibilities requires multivariable modelling, which the existing NCL literature has largely not performed. Second, NCL has been proposed to shape the tumour immune microenvironment. Nucleolin targeting reduces immunosuppression in pancreatic cancer models, and the MDK–NCL axis has been implicated in an immunosuppressive niche in endometrial carcinoma.^14,19^ Whether a relationship between NCL and immune phenotype exists consistently across cancers, in which direction, and whether it can be separated from the generic association between proliferation and immune content, has not been established.

A further consideration is methodological and specific to immune analyses. Transcriptome deconvolution algorithms differ in their reference signatures, in whether their outputs are comparable between cell types or only between samples, and in their sensitivity to tumour purity.^20–24^ Estimates for nominally the same cell type can therefore disagree, including in sign. When a single algorithm is selected for each result, the reported associations may reflect that choice rather than the underlying biology. Reporting several algorithms and quantifying their concordance is a more conservative and more interpretable approach.

The present study addresses these questions systematically. Using uniformly reprocessed transcriptomes so that tumour and normal tissue are directly comparable, we characterise NCL expression across 33 cancer types against both GTEx and adjacent normal references; we test stage association with a test appropriate to ordered categories; we model three survival endpoints with multivariable Cox regression including formal proportional-hazards diagnostics; we quantify immune associations across seven deconvolution algorithms with purity adjustment and report their concordance; we test the immune-checkpoint associations directly, with adjustment for proliferation as a potential confounder; and we seek independent confirmation in CPTAC proteomic cohorts. Every family of tests is corrected for multiple comparisons and every association is reported with an effect size and confidence interval.

This work is a substantially revised and expanded reanalysis. Several conclusions differ from those we previously drew, and we identify and correct specific errors in our earlier analysis; these are set out explicitly in Section 4.5 rather than left implicit.

## 2. Materials and Methods

All analyses were performed in Python 3.14. Analysis code, intermediate result tables and the exact software environment are openly available (see Data availability), and the workflow is documented step by step in the accompanying repository.

### 2.1 Expression data and cohort definition

Transcriptome data were obtained from the UCSC Xena Toil recompute of TCGA, TARGET and GTEx (dataset `TcgaTargetGtex_rsem_gene_tpm`), comprising 60,498 genes across 19,131 samples quantified as log₂(TPM + 0.001).^3,25^ This resource was chosen deliberately: TCGA and GTEx were originally processed with different pipelines, and comparing them directly introduces batch effects that can exceed the biological differences of interest. The Toil recompute applies a single alignment and quantification pipeline to both, which is what makes a tumour-versus-GTEx comparison defensible.

Samples were assigned to 33 TCGA cohorts using the Xena phenotype annotation. Tumour samples were those designated "Primary Tumor" or, for LAML, "Primary Blood Derived Cancer – Peripheral Blood". Two classes of normal comparator were defined: TCGA adjacent normal tissue ("Solid Tissue Normal"), and GTEx normal tissue matched to each cancer's tissue of origin. Each TCGA–GTEx pairing was annotated with a match quality of *good*, *approximate* or *none*; approximate pairings (CHOL–liver, HNSC–salivary gland, READ–colon, SARC–adipose, DLBC–spleen, LAML–bone marrow) are reported but flagged, and cancers with no acceptable GTEx counterpart (MESO, THYM, UVM) were analysed against adjacent normals only. The full mapping with quality annotations is given in Supplementary Table S1 and in `scripts/cohorts.py`.

Gene identifiers were mapped from GENCODE v23 to gene symbols. Where a symbol mapped to several identifiers, the row with the highest mean expression was retained for single-gene analyses and the most variable row for correlation analyses.

### 2.2 Clinical, genomic and immune annotation

Clinical data, survival endpoints, tumour mutational burden and microsatellite instability scores for all 32 TCGA PanCancer Atlas studies were retrieved programmatically from cBioPortal.^26,27^ Survival endpoints are the curated definitions of the TCGA Clinical Data Resource — overall survival (OS), disease-specific survival (DSS) and progression-free survival (PFS) — which harmonise endpoint definitions across studies and are recommended in preference to raw TCGA follow-up fields.^28^ Covariates extracted were age at diagnosis, sex, AJCC pathological stage and histological grade. Tumour mutational burden (non-synonymous), MANTIS microsatellite-instability score, aneuploidy score and fraction of genome altered were also retrieved.^29^

Immune infiltration estimates for TCGA samples were downloaded from TIMER2.0, which provides pre-computed estimates from seven algorithms: TIMER, CIBERSORT, CIBERSORT-ABS, quanTIseq, MCP-counter, xCell and EPIC.^20–24,30^ xCell immune, stromal and microenvironment scores were taken from the same resource. Tumour purity was represented by the EPIC "uncharacterised cell" fraction, which is EPIC's explicit estimate of the compartment that is neither immune nor stromal.

### 2.3 Differential expression and stage association

For each cancer, NCL expression in tumours was compared with each normal comparator using two-sided Wilcoxon rank-sum tests. Effect sizes are reported as Cliff's delta with 95% percentile bootstrap confidence intervals (2,000 resamples), and additionally as Hedges' *g*. A non-parametric effect size was preferred because expression distributions are skewed and several cohorts have very small numbers of adjacent normals. Bootstrap intervals were used in preference to Cliff's asymptotic variance because that approximation is unreliable at the sample sizes some cohorts provide.

Association with pathological stage was assessed with the Jonckheere–Terpstra test for monotonic trend across ordered stages I–IV, using the tie-corrected normal approximation, alongside a Kruskal–Wallis test for any difference among stages. The distinction matters: a Kruskal–Wallis result establishes only that stages differ, not that expression rises with stage, and comparisons of each stage against normal tissue — as displayed by several web tools — do not test a stage-ordered trend at all. Cancers were included when at least three stage groups contained ≥5 patients and the cohort totalled ≥40 staged patients.

### 2.4 Survival analysis

For each cancer and each endpoint (OS, DSS, PFS), NCL expression was z-scored within cohort so that hazard ratios are expressed per standard deviation and are comparable across cancers. Three models were fitted: a log-rank test comparing patients above and below the cohort median (the univariate analysis reported in most prior work); a univariate Cox model on continuous NCL; and a multivariable Cox model including age, sex, stage and grade. Covariates were retained only where available for >80% of the cohort and showing more than one distinct value, and the covariates actually used are recorded for every model. Cancers were modelled when ≥30 patients and ≥10 events were available.

The proportional-hazards assumption was tested for every multivariable model using Schoenfeld residuals with rank-transformed time.^31^ Models in which the NCL term violated the assumption (p<0.05) are identified explicitly in Results and in Supplementary Table S3, because a hazard ratio from such a model summarises a time-varying effect and should not be read as a single constant risk multiplier.

### 2.5 Immune microenvironment analysis

For every combination of cancer, cell type and algorithm, the association between NCL expression and estimated infiltration was quantified as a Spearman correlation and as a partial Spearman correlation adjusting for tumour purity, implemented by regressing the ranks of both variables on the ranks of purity and correlating the residuals (t-test on n−3 degrees of freedom).

To assess robustness to algorithm choice, ten cell types resolvable by more than one algorithm were designated canonical. For each cancer and canonical cell type we recorded how many algorithms yielded a significant purity-adjusted association after FDR correction and whether they agreed in sign. A combination was called *concordant* when at least two algorithms were significant and none disagreed in direction, and *conflicting* when significant estimates of opposite sign were both present.

Associations between NCL and immune, stromal and microenvironment scores, tumour mutational burden, MANTIS MSI score, aneuploidy score and fraction of genome altered were computed as Spearman correlations with Fisher-*z* confidence intervals.

### 2.6 Immune checkpoint analysis

Sixteen immune checkpoint and immunomodulatory genes were examined: *CD274* (PD-L1), *PDCD1LG2* (PD-L2), *PDCD1* (PD-1), *CTLA4*, *HAVCR2* (TIM-3), *LAG3*, *TIGIT*, *IDO1*, *BTLA*, *VSIR* (VISTA), *SIGLEC15*, *CD276* (B7-H3), *IL10*, *TGFB1*, *ENTPD1* (CD39) and *NT5E* (CD73). For each cancer, the Spearman correlation between NCL and each gene was computed and reported in three forms: unadjusted; adjusted for tumour purity; and adjusted for proliferation.

The proliferation adjustment is essential to interpretation. NCL is a ribosome-biogenesis gene whose expression tracks proliferative rate, and proliferative tumours differ systematically in immune composition. A raw correlation between NCL and a checkpoint gene may therefore reflect nothing more specific than shared covariation with proliferation. The proliferation score was defined as the mean within-cohort z-score of ten canonical markers (*MKI67*, *PCNA*, *TOP2A*, *CCNB1*, *BUB1*, *AURKA*, *CDK1*, *TYMS*, *RRM2*, *TK1*). Associations described below as robust are those significant after both purity and proliferation adjustment.

### 2.7 Gene set enrichment analysis

For each cancer, all genes were ranked by their within-cancer Spearman correlation with NCL and analysed by pre-ranked GSEA against the MSigDB Hallmark and Reactome collections (v2024.1.Hs), with 1,000 permutations and gene sets restricted to 15–500 genes.^32–35^ Correlations were computed within cancers rather than across pooled samples: pooled pan-cancer correlations are dominated by tissue-of-origin differences in composition and do not describe any within-tumour relationship.

Ranked lists were restricted to genes expressed in the cancer concerned, defined as TPM > 1 in at least 25% of that cohort's tumours. Without this restriction approximately 10% of genes carried tied ranking statistics, because genes that are essentially undetected yield degenerate correlations that GSEA then orders arbitrarily; the restriction reduced ties to below 0.2%.

### 2.8 Independent proteomic validation

Protein-level validation used CPTAC cohorts accessed through the `cptac` Python package: BRCA, ccRCC, COAD, GBM, HNSCC, LSCC, LUAD, OV, PDAC and UCEC.^5^ CPTAC provides an independent patient series, measured by mass spectrometry rather than sequencing, at the protein rather than the transcript level, with adjacent normal tissue from the same patients — making it the most stringent validation available without new experimental work. Tumour and normal samples were compared by Wilcoxon rank-sum test with Cliff's delta, and by Wilcoxon signed-rank test where the same patient contributed both.

### 2.9 Multiple testing and reporting

Benjamini–Hochberg false discovery rate correction was applied within each family of tests: across cancers within each differential-expression comparator; across cancers within each survival endpoint and model type; across all cancer × cell type × algorithm infiltration tests; across all cancer × gene checkpoint tests; across all cancer × measure genomic-score tests; and across all GSEA tests. Adjusted values are reported as *q*. Significance was defined as q<0.05 throughout. Every association is reported with a point estimate, a 95% confidence interval and the sample size on which it rests.

### 2.10 Reproducibility

Analyses are deterministic given fixed random seeds (bootstrap resampling, GSEA permutation). Statistical utility functions were verified against brute-force implementations and against `statsmodels` before use, including a null-calibration check confirming a 5.7% type-I error rate for the trend test under the null and confirmation that a non-monotonic pattern is not reported as a trend. The complete workflow, environment specification and verification suite are provided in the accompanying repository.

## 3. Results

### 3.1 NCL is overexpressed across most cancers, and the choice of normal comparator matters

Across 29 cancers with an acceptable GTEx comparator, NCL was significantly overexpressed in tumours in 24 (all q<0.05), showing no significant difference in four (ACC, CESC, BLCA, KICH) and significantly reduced expression in one, ovarian serous cystadenocarcinoma (Cliff's delta −0.27, 95% CI −0.36 to −0.17, q=1.1×10⁻⁴) (Figure 1a; Table 1). Effect sizes were large in several cancers: LGG (delta +0.93, 95% CI +0.91 to +0.94, q=1.3×10⁻¹⁹⁹), CHOL (+0.92, +0.85 to +0.97), TGCT (+0.91, +0.84 to +0.96), READ (+0.90, +0.82 to +0.96), PAAD (+0.86, +0.80 to +0.91) and COAD (+0.86, +0.81 to +0.90).

An important qualification emerged when the same comparison was repeated against TCGA adjacent normal tissue (Figure 1b). Twelve cancers were significant against both comparators; ten agreed in direction, but two did not. In KIRP, NCL was higher than GTEx kidney (delta +0.63, q=5.1×10⁻⁸) yet lower than adjacent normal kidney (delta −0.29, q=1.1×10⁻²). The same reversal occurred in THCA (+0.37 versus −0.27), and the discordance extended to KICH (+0.19, not significant, versus −0.80, q=1.6×10⁻⁸) and KIRC (+0.78 versus −0.07, not significant). The kidney and thyroid are therefore cancers in which the reported direction of NCL dysregulation is determined by the reference tissue chosen rather than by the tumour. Adjacent normal kidney is not histologically inert — it is subject to field effects, inflammation and, in nephrectomy specimens, ischaemic injury — whereas GTEx tissue is post-mortem, and ribosome-biogenesis genes such as NCL are sensitive to both. We therefore report both comparators throughout and restrict confident claims of overexpression to cancers where they agree.

This finding also resolves a discrepancy with our previously reported analysis, which described NCL as reduced in KIRP, READ and cutaneous melanoma. In the present analysis READ (delta +0.90, q=2.2×10⁻³⁸) and SKCM (+0.65, q=4.3×10⁻²⁵) are unambiguously elevated against both available comparators, and only KIRP reproduces a reduction, and only against adjacent normals (Section 4.5).

### 3.2 A stage-ordered increase in NCL is the exception, not the rule

Across 17 cancers with sufficient staging data, a significant monotonic trend across stages I–IV was present in only three, and in one of these the direction was negative: LIHC (Jonckheere–Terpstra z=+3.47, q=0.005), LUAD (z=+3.45, q=0.005) and, decreasing, KIRC (z=−2.95, q=0.018) (Figure 1c; Table 2). Fourteen cancers showed no stage-ordered trend.

This directly revises a claim in our earlier report that NCL "progressively increased as the cancer advanced through stages I–IV" in BRCA, COAD, ESCA, HNSC, KICH, LIHC, STAD and READ. Of those eight, only LIHC is supported by a formal trend test. The earlier claim rested on visual inspection of plots in which the significance annotations compare each stage against normal tissue, not against the preceding stage; such annotations can be uniformly significant while expression is flat or falling across stages. We note additionally that even in LIHC the stage IV group contains only six patients, so that the trend is driven by stages I–III and the late-versus-early effect size is not distinguishable from zero (delta −0.01, 95% CI −0.53 to +0.53).

### 3.3 NCL is rarely an independent prognostic factor

Univariate analyses reproduced the pattern that has led NCL to be described as a prognostic biomarker: NCL was associated with overall survival in 6 of 24 cancers, with disease-specific survival in 7 of 21, and with progression-free survival in 8 of 29 (q<0.05).

After adjustment for age, sex, stage and grade, and FDR correction, this largely did not survive (Figure 2; Table 3). NCL remained independently associated with outcome in a single cancer per endpoint: KIRP for overall survival (HR 2.12 per SD, 95% CI 1.37–3.29, q=0.019, 38 events) and disease-specific survival (HR 2.90, 1.60–5.27, q=0.0099, 25 events), and ACC for progression-free survival (HR 2.41, 1.40–4.16, q=0.044, 38 events). The proportional-hazards assumption held for the NCL term in both, supporting interpretation of these hazard ratios as constant effects. Across all endpoints the assumption was violated in 10 of 74 models — consistently in KIRC and LUAD — and those hazard ratios are flagged in Table 3 and should not be read as time-invariant.

Re-examining the specific survival claims of our earlier report is instructive. Of four cancers previously reported as showing worse overall survival with high NCL, one is confirmed (KIRP, adjusted HR 2.12, q=0.019); one is attenuated to non-significance after adjustment despite a strong univariate signal (LIHC, univariate q=0.006, adjusted HR 1.32, 95% CI 1.07–1.62, q=0.079); one is null (HNSC, adjusted HR 1.06, 95% CI 0.91–1.25, q=0.72); and one could not be modelled at all because TCGA KICH contains too few events. The general conclusion is that NCL's apparent prognostic value in most cancers is largely explained by variables already used clinically, and that treating it as a broadly applicable independent prognostic biomarker is not supported.

### 3.4 NCL marks an immune-excluded tumour phenotype

Across 3,910 cancer × cell type × algorithm tests, 1,229 purity-adjusted associations were significant at q<0.05. The dominant pattern was negative. NCL correlated inversely with the xCell microenvironment score in all 20 cancers where the association was significant, with the immune score in 18 of 20, and with the stromal score in 14 of 15 (Figure 3a; Table 4). The strongest inverse immune-score associations were in GBM (rho −0.38, 95% CI −0.52 to −0.22), SKCM (−0.36, −0.53 to −0.17), LUSC (−0.35, −0.43 to −0.27) and ESCA (−0.34, −0.47 to −0.19). Only PRAD and UVM showed positive associations.

Cell-type-level results were consistent with this: among the strongest concordant associations were reduced CD8⁺ T cells in SKCM (median rho −0.32 across 5 of 7 algorithms) and ACC (−0.31, 4 of 7), reduced macrophages in GBM (−0.46), KIRP (−0.39), LUSC (−0.33) and OV (−0.32), and reduced cancer-associated fibroblasts in PRAD (−0.46). NCL-high tumours are thus, in general, immune- and stroma-poor — a phenotype consistent with a proliferative, tumour-cell-rich compartment.

The one prominent positive association in our earlier report, between NCL and neutrophil infiltration in thyroid carcinoma, is directionally reproduced but substantially smaller than previously stated: median rho +0.30 across 4 of 6 algorithms after purity adjustment, compared with the rho=0.539 obtained from MCP-counter alone without adjustment.

### 3.5 Deconvolution algorithms frequently disagree

Of 330 cancer × canonical cell type combinations, only 83 (25%) were concordant, meaning at least two algorithms gave a significant purity-adjusted association with no disagreement in sign. In 109 combinations (33%) algorithms produced significant estimates of *opposite* direction for nominally the same cell type in the same cancer (Figure 3b).

This has a direct methodological implication. A pan-cancer immune association supported by a single deconvolution algorithm has roughly a one-in-three chance of being contradicted by another algorithm applied to the same data. Results in this section are therefore reported only where multiple algorithms agree, and we regard single-algorithm immune associations — including several in our earlier report, which used a different algorithm for each panel — as insufficiently robust to support biological interpretation.

### 3.6 NCL correlates with tumour-intrinsic immunosuppressive ligands, not T-cell checkpoints

Testing 495 cancer × gene associations produced a clearly structured result (Figure 4; Table 5). After adjustment for both tumour purity and proliferation, the most consistent association by a wide margin was with *CD276* (B7-H3): positive and significant in 21 of 33 cancers, with no cancer showing a significant negative association (median rho +0.34; ACC +0.58, PRAD +0.50, GBM +0.47, THCA +0.47, KIRC +0.46). Associations with the two principal components of the adenosine axis followed the same pattern: *NT5E* (CD73) robust in 17 of 33 cancers (16 positive) and *ENTPD1* (CD39) in 17 of 33 (16 positive).

By contrast, the T-cell checkpoint receptors were weakly and inconsistently associated: PD-1 robust in 3 of 33 cancers, CTLA-4 in 4 of 33 (median rho +0.02), TIGIT in 5, LAG-3 in 6 (all negative), and TIM-3 in 9 of 33 (median rho +0.05). PD-L1 was robust in 7 of 33 (median rho +0.17) and TGF-β1 in 6 of 33.

This distinction is biologically coherent. B7-H3, CD73 and CD39 are expressed by tumour cells themselves, whereas PD-1, CTLA-4, LAG-3 and TIGIT are expressed principally by infiltrating lymphocytes. A tumour-cell-expressed gene such as NCL would be expected to covary with the former and not with the latter — and, given the inverse relationship between NCL and immune content described above, to correlate negatively rather than positively with lymphocyte-restricted transcripts. That the observed pattern matches this expectation, and persists after adjustment for proliferation, argues that it is not merely an artefact of tumour cellularity.

These results also require us to withdraw a specific claim from our earlier report. That report's title and abstract asserted correlations between NCL and PD-L1, CTLA-4, TIM-3, IL-10 and TGF-β. No such analysis was in fact performed in that work. Having now performed it, we find these particular associations to be weak or absent, while a robust association with B7-H3 — not previously proposed — emerges in their place.

### 3.7 NCL is largely unrelated to tumour mutational burden and microsatellite instability

NCL showed little relationship with genomic instability measures: tumour mutational burden was significantly associated in 7 of 33 cancers (median rho +0.06), MANTIS MSI score in 10 of 32 (median rho +0.04, with equal numbers of positive and negative associations), and aneuploidy score in 7 of 33 (median rho +0.04). Fraction of genome altered was slightly more often associated (11 of 33, median rho +0.12). NCL therefore does not track the mutational processes that predict immune-checkpoint blockade response, which further distinguishes its immune association from a conventional immunogenicity signal.

### 3.8 Independent proteomic validation

In nine CPTAC cohorts with adjacent normal tissue, NCL protein was significantly more abundant in tumour in seven (Figure 5; Table 6): LUAD (Cliff's delta +0.99, q=5.1×10⁻³⁵, 102 matched pairs), LSCC (+0.97), GBM (+0.99), OV (+0.92), COAD (+0.89), ccRCC (+0.65) and HNSCC (+0.53). Paired analyses within patients gave the same conclusion wherever pairing was available.

Two cohorts diverged from the transcript-level results and are reported as such. In UCEC, NCL protein did not differ between tumour and normal (delta −0.04, q=0.69) despite modest transcript elevation. In PDAC, NCL protein was significantly *lower* in tumour (delta −0.28, 95% CI −0.43 to −0.13, q=3.2×10⁻⁴) whereas the transcript was strongly elevated in TCGA PAAD (delta +0.86). Pancreatic tumours have a high stromal content and low epithelial cellularity, and bulk transcript and bulk protein measurements weight these compartments differently; discordance of this kind is a known feature of pancreatic proteogenomic data. We therefore do not claim protein-level overexpression of NCL in pancreatic cancer.

Notably, the ccRCC proteomic result (+0.65) supports the GTEx-based direction in kidney cancer rather than the adjacent-normal direction, providing independent evidence relevant to the comparator discordance described in Section 3.1.

### 3.9 NCL is coupled to proliferation and RNA metabolism in every cancer

Pre-ranked GSEA of genes ranked by their within-cancer correlation with NCL produced 32,491 tests across 32 cancers, of which 9,349 were significant at q<0.05 (Figure 6; Supplementary Table S7). The result was strikingly uniform.

In the Hallmark collection, the G2M checkpoint signature was positively enriched in **all 32 cancers** (median NES +1.91), MYC targets V1 in 31 of 32 (+1.83), the mitotic spindle signature in 31 of 32 (+1.60), E2F targets in 30 of 32 (+1.91) and MYC targets V2 in 26 of 32 (+1.78). Further positively enriched sets included mTORC1 signalling (21 of 32), PI3K–AKT–mTOR signalling (20 of 32), the unfolded protein response (22 of 32) and DNA repair (14 of 32). Only one signature was consistently depleted, and in a minority of cancers (coagulation, 3 of 32).

Reactome results were concordant and dominated by nuclear RNA metabolism: transport of mature transcripts to the cytoplasm, processing of capped intron-containing pre-mRNA, snRNP assembly, nuclear tRNA processing, gene silencing by RNA, SUMOylation of RNA-binding and DNA-repair proteins, chromatin-modifying enzymes and histone acetylation were each positively enriched in 31 of 32 cancers (median NES +1.7 to +2.0).

Two features of this result matter for interpretation. First, it recovers NCL's established biology — ribosome biogenesis, RNA processing, chromatin function — from an unsupervised, correlation-ranked analysis, which supports the validity of the ranking statistic. Second, and more importantly, **no immune or inflammatory signature appears among the consistently enriched pathways in either collection.** The pathways associated with NCL across cancers are those of a proliferating, biosynthetically active cell. This is the strongest single piece of evidence that NCL's immune associations (Sections 3.4–3.6) follow from the proliferative state of the tumour compartment rather than from a dedicated immunoregulatory programme, and it is why the checkpoint analysis in Section 3.6 was adjusted for proliferation before any association was described as robust.

We note that TGF-β signalling was positively enriched in 22 of 32 cancers, which is of interest given that our earlier report asserted an NCL–TGF-β association. At the level of the individual gene, however, *TGFB1* was robustly correlated with NCL in only 6 of 33 cancers (Section 3.6); the pathway-level and gene-level results should not be conflated.

## 4. Discussion

### 4.1 What this analysis establishes

Three conclusions follow from this reanalysis. First, NCL overexpression in cancer is real, widespread and reproducible: it is elevated in 24 of 29 evaluable cancers at the transcript level, with large effect sizes, and confirmed at the protein level in an independent patient series in seven of nine cohorts. This is the best-supported claim in the NCL literature and it survives rigorous testing.

Second, that overexpression translates into independent prognostic information far less often than the literature implies. In 74 multivariable models spanning three endpoints, NCL remained significant after covariate adjustment and FDR correction in two cancers. This is not a null result — the KIRP association is substantial and consistent across two endpoints — but it is a considerably narrower claim than "pan-cancer prognostic biomarker". The discrepancy with univariate analyses is precisely what one expects of a proliferation-associated gene: proliferation is embedded in grade and stage, so a marker of proliferation will appear prognostic until those variables are included, and then will often cease to.

Third, the relationship between NCL and the tumour immune microenvironment is real but is close to the opposite of what has been proposed. NCL-high tumours are immune- and stroma-poor. Where NCL does correlate positively with an immune-relevant molecule, that molecule is characteristically expressed by tumour cells rather than by lymphocytes.

### 4.2 B7-H3 as the principal immune association

The most robust immune finding is the association between NCL and B7-H3 (CD276), positive in 21 of 33 cancers after adjustment for both purity and proliferation, and negative in none. B7-H3 is an attractive comparison for NCL: both are broadly overexpressed across carcinomas, both reach the tumour-cell surface, and both are being pursued as surface-targeting opportunities. B7-H3 is currently the subject of extensive clinical development, including antibody–drug conjugates and CAR-T approaches, and is notable for being expressed in tumours that are not inflamed — the immune-excluded setting in which conventional PD-1/PD-L1 blockade performs poorly.

The co-occurrence of NCL and B7-H3 in an immune-excluded context, together with the parallel associations with CD73 and CD39, suggests that NCL-high tumours may be characterised less by classical T-cell exhaustion than by an adenosine-mediated, tumour-cell-intrinsic mode of immune evasion. We emphasise that this is an association-level observation and a hypothesis for testing, not a demonstrated mechanism; nothing in these data establishes that NCL regulates B7-H3 or that either is causal for the immune phenotype. The relevant experiment — perturbing NCL and measuring B7-H3 and immune composition — has not been performed here.

### 4.3 Pathway context

The enrichment results are unusually consistent and place the immune findings in context. Across 32 cancers, the programmes associated with NCL are those of proliferation and nuclear RNA metabolism: the G2M checkpoint signature in every cancer without exception, MYC and E2F target sets in almost all, and Reactome sets covering mRNA processing and export, snRNP assembly, tRNA processing and chromatin modification in 31 of 32.

That this analysis recovers NCL's canonical biology from an unsupervised, correlation-ranked procedure is reassuring about the method. What it adds is the absence: no immune or inflammatory programme is consistently enriched in either collection. Taken with the negative correlations against immune, stromal and microenvironment scores, this supports reading NCL as a marker of the biosynthetic and proliferative state of the tumour compartment, whose immune associations follow from that state rather than from a dedicated immunoregulatory function.

This is also why we regard the B7-H3 result as the most interesting finding rather than a corollary of proliferation. Proliferation is the dominant axis of NCL covariation, so an immune association that persists after explicit adjustment for a proliferation score — as the B7-H3, CD73 and CD39 associations do, while most T-cell checkpoint associations do not — is not simply a restatement of that axis.

We deliberately did not repeat the miRNA network analysis of our earlier report. That analysis described a network in which NCL was one node among many as a set of "NCL-interacting miRNAs", assigned biological roles to individual miRNAs that were never tested, and drew on a prediction resource named inconsistently between Methods and Results. We do not consider it to have supported any conclusion, and no replacement analysis is offered in its place.

### 4.4 Methodological implications

Two methodological observations from this work generalise beyond NCL.

The first concerns deconvolution. One third of the cancer × cell type combinations we examined produced significant estimates of *opposite* sign depending on the algorithm used. Since a typical pan-cancer immune analysis reports one algorithm, and since algorithms are often selected after inspecting results, the field's immune-association literature is likely to contain a substantial number of associations that would reverse under a different but equally defensible choice. Reporting concordance across algorithms costs little and is considerably more informative than a single point estimate.

The second concerns the choice of normal comparator. In kidney and thyroid cancers the direction of NCL dysregulation depends on whether GTEx or adjacent normal tissue is used. Neither is simply correct: adjacent normal tissue shares the patient's genetic background and exposures but is subject to field effects and surgical ischaemia, while GTEx tissue is unaffected by the tumour but post-mortem and from different individuals. For genes sensitive to metabolic and ischaemic state — which includes most of the ribosome-biogenesis machinery — this choice is consequential and should be reported rather than made silently.

### 4.5 Relationship to our earlier report and corrections

Because this manuscript revises conclusions we previously published in preprint form, we state the specific corrections explicitly.

1. **Immune checkpoint claims.** The earlier title and abstract asserted correlations between NCL and PD-L1, CTLA-4, TIM-3, IL-10 and TGF-β. No checkpoint analysis was performed in that work. That analysis has now been done; those specific associations are weak or absent, and B7-H3 emerges as the robust association instead. The title of the present manuscript has been changed accordingly, and uses "association" rather than "regulation" because the data are correlational.

2. **CancerSEA and cBioPortal.** The earlier report referred in its title and text to CancerSEA functional-state analysis and to cBioPortal genetic-alteration analysis. Neither was performed. References to CancerSEA have been removed; cBioPortal is now genuinely used, as the source of curated survival endpoints and clinical covariates.

3. **Direction of dysregulation.** The earlier report described NCL as reduced in KIRP, READ and cutaneous melanoma. READ and SKCM are elevated against every available comparator. KIRP is reduced only against adjacent normals and elevated against GTEx.

4. **Stage progression.** The earlier claim of a progressive stage I–IV increase across eight cancers is supported by a formal trend test in one.

5. **Survival panels.** Three of the eight Kaplan–Meier panels in the earlier report carried cohort sizes incompatible with the cancer they were labelled with — most clearly a "KICH" panel showing 877 patients, whereas TCGA KICH comprises 65 primary tumours. Those panels derived from pooled Human Protein Atlas cohorts (renal, lung) rather than from the individual TCGA projects named. All survival analyses here are computed directly from TCGA data with stated sample sizes.

6. **Cancer type nomenclature.** The earlier Table 1, presented as the TCGA cancer list, contained 39 entries including entities that are not TCGA projects, and defined ACC as adenoid cystic carcinoma whereas TCGA ACC is adrenocortical carcinoma. The corrected list of 33 TCGA projects is given in Supplementary Table S1.

7. **Differential expression method.** The earlier Methods described GEPIA2 with GTEx normals, but the figure presented was generated by TIMER2's TCGA-only module; cancers such as PAAD therefore rested on four adjacent normals. The present analysis uses 167 GTEx pancreas samples for that comparison.

We report these corrections in full because the conclusions of the earlier analysis have been cited, and because the errors are of a kind — tool misattribution, figure mislabelling, claims outrunning analyses — that is difficult for readers to detect from the published record alone.

## 5. Limitations

The central limitation is that this study generates no new experimental data. Reviewers of an earlier version reasonably asked for validation in patient cohorts, tissue samples or functional in vitro experiments. We have addressed this as far as is possible computationally, by validating in an independent patient series measured on an orthogonal platform at the protein level with paired adjacent normal tissue (CPTAC). This is a genuine and reasonably stringent form of external validation, and it substantiates the overexpression claim. It is not, however, a substitute for functional experimentation, and none of the mechanistic hypotheses advanced in Section 4.2 is tested by it. The association between NCL and B7-H3 in particular requires perturbation experiments — NCL knockdown or overexpression with measurement of B7-H3 and of immune composition — before any regulatory relationship can be claimed. We state plainly that this work does not provide that evidence.

Further limitations follow from the data sources. Transcriptome deconvolution infers cell composition from bulk expression and is not a measurement of infiltration; the poor cross-algorithm concordance documented in Section 3.5 is itself evidence of this limitation, and single-cell or spatial data would be required to resolve the cell-type associations convincingly. Six TCGA–GTEx tissue pairings are approximate and are flagged as such; three cancers have no acceptable GTEx counterpart. Several cohorts are small — CHOL (36 tumours), KICH (65), UCS (57) — so that the absence of a significant association in those cancers reflects limited power as much as absence of effect, and KICH in particular could not be modelled for overall survival. Stage and grade are unavailable for some cancers, so the covariate set differs between models; the covariates used are recorded for every model. The proportional-hazards assumption was violated for the NCL term in 10 of 74 models, and those hazard ratios summarise time-varying effects. Finally, TCGA is predominantly composed of treatment-naive primary resections from a limited set of populations, and prognostic associations derived from it need not transfer to treated, metastatic or differently ascertained populations.

## 6. Conclusion

Nucleolin is broadly and reproducibly overexpressed across human cancers at both transcript and protein level, but it is only rarely an independent prognostic factor once the clinical variables already in routine use are accounted for; on the evidence presented here, kidney renal papillary cell carcinoma is the clearest exception. Its relationship to the tumour immune microenvironment is consistent but is not the one previously proposed: NCL marks an immune-excluded phenotype and co-expresses with tumour-cell-intrinsic immunosuppressive ligands, most consistently B7-H3 and the adenosine-pathway ectoenzymes CD73 and CD39, rather than with T-cell checkpoint receptors. These results argue for a more circumscribed view of NCL as a prognostic biomarker, and identify the NCL–B7-H3 relationship in immune-excluded tumours as the hypothesis most deserving of experimental test.

## Declarations

**Ethics approval and consent to participate:** Not applicable. This study analysed only publicly available, de-identified data from TCGA, GTEx and CPTAC, obtained under the data use policies of those consortia.

**Clinical trial number:** Not applicable.

**Consent for publication:** Not applicable.

**Availability of data and materials:** All primary data are publicly available. Transcriptome data: UCSC Xena Toil recompute (https://xenabrowser.net/datapages/, dataset `TcgaTargetGtex_rsem_gene_tpm`). Clinical, survival, TMB and MSI data: cBioPortal (https://www.cbioportal.org), TCGA PanCancer Atlas studies. Immune infiltration estimates: TIMER2.0 (http://timer.cistrome.org). Gene sets: MSigDB v2024.1.Hs (https://www.gsea-msigdb.org). Proteomic data: CPTAC via the `cptac` Python package. All analysis code, the complete computational environment specification, intermediate result tables and a step-by-step methodology document are available at [REPOSITORY URL]. No new data were generated.

**Competing interests:** The authors declare that the research was conducted without any commercial or financial relationships that could be construed as a potential conflict of interest.

**Funding:** This work was supported by the National Research Foundation of Korea (NRF) grant funded by the Korea government (MSIT) (NRF-2022R1A2C2005057).

**Authors' contributions:** [To be confirmed by all authors for the revised manuscript, whose conclusions differ substantively from the previous version.]

**Use of AI and AI-assisted technologies:** An AI assistant (Anthropic Claude) was used during preparation of this revision for analysis code development, statistical implementation, execution of the computational workflow, and drafting and editing of manuscript text. All analyses were specified, inspected and verified by the authors; all statistical implementations were validated against reference implementations prior to use; and the authors take full responsibility for the integrity of the data, the accuracy of the analyses and the content of the manuscript. No AI system is an author, and no AI system was used to generate or fabricate data.

**Acknowledgements:** The authors thank the reviewers of the previous version of this manuscript, whose criticisms prompted the reanalysis reported here and directly identified several of the errors corrected in Section 4.5.

## Figure legends

**Figure 1. NCL expression across cancers and pathological stages.** (a) Cliff's delta with 95% bootstrap confidence intervals for NCL expression in tumours versus GTEx normal tissue, ordered by effect size; point size indicates cohort size, and colour indicates TCGA–GTEx tissue match quality. Filled markers denote q<0.05 (Benjamini–Hochberg). (b) Comparison of effect sizes obtained against GTEx normals and against TCGA adjacent normals; points off the diagonal indicate comparator-dependent results, with KIRP, KICH, KIRC and THCA labelled. (c) Jonckheere–Terpstra standardised trend statistic for NCL across pathological stages I–IV; positive values indicate increasing expression with stage. Only cancers with ≥3 stage groups of ≥5 patients are shown. Sample sizes are given for all groups.

**Figure 2. Survival analysis.** (a) Forest plot of hazard ratios per standard deviation of NCL from multivariable Cox models adjusted for age, sex, stage and grade, for overall survival; covariates included and numbers of patients and events are annotated. (b) The corresponding disease-specific and progression-free survival results. (c) Comparison of univariate and covariate-adjusted hazard ratios, illustrating attenuation after adjustment. (d) Kaplan–Meier curves for the two cancers retaining significance after adjustment (KIRP, ACC), stratified at the cohort median. Models in which the proportional-hazards assumption was violated for the NCL term are marked.

**Figure 3. Immune infiltration and cross-algorithm concordance.** (a) Purity-adjusted partial Spearman correlations between NCL expression and immune cell infiltration across cancers, shown for each of seven deconvolution algorithms; only associations with q<0.05 are coloured. (b) Concordance across algorithms for ten canonical cell types: the number of algorithms yielding significant associations and whether they agree in direction. Combinations in which algorithms disagree in sign are highlighted. (c) Correlations between NCL and xCell immune, stromal and microenvironment scores across cancers.

**Figure 4. Immune checkpoint associations.** (a) Spearman correlations between NCL and 16 immune checkpoint and immunomodulatory genes across 33 cancers, before adjustment, after purity adjustment, and after proliferation adjustment. (b) Number of cancers in which each molecule shows a robust association, defined as significant after both adjustments, separated by direction. (c) Per-cancer detail for CD276 (B7-H3), NT5E (CD73) and ENTPD1 (CD39), with 95% confidence intervals.

**Figure 5. Independent proteomic validation in CPTAC.** Cliff's delta with 95% bootstrap confidence intervals for NCL protein abundance in tumour versus adjacent normal tissue across nine CPTAC cohorts, with numbers of tumours, normals and matched pairs annotated. Cohorts discordant with the transcript-level result (PDAC, UCEC) are indicated.

**Figure 6. Pathway associations.** Normalised enrichment scores from pre-ranked GSEA of genes ranked by within-cancer correlation with NCL, for Hallmark and Reactome gene sets showing consistent direction across cancers. Only gene sets significant at q<0.05 in at least three cancers are shown.

## Tables

**Table 1.** NCL differential expression by cancer type, against GTEx and TCGA adjacent normal comparators: sample sizes, median expression, Cliff's delta with 95% CI, Hedges' *g*, and FDR-adjusted p-values.

**Table 2.** Association between NCL expression and pathological stage: per-stage sample sizes, Kruskal–Wallis and Jonckheere–Terpstra statistics, and late-versus-early effect sizes.

**Table 3.** Multivariable Cox regression of NCL against overall, disease-specific and progression-free survival: hazard ratios per standard deviation with 95% CI, covariates included, numbers of patients and events, FDR-adjusted p-values, and proportional-hazards test results.

**Table 4.** NCL correlations with immune, stromal and microenvironment scores, tumour mutational burden, MANTIS MSI score, aneuploidy score and fraction of genome altered.

**Table 5.** NCL correlations with immune checkpoint and immunomodulatory genes, unadjusted and adjusted for tumour purity and proliferation.

**Table 6.** Independent validation of NCL protein abundance in nine CPTAC cohorts.

**Supplementary Table S1.** Corrected list of the 33 TCGA study abbreviations with full names, and the TCGA–GTEx normal tissue mapping with match quality annotations.

**Supplementary Tables S2–S7.** Complete per-test result tables for differential expression, stage association, survival, immune infiltration across all algorithms, checkpoint associations, and GSEA.
