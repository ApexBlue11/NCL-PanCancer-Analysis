# Run manifest

Generated 2026-07-27T19:32:58Z on Windows-11-10.0.26200-SP0, Python 3.14.3.

Regenerate with `python scripts/11_manifest.py`. Diff two manifests to localise any disagreement in results to an input file, a package version or the code.

## Package versions

| package | version |
|---|---|
| numpy | 2.4.3 |
| scipy | 1.17.1 |
| pandas | 2.3.3 |
| statsmodels | 0.14.6 |
| lifelines | 0.30.3 |
| gseapy | 1.3.1 |
| cptac | 1.5.14 |
| matplotlib | 3.10.8 |
| requests | 2.32.5 |

## Random seeds

| procedure | seed | draws | location |
|---|---|---|---|
| cliffs_delta_bootstrap | 0 | 2000 | `statsutil.cliffs_delta` |
| gsea_permutation | 0 | 1000 | `07_gsea.py -> gseapy.prerank` |
| jonckheere_permutation | 0 | 0 | `statsutil.jonckheere_terpstra (asymptotic p reported; permutation optional)` |

## Raw inputs

| file | bytes | sha256 (first 16) |
|---|---|---|
| `TcgaTargetGtex_rsem_gene_tpm.gz` | 1,323,254,426 | `a8c36cb16ef82ecc` |
| `TcgaTargetGTEX_phenotype.txt.gz` | 135,753 | `ba4d4461cff0fe5e` |
| `gencode.v23.annotation.gene.probemap` | 3,244,244 | `6783ea58791ae876` |
| `infiltration_estimation_for_tcga.csv.gz` | 8,743,177 | `581b2969642c20f5` |
| `h.all.v2024.1.Hs.symbols.gmt` | 48,690 | `ee2463540042078b` |
| `c2.cp.reactome.v2024.1.Hs.symbols.gmt` | 852,693 | `9ea1b5e656597daf` |

## Result tables

| table | rows | cols | sha256 (first 16) |
|---|---|---|---|
| `S1_cohorts_and_tissue_mapping.tsv` | 33 | 7 | `0fb3c0cc1fd368eb` |
| `SUPPLEMENTARY_INDEX.md` | 30 | 1 | `7f96db381919fa5f` |
| `T10_gsea_per_cancer.tsv.gz` | 32491 | 13 | `313a80daace3dd82` |
| `T1_differential_expression.tsv` | 66 | 19 | `c68c5e4fb6f28c30` |
| `T2_stage_association.tsv` | 17 | 17 | `e89ca282bad54263` |
| `T3_survival.tsv` | 84 | 24 | `c18c9de3bedd935a` |
| `T4_immune_infiltration.tsv` | 3910 | 12 | `dac63ea8926f17c9` |
| `T5_algorithm_concordance.tsv` | 330 | 9 | `8e5e2029c0783f98` |
| `T6_checkpoints.tsv` | 495 | 15 | `99bb890b190665d3` |
| `T7_cptac_validation.tsv` | 10 | 19 | `d0a0e4e53951ac53` |
| `T8_genomic_scores.tsv` | 230 | 8 | `a4939e4a197e2deb` |
| `T9_ncl_gene_correlations.tsv.gz` | 58581 | 33 | `4c5acec60db282d9` |

## Live services (not checksummable at source)

- **cBioPortal** — Live API with no version pin. Survival endpoints are the TCGA-CDR definitions; clinical values are stable in practice but the service is not archival. The retrieved table is checksummed below as clinical_pancan.tsv.
- **CPTAC** — Datasets are downloaded at run time by the package; the package version below determines the release retrieved.
