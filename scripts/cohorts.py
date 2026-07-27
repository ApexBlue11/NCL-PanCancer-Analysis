"""TCGA cohort definitions and TCGA<->GTEx normal-tissue mapping.

`match` records how well the GTEx tissue substitutes for the tumour's true
tissue of origin. Only 'good' pairings carry the primary tumour-vs-normal
claim; 'approximate' ones are reported but flagged in the manuscript, and
'none' are restricted to TCGA adjacent normals.

This is the correction for the original manuscript, which described GEPIA2 +
GTEx in Methods while the figure it presented was TIMER2's TCGA-only module --
leaving cancers such as PAAD resting on 4 adjacent normals.
"""

# TCGA code -> (Xena detailed_category, GTEx _primary_site or None, match quality)
COHORTS = {
    "ACC":  ("Adrenocortical Cancer",                    "Adrenal Gland", "good"),
    "BLCA": ("Bladder Urothelial Carcinoma",             "Bladder",       "good"),
    "BRCA": ("Breast Invasive Carcinoma",                "Breast",        "good"),
    "CESC": ("Cervical & Endocervical Cancer",           "Cervix Uteri",  "good"),
    "CHOL": ("Cholangiocarcinoma",                       "Liver",         "approximate"),
    "COAD": ("Colon Adenocarcinoma",                     "Colon",         "good"),
    "DLBC": ("Diffuse Large B-Cell Lymphoma",            "Spleen",        "approximate"),
    "ESCA": ("Esophageal Carcinoma",                     "Esophagus",     "good"),
    "GBM":  ("Glioblastoma Multiforme",                  "Brain",         "good"),
    "HNSC": ("Head & Neck Squamous Cell Carcinoma",      "Salivary Gland","approximate"),
    "KICH": ("Kidney Chromophobe",                       "Kidney",        "good"),
    "KIRC": ("Kidney Clear Cell Carcinoma",              "Kidney",        "good"),
    "KIRP": ("Kidney Papillary Cell Carcinoma",          "Kidney",        "good"),
    "LAML": ("Acute Myeloid Leukemia",                   "Bone Marrow",   "approximate"),
    "LGG":  ("Brain Lower Grade Glioma",                 "Brain",         "good"),
    "LIHC": ("Liver Hepatocellular Carcinoma",           "Liver",         "good"),
    "LUAD": ("Lung Adenocarcinoma",                      "Lung",          "good"),
    "LUSC": ("Lung Squamous Cell Carcinoma",             "Lung",          "good"),
    "MESO": ("Mesothelioma",                             None,            "none"),
    "OV":   ("Ovarian Serous Cystadenocarcinoma",        "Ovary",         "good"),
    "PAAD": ("Pancreatic Adenocarcinoma",                "Pancreas",      "good"),
    "PCPG": ("Pheochromocytoma & Paraganglioma",         "Adrenal Gland", "good"),
    "PRAD": ("Prostate Adenocarcinoma",                  "Prostate",      "good"),
    "READ": ("Rectum Adenocarcinoma",                    "Colon",         "approximate"),
    "SARC": ("Sarcoma",                                  "Adipose Tissue","approximate"),
    "SKCM": ("Skin Cutaneous Melanoma",                  "Skin",          "good"),
    "STAD": ("Stomach Adenocarcinoma",                   "Stomach",       "good"),
    "TGCT": ("Testicular Germ Cell Tumor",               "Testis",        "good"),
    "THCA": ("Thyroid Carcinoma",                        "Thyroid",       "good"),
    "THYM": ("Thymoma",                                  None,            "none"),
    "UCEC": ("Uterine Corpus Endometrioid Carcinoma",    "Uterus",        "good"),
    "UCS":  ("Uterine Carcinosarcoma",                   "Uterus",        "good"),
    "UVM":  ("Uveal Melanoma",                           None,            "none"),
}

# Full names, corrected against the official TCGA study abbreviations.
# The submitted Table 1 was a CancerSEA list mislabelled "all TCGA cancers":
# it defined ACC as adenoid cystic carcinoma (TCGA ACC is adrenocortical),
# used non-TCGA codes (HCC/RCC/CM/CCC/UEC/AML) and non-TCGA entities
# (AST, CML, HGG, NSCLC, ODG, MEL), and listed 39 rows for 33 projects.
FULL_NAME = {
    "ACC":  "Adrenocortical carcinoma",
    "BLCA": "Bladder urothelial carcinoma",
    "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "CHOL": "Cholangiocarcinoma",
    "COAD": "Colon adenocarcinoma",
    "DLBC": "Lymphoid neoplasm diffuse large B-cell lymphoma",
    "ESCA": "Esophageal carcinoma",
    "GBM":  "Glioblastoma multiforme",
    "HNSC": "Head and neck squamous cell carcinoma",
    "KICH": "Kidney chromophobe",
    "KIRC": "Kidney renal clear cell carcinoma",
    "KIRP": "Kidney renal papillary cell carcinoma",
    "LAML": "Acute myeloid leukemia",
    "LGG":  "Brain lower grade glioma",
    "LIHC": "Liver hepatocellular carcinoma",
    "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous cell carcinoma",
    "MESO": "Mesothelioma",
    "OV":   "Ovarian serous cystadenocarcinoma",
    "PAAD": "Pancreatic adenocarcinoma",
    "PCPG": "Pheochromocytoma and paraganglioma",
    "PRAD": "Prostate adenocarcinoma",
    "READ": "Rectum adenocarcinoma",
    "SARC": "Sarcoma",
    "SKCM": "Skin cutaneous melanoma",
    "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular germ cell tumors",
    "THCA": "Thyroid carcinoma",
    "THYM": "Thymoma",
    "UCEC": "Uterine corpus endometrial carcinoma",
    "UCS":  "Uterine carcinosarcoma",
    "UVM":  "Uveal melanoma",
}

# Immune checkpoint and immunomodulatory genes. The submitted abstract asserted
# correlations with PD-L1, CTLA-4, TIM-3, IL-10 and TGF-beta without analysing
# any of them; these are the genes that claim actually requires.
CHECKPOINTS = {
    "CD274":    "PD-L1",
    "PDCD1LG2": "PD-L2",
    "PDCD1":    "PD-1",
    "CTLA4":    "CTLA-4",
    "HAVCR2":   "TIM-3",
    "LAG3":     "LAG-3",
    "TIGIT":    "TIGIT",
    "IDO1":     "IDO1",
    "BTLA":     "BTLA",
    "VSIR":     "VISTA",
    "SIGLEC15": "Siglec-15",
    "CD276":    "B7-H3",
    "IL10":     "IL-10",
    "TGFB1":    "TGF-beta1",
    "ENTPD1":   "CD39",
    "NT5E":     "CD73",
}

# Proliferation markers. NCL is a ribosome-biogenesis/proliferation gene, so any
# checkpoint correlation must be tested against proliferation as a confounder
# before it can be described as immune-specific.
PROLIFERATION = ["MKI67", "PCNA", "TOP2A", "CCNB1", "BUB1", "AURKA",
                 "CDK1", "TYMS", "RRM2", "TK1"]

TUMOUR_TYPES = ["Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"]
TCGA_NORMAL = "Solid Tissue Normal"
GTEX_NORMAL = "Normal Tissue"
