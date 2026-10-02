#!/usr/bin/env python3
# Purpose: Build the terminology for the 33 TCGA cancer categories.
# Input: No external input file
# Output: TCGA_Cancer_Vocabulary.csv; TCGA_Cancer_PubMed_Terms_Long.csv; TCGA_Cancer_UI_Aliases_Long.csv

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import re
import sys
from collections import defaultdict

import pandas as pd

OUT_MAIN = "TCGA_Cancer_Vocabulary.csv"
OUT_PUBMED = "TCGA_Cancer_PubMed_Terms_Long.csv"
OUT_UI = "TCGA_Cancer_UI_Aliases_Long.csv"
OUT_REPORT = "TCGA_Cancer_Vocabulary_Integrity_Report.txt"

EXPECTED_TCGA_CODES = {
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS",
    "UVM",
}

EXPECTED_COUNT = 33

SOURCE_TCGA = (
    "https://www.cancer.gov/ccg/research/genome-sequencing/tcga"
)
SOURCE_CODES = (
    "https://gdc.cancer.gov/resources-tcga-users/"
    "tcga-code-tables/tcga-study-abbreviations"
)
SOURCE_PANCAN = (
    "https://gdc.cancer.gov/about-data/publications/pancanatlas"
)

CANCERS = [
    {
        "code": "LAML",
        "tcga_name": "Acute Myeloid Leukemia",
        "display": "Acute myeloid leukemia",
        "pubmed": [
            "acute myeloid leukemia",
            "acute myelogenous leukemia",
        ],
        "ui": [
            "AML", "leukemia", "myeloid leukemia",
        ],
    },
    {
        "code": "ACC",
        "tcga_name": "Adrenocortical carcinoma",
        "display": "Adrenocortical carcinoma",
        "pubmed": [
            "adrenocortical carcinoma",
            "adrenal cortical carcinoma",
            "adrenocortical cancer",
        ],
        "ui": [
            "adrenal", "adrenal cortex", "adrenocortical", "carcinoma",
        ],
    },
    {
        "code": "BLCA",
        "tcga_name": "Bladder Urothelial Carcinoma",
        "display": "Bladder cancer",
        "pubmed": [
            "bladder urothelial carcinoma",
            "urothelial carcinoma of the bladder",
            "bladder cancer",
            "bladder carcinoma",
        ],
        "ui": [
            "bladder", "urothelial", "urothelial carcinoma", "carcinoma",
        ],
    },
    {
        "code": "LGG",
        "tcga_name": "Brain Lower Grade Glioma",
        "display": "Lower grade glioma",
        "pubmed": [
            "lower grade glioma",
            "lower-grade glioma",
            "low grade glioma",
            "low-grade glioma",
        ],
        "ui": [
            "glioma", "brain glioma", "low grade", "lower grade", "brain",
        ],
    },
    {
        "code": "BRCA",
        "tcga_name": "Breast invasive carcinoma",
        "display": "Breast cancer",
        "pubmed": [
            "breast cancer",
            "breast carcinoma",
            "breast adenocarcinoma",
        ],
        "ui": [
            "breast", "mammary", "carcinoma",
        ],
    },
    {
        "code": "CESC",
        "tcga_name": (
            "Cervical squamous cell carcinoma and endocervical adenocarcinoma"
        ),
        "display": "Cervical cancer",
        "pubmed": [
            "cervical cancer",
            "cervical carcinoma",
            "cervical squamous cell carcinoma",
            "endocervical adenocarcinoma",
            "cancer of the cervix",
        ],
        "ui": [
            "cervical", "cervix", "endocervical", "carcinoma",
            "adenocarcinoma",
        ],
    },
    {
        "code": "CHOL",
        "tcga_name": "Cholangiocarcinoma",
        "display": "Cholangiocarcinoma",
        "pubmed": [
            "cholangiocarcinoma",
            "bile duct cancer",
            "bile duct carcinoma",
        ],
        "ui": [
            "bile duct", "biliary", "cholangio", "carcinoma",
        ],
    },
    {
        "code": "COAD",
        "tcga_name": "Colon adenocarcinoma",
        "display": "Colon cancer",
        "pubmed": [
            "colon cancer",
            "colon carcinoma",
            "colon adenocarcinoma",
            "colonic adenocarcinoma",
        ],
        "ui": [
            "colon", "colonic", "colorectal", "CRC", "carcinoma",
            "adenocarcinoma",
        ],
    },
    {
        "code": "ESCA",
        "tcga_name": "Esophageal carcinoma",
        "display": "Esophageal cancer",
        "pubmed": [
            "esophageal cancer",
            "esophageal carcinoma",
            "oesophageal cancer",
            "oesophageal carcinoma",
            "esophageal adenocarcinoma",
            "oesophageal adenocarcinoma",
            "esophageal squamous cell carcinoma",
            "oesophageal squamous cell carcinoma",
        ],
        "ui": [
            "esophagus", "oesophagus", "esophageal", "oesophageal",
            "carcinoma", "adenocarcinoma",
        ],
    },
    {
        "code": "GBM",
        "tcga_name": "Glioblastoma multiforme",
        "display": "Glioblastoma",
        "pubmed": [
            "glioblastoma",
            "glioblastoma multiforme",
        ],
        "ui": [
            "GBM", "brain", "glioma", "glioblastoma",
        ],
    },
    {
        "code": "HNSC",
        "tcga_name": "Head and Neck squamous cell carcinoma",
        "display": "Head and neck squamous cell carcinoma",
        "pubmed": [
            "head and neck squamous cell carcinoma",
            "head and neck cancer",
            "HNSCC",
        ],
        "ui": [
            "head and neck", "HNSC", "HNSCC", "squamous",
            "carcinoma",
        ],
    },
    {
        "code": "KICH",
        "tcga_name": "Kidney Chromophobe",
        "display": "Kidney chromophobe carcinoma",
        "pubmed": [
            "kidney chromophobe",
            "chromophobe renal cell carcinoma",
            "renal chromophobe carcinoma",
            "chromophobe kidney carcinoma",
        ],
        "ui": [
            "kidney", "renal", "chromophobe", "carcinoma",
        ],
    },
    {
        "code": "KIRC",
        "tcga_name": "Kidney renal clear cell carcinoma",
        "display": "Clear cell renal cell carcinoma",
        "pubmed": [
            "kidney renal clear cell carcinoma",
            "clear cell renal cell carcinoma",
            "renal clear cell carcinoma",
            "clear cell kidney carcinoma",
        ],
        "ui": [
            "kidney", "renal", "clear cell", "ccRCC", "carcinoma",
        ],
    },
    {
        "code": "KIRP",
        "tcga_name": "Kidney renal papillary cell carcinoma",
        "display": "Papillary renal cell carcinoma",
        "pubmed": [
            "kidney renal papillary cell carcinoma",
            "papillary renal cell carcinoma",
            "renal papillary cell carcinoma",
            "papillary kidney carcinoma",
        ],
        "ui": [
            "kidney", "renal", "papillary", "pRCC", "carcinoma",
        ],
    },
    {
        "code": "LIHC",
        "tcga_name": "Liver hepatocellular carcinoma",
        "display": "Liver hepatocellular carcinoma",
        "pubmed": [
            "liver hepatocellular carcinoma",
            "hepatocellular carcinoma",
            "hepatoma",
        ],
        "ui": [
            "liver", "liver cancer", "hepatic", "hepatocellular",
            "HCC", "carcinoma",
        ],
    },
    {
        "code": "LUAD",
        "tcga_name": "Lung adenocarcinoma",
        "display": "Lung adenocarcinoma",
        "pubmed": [
            "lung adenocarcinoma",
            "pulmonary adenocarcinoma",
        ],
        "ui": [
            "lung", "pulmonary", "adenocarcinoma",
        ],
    },
    {
        "code": "LUSC",
        "tcga_name": "Lung squamous cell carcinoma",
        "display": "Lung squamous cell carcinoma",
        "pubmed": [
            "lung squamous cell carcinoma",
            "squamous cell carcinoma of the lung",
            "pulmonary squamous cell carcinoma",
        ],
        "ui": [
            "lung", "pulmonary", "squamous", "carcinoma",
        ],
    },
    {
        "code": "DLBC",
        "tcga_name": "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
        "display": "Diffuse large B-cell lymphoma",
        "pubmed": [
            "diffuse large B-cell lymphoma",
            "diffuse large B cell lymphoma",
            "DLBCL",
        ],
        "ui": [
            "lymphoma", "B-cell lymphoma", "B cell lymphoma",
            "diffuse large B-cell", "DLBCL",
        ],
    },
    {
        "code": "MESO",
        "tcga_name": "Mesothelioma",
        "display": "Mesothelioma",
        "pubmed": [
            "mesothelioma",
            "malignant mesothelioma",
            "pleural mesothelioma",
        ],
        "ui": [
            "mesothelioma", "pleural",
        ],
    },
    {
        "code": "OV",
        "tcga_name": "Ovarian serous cystadenocarcinoma",
        "display": "Ovarian cancer",
        "pubmed": [
            "ovarian serous cystadenocarcinoma",
            "ovarian serous carcinoma",
            "high grade serous ovarian cancer",
            "high-grade serous ovarian cancer",
            "ovarian cancer",
        ],
        "ui": [
            "ovarian", "ovary", "serous ovarian", "HGSOC",
            "carcinoma",
        ],
    },
    {
        "code": "PAAD",
        "tcga_name": "Pancreatic adenocarcinoma",
        "display": "Pancreatic cancer",
        "pubmed": [
            "pancreatic adenocarcinoma",
            "pancreatic ductal adenocarcinoma",
            "pancreatic cancer",
            "PDAC",
        ],
        "ui": [
            "pancreas", "pancreatic", "PDAC", "adenocarcinoma",
        ],
    },
    {
        "code": "PCPG",
        "tcga_name": "Pheochromocytoma and Paraganglioma",
        "display": "Pheochromocytoma and paraganglioma",
        "pubmed": [
            "pheochromocytoma",
            "phaeochromocytoma",
            "paraganglioma",
            "pheochromocytoma and paraganglioma",
            "phaeochromocytoma and paraganglioma",
        ],
        "ui": [
            "pheochromocytoma", "phaeochromocytoma", "pheo",
            "paraganglioma",
        ],
    },
    {
        "code": "PRAD",
        "tcga_name": "Prostate adenocarcinoma",
        "display": "Prostate cancer",
        "pubmed": [
            "prostate adenocarcinoma",
            "prostate cancer",
            "prostatic adenocarcinoma",
            "prostatic cancer",
        ],
        "ui": [
            "prostate", "prostatic", "adenocarcinoma",
        ],
    },
    {
        "code": "READ",
        "tcga_name": "Rectum adenocarcinoma",
        "display": "Rectal cancer",
        "pubmed": [
            "rectal cancer",
            "rectal carcinoma",
            "rectal adenocarcinoma",
            "rectum adenocarcinoma",
        ],
        "ui": [
            "rectum", "rectal", "colorectal", "CRC", "carcinoma",
            "adenocarcinoma",
        ],
    },
    {
        "code": "SARC",
        "tcga_name": "Sarcoma",
        "display": "Sarcoma",
        "pubmed": [
            "sarcoma",
            "soft tissue sarcoma",
        ],
        "ui": [
            "sarcoma", "soft tissue",
        ],
    },
    {
        "code": "SKCM",
        "tcga_name": "Skin Cutaneous Melanoma",
        "display": "Cutaneous melanoma",
        "pubmed": [
            "skin cutaneous melanoma",
            "cutaneous melanoma",
            "skin melanoma",
        ],
        "ui": [
            "skin", "cutaneous", "melanoma",
        ],
    },
    {
        "code": "STAD",
        "tcga_name": "Stomach adenocarcinoma",
        "display": "Stomach cancer",
        "pubmed": [
            "stomach adenocarcinoma",
            "gastric adenocarcinoma",
            "stomach cancer",
            "gastric cancer",
            "gastric carcinoma",
        ],
        "ui": [
            "stomach", "gastric", "adenocarcinoma", "carcinoma",
        ],
    },
    {
        "code": "TGCT",
        "tcga_name": "Testicular Germ Cell Tumors",
        "display": "Testicular germ cell cancer",
        "pubmed": [
            "testicular germ cell tumor",
            "testicular germ cell tumors",
            "testicular germ cell cancer",
            "testicular cancer",
        ],
        "ui": [
            "testicular", "testis", "germ cell",
        ],
    },
    {
        "code": "THYM",
        "tcga_name": "Thymoma",
        "display": "Thymoma",
        "pubmed": [
            "thymoma",
            "thymic epithelial tumor",
            "thymic epithelial tumour",
        ],
        "ui": [
            "thymoma", "thymic",
        ],
    },
    {
        "code": "THCA",
        "tcga_name": "Thyroid carcinoma",
        "display": "Thyroid cancer",
        "pubmed": [
            "thyroid carcinoma",
            "thyroid cancer",
            "papillary thyroid carcinoma",
            "papillary thyroid cancer",
        ],
        "ui": [
            "thyroid", "papillary thyroid", "carcinoma",
        ],
    },
    {
        "code": "UCS",
        "tcga_name": "Uterine Carcinosarcoma",
        "display": "Uterine carcinosarcoma",
        "pubmed": [
            "uterine carcinosarcoma",
            "uterine malignant mixed mullerian tumor",
            "uterine malignant mixed müllerian tumor",
            "malignant mixed mullerian tumor of the uterus",
            "malignant mixed müllerian tumor of the uterus",
        ],
        "ui": [
            "uterus", "uterine", "carcinosarcoma", "MMMT",
            "malignant mixed mullerian", "malignant mixed müllerian",
            "carcinoma", "sarcoma",
        ],
    },
    {
        "code": "UCEC",
        "tcga_name": "Uterine Corpus Endometrial Carcinoma",
        "display": "Endometrial cancer",
        "pubmed": [
            "uterine corpus endometrial carcinoma",
            "endometrial cancer",
            "endometrial carcinoma",
            "uterine endometrial carcinoma",
        ],
        "ui": [
            "uterus", "uterine", "endometrial", "endometrium",
            "carcinoma",
        ],
    },
    {
        "code": "UVM",
        "tcga_name": "Uveal Melanoma",
        "display": "Uveal melanoma",
        "pubmed": [
            "uveal melanoma",
            "ocular melanoma",
            "choroidal melanoma",
        ],
        "ui": [
            "uveal", "ocular", "eye", "choroidal", "melanoma",
        ],
    },
]

# Standardize whitespace and missing values in text.
def clean(value):
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Return dedupe for the current record.
def dedupe(values):
    seen = set()
    out = []

    for value in values:
        value = clean(value)
        key = norm(value)

        if value and key not in seen:
            seen.add(key)
            out.append(value)

    return out

# Raise a descriptive validation error.
def fail(message):
    raise ValueError(message)

# Reject blank or missing terminology values.
def check_no_blank(label, values):
    blanks = [i for i, x in enumerate(values) if not clean(x)]
    if blanks:
        fail(f"Blank {label} value(s) found at indices: {blanks}")

# Build main table from the input records.
def build_main_table():
    rows = []

    for cancer in CANCERS:
        pubmed_terms = dedupe(cancer["pubmed"])

        ui_terms = dedupe(
            [cancer["code"], cancer["tcga_name"], cancer["display"]]
            + cancer["ui"]
            + pubmed_terms
        )

        rows.append({
            "TCGA_Code": cancer["code"],
            "TCGA_Study_Name": cancer["tcga_name"],
            "Cancer": cancer["display"],
            "PubMed_Search_Terms": ";".join(pubmed_terms),
            "UI_Search_Aliases": ";".join(ui_terms),
            "PubMed_Term_Count": len(pubmed_terms),
            "UI_Alias_Count": len(ui_terms),
        })

    return pd.DataFrame(rows)

# Build long pubmed from the input records.
def build_long_pubmed(main):
    rows = []

    for _, row in main.iterrows():
        for term in row["PubMed_Search_Terms"].split(";"):
            rows.append({
                "TCGA_Code": row["TCGA_Code"],
                "TCGA_Study_Name": row["TCGA_Study_Name"],
                "Cancer": row["Cancer"],
                "PubMed_Search_Term": term,
            })

    return pd.DataFrame(rows)

# Build long ui from the input records.
def build_long_ui(main):
    rows = []

    for _, row in main.iterrows():
        for term in row["UI_Search_Aliases"].split(";"):
            rows.append({
                "TCGA_Code": row["TCGA_Code"],
                "Cancer": row["Cancer"],
                "UI_Search_Alias": term,
            })

    return pd.DataFrame(rows)

# Find PubMed terms shared across cancer types.
def find_cross_cancer_pubmed_terms(pubmed):
    owners = defaultdict(set)
    spellings = defaultdict(set)

    for _, row in pubmed.iterrows():
        key = norm(row["PubMed_Search_Term"])
        owners[key].add(row["TCGA_Code"])
        spellings[key].add(row["PubMed_Search_Term"])

    rows = []

    for key, codes in owners.items():
        if len(codes) > 1:
            rows.append({
                "PubMed_Search_Term": ";".join(sorted(spellings[key])),
                "TCGA_Codes": ";".join(sorted(codes)),
            })

    return pd.DataFrame(rows)

# Find display aliases shared by cancer types.
def find_ui_collisions(ui):
    owners = defaultdict(set)

    for _, row in ui.iterrows():
        owners[norm(row["UI_Search_Alias"])].add(row["TCGA_Code"])

    collisions = {
        term: sorted(codes)
        for term, codes in owners.items()
        if len(codes) > 1
    }

    return collisions

# Validate  for the current record.
def validate(main, pubmed, ui):
    errors = []

    codes = set(main["TCGA_Code"])

    if len(main) != EXPECTED_COUNT:
        errors.append(
            f"Expected {EXPECTED_COUNT} cancer rows; found {len(main)}."
        )

    if codes != EXPECTED_TCGA_CODES:
        missing = sorted(EXPECTED_TCGA_CODES - codes)
        extra = sorted(codes - EXPECTED_TCGA_CODES)
        errors.append(
            f"TCGA code set mismatch. Missing={missing}; Extra={extra}"
        )

    if main["TCGA_Code"].duplicated().any():
        errors.append("Duplicate TCGA codes found.")

    if main["Cancer"].map(norm).duplicated().any():
        errors.append("Duplicate Cancer display names found.")

    if (main["PubMed_Term_Count"] <= 0).any():
        errors.append("At least one cancer has no PubMed search terms.")

    if (main["UI_Alias_Count"] <= 0).any():
        errors.append("At least one cancer has no UI search aliases.")

    pubmed_collisions = find_cross_cancer_pubmed_terms(pubmed)

    if not pubmed_collisions.empty:
        errors.append(
            "One or more exact PubMed terms are assigned to multiple "
            "TCGA cancer categories."
        )

    forbidden_pubmed_terms = {
        "cancer", "carcinoma", "adenocarcinoma", "tumor", "tumour",
        "neoplasm", "malignancy", "kidney cancer", "renal cancer",
        "lung cancer", "colorectal cancer", "melanoma", "uterine cancer",
    }

    dangerous = sorted(
        {
            term
            for term in pubmed["PubMed_Search_Term"]
            if norm(term) in forbidden_pubmed_terms
        },
        key=str.casefold,
    )

    if dangerous:
        errors.append(
            "Forbidden broad PubMed term(s) found: "
            + "; ".join(dangerous)
        )

    return errors, pubmed_collisions

# Write report for the current record.
def write_report(main, pubmed, ui, ui_collisions):
    with open(OUT_REPORT, "w", encoding="utf-8") as fh:
        fh.write("GPCR Cancer Pipeline - Step 05 Integrity Report\n")
        fh.write("=" * 64 + "\n\n")

        fh.write("SCOPE\n")
        fh.write("-" * 64 + "\n")
        fh.write("Cancer universe: TCGA PanCancer tumor types\n")
        fh.write(f"Expected cancer types : {EXPECTED_COUNT}\n")
        fh.write(f"Observed cancer types : {len(main)}\n")
        fh.write(f"PubMed search terms   : {len(pubmed)}\n")
        fh.write(f"UI search aliases     : {len(ui)}\n\n")

        fh.write("INTEGRITY CHECKS\n")
        fh.write("-" * 64 + "\n")
        fh.write(
            "Exactly 33 TCGA codes: "
            + ("PASS" if len(main) == EXPECTED_COUNT else "FAIL")
            + "\n"
        )
        fh.write(
            "Exact expected PanCancer code set: "
            + (
                "PASS"
                if set(main["TCGA_Code"]) == EXPECTED_TCGA_CODES
                else "FAIL"
            )
            + "\n"
        )
        fh.write(
            "Unique TCGA codes: "
            + ("PASS" if not main["TCGA_Code"].duplicated().any() else "FAIL")
            + "\n"
        )
        fh.write(
            "Every cancer has PubMed terms: "
            + ("PASS" if (main["PubMed_Term_Count"] > 0).all() else "FAIL")
            + "\n"
        )
        fh.write(
            "Every cancer has UI aliases: "
            + ("PASS" if (main["UI_Alias_Count"] > 0).all() else "FAIL")
            + "\n"
        )
        fh.write(
            "No exact PubMed term shared across cancers: PASS\n"
        )

        fh.write("\nDESIGN NOTE\n")
        fh.write("-" * 64 + "\n")
        fh.write(
            "UI aliases are allowed to overlap across cancer categories. "
            "For example, 'kidney' should return KICH, KIRC, and KIRP; "
            "'lung' should return LUAD and LUSC; 'colorectal' should "
            "return COAD and READ; and 'melanoma' should return SKCM "
            "and UVM. These broad UI aliases are never used as PubMed "
            "disease queries.\n"
        )

        fh.write("\nUI ALIASES SHARED BY MULTIPLE CANCERS\n")
        fh.write("-" * 64 + "\n")

        if not ui_collisions:
            fh.write("(none)\n")
        else:
            for term in sorted(ui_collisions):
                fh.write(
                    f"{term}: {', '.join(ui_collisions[term])}\n"
                )

        fh.write("\nOFFICIAL SOURCES\n")
        fh.write("-" * 64 + "\n")
        fh.write(SOURCE_TCGA + "\n")
        fh.write(SOURCE_CODES + "\n")
        fh.write(SOURCE_PANCAN + "\n")

# Run this stage and write its output files.
def main():
    print("=" * 72)
    print("GPCR Cancer Pipeline - Step 05 Build TCGA Cancer Vocabulary")
    print("=" * 72)

    try:
        main_df = build_main_table()
        pubmed_df = build_long_pubmed(main_df)
        ui_df = build_long_ui(main_df)

        errors, pubmed_collisions = validate(
            main_df,
            pubmed_df,
            ui_df,
        )

        if errors:
            print()
            print("STEP 05 FAILED")
            print("-" * 72)

            for error in errors:
                print(f"ERROR: {error}")

            if not pubmed_collisions.empty:
                print()
                print("Cross-cancer PubMed term collisions:")
                print(pubmed_collisions.to_string(index=False))

            sys.exit(1)

        ui_collisions = find_ui_collisions(ui_df)

        main_df.to_csv(OUT_MAIN, index=False, na_rep="")
        pubmed_df.to_csv(OUT_PUBMED, index=False, na_rep="")
        ui_df.to_csv(OUT_UI, index=False, na_rep="")

        write_report(
            main_df,
            pubmed_df,
            ui_df,
            ui_collisions,
        )

        print()
        print("STEP 05 COMPLETE - ALL INTEGRITY CHECKS PASSED")
        print("=" * 72)
        print(f"TCGA cancer types     : {len(main_df)}")
        print(f"PubMed search terms   : {len(pubmed_df)}")
        print(f"UI search aliases     : {len(ui_df)}")
        print(f"Shared PubMed terms   : 0")
        print(f"Shared UI aliases     : {len(ui_collisions)}")
        print()
        print("Saved:")
        print(f"  {OUT_MAIN}")
        print(f"  {OUT_PUBMED}")
        print(f"  {OUT_UI}")
        print(f"  {OUT_REPORT}")
        print("=" * 72)

    except Exception as exc:
        print()
        print("=" * 72)
        print("STEP 05 FAILED")
        print("=" * 72)
        print(f"ERROR: {exc}")
        sys.exit(1)

if __name__ == "__main__":
    main()
