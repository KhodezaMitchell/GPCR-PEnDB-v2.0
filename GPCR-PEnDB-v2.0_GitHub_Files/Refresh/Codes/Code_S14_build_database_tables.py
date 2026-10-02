#!/usr/bin/env python3
# Purpose: Build the final GPCR-centered, cancer-centered, and website tables.
# Input: CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv; TCGA cancer terminology files
# Output: CORRECTED_FINAL16_DB_*.csv

import os
import re
import sys
from pathlib import Path

import pandas as pd

CLASSIFIED_FILE = "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv"
VOCABULARY_FILE = "TCGA_Cancer_Vocabulary.csv"
UI_ALIASES_FILE = "TCGA_Cancer_UI_Aliases_Long.csv"

OUT_PAPERS = "CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv"
OUT_ASSOCIATIONS = "CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv"
OUT_GPCR_PAGES = "CORRECTED_FINAL16_DB_Page_GPCR_to_Cancers.csv"
OUT_CANCER_PAGES = "CORRECTED_FINAL16_DB_Page_Cancer_to_GPCRs.csv"
OUT_CANCER_ALIASES = "CORRECTED_FINAL16_DB_Cancer_Search_Aliases.csv"
OUT_SUMMARY = "CORRECTED_FINAL16_DB_Summary_Counts.csv"
OUT_REPORT = "CORRECTED_FINAL16_GPCR_Cancer_Database_Integrity_Report.txt"

EXPECTED_CANCERS = 33

ASSOCIATION_CATEGORIES = [
    "Proliferation/Apoptosis",
    "Migration/Invasion/Metastasis",
    "Tumor Growth/Progression",
    "Biomarker/Diagnosis",
    "Prognosis",
    "Survival",
    "Therapeutic Target",
    "Expression",
    "Signaling",
    "Genetic Association",
    "Cancer Risk/Susceptibility",
    "Drug Response/Resistance",
    "Immune Microenvironment",
    "Imaging/Theranostics",
    "Computational/Network Evidence",
    "Other/Unspecified",
]

# Standardize whitespace and missing values in text.
def clean(value):
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Split a semicolon-delimited term field into names.
def split_terms(value):
    return [
        clean(x)
        for x in clean(value).split(";")
        if clean(x)
    ]

# Return unique terms while preserving the chosen order.
def unique_terms(values):
    seen = set()
    out = []

    for value in values:
        for item in split_terms(value):
            key = norm(item)
            if key and key not in seen:
                seen.add(key)
                out.append(item)

    return out

# Join unique for the current record.
def join_unique(values, sep="; "):
    return sep.join(unique_terms(values))

# Count distinct normalized terms.
def unique_term_count(values):
    return len({
        norm(x)
        for value in values
        for x in split_terms(value)
        if norm(x)
    })

# Join text for the current record.
def join_text(values, sep=" || "):
    seen = set()
    out = []

    for value in values:
        text = clean(value)
        if text and text not in seen:
            seen.add(text)
            out.append(text)

    return sep.join(out)

# Return a numeric key for sorting PMIDs.
def pmid_sort_key(value):
    text = clean(value)
    if text.isdigit():
        return (0, int(text))
    return (1, text)

# Join unique PMIDs in numeric order.
def join_pmids(values):
    pmids = {
        x
        for value in values
        for x in split_terms(value)
        if x
    }
    return ";".join(sorted(pmids, key=pmid_sort_key))

# Count distinct PMIDs in a field.
def unique_pmid_count(values):
    return len({
        x
        for value in values
        for x in split_terms(value)
        if x
    })

# Join distinct publication years in order.
def join_years(values):
    years = {
        clean(x)
        for x in values
        if clean(x)
    }
    return ";".join(sorted(years))

# Check that required columns are present.
def require_columns(df, filename, required):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(
            f"{filename} is missing required column(s): "
            + ", ".join(missing)
        )

# Return fail if for the current record.
def fail_if(condition, message):
    if condition:
        raise ValueError(message)

# Load inputs from the input records.
def load_inputs():
    for filename in [
        CLASSIFIED_FILE,
        VOCABULARY_FILE,
        UI_ALIASES_FILE,
    ]:
        if not os.path.exists(filename):
            raise FileNotFoundError(
                f"Missing required input: {filename}"
            )

    classified = pd.read_csv(
        CLASSIFIED_FILE,
        dtype={"PMID": str},
        low_memory=False,
    ).fillna("")

    vocabulary = pd.read_csv(
        VOCABULARY_FILE,
        dtype=str,
    ).fillna("")

    ui_aliases = pd.read_csv(
        UI_ALIASES_FILE,
        dtype=str,
    ).fillna("")

    require_columns(
        classified,
        CLASSIFIED_FILE,
        [
            "UniProt_ID",
            "Gene",
            "Official_Gene",
            "Original_Gene",
            "Protein_Name",
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PMID",
            "Title",
            "Abstract",
            "Journal",
            "Year",
            "Authors",
            "Search_ID",
            "Validated_Search_Terms",
            "Matched_Aliases",
            "Matched_Alias_Count",
            "Supporting_Sentence_Count",
            "Supporting_Sentences",
            "Association_Categories",
            "Association_Category_Count",
            "Primary_Association_Category",
            "Matched_Keywords",
            "Association_Status",
            "Review_Reason",
            "PubMed_Query",
        ],
    )

    require_columns(
        vocabulary,
        VOCABULARY_FILE,
        [
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PubMed_Search_Terms",
            "UI_Search_Aliases",
        ],
    )

    require_columns(
        ui_aliases,
        UI_ALIASES_FILE,
        [
            "TCGA_Code",
            "Cancer",
            "UI_Search_Alias",
        ],
    )

    fail_if(
        len(vocabulary) != EXPECTED_CANCERS,
        f"Expected {EXPECTED_CANCERS} TCGA cancer rows; "
        f"found {len(vocabulary)}.",
    )

    fail_if(
        vocabulary["TCGA_Code"].map(norm).duplicated().any(),
        "Duplicate TCGA codes in cancer vocabulary.",
    )

    fail_if(
        classified.duplicated(
            ["UniProt_ID", "TCGA_Code", "PMID"]
        ).any(),
        "Duplicate GPCR-cancer-PMID rows in classified input.",
    )

    statuses = set(
        classified["Association_Status"].map(clean)
    )
    fail_if(
        not statuses or not statuses.issubset({"Retain", "Review"}),
        f"Unexpected Association_Status values: {sorted(statuses)}",
    )

    vocab_codes = set(
        vocabulary["TCGA_Code"].map(clean)
    )
    classified_codes = set(
        classified["TCGA_Code"].map(clean)
    )

    fail_if(
        not classified_codes.issubset(vocab_codes),
        "Classified input contains TCGA codes absent from vocabulary.",
    )

    alias_codes = set(
        ui_aliases["TCGA_Code"].map(clean)
    )
    fail_if(
        not alias_codes.issubset(vocab_codes),
        "UI alias file contains TCGA codes absent from vocabulary.",
    )

    duplicate_alias_pairs = pd.DataFrame({
        "TCGA_Code": ui_aliases["TCGA_Code"].map(norm),
        "Alias": ui_aliases["UI_Search_Alias"].map(norm),
    }).duplicated()

    fail_if(
        duplicate_alias_pairs.any(),
        "Duplicate case-insensitive TCGA-code/UI-alias rows found.",
    )

    return classified, vocabulary, ui_aliases

# Validate retained for the current record.
def validate_retained(retained):
    fail_if(
        retained.empty,
        "No retained classification rows."
    )

    fail_if(
        retained["Matched_Alias_Count"].astype(int).le(0).any(),
        "Retained row without matched validated GPCR alias.",
    )

    fail_if(
        retained["Supporting_Sentence_Count"].astype(int).le(0).any(),
        "Retained row without supporting sentence.",
    )

    fail_if(
        retained["Association_Category_Count"].astype(int).le(0).any(),
        "Retained row without association category.",
    )

    observed = {
        item
        for value in retained["Association_Categories"]
        for item in split_terms(value)
    }

    fail_if(
        observed != set(ASSOCIATION_CATEGORIES),
        "Observed retained category set does not exactly match "
        "the frozen 16-category vocabulary, including Other/Unspecified."
    )

# Build papers from the input records.
def build_papers(retained):
    columns = [
        "UniProt_ID",
        "Gene",
        "Official_Gene",
        "Original_Gene",
        "Protein_Name",
        "TCGA_Code",
        "TCGA_Study_Name",
        "Cancer",
        "PMID",
        "Title",
        "Abstract",
        "Journal",
        "Year",
        "Authors",
        "Validated_Search_Terms",
        "Matched_Aliases",
        "Matched_Alias_Count",
        "Supporting_Sentence_Count",
        "Supporting_Sentences",
        "Association_Categories",
        "Association_Category_Count",
        "Primary_Association_Category",
        "Matched_Keywords",
        "PubMed_Query",
        "Search_ID",
    ]

    return (
        retained[columns]
        .copy()
        .sort_values(
            ["Official_Gene", "TCGA_Code", "Year", "PMID"],
            kind="stable",
        )
        .reset_index(drop=True)
    )

# Build associations from the input records.
def build_associations(papers):
    return (
        papers.groupby(
            [
                "UniProt_ID",
                "Gene",
                "Official_Gene",
                "Original_Gene",
                "Protein_Name",
                "TCGA_Code",
                "TCGA_Study_Name",
                "Cancer",
            ],
            dropna=False,
        )
        .agg(
            Paper_Count=("PMID", "nunique"),
            PMIDs=("PMID", join_pmids),
            Years=("Year", join_years),
            Journals=("Journal", join_unique),
            Validated_Search_Terms=(
                "Validated_Search_Terms",
                join_unique,
            ),
            Matched_Aliases=(
                "Matched_Aliases",
                join_unique,
            ),
            Matched_Alias_Count=(
                "Matched_Aliases",
                unique_term_count,
            ),
            Association_Categories=(
                "Association_Categories",
                join_unique,
            ),
            Primary_Association_Categories=(
                "Primary_Association_Category",
                join_unique,
            ),
            Matched_Keywords=(
                "Matched_Keywords",
                join_unique,
            ),
            Supporting_Sentences=(
                "Supporting_Sentences",
                join_text,
            ),
        )
        .reset_index()
        .sort_values(
            ["Paper_Count", "Official_Gene", "TCGA_Code"],
            ascending=[False, True, True],
            kind="stable",
        )
        .reset_index(drop=True)
    )

# Build gpcr pages from the input records.
def build_gpcr_pages(associations):
    return (
        associations.groupby(
            [
                "UniProt_ID",
                "Gene",
                "Official_Gene",
                "Original_Gene",
                "Protein_Name",
            ],
            dropna=False,
        )
        .agg(
            Cancer_Count=("TCGA_Code", "nunique"),
            Paper_Count=("Paper_Count", "sum"),
            Unique_PMID_Count=("PMIDs", unique_pmid_count),
            TCGA_Codes=("TCGA_Code", join_unique),
            Cancers=("Cancer", join_unique),
            PMIDs=("PMIDs", join_pmids),
            Validated_Search_Terms=(
                "Validated_Search_Terms",
                join_unique,
            ),
            Matched_Aliases=(
                "Matched_Aliases",
                join_unique,
            ),
            Association_Categories=(
                "Association_Categories",
                join_unique,
            ),
            Primary_Association_Categories=(
                "Primary_Association_Categories",
                join_unique,
            ),
        )
        .reset_index()
        .sort_values(
            ["Cancer_Count", "Paper_Count", "Official_Gene"],
            ascending=[False, False, True],
            kind="stable",
        )
        .reset_index(drop=True)
    )

# Build cancer pages from the input records.
def build_cancer_pages(associations, vocabulary):
    agg = (
        associations.groupby(
            ["TCGA_Code", "TCGA_Study_Name", "Cancer"],
            dropna=False,
        )
        .agg(
            GPCR_Count=("UniProt_ID", "nunique"),
            Paper_Count=("Paper_Count", "sum"),
            Unique_PMID_Count=("PMIDs", unique_pmid_count),
            GPCRs=("Official_Gene", join_unique),
            UniProt_IDs=("UniProt_ID", join_unique),
            PMIDs=("PMIDs", join_pmids),
            Association_Categories=(
                "Association_Categories",
                join_unique,
            ),
            Primary_Association_Categories=(
                "Primary_Association_Categories",
                join_unique,
            ),
        )
        .reset_index()
    )

    base = vocabulary[
        ["TCGA_Code", "TCGA_Study_Name", "Cancer"]
    ].copy()

    pages = base.merge(
        agg,
        on=["TCGA_Code", "TCGA_Study_Name", "Cancer"],
        how="left",
        validate="one_to_one",
    )

    numeric_cols = [
        "GPCR_Count",
        "Paper_Count",
        "Unique_PMID_Count",
    ]

    for col in numeric_cols:
        pages[col] = (
            pd.to_numeric(pages[col], errors="coerce")
            .fillna(0)
            .astype(int)
        )

    for col in [
        "GPCRs",
        "UniProt_IDs",
        "PMIDs",
        "Association_Categories",
        "Primary_Association_Categories",
    ]:
        pages[col] = pages[col].fillna("")

    return (
        pages.sort_values(
            ["GPCR_Count", "Paper_Count", "TCGA_Code"],
            ascending=[False, False, True],
            kind="stable",
        )
        .reset_index(drop=True)
    )

# Build cancer aliases from the input records.
def build_cancer_aliases(vocabulary, ui_aliases):
    studies = vocabulary[
        ["TCGA_Code", "TCGA_Study_Name", "Cancer"]
    ].copy()

    aliases = ui_aliases.merge(
        studies,
        on=["TCGA_Code", "Cancer"],
        how="left",
        validate="many_to_one",
    )

    fail_if(
        aliases["TCGA_Study_Name"].map(clean).eq("").any(),
        "UI alias row failed to map to TCGA vocabulary.",
    )

    aliases["UI_Search_Alias_Normalized"] = (
        aliases["UI_Search_Alias"].map(norm)
    )

    owners = (
        aliases.groupby("UI_Search_Alias_Normalized")
        .agg(
            Matched_Cancer_Count=(
                "TCGA_Code",
                "nunique",
            ),
            Matched_TCGA_Codes=(
                "TCGA_Code",
                join_unique,
            ),
            Matched_Cancers=(
                "Cancer",
                join_unique,
            ),
        )
        .reset_index()
    )

    aliases = aliases.merge(
        owners,
        on="UI_Search_Alias_Normalized",
        how="left",
        validate="many_to_one",
    )

    aliases["Shared_UI_Alias"] = (
        aliases["Matched_Cancer_Count"].astype(int) > 1
    )

    return aliases[
        [
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "UI_Search_Alias",
            "UI_Search_Alias_Normalized",
            "Matched_Cancer_Count",
            "Shared_UI_Alias",
            "Matched_TCGA_Codes",
            "Matched_Cancers",
        ]
    ].sort_values(
        ["UI_Search_Alias_Normalized", "TCGA_Code"],
        kind="stable",
    ).reset_index(drop=True)

# Build summary from the input records.
def build_summary(
    classified,
    papers,
    associations,
    gpcr_pages,
    cancer_pages,
    cancer_aliases,
):
    retained = int(
        classified["Association_Status"].eq("Retain").sum()
    )
    review = int(
        classified["Association_Status"].eq("Review").sum()
    )

    shared_aliases = int(
        cancer_aliases.loc[
            cancer_aliases["Shared_UI_Alias"],
            "UI_Search_Alias_Normalized",
        ].nunique()
    )

    return pd.DataFrame([{
        "Classified_Rows": len(classified),
        "Retained_Rows": retained,
        "Review_Rows": review,
        "Paper_Level_Rows": len(papers),
        "Unique_GPCRs": len(gpcr_pages),
        "Unique_Cancers": len(cancer_pages),
        "Unique_GPCR_Cancer_Associations": len(associations),
        "Unique_PMIDs": papers["PMID"].nunique(),
        "Association_Categories": len(ASSOCIATION_CATEGORIES),
        "UI_Alias_Rows": len(cancer_aliases),
        "Shared_UI_Aliases": shared_aliases,
    }])

# Write report for the current record.
def write_report(
    classified,
    retained,
    papers,
    associations,
    gpcr_pages,
    cancer_pages,
    cancer_aliases,
):
    review = classified[
        classified["Association_Status"].eq("Review")
    ]

    association_count_expected = (
        retained[
            ["UniProt_ID", "TCGA_Code"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    gpcr_count_expected = retained[
        "UniProt_ID"
    ].nunique()

    unique_pmids_expected = retained[
        "PMID"
    ].nunique()

    checks = {
        "Retained paper rows preserved": (
            len(papers) == len(retained)
        ),
        "No duplicate paper associations": (
            not papers.duplicated(
                ["UniProt_ID", "TCGA_Code", "PMID"]
            ).any()
        ),
        "Association count matches retained input": (
            len(associations) == association_count_expected
        ),
        "GPCR page count matches retained GPCRs": (
            len(gpcr_pages) == gpcr_count_expected
        ),
        "All 33 cancer page rows present": (
            len(cancer_pages) == EXPECTED_CANCERS
        ),
        "Unique retained PMID count preserved": (
            papers["PMID"].nunique() == unique_pmids_expected
        ),
        "All 16 categories preserved": (
            {
                x
                for value in papers["Association_Categories"]
                for x in split_terms(value)
            }
            == set(ASSOCIATION_CATEGORIES)
        ),
        "Every UI alias maps to known cancer": (
            cancer_aliases["TCGA_Code"].map(clean).ne("").all()
        ),
    }

    failed = [
        label
        for label, passed in checks.items()
        if not passed
    ]

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 11 Integrity Report",
        "=" * 78,
        "",
        "INPUT",
        "-" * 78,
        f"Classified rows                    : {len(classified)}",
        f"Retained rows                      : {len(retained)}",
        f"Review rows                        : {len(review)}",
        "",
        "FINAL DATABASE TABLES",
        "-" * 78,
        f"Paper-level rows                   : {len(papers)}",
        f"Unique retained PMIDs              : {papers['PMID'].nunique()}",
        f"GPCR-cancer associations           : {len(associations)}",
        f"GPCR page rows                     : {len(gpcr_pages)}",
        f"Cancer page rows                   : {len(cancer_pages)}",
        f"Association categories             : {len(ASSOCIATION_CATEGORIES)}",
        f"Website cancer alias rows          : {len(cancer_aliases)}",
        f"Shared website aliases             : "
        f"{cancer_aliases.loc[cancer_aliases['Shared_UI_Alias'], 'UI_Search_Alias_Normalized'].nunique()}",
        "",
        "INTEGRITY CHECKS",
        "-" * 78,
    ]

    for label, passed in checks.items():
        lines.append(
            f"{label}: {'PASS' if passed else 'FAIL'}"
        )

    lines += [
        "",
        "METHOD NOTE",
        "-" * 78,
        (
            "Step 11 performs no scientific reclassification. "
            "Only Step 10 rows labeled Association_Status=Retain are "
            "aggregated into database and website tables. The complete "
            "15-category assignments, supporting sentences, matched "
            "validated terms, and PubMed identifiers are preserved."
        ),
        "",
        f"Integrity failures                 : {len(failed)}",
    ]

    if failed:
        lines += [
            "",
            "FAILED CHECKS",
            "-" * 78,
        ]
        lines.extend(f"- {x}" for x in failed)
    else:
        lines += [
            "",
            "ALL STEP 11 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return failed

# Run this stage and write its output files.
def main():
    print("=" * 78)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 11")
    print("Database and Website Table Builder")
    print("=" * 78)

    try:
        classified, vocabulary, ui_aliases = load_inputs()

        retained = classified[
            classified["Association_Status"].eq("Retain")
        ].copy()

        validate_retained(retained)

        papers = build_papers(retained)
        associations = build_associations(papers)
        gpcr_pages = build_gpcr_pages(associations)
        cancer_pages = build_cancer_pages(
            associations,
            vocabulary,
        )
        cancer_aliases = build_cancer_aliases(
            vocabulary,
            ui_aliases,
        )
        summary = build_summary(
            classified,
            papers,
            associations,
            gpcr_pages,
            cancer_pages,
            cancer_aliases,
        )

        failed = write_report(
            classified,
            retained,
            papers,
            associations,
            gpcr_pages,
            cancer_pages,
            cancer_aliases,
        )

    except Exception as exc:
        print(f"STEP 11 VALIDATION FAILED:\n{exc}")
        sys.exit(1)

    papers.to_csv(
        OUT_PAPERS,
        index=False,
        na_rep="",
    )
    associations.to_csv(
        OUT_ASSOCIATIONS,
        index=False,
        na_rep="",
    )
    gpcr_pages.to_csv(
        OUT_GPCR_PAGES,
        index=False,
        na_rep="",
    )
    cancer_pages.to_csv(
        OUT_CANCER_PAGES,
        index=False,
        na_rep="",
    )
    cancer_aliases.to_csv(
        OUT_CANCER_ALIASES,
        index=False,
        na_rep="",
    )
    summary.to_csv(
        OUT_SUMMARY,
        index=False,
        na_rep="",
    )

    print()
    print("=" * 78)
    print(
        "STEP 11 COMPLETE"
        if not failed
        else "STEP 11 COMPLETE WITH INTEGRITY ISSUES"
    )
    print("=" * 78)
    print(f"Classified rows        : {len(classified):,}")
    print(f"Retained paper rows    : {len(papers):,}")
    print(f"GPCR-cancer associations: {len(associations):,}")
    print(f"GPCR page rows         : {len(gpcr_pages):,}")
    print(f"Cancer page rows       : {len(cancer_pages):,}")
    print(f"Unique retained PMIDs  : {papers['PMID'].nunique():,}")
    print(f"Association categories : {len(ASSOCIATION_CATEGORIES)}")
    print(f"Website alias rows     : {len(cancer_aliases):,}")
    print(f"Integrity errors       : {len(failed)}")
    print()
    print("Saved:")
    print(f"  {OUT_PAPERS}")
    print(f"  {OUT_ASSOCIATIONS}")
    print(f"  {OUT_GPCR_PAGES}")
    print(f"  {OUT_CANCER_PAGES}")
    print(f"  {OUT_CANCER_ALIASES}")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_REPORT}")
    print("=" * 78)

    if failed:
        sys.exit(2)

if __name__ == "__main__":
    main()
