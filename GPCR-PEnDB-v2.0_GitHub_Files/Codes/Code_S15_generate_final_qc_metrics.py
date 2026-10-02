#!/usr/bin/env python3
# Purpose: Reconcile final tables and generate manuscript summary metrics.
# Input: CORRECTED_FINAL16_DB_*.csv
# Output: CORRECTED_FINAL12_*.csv; CORRECTED_FINAL12_Manuscript_Metrics.txt; CORRECTED_FINAL12_Integrity_Report.txt

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import os
import re
import statistics
import sys
from pathlib import Path

import pandas as pd

PAPERS_FILE = "CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv"
ASSOCIATIONS_FILE = "CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv"
GPCR_PAGES_FILE = "CORRECTED_FINAL16_DB_Page_GPCR_to_Cancers.csv"
CANCER_PAGES_FILE = "CORRECTED_FINAL16_DB_Page_Cancer_to_GPCRs.csv"
CANCER_ALIASES_FILE = "CORRECTED_FINAL16_DB_Cancer_Search_Aliases.csv"
SUMMARY_FILE = "CORRECTED_FINAL16_DB_Summary_Counts.csv"

OUT_OVERALL = "CORRECTED_FINAL12_Overall_Summary.csv"
OUT_CATEGORY = "CORRECTED_FINAL12_Category_Summary.csv"
OUT_GPCR = "CORRECTED_FINAL12_GPCR_Summary.csv"
OUT_CANCER = "CORRECTED_FINAL12_Cancer_Summary.csv"
OUT_GPCR_DISTRIBUTION = "CORRECTED_FINAL12_GPCR_Cancer_Count_Distribution.csv"
OUT_ASSOC_DISTRIBUTION = "CORRECTED_FINAL12_Association_Paper_Count_Distribution.csv"
OUT_YEAR = "CORRECTED_FINAL12_Year_Summary.csv"
OUT_MANUSCRIPT = "CORRECTED_FINAL12_Manuscript_Metrics.txt"
OUT_REPORT = "CORRECTED_FINAL12_Integrity_Report.txt"

EXPECTED_CANCERS = 33
EXPECTED_CATEGORIES = 16

CATEGORY_ORDER = [
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

# Check that required columns are present.
def require_columns(df, filename, required):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(
            f"{filename} is missing required column(s): "
            + ", ".join(missing)
        )

# Extract or validate a publication year.
def extract_year(value):
    match = re.search(r"\b(18|19|20|21)\d{2}\b", clean(value))
    if not match:
        return None
    return int(match.group(0))

# Calculate a mean while handling empty input.
def safe_mean(values):
    values = list(values)
    return float(statistics.mean(values)) if values else 0.0

# Calculate a median while handling empty input.
def safe_median(values):
    values = list(values)
    return float(statistics.median(values)) if values else 0.0

# Format a count for the manuscript summary.
def fmt_num(value, digits=2):
    if isinstance(value, int):
        return f"{value:,}"
    return f"{value:,.{digits}f}"

# Collect distinct category labels from records.
def category_set_from_series(series):
    out = set()
    for value in series:
        out.update(split_terms(value))
    return out

# Select records assigned a category.
def rows_with_category(df, column, category):
    return df[column].fillna("").apply(
        lambda value: category in split_terms(value)
    )

# Load inputs from the input records.
def load_inputs():
    files = [
        PAPERS_FILE,
        ASSOCIATIONS_FILE,
        GPCR_PAGES_FILE,
        CANCER_PAGES_FILE,
        CANCER_ALIASES_FILE,
        SUMMARY_FILE,
    ]

    for filename in files:
        if not os.path.exists(filename):
            raise FileNotFoundError(
                f"Missing required input: {filename}"
            )

    papers = pd.read_csv(
        PAPERS_FILE,
        dtype={"PMID": str},
        low_memory=False,
    ).fillna("")

    associations = pd.read_csv(
        ASSOCIATIONS_FILE,
        dtype=str,
        low_memory=False,
    ).fillna("")

    gpcr_pages = pd.read_csv(
        GPCR_PAGES_FILE,
        dtype=str,
        low_memory=False,
    ).fillna("")

    cancer_pages = pd.read_csv(
        CANCER_PAGES_FILE,
        dtype=str,
        low_memory=False,
    ).fillna("")

    cancer_aliases = pd.read_csv(
        CANCER_ALIASES_FILE,
        dtype=str,
        low_memory=False,
    ).fillna("")

    summary = pd.read_csv(
        SUMMARY_FILE,
        dtype=str,
        low_memory=False,
    ).fillna("")

    require_columns(
        papers,
        PAPERS_FILE,
        [
            "UniProt_ID",
            "Official_Gene",
            "TCGA_Code",
            "Cancer",
            "PMID",
            "Year",
            "Association_Categories",
            "Primary_Association_Category",
        ],
    )

    require_columns(
        associations,
        ASSOCIATIONS_FILE,
        [
            "UniProt_ID",
            "Official_Gene",
            "TCGA_Code",
            "Cancer",
            "Paper_Count",
            "PMIDs",
            "Association_Categories",
        ],
    )

    require_columns(
        gpcr_pages,
        GPCR_PAGES_FILE,
        [
            "UniProt_ID",
            "Official_Gene",
            "Cancer_Count",
            "Paper_Count",
            "Unique_PMID_Count",
        ],
    )

    require_columns(
        cancer_pages,
        CANCER_PAGES_FILE,
        [
            "TCGA_Code",
            "Cancer",
            "GPCR_Count",
            "Paper_Count",
            "Unique_PMID_Count",
        ],
    )

    require_columns(
        cancer_aliases,
        CANCER_ALIASES_FILE,
        [
            "TCGA_Code",
            "Cancer",
            "UI_Search_Alias",
            "Shared_UI_Alias",
        ],
    )

    return (
        papers,
        associations,
        gpcr_pages,
        cancer_pages,
        cancer_aliases,
        summary,
    )

# Build category summary from the input records.
def build_category_summary(papers, associations):
    rows = []

    paper_total = len(papers)
    assoc_total = len(associations)

    for category in CATEGORY_ORDER:
        paper_mask = rows_with_category(
            papers,
            "Association_Categories",
            category,
        )
        assoc_mask = rows_with_category(
            associations,
            "Association_Categories",
            category,
        )

        primary_count = int(
            papers[
                "Primary_Association_Category"
            ].eq(category).sum()
        )

        rows.append({
            "Category": category,
            "Paper_Level_Count": int(paper_mask.sum()),
            "Paper_Level_Percent": (
                100.0 * paper_mask.sum() / paper_total
                if paper_total
                else 0.0
            ),
            "Association_Level_Count": int(assoc_mask.sum()),
            "Association_Level_Percent": (
                100.0 * assoc_mask.sum() / assoc_total
                if assoc_total
                else 0.0
            ),
            "Primary_Category_Paper_Count": primary_count,
            "GPCR_Count": papers.loc[
                paper_mask,
                "UniProt_ID",
            ].nunique(),
            "Cancer_Count": papers.loc[
                paper_mask,
                "TCGA_Code",
            ].nunique(),
        })

    return pd.DataFrame(rows).sort_values(
        [
            "Association_Level_Count",
            "Paper_Level_Count",
            "Category",
        ],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)

# Build gpcr summary from the input records.
def build_gpcr_summary(gpcr_pages):
    out = gpcr_pages.copy()

    for column in [
        "Cancer_Count",
        "Paper_Count",
        "Unique_PMID_Count",
    ]:
        out[column] = pd.to_numeric(
            out[column],
            errors="raise",
        ).astype(int)

    return out.sort_values(
        [
            "Cancer_Count",
            "Paper_Count",
            "Official_Gene",
        ],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)

# Build cancer summary from the input records.
def build_cancer_summary(cancer_pages):
    out = cancer_pages.copy()

    for column in [
        "GPCR_Count",
        "Paper_Count",
        "Unique_PMID_Count",
    ]:
        out[column] = pd.to_numeric(
            out[column],
            errors="raise",
        ).astype(int)

    return out.sort_values(
        [
            "GPCR_Count",
            "Paper_Count",
            "TCGA_Code",
        ],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)

# Build gpcr distribution from the input records.
def build_gpcr_distribution(gpcr_summary):
    dist = (
        gpcr_summary[
            "Cancer_Count"
        ]
        .value_counts()
        .sort_index()
        .rename_axis("Cancer_Count")
        .reset_index(name="GPCR_Count")
    )

    dist["Percent_of_GPCRs"] = (
        100.0 * dist["GPCR_Count"] / len(gpcr_summary)
    )

    return dist

# Build association distribution from the input records.
def build_association_distribution(associations):
    counts = pd.to_numeric(
        associations["Paper_Count"],
        errors="raise",
    ).astype(int)

    dist = (
        counts
        .value_counts()
        .sort_index()
        .rename_axis("Paper_Count")
        .reset_index(name="Association_Count")
    )

    dist["Percent_of_Associations"] = (
        100.0 * dist["Association_Count"] / len(associations)
    )

    return dist

# Build year summary from the input records.
def build_year_summary(papers):
    years = papers["Year"].map(extract_year)
    valid = papers.assign(
        Parsed_Year=years
    ).dropna(subset=["Parsed_Year"]).copy()

    valid["Parsed_Year"] = (
        valid["Parsed_Year"].astype(int)
    )

    annual = (
        valid.groupby("Parsed_Year")
        .agg(
            GPCR_Cancer_Paper_Rows=("PMID", "size"),
            Unique_PMIDs=("PMID", "nunique"),
            GPCRs=("UniProt_ID", "nunique"),
            Cancers=("TCGA_Code", "nunique"),
        )
        .reset_index()
        .rename(columns={"Parsed_Year": "Year"})
        .sort_values("Year")
        .reset_index(drop=True)
    )

    return annual

# Build overall summary from the input records.
def build_overall_summary(
    papers,
    associations,
    gpcr_summary,
    cancer_summary,
    cancer_aliases,
    category_summary,
):
    cancer_counts = gpcr_summary[
        "Cancer_Count"
    ].tolist()

    gpcr_counts_by_cancer = cancer_summary[
        "GPCR_Count"
    ].tolist()

    assoc_paper_counts = pd.to_numeric(
        associations["Paper_Count"],
        errors="raise",
    ).astype(int).tolist()

    parsed_years = [
        year
        for year in papers["Year"].map(extract_year)
        if year is not None
    ]

    other_paper_mask = rows_with_category(
        papers,
        "Association_Categories",
        "Other/Unspecified",
    )

    only_other_paper_mask = (
        papers["Association_Categories"]
        .map(clean)
        .eq("Other/Unspecified")
    )

    other_assoc_mask = rows_with_category(
        associations,
        "Association_Categories",
        "Other/Unspecified",
    )

    only_other_assoc_mask = (
        associations[
            "Association_Categories"
        ]
        .map(clean)
        .eq("Other/Unspecified")
    )

    shared_alias_count = (
        cancer_aliases.loc[
            cancer_aliases[
                "Shared_UI_Alias"
            ].astype(str).str.lower().isin(
                ["true", "1", "yes"]
            ),
            "UI_Search_Alias",
        ]
        .map(norm)
        .nunique()
    )

    top_category = category_summary.iloc[0]

    row = {
        "Human_GPCRs_In_Final_Database":
            int(gpcr_summary["UniProt_ID"].nunique()),
        "Cancer_Types_Represented":
            int(cancer_summary["TCGA_Code"].nunique()),
        "GPCR_Cancer_Associations":
            len(associations),
        "GPCR_Cancer_PMID_Rows":
            len(papers),
        "Unique_PMIDs":
            int(papers["PMID"].nunique()),
        "Association_Categories":
            EXPECTED_CATEGORIES,
        "Earliest_Publication_Year":
            min(parsed_years) if parsed_years else "",
        "Latest_Publication_Year":
            max(parsed_years) if parsed_years else "",
        "Min_Cancers_Per_GPCR":
            min(cancer_counts),
        "Max_Cancers_Per_GPCR":
            max(cancer_counts),
        "Median_Cancers_Per_GPCR":
            safe_median(cancer_counts),
        "Mean_Cancers_Per_GPCR":
            safe_mean(cancer_counts),
        "Min_GPCRs_Per_Cancer":
            min(gpcr_counts_by_cancer),
        "Max_GPCRs_Per_Cancer":
            max(gpcr_counts_by_cancer),
        "Median_GPCRs_Per_Cancer":
            safe_median(gpcr_counts_by_cancer),
        "Mean_GPCRs_Per_Cancer":
            safe_mean(gpcr_counts_by_cancer),
        "Min_Papers_Per_Association":
            min(assoc_paper_counts),
        "Max_Papers_Per_Association":
            max(assoc_paper_counts),
        "Median_Papers_Per_Association":
            safe_median(assoc_paper_counts),
        "Mean_Papers_Per_Association":
            safe_mean(assoc_paper_counts),
        "Other_Paper_Rows":
            int(other_paper_mask.sum()),
        "GPCRs_With_Any_Other_Paper":
            int(
                papers.loc[
                    other_paper_mask,
                    "UniProt_ID",
                ].nunique()
            ),
        "GPCRs_With_Only_Other_Papers":
            int(
                set(
                    papers.loc[
                        only_other_paper_mask,
                        "UniProt_ID",
                    ]
                )
                .difference(
                    set(
                        papers.loc[
                            ~only_other_paper_mask,
                            "UniProt_ID",
                        ]
                    )
                )
                .__len__()
            ),
        "Associations_With_Any_Other":
            int(other_assoc_mask.sum()),
        "Associations_Only_Other":
            int(only_other_assoc_mask.sum()),
        "Website_Cancer_Alias_Rows":
            len(cancer_aliases),
        "Shared_Website_Aliases":
            int(shared_alias_count),
        "Most_Frequent_Association_Category":
            clean(top_category["Category"]),
        "Most_Frequent_Category_Associations":
            int(top_category["Association_Level_Count"]),
        "Most_Frequent_Category_Association_Percent":
            float(top_category["Association_Level_Percent"]),
    }

    return pd.DataFrame([row])

# Return integrity checks for the current record.
def integrity_checks(
    papers,
    associations,
    gpcr_summary,
    cancer_summary,
    cancer_aliases,
    category_summary,
):
    paper_keys = (
        papers[
            ["UniProt_ID", "TCGA_Code", "PMID"]
        ]
        .drop_duplicates()
    )

    association_keys_from_papers = (
        papers[
            ["UniProt_ID", "TCGA_Code"]
        ]
        .drop_duplicates()
    )

    association_keys = (
        associations[
            ["UniProt_ID", "TCGA_Code"]
        ]
        .drop_duplicates()
    )

    paper_count_sum = int(
        pd.to_numeric(
            associations["Paper_Count"],
            errors="raise",
        ).sum()
    )

    gpcr_cancer_count_sum = int(
        gpcr_summary["Cancer_Count"].sum()
    )

    cancer_gpcr_count_sum = int(
        cancer_summary["GPCR_Count"].sum()
    )

    gpcr_paper_count_sum = int(
        gpcr_summary["Paper_Count"].sum()
    )

    cancer_paper_count_sum = int(
        cancer_summary["Paper_Count"].sum()
    )

    observed_paper_categories = category_set_from_series(
        papers["Association_Categories"]
    )

    observed_assoc_categories = category_set_from_series(
        associations["Association_Categories"]
    )

    checks = {
        "No duplicate GPCR-cancer-PMID paper rows":
            len(paper_keys) == len(papers),

        "Association table has unique GPCR-cancer rows":
            len(association_keys) == len(associations),

        "Association keys exactly match paper table":
            set(map(tuple, association_keys_from_papers.values))
            == set(map(tuple, association_keys.values)),

        "Association Paper_Count sums to paper rows":
            paper_count_sum == len(papers),

        "GPCR Cancer_Count sums to association rows":
            gpcr_cancer_count_sum == len(associations),

        "Cancer GPCR_Count sums to association rows":
            cancer_gpcr_count_sum == len(associations),

        "GPCR Paper_Count sums to paper rows":
            gpcr_paper_count_sum == len(papers),

        "Cancer Paper_Count sums to paper rows":
            cancer_paper_count_sum == len(papers),

        "All 33 cancers retained":
            len(cancer_summary) == EXPECTED_CANCERS,

        "GPCR page rows match unique paper-table GPCRs":
            len(gpcr_summary)
            == papers["UniProt_ID"].nunique(),

        "Paper categories exactly match frozen 16":
            observed_paper_categories == set(CATEGORY_ORDER),

        "Association categories exactly match frozen 16":
            observed_assoc_categories == set(CATEGORY_ORDER),

        "Category summary contains all 16 categories":
            set(category_summary["Category"])
            == set(CATEGORY_ORDER),

        "All PMIDs are nonblank":
            papers["PMID"].map(clean).ne("").all(),

        "All GPCR page Cancer_Count values >=1":
            gpcr_summary["Cancer_Count"].ge(1).all(),

        "All cancer page rows have >=1 GPCR":
            cancer_summary["GPCR_Count"].ge(1).all(),

        "Cancer alias table maps only to the 33 cancers":
            set(cancer_aliases["TCGA_Code"])
            .issubset(set(cancer_summary["TCGA_Code"])),
    }

    return checks

# Return manuscript metrics for the current record.
def manuscript_metrics(
    overall,
    category_summary,
    gpcr_summary,
    cancer_summary,
    associations,
):
    row = overall.iloc[0]

    top_categories = category_summary.head(8)

    top_gpcrs = gpcr_summary.head(10)

    top_cancers = cancer_summary.head(10)

    assoc_counts = pd.to_numeric(
        associations["Paper_Count"],
        errors="raise",
    ).astype(int)

    one_paper = int((assoc_counts == 1).sum())
    multi_paper = int((assoc_counts > 1).sum())

    lines = [
        "GPCR Cancer Pipeline - Final Manuscript Metrics",
        "=" * 78,
        "",
        "CORE COUNTS",
        "-" * 78,
        f"Human GPCRs searched upstream                      : 832",
        f"GPCRs in final validated/categorized database      : "
        f"{int(row['Human_GPCRs_In_Final_Database'])}",
        f"TCGA cancer types represented                      : "
        f"{int(row['Cancer_Types_Represented'])}/33",
        f"Unique GPCR-cancer associations                    : "
        f"{int(row['GPCR_Cancer_Associations'])}",
        f"GPCR-cancer-PMID rows                              : "
        f"{int(row['GPCR_Cancer_PMID_Rows'])}",
        f"Unique PMIDs                                       : "
        f"{int(row['Unique_PMIDs'])}",
        f"Association categories                             : "
        f"{int(row['Association_Categories'])}",
        f"Publication years                                  : "
        f"{row['Earliest_Publication_Year']}-"
        f"{row['Latest_Publication_Year']}",
        "",
        "GPCR DISTRIBUTION",
        "-" * 78,
        f"Cancers per GPCR                                   : "
        f"range {int(row['Min_Cancers_Per_GPCR'])}-"
        f"{int(row['Max_Cancers_Per_GPCR'])}; "
        f"median {row['Median_Cancers_Per_GPCR']:.1f}; "
        f"mean {row['Mean_Cancers_Per_GPCR']:.2f}",
        "",
        "CANCER DISTRIBUTION",
        "-" * 78,
        f"GPCRs per cancer                                   : "
        f"range {int(row['Min_GPCRs_Per_Cancer'])}-"
        f"{int(row['Max_GPCRs_Per_Cancer'])}; "
        f"median {row['Median_GPCRs_Per_Cancer']:.1f}; "
        f"mean {row['Mean_GPCRs_Per_Cancer']:.2f}",
        "",
        "EVIDENCE DENSITY",
        "-" * 78,
        f"Papers per GPCR-cancer association                 : "
        f"range {int(row['Min_Papers_Per_Association'])}-"
        f"{int(row['Max_Papers_Per_Association'])}; "
        f"median {row['Median_Papers_Per_Association']:.1f}; "
        f"mean {row['Mean_Papers_Per_Association']:.2f}",
        f"Associations represented by one paper              : "
        f"{one_paper}",
        f"Associations represented by >1 paper               : "
        f"{multi_paper}",
        "",
        "OTHER / UNSPECIFIED",
        "-" * 78,
        f"Paper rows assigned Other/Unspecified              : "
        f"{int(row['Other_Paper_Rows'])}",
        f"GPCRs with >=1 Other/Unspecified paper             : "
        f"{int(row['GPCRs_With_Any_Other_Paper'])}",
        f"GPCRs represented exclusively by Other papers      : "
        f"{int(row['GPCRs_With_Only_Other_Papers'])}",
        f"GPCR-cancer associations containing Other          : "
        f"{int(row['Associations_With_Any_Other'])}",
        f"GPCR-cancer associations exclusively Other         : "
        f"{int(row['Associations_Only_Other'])}",
        "",
        "CATEGORY COUNTS AT GPCR-CANCER ASSOCIATION LEVEL",
        "-" * 78,
        (
            "NOTE: multi-label categories are permitted; percentages therefore "
            "do not sum to 100%."
        ),
    ]

    for _, cat in top_categories.iterrows():
        lines.append(
            f"{clean(cat['Category']):36s} : "
            f"{int(cat['Association_Level_Count']):5d} "
            f"({cat['Association_Level_Percent']:.1f}%)"
        )

    lines += [
        "",
        "TOP GPCRs BY NUMBER OF CANCER TYPES",
        "-" * 78,
    ]

    for _, item in top_gpcrs.iterrows():
        lines.append(
            f"{clean(item['Official_Gene']):16s} : "
            f"{int(item['Cancer_Count']):2d} cancers; "
            f"{int(item['Paper_Count']):4d} GPCR-cancer-PMID rows"
        )

    lines += [
        "",
        "TOP CANCERS BY NUMBER OF GPCRs",
        "-" * 78,
    ]

    for _, item in top_cancers.iterrows():
        lines.append(
            f"{clean(item['TCGA_Code']):5s} "
            f"{clean(item['Cancer']):28s} : "
            f"{int(item['GPCR_Count']):3d} GPCRs; "
            f"{int(item['Paper_Count']):5d} GPCR-cancer-PMID rows"
        )

    return "\n".join(lines)

# Run this stage and write its output files.
def main():
    print("=" * 78)
    print("GPCR Cancer Pipeline - Code S15")
    print("Final QC, Summary Statistics, and Manuscript Metrics")
    print("=" * 78)

    try:
        (
            papers,
            associations,
            gpcr_pages,
            cancer_pages,
            cancer_aliases,
            summary_input,
        ) = load_inputs()

        gpcr_summary = build_gpcr_summary(
            gpcr_pages
        )
        cancer_summary = build_cancer_summary(
            cancer_pages
        )
        category_summary = build_category_summary(
            papers,
            associations,
        )
        gpcr_distribution = build_gpcr_distribution(
            gpcr_summary
        )
        assoc_distribution = (
            build_association_distribution(
                associations
            )
        )
        year_summary = build_year_summary(
            papers
        )
        overall = build_overall_summary(
            papers,
            associations,
            gpcr_summary,
            cancer_summary,
            cancer_aliases,
            category_summary,
        )

        checks = integrity_checks(
            papers,
            associations,
            gpcr_summary,
            cancer_summary,
            cancer_aliases,
            category_summary,
        )

    except Exception as exc:
        print(f"CODE S15 FAILED:\n{exc}")
        sys.exit(1)

    failed = [
        label
        for label, passed in checks.items()
        if not passed
    ]

    overall.to_csv(
        OUT_OVERALL,
        index=False,
        na_rep="",
    )
    category_summary.to_csv(
        OUT_CATEGORY,
        index=False,
        na_rep="",
    )
    gpcr_summary.to_csv(
        OUT_GPCR,
        index=False,
        na_rep="",
    )
    cancer_summary.to_csv(
        OUT_CANCER,
        index=False,
        na_rep="",
    )
    gpcr_distribution.to_csv(
        OUT_GPCR_DISTRIBUTION,
        index=False,
        na_rep="",
    )
    assoc_distribution.to_csv(
        OUT_ASSOC_DISTRIBUTION,
        index=False,
        na_rep="",
    )
    year_summary.to_csv(
        OUT_YEAR,
        index=False,
        na_rep="",
    )

    manuscript = manuscript_metrics(
        overall,
        category_summary,
        gpcr_summary,
        cancer_summary,
        associations,
    )
    Path(OUT_MANUSCRIPT).write_text(
        manuscript + "\n",
        encoding="utf-8",
    )

    row = overall.iloc[0]

    report_lines = [
        "GPCR Cancer Pipeline - Code S15 Integrity Report",
        "=" * 78,
        "",
        "FINAL DATABASE COUNTS",
        "-" * 78,
        f"GPCRs                              : "
        f"{int(row['Human_GPCRs_In_Final_Database'])}",
        f"Cancer types                       : "
        f"{int(row['Cancer_Types_Represented'])}",
        f"GPCR-cancer associations           : "
        f"{int(row['GPCR_Cancer_Associations'])}",
        f"GPCR-cancer-PMID rows              : "
        f"{int(row['GPCR_Cancer_PMID_Rows'])}",
        f"Unique PMIDs                       : "
        f"{int(row['Unique_PMIDs'])}",
        f"Categories                         : "
        f"{int(row['Association_Categories'])}",
        "",
        "INTEGRITY CHECKS",
        "-" * 78,
    ]

    for label, passed in checks.items():
        report_lines.append(
            f"{label}: {'PASS' if passed else 'FAIL'}"
        )

    report_lines += [
        "",
        "METHOD NOTE",
        "-" * 78,
        (
            "Code S15 performs no scientific filtering or reclassification. "
            "It reconciles the completed final database tables and calculates "
            "descriptive statistics for reporting. Category percentages are "
            "computed against the total number of unique GPCR-cancer "
            "associations; categories are multi-label and therefore their "
            "percentages are not expected to sum to 100%."
        ),
        "",
        f"Integrity failures                 : {len(failed)}",
    ]

    if failed:
        report_lines += [
            "",
            "FAILED CHECKS",
            "-" * 78,
        ]
        report_lines.extend(
            f"- {item}"
            for item in failed
        )
    else:
        report_lines += [
            "",
            "ALL CODE S15 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(report_lines) + "\n",
        encoding="utf-8",
    )

    print()
    print("=" * 78)
    print(
        "CODE S15 COMPLETE"
        if not failed
        else "CODE S15 COMPLETE WITH INTEGRITY ISSUES"
    )
    print("=" * 78)
    print(
        f"GPCRs             : "
        f"{int(row['Human_GPCRs_In_Final_Database']):,}"
    )
    print(
        f"Cancer types      : "
        f"{int(row['Cancer_Types_Represented']):,}"
    )
    print(
        f"Associations      : "
        f"{int(row['GPCR_Cancer_Associations']):,}"
    )
    print(
        f"Paper rows        : "
        f"{int(row['GPCR_Cancer_PMID_Rows']):,}"
    )
    print(
        f"Unique PMIDs      : "
        f"{int(row['Unique_PMIDs']):,}"
    )
    print(
        f"Categories        : "
        f"{int(row['Association_Categories']):,}"
    )
    print(
        f"Integrity errors  : "
        f"{len(failed)}"
    )
    print()
    print("Saved:")
    print(f"  {OUT_OVERALL}")
    print(f"  {OUT_CATEGORY}")
    print(f"  {OUT_GPCR}")
    print(f"  {OUT_CANCER}")
    print(f"  {OUT_GPCR_DISTRIBUTION}")
    print(f"  {OUT_ASSOC_DISTRIBUTION}")
    print(f"  {OUT_YEAR}")
    print(f"  {OUT_MANUSCRIPT}")
    print(f"  {OUT_REPORT}")

    if failed:
        sys.exit(2)

if __name__ == "__main__":
    main()
