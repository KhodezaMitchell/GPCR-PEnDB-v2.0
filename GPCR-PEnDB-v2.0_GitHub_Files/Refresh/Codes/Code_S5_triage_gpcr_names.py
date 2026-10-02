#!/usr/bin/env python3
# Purpose: Apply reproducible triage rules to candidate GPCR names.
# Input: Alias_Audit_Recommendations.csv; Alias_Context_Review.csv
# Output: Alias_Systematic_Triage.csv; Alias_Manual_Review_Queue.csv

import argparse
import os
import re
from pathlib import Path

import pandas as pd

DEFAULT_AUDIT = "Alias_Audit_Recommendations.csv"
DEFAULT_CONTEXT = "Alias_Context_Review.csv"

OUT_TRIAGE = "Alias_Systematic_Triage.csv"
OUT_MANUAL = "Alias_Manual_Review_Queue.csv"
OUT_SUMMARY = "Alias_Systematic_Triage_Summary.txt"
OUT_INTEGRITY = "Alias_Systematic_Triage_Integrity_Report.txt"

ALLOWED_TRIAGE = {
    "AUTO_KEEP",
    "CONTEXT_REQUIRED",
    "NOT_VALIDATED",
    "MANUAL_REVIEW",
}

# Standardize whitespace and missing values in text.
def clean(value):
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Interpret a text value as a Boolean flag.
def yes(value):
    return clean(value).upper() in {"TRUE", "YES", "1"}

# Convert a value to an integer when possible.
def as_int(value):
    try:
        return int(float(clean(value) or 0))
    except Exception:
        return 0

# Split types for the current record.
def split_types(value):
    return {clean(x) for x in clean(value).split(";") if clean(x)}

# Build normalized keys for candidate rows.
def key_series(df):
    return (
        df["UniProt_ID"].map(clean)
        + "\t"
        + df["Alias"].map(norm)
    )

# Check that required columns are present.
def require_columns(df, required, filename):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(
            f"ERROR: {filename} is missing required column(s): "
            + ", ".join(missing)
        )

AUDIT_REQUIRED = {
    "UniProt_ID",
    "Original_Gene",
    "Official_Gene",
    "Protein_Name",
    "HGNC_Approved_Symbol",
    "Alias",
    "Source_Provenance",
    "Term_Types",
    "Shared_By_UniProt_HGNC",
    "Official_Symbol_Disagreement",
    "Is_Official_Gene",
    "Is_Original_Gene",
    "Is_Protein_Name",
    "Is_Short_Term",
    "Is_Acronym_Risk",
    "Is_Generic_Numbered_Receptor",
    "Cross_GPCR_Collision",
    "Collision_Type",
    "Collision_Genes",
    "Flag_Z_Reasons",
    "Retrieved_PMIDs",
    "PMID_List_Truncated",
    "Independent_Reference_PMIDs",
    "Reference_PMID_List_Truncated",
    "Overlap",
    "Unique_PMIDs",
    "Preliminary_Status",
    "Preliminary_Reason",
}

CONTEXT_REQUIRED = {
    "UniProt_ID",
    "Alias",
    "PMID",
    "Context_Match_Mode",
    "Independent_Support_In_Record",
    "Metadata_Retrieved",
}

# Aggregate context from the input records.
def aggregate_context(context):
    if context.empty:
        return pd.DataFrame(columns=[
            "_key",
            "Sampled_Context_Rows",
            "Sampled_Unique_PMIDs",
            "Exact_Context_Rows",
            "Flexible_Context_Rows",
            "No_Match_Context_Rows",
            "Independent_Support_Rows",
            "Metadata_Missing_Rows",
        ])

    c = context.copy()
    c["_key"] = key_series(c)

    rows = []

    for key, group in c.groupby("_key", sort=False):
        rows.append({
            "_key": key,
            "Sampled_Context_Rows": len(group),
            "Sampled_Unique_PMIDs": group["PMID"].map(clean).nunique(),
            "Exact_Context_Rows": int(
                group["Context_Match_Mode"].map(clean).eq("EXACT").sum()
            ),
            "Flexible_Context_Rows": int(
                group["Context_Match_Mode"].map(clean).eq(
                    "FLEXIBLE_VARIANT"
                ).sum()
            ),
            "No_Match_Context_Rows": int(
                group["Context_Match_Mode"].map(clean).eq("NO_MATCH").sum()
            ),
            "Independent_Support_Rows": int(
                group["Independent_Support_In_Record"].map(yes).sum()
            ),
            "Metadata_Missing_Rows": int(
                (~group["Metadata_Retrieved"].map(yes)).sum()
            ),
        })

    return pd.DataFrame(rows)

HGNC_NAME_ONLY_TYPES = {
    "HGNC_APPROVED_NAME",
    "HGNC_ALIAS_NAME",
    "HGNC_PREVIOUS_NAME",
}

# Classify  for the current record.
def classify(row):
    """
    Apply rules in strict priority order.

    Returns
    -------
    status, reason, evidence_note
    """

    alias = clean(row["Alias"])
    hgnc_symbol = clean(row["HGNC_Approved_Symbol"])

    retrieved = as_int(row["Retrieved_PMIDs"])
    reference = as_int(row["Independent_Reference_PMIDs"])
    overlap = as_int(row["Overlap"])
    context_support = as_int(row["Independent_Support_Rows"])

    source_disagreement = yes(row["Official_Symbol_Disagreement"])
    collision = yes(row["Cross_GPCR_Collision"])

    is_official_gene = yes(row["Is_Official_Gene"])
    is_original_gene = yes(row["Is_Original_Gene"])
    is_protein_name = yes(row["Is_Protein_Name"])

    is_hgnc_approved_symbol = (
        bool(alias)
        and bool(hgnc_symbol)
        and norm(alias) == norm(hgnc_symbol)
    )

    is_current_symbol = is_official_gene or is_hgnc_approved_symbol

    shared_sources = yes(row["Shared_By_UniProt_HGNC"])

    is_short = yes(row["Is_Short_Term"])
    is_acronym = yes(row["Is_Acronym_Risk"])
    is_generic_numbered = yes(row["Is_Generic_Numbered_Receptor"])

    z_reasons = clean(row["Flag_Z_Reasons"])
    alias_truncated = yes(row["PMID_List_Truncated"])
    reference_truncated = yes(row["Reference_PMID_List_Truncated"])

    truncation = (
        alias_truncated
        or reference_truncated
        or "TRUNCATED" in z_reasons
    )

    high_unique = (
        "HIGH_UNIQUE_RETRIEVAL" in z_reasons
        or "VERY_HIGH_UNIQUE_RETRIEVAL" in z_reasons
    )

    independent_support = overlap > 0 or context_support > 0

    term_types = split_types(row["Term_Types"])
    hgnc_name_only = (
        bool(term_types)
        and term_types.issubset(HGNC_NAME_ONLY_TYPES)
    )

    if source_disagreement:
        return (
            "MANUAL_REVIEW",
            "OFFICIAL_SYMBOL_DISAGREEMENT",
            "UniProt/HGNC current-symbol information disagrees.",
        )

    if is_current_symbol:
        return (
            "AUTO_KEEP",
            "CURRENT_OFFICIAL_SYMBOL",
            "Current UniProt/HGNC symbol for this GPCR.",
        )

    if is_protein_name and not collision:
        return (
            "AUTO_KEEP",
            "PRIMARY_PROTEIN_NAME",
            "Primary protein name for this GPCR.",
        )

    if retrieved == 0:
        if is_original_gene and not collision:
            return (
                "AUTO_KEEP",
                "ORIGINAL_GENE_NO_RETRIEVAL",
                "Original GPCR-PEnDB gene term retained as identity terminology.",
            )

        return (
            "NOT_VALIDATED",
            "NO_PUBMED_RETRIEVAL",
            "No PubMed Title/Abstract retrieval; candidate remains in audit trail.",
        )

    if collision:
        return (
            "MANUAL_REVIEW",
            "CROSS_GPCR_COLLISION",
            f"Same normalized term is associated with: "
            f"{clean(row['Collision_Genes']) or 'another GPCR'}.",
        )

    if truncation:
        return (
            "MANUAL_REVIEW",
            "PUBMED_TRUNCATION",
            "Alias/reference PMID retrieval was truncated.",
        )

    if is_generic_numbered:
        if independent_support or shared_sources:
            return (
                "CONTEXT_REQUIRED",
                "GENERIC_NUMBERED_RECEPTOR_WITH_SUPPORT",
                "Generic receptor name has supporting evidence but is unsafe standalone.",
            )

        return (
            "MANUAL_REVIEW",
            "GENERIC_NUMBERED_RECEPTOR_WITHOUT_SUPPORT",
            "Generic receptor name retrieves literature without independent support.",
        )

    if is_short or is_acronym:
        if independent_support:
            return (
                "CONTEXT_REQUIRED",
                "AMBIGUOUS_TERM_WITH_INDEPENDENT_SUPPORT",
                "Short/acronym term has independent GPCR support but may be ambiguous.",
            )

        if shared_sources:
            return (
                "CONTEXT_REQUIRED",
                "AMBIGUOUS_TERM_SHARED_BY_HGNC_UNIPROT",
                "Short/acronym term is independently documented by HGNC and UniProt.",
            )

        return (
            "MANUAL_REVIEW",
            "AMBIGUOUS_TERM_WITHOUT_SUPPORT",
            "Short/acronym term retrieves papers without independent GPCR support.",
        )

    if high_unique:
        if independent_support:
            return (
                "CONTEXT_REQUIRED",
                "HIGH_UNIQUE_RETRIEVAL_WITH_SUPPORT",
                "Term has GPCR support but a large unique retrieval component.",
            )

        return (
            "MANUAL_REVIEW",
            "HIGH_UNIQUE_RETRIEVAL_WITHOUT_SUPPORT",
            "Large unique retrieval component without independent GPCR support.",
        )

    if independent_support:
        return (
            "AUTO_KEEP",
            "NONAMBIGUOUS_WITH_INDEPENDENT_SUPPORT",
            "Non-ambiguous term has independent GPCR PMID support.",
        )

    if shared_sources:
        return (
            "AUTO_KEEP",
            "NONAMBIGUOUS_SHARED_BY_HGNC_UNIPROT",
            "Non-ambiguous term is independently documented by HGNC and UniProt.",
        )

    if hgnc_name_only:
        return (
            "NOT_VALIDATED",
            "HGNC_NAME_WITHOUT_INDEPENDENT_SUPPORT",
            "HGNC descriptive name retrieved PubMed records but lacked independent GPCR support.",
        )

    if is_original_gene:
        return (
            "MANUAL_REVIEW",
            "ORIGINAL_GENE_WITHOUT_SUPPORT",
            "Historical/original gene term retrieves papers without independent support.",
        )

    return (
        "MANUAL_REVIEW",
        "RETRIEVAL_WITHOUT_INDEPENDENT_SUPPORT",
        "Candidate retrieves papers but could not be independently resolved.",
    )

# Write summary for the current record.
def write_summary(df, manual, path):
    counts = df["Triage_Status"].value_counts().to_dict()
    reasons = df["Triage_Reason"].value_counts().to_dict()

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 03 Systematic Triage",
        "=" * 76,
        f"Total GPCR-term candidates          : {len(df)}",
        f"AUTO_KEEP                           : {counts.get('AUTO_KEEP', 0)}",
        f"CONTEXT_REQUIRED                    : {counts.get('CONTEXT_REQUIRED', 0)}",
        f"NOT_VALIDATED                       : {counts.get('NOT_VALIDATED', 0)}",
        f"MANUAL_REVIEW                       : {counts.get('MANUAL_REVIEW', 0)}",
        "",
        "MANUAL REVIEW REASONS",
        "-" * 76,
    ]

    manual_reasons = (
        manual["Triage_Reason"].value_counts()
        if not manual.empty
        else pd.Series(dtype=int)
    )

    for reason, count in manual_reasons.items():
        lines.append(f"{reason:<45} : {count}")

    lines += [
        "",
        "ALL TRIAGE REASONS",
        "-" * 76,
    ]

    for reason, count in sorted(
        reasons.items(),
        key=lambda x: (-x[1], x[0]),
    ):
        lines.append(f"{reason:<45} : {count}")

    lines += [
        "",
        "Code S5 lists unresolved aliases for review; Code S6 applies the reviewed decisions.",
        "The final alias table is built in Code S6.",
    ]

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

# Write integrity for the current record.
def write_integrity(audit, triage, manual, path):
    errors = []

    if len(triage) != len(audit):
        errors.append(
            f"Row count changed: audit={len(audit)}, triage={len(triage)}"
        )

    if key_series(triage).duplicated().any():
        errors.append("Duplicate UniProt_ID + Alias rows detected")

    statuses = set(triage["Triage_Status"].map(clean))
    invalid = statuses - ALLOWED_TRIAGE
    if invalid:
        errors.append(
            "Invalid triage status(es): " + ", ".join(sorted(invalid))
        )

    if (triage["Triage_Status"].map(clean) == "IGNORE").any():
        errors.append("Step 03 must not assign IGNORE automatically")

    expected_manual = triage[
        triage["Triage_Status"].eq("MANUAL_REVIEW")
    ]

    if len(expected_manual) != len(manual):
        errors.append(
            "Manual queue size does not match MANUAL_REVIEW subset"
        )

    if set(key_series(expected_manual)) != set(key_series(manual)):
        errors.append(
            "Manual queue contents do not exactly match MANUAL_REVIEW rows"
        )

    for _, row in triage.iterrows():
        alias = clean(row["Alias"])
        hgnc_symbol = clean(row["HGNC_Approved_Symbol"])
        current = (
            yes(row["Is_Official_Gene"])
            or (
                alias
                and hgnc_symbol
                and norm(alias) == norm(hgnc_symbol)
            )
        )

        if (
            current
            and not yes(row["Official_Symbol_Disagreement"])
            and clean(row["Triage_Status"]) != "AUTO_KEEP"
        ):
            errors.append(
                f"Current official symbol not AUTO_KEEP: "
                f"{clean(row['UniProt_ID'])} / {alias}"
            )
            break

    for _, row in triage.iterrows():
        alias = clean(row["Alias"])
        hgnc_symbol = clean(row["HGNC_Approved_Symbol"])
        current = (
            yes(row["Is_Official_Gene"])
            or (
                alias
                and hgnc_symbol
                and norm(alias) == norm(hgnc_symbol)
            )
        )

        if (
            yes(row["Cross_GPCR_Collision"])
            and not current
            and clean(row["Triage_Status"]) == "AUTO_KEEP"
        ):
            errors.append(
                f"Non-official collision was AUTO_KEEP: "
                f"{clean(row['UniProt_ID'])} / {alias}"
            )
            break

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 03 Integrity Report",
        "=" * 76,
        f"Input candidates                    : {len(audit)}",
        f"Triage rows                         : {len(triage)}",
        f"Manual review rows                  : {len(manual)}",
        f"Integrity errors                    : {len(errors)}",
        "",
    ]

    if errors:
        lines += ["ERRORS", "-" * 76]
        lines.extend(f"- {x}" for x in errors)
    else:
        lines.append("ALL STEP 03 INTEGRITY CHECKS PASSED")

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")

    return errors

# Run this stage and write its output files.
def main():
    parser = argparse.ArgumentParser(
        description="Systematic triage of GPCR alias candidates"
    )
    parser.add_argument(
        "--audit",
        default=DEFAULT_AUDIT,
        help=f"Step 01 recommendations CSV (default: {DEFAULT_AUDIT})",
    )
    parser.add_argument(
        "--context",
        default=DEFAULT_CONTEXT,
        help=f"Step 02 context CSV (default: {DEFAULT_CONTEXT})",
    )
    args = parser.parse_args()

    print("=" * 76)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 03")
    print("Systematic Alias Triage")
    print("=" * 76)

    if not os.path.exists(args.audit):
        raise SystemExit(f"ERROR: missing {args.audit}")

    if not os.path.exists(args.context):
        raise SystemExit(f"ERROR: missing {args.context}")

    audit = pd.read_csv(args.audit, dtype=str).fillna("")
    context = pd.read_csv(args.context, dtype=str).fillna("")

    require_columns(audit, AUDIT_REQUIRED, args.audit)
    require_columns(context, CONTEXT_REQUIRED, args.context)

    print(f"Step 01 candidates       : {len(audit):,}")
    print(f"Step 02 context rows     : {len(context):,}")

    context_agg = aggregate_context(context)

    triage = audit.copy()
    triage["_key"] = key_series(triage)

    triage = triage.merge(
        context_agg,
        on="_key",
        how="left",
        validate="one_to_one",
    )

    context_cols = [
        "Sampled_Context_Rows",
        "Sampled_Unique_PMIDs",
        "Exact_Context_Rows",
        "Flexible_Context_Rows",
        "No_Match_Context_Rows",
        "Independent_Support_Rows",
        "Metadata_Missing_Rows",
    ]

    for col in context_cols:
        triage[col] = triage[col].fillna(0).astype(int)

    results = triage.apply(classify, axis=1, result_type="expand")
    results.columns = [
        "Triage_Status",
        "Triage_Reason",
        "Triage_Evidence_Note",
    ]

    triage = pd.concat([triage, results], axis=1)
    triage = triage.drop(columns=["_key"])

    manual = triage[
        triage["Triage_Status"].eq("MANUAL_REVIEW")
    ].copy()

    triage.to_csv(OUT_TRIAGE, index=False)
    manual.to_csv(OUT_MANUAL, index=False)

    write_summary(triage, manual, OUT_SUMMARY)
    errors = write_integrity(
        audit,
        triage,
        manual,
        OUT_INTEGRITY,
    )

    counts = triage["Triage_Status"].value_counts().to_dict()

    print("\n" + "=" * 76)
    print("STEP 03 COMPLETE")
    print("=" * 76)
    print(f"AUTO_KEEP               : {counts.get('AUTO_KEEP', 0):,}")
    print(f"CONTEXT_REQUIRED        : {counts.get('CONTEXT_REQUIRED', 0):,}")
    print(f"NOT_VALIDATED           : {counts.get('NOT_VALIDATED', 0):,}")
    print(f"MANUAL_REVIEW           : {counts.get('MANUAL_REVIEW', 0):,}")
    print(f"Integrity errors        : {len(errors)}")
    print()
    print(f"Triage table            : {OUT_TRIAGE}")
    print(f"Manual review queue     : {OUT_MANUAL}")
    print(f"Summary                 : {OUT_SUMMARY}")
    print(f"Integrity report        : {OUT_INTEGRITY}")
    print()
    print("Code S5 lists unresolved aliases for review; Code S6 applies the reviewed decisions.")
    print("The final alias table is built in Code S6.")

    if errors:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
