#!/usr/bin/env python3
# Purpose: Apply the frozen collision rule and finalize context resolution.
# Input: Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv; GPCR_Cancer_Context_Manual_Review.csv
# Output: Human_GPCR_Cancer_PubMed_Papers_FINAL.csv

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import os
import re
import sys
from pathlib import Path

import pandas as pd

INPUT_VALIDATED = "Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv"
INPUT_MANUAL = "GPCR_Cancer_Context_Manual_Review.csv"

OUT_FINAL = "Human_GPCR_Cancer_PubMed_Papers_FINAL.csv"
OUT_MANUAL_DECISIONS = "GPCR_Cancer_Context_Manual_Resolution.csv"
OUT_REPORT = "GPCR_Cancer_Context_Finalization_Integrity_Report.txt"

ALLOWED_COLLISION_TERMS = {
    "ccr9",
    "ccr10",
    "chemokine receptor ccr-9",
}

ALLOWED_COLLISION_GENES = {
    "ackr2",
    "ccr9",
    "ccr10",
}

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
    seen = set()
    out = []

    for item in clean(value).split(";"):
        term = clean(item)
        key = norm(term)

        if term and key not in seen:
            seen.add(key)
            out.append(term)

    return out

# Check that required columns are present.
def require_columns(df, required, filename):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(
            f"{filename} is missing required column(s): "
            + ", ".join(missing)
        )

# Run this stage and write its output files.
def main():
    for filename in [INPUT_VALIDATED, INPUT_MANUAL]:
        if not os.path.exists(filename):
            raise SystemExit(f"ERROR: missing {filename}")

    validated = pd.read_csv(
        INPUT_VALIDATED,
        dtype={"PMID": str},
    ).fillna("")

    manual = pd.read_csv(
        INPUT_MANUAL,
        dtype={"PMID": str},
    ).fillna("")

    required = {
        "UniProt_ID",
        "Official_Gene",
        "TCGA_Code",
        "PMID",
        "Matched_KEEP_GPCR_Terms",
        "Matched_CONTEXT_REQUIRED_GPCR_Terms",
        "Context_Resolution",
        "Context_Resolution_Reason",
        "Strong_Context_Terms",
        "Active_Collision_Terms",
        "Active_Collision_Genes",
    }

    require_columns(manual, required, INPUT_MANUAL)

    errors = []
    decision_rows = []

    if manual.duplicated(
        subset=["UniProt_ID", "TCGA_Code", "PMID"]
    ).any():
        errors.append("Duplicate manual-review association rows detected")

    for _, row in manual.iterrows():
        gene = norm(row["Official_Gene"])
        matched_keep = split_terms(row["Matched_KEEP_GPCR_Terms"])
        matched_context = {
            norm(x)
            for x in split_terms(
                row["Matched_CONTEXT_REQUIRED_GPCR_Terms"]
            )
        }
        collision_terms = {
            norm(x)
            for x in split_terms(row["Active_Collision_Terms"])
        }
        collision_genes = {
            norm(x)
            for x in split_terms(row["Active_Collision_Genes"])
        }

        reasons = []

        if gene != "ackr2":
            reasons.append("INTENDED_GPCR_NOT_ACKR2")

        if matched_keep:
            reasons.append("FINAL_KEEP_TERM_PRESENT")

        if not matched_context:
            reasons.append("NO_CONTEXT_TERM_PRESENT")

        if not matched_context.issubset(ALLOWED_COLLISION_TERMS):
            reasons.append("UNEXPECTED_CONTEXT_TERM")

        if not collision_terms:
            reasons.append("NO_ACTIVE_COLLISION_TERM")

        if not collision_terms.issubset(ALLOWED_COLLISION_TERMS):
            reasons.append("UNEXPECTED_COLLISION_TERM")

        if not collision_genes:
            reasons.append("NO_COLLISION_GENES")

        if not collision_genes.issubset(ALLOWED_COLLISION_GENES):
            reasons.append("UNEXPECTED_COLLISION_GENE")

        if not ({"ccr9", "ccr10"} & collision_genes):
            reasons.append("NO_CURRENT_CCR9_OR_CCR10_OWNER")

        if clean(row["Context_Resolution"]) != "MANUAL_REVIEW":
            reasons.append("ROW_NOT_MANUAL_REVIEW")

        if reasons:
            decision = "UNRESOLVED"
            comment = ";".join(reasons)
            errors.append(
                f"{clean(row['Official_Gene'])}/"
                f"{clean(row['TCGA_Code'])}/"
                f"{clean(row['PMID'])}: "
                + ";".join(reasons)
            )
        else:
            decision = "REJECT"
            comment = (
                "ACKR2 assignment is supported only by active CCR9/CCR10 "
                "collision terminology, with no final KEEP ACKR2 term in "
                "the article; reject this paper-to-ACKR2 assignment."
            )

        out = dict(row)
        out["Manual_Final_Decision"] = decision
        out["Manual_Final_Reason"] = (
            "ACKR2_CCR9_CCR10_COLLISION"
            if decision == "REJECT"
            else "UNRESOLVED_RULE_MISMATCH"
        )
        out["Manual_Final_Comment"] = comment
        decision_rows.append(out)

    manual_decisions = pd.DataFrame(decision_rows)

    if errors:
        manual_decisions.to_csv(
            OUT_MANUAL_DECISIONS,
            index=False,
            na_rep="",
        )

        report = [
            "GPCR Cancer Pipeline - Step 09 Finalization Integrity Report",
            "=" * 78,
            f"Automatic/context-resolved validated rows : {len(validated)}",
            f"Manual-review rows                        : {len(manual)}",
            f"Resolved as REJECT                        : "
            f"{int((manual_decisions['Manual_Final_Decision'] == 'REJECT').sum())}",
            f"Unresolved rule mismatches                : {len(errors)}",
            "",
            "FAILED - manual queue contained rows outside the frozen rule.",
            "",
        ]
        report.extend(f"- {x}" for x in errors)

        Path(OUT_REPORT).write_text(
            "\n".join(report) + "\n",
            encoding="utf-8",
        )
        raise SystemExit(1)

    final = validated.copy()

    final["Final_Literature_Status"] = "KEEP"

    # Identify the evidence source used for validation.
    def validation_source(row):
        context_resolution = clean(row.get("Context_Resolution", ""))

        if context_resolution == "AUTO_KEEP":
            return "STEP08_CONTEXT_RESOLVED"

        if context_resolution == "NOT_REQUIRED":
            return "STEP07_FINAL_KEEP_TERM"

        return "VALIDATED_PRIOR_STEP"

    final["Final_Validation_Source"] = final.apply(
        validation_source,
        axis=1,
    )

    # Collect receptor terms that support a paper.
    def supporting_terms(row):
        if clean(row.get("Context_Resolution", "")) == "AUTO_KEEP":
            return clean(row.get("Strong_Context_Terms", ""))

        return clean(row.get("Matched_KEEP_GPCR_Terms", ""))

    final["Validated_Search_Terms"] = final.apply(
        supporting_terms,
        axis=1,
    )

    if final["Validated_Search_Terms"].map(clean).eq("").any():
        bad = final[
            final["Validated_Search_Terms"].map(clean).eq("")
        ][
            ["Official_Gene", "TCGA_Code", "PMID", "Context_Resolution"]
        ]
        raise SystemExit(
            "ERROR: final validated row has no receptor-supporting term:\n"
            + bad.head(30).to_string(index=False)
        )

    if final.duplicated(
        subset=["UniProt_ID", "TCGA_Code", "PMID"]
    ).any():
        raise SystemExit(
            "ERROR: duplicate GPCR-cancer-PMID rows in final literature set"
        )

    manual_decisions.to_csv(
        OUT_MANUAL_DECISIONS,
        index=False,
        na_rep="",
    )

    final.to_csv(
        OUT_FINAL,
        index=False,
        na_rep="",
    )

    source_counts = (
        final["Final_Validation_Source"]
        .value_counts()
        .to_dict()
    )

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 09 Integrity Report",
        "=" * 78,
        "",
        "MANUAL CONTEXT RESOLUTION",
        "-" * 78,
        f"Manual-review input rows               : {len(manual)}",
        f"Rejected by frozen collision rule      : {len(manual_decisions)}",
        f"Unresolved rows                        : 0",
        "",
        "FINAL VALIDATED LITERATURE",
        "-" * 78,
        f"Final GPCR-cancer-PMID rows            : {len(final)}",
        f"Unique PMIDs                           : {final['PMID'].nunique()}",
        f"GPCRs with >=1 validated paper         : {final['Official_Gene'].nunique()}",
        f"GPCR-cancer associations before category classification: "
        f"{final[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]}",
        f"Step07 final-KEEP-supported rows       : "
        f"{source_counts.get('STEP07_FINAL_KEEP_TERM', 0)}",
        f"Step08 context-resolved rows           : "
        f"{source_counts.get('STEP08_CONTEXT_RESOLVED', 0)}",
        "",
        "INTEGRITY CHECKS",
        "-" * 78,
        "All manual rows matched frozen rule: PASS",
        "All manual rows resolved: PASS",
        "No duplicate final associations: PASS",
        "Every final row has supporting GPCR term(s): PASS",
        "",
        "METHOD NOTE",
        "-" * 78,
        (
            "The Step 08 manual-review remainder consisted exclusively of "
            "ACKR2 records supported by active CCR9/CCR10 collision terms. "
            "Because no final KEEP ACKR2 term occurred in these records, "
            "the paper-to-ACKR2 assignments were rejected as a single "
            "predefined collision-resolution case. No alias was deleted "
            "from the terminology audit by this step."
        ),
        "",
        "ALL STEP 09 INTEGRITY CHECKS PASSED",
    ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(lines))
    print()
    print("Saved:")
    print(f"  {OUT_FINAL}")
    print(f"  {OUT_MANUAL_DECISIONS}")
    print(f"  {OUT_REPORT}")

if __name__ == "__main__":
    main()
