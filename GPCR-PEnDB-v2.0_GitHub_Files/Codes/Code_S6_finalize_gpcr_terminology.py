#!/usr/bin/env python3
# Purpose: Combine triage and frozen review decisions into the final terminology table.
# Input: Alias_Systematic_Triage.csv; Alias_Manual_Review_Assistant_Decisions.csv
# Output: GPCR_Alias_Table_FINAL_VALIDATED.csv; Alias_Final_Decisions.csv

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

TRIAGE_FILE = "Alias_Systematic_Triage.csv"
MANUAL_FILE = "Alias_Manual_Review_Assistant_Decisions.csv"

OUT_VALIDATED = "GPCR_Alias_Table_FINAL_VALIDATED.csv"
OUT_DECISIONS = "Alias_Final_Decisions.csv"
OUT_COLLISIONS = "Alias_Final_Active_Collisions.csv"
OUT_REPORT = "Alias_Finalization_Integrity_Report.txt"

EXPECTED_GPCR_COUNT = 832

ALLOWED_TRIAGE = {
    "AUTO_KEEP",
    "CONTEXT_REQUIRED",
    "NOT_VALIDATED",
    "MANUAL_REVIEW",
}

ALLOWED_MANUAL = {
    "KEEP",
    "IGNORE",
    "CONTEXT_REQUIRED",
}

ALLOWED_FINAL = {
    "KEEP",
    "IGNORE",
    "CONTEXT_REQUIRED",
    "NOT_VALIDATED",
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

# Return unique terms while preserving the chosen order.
def unique_terms(terms):
    seen = set()
    out = []
    for term in terms:
        term = clean(term)
        key = norm(term)
        if term and key not in seen:
            seen.add(key)
            out.append(term)
    return out

# Check that required columns are present.
def require_columns(df, columns, filename):
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(
            f"{filename} is missing required column(s): "
            + ", ".join(missing)
        )

# Identify a GPCR–alias row for the decision merge.
def row_key(row):
    return (norm(row["UniProt_ID"]), norm(row["Alias"]))

# List active aliases shared between receptors.
def active_collision_table(decisions):
    """
    Report terms active for more than one GPCR.

    KEEP/KEEP collisions should normally not remain.
    CONTEXT_REQUIRED collisions are allowed and explicitly reported because
    their ambiguity is the reason they require downstream context validation.
    """
    active = decisions[
        decisions["Final_Decision"].isin(["KEEP", "CONTEXT_REQUIRED"])
    ].copy()

    owners = defaultdict(list)

    for _, row in active.iterrows():
        owners[norm(row["Alias"])].append({
            "Alias": clean(row["Alias"]),
            "Gene": clean(row["Official_Gene"]),
            "UniProt_ID": clean(row["UniProt_ID"]),
            "Decision": clean(row["Final_Decision"]),
            "Source": clean(row["Source_Provenance"]),
        })

    rows = []

    for records in owners.values():
        genes = sorted({r["Gene"] for r in records}, key=str.casefold)
        if len(genes) <= 1:
            continue

        decisions_here = sorted(
            {r["Decision"] for r in records},
            key=str.casefold,
        )

        rows.append({
            "Term": records[0]["Alias"],
            "Genes": ";".join(genes),
            "UniProt_IDs": ";".join(
                sorted(
                    {r["UniProt_ID"] for r in records},
                    key=str.casefold,
                )
            ),
            "Active_Decisions": ";".join(decisions_here),
            "Owner_Count": len(genes),
            "Contains_KEEP": "KEEP" in decisions_here,
            "Contains_CONTEXT_REQUIRED": (
                "CONTEXT_REQUIRED" in decisions_here
            ),
        })

    return pd.DataFrame(
        rows,
        columns=[
            "Term",
            "Genes",
            "UniProt_IDs",
            "Active_Decisions",
            "Owner_Count",
            "Contains_KEEP",
            "Contains_CONTEXT_REQUIRED",
        ],
    )

# Run this stage and write its output files.
def main():
    if not Path(TRIAGE_FILE).exists():
        raise SystemExit(f"ERROR: missing {TRIAGE_FILE}")
    if not Path(MANUAL_FILE).exists():
        raise SystemExit(f"ERROR: missing {MANUAL_FILE}")

    triage = pd.read_csv(TRIAGE_FILE, dtype=str).fillna("")
    manual = pd.read_csv(MANUAL_FILE, dtype=str).fillna("")

    require_columns(
        triage,
        [
            "UniProt_ID",
            "Original_Gene",
            "Official_Gene",
            "Protein_Name",
            "Alias",
            "Source_Provenance",
            "Term_Types",
            "Triage_Status",
            "Triage_Reason",
            "Triage_Evidence_Note",
        ],
        TRIAGE_FILE,
    )

    require_columns(
        manual,
        [
            "UniProt_ID",
            "Alias",
            "Final_Decision",
            "Curator_Comment",
        ],
        MANUAL_FILE,
    )

    errors = []

    if len(triage) == 0:
        errors.append("Triage table is empty")

    invalid_triage = sorted(
        set(triage["Triage_Status"].map(clean)) - ALLOWED_TRIAGE
    )
    if invalid_triage:
        errors.append(
            "Invalid triage status(es): " + ";".join(invalid_triage)
        )

    if triage.apply(row_key, axis=1).duplicated().any():
        errors.append("Duplicate UniProt_ID + Alias rows in triage table")

    if manual.apply(row_key, axis=1).duplicated().any():
        errors.append("Duplicate UniProt_ID + Alias rows in manual file")

    manual_keys_expected = {
        row_key(row)
        for _, row in triage[
            triage["Triage_Status"].eq("MANUAL_REVIEW")
        ].iterrows()
    }

    manual_keys_found = {
        row_key(row)
        for _, row in manual.iterrows()
    }

    if manual_keys_found != manual_keys_expected:
        missing = manual_keys_expected - manual_keys_found
        extra = manual_keys_found - manual_keys_expected
        if missing:
            errors.append(
                f"{len(missing)} MANUAL_REVIEW rows are missing decisions"
            )
        if extra:
            errors.append(
                f"{len(extra)} manual rows are not in MANUAL_REVIEW"
            )

    manual_map = {}
    for _, row in manual.iterrows():
        key = row_key(row)
        decision = clean(row["Final_Decision"]).upper()
        comment = clean(row["Curator_Comment"])

        if decision not in ALLOWED_MANUAL:
            errors.append(
                f"Invalid manual decision for "
                f"{clean(row['UniProt_ID'])}/{clean(row['Alias'])}: "
                f"{decision or '<blank>'}"
            )
            continue

        if not comment:
            errors.append(
                f"Blank curator comment for "
                f"{clean(row['UniProt_ID'])}/{clean(row['Alias'])}"
            )

        manual_map[key] = row

    decision_rows = []

    for _, row in triage.iterrows():
        status = clean(row["Triage_Status"])
        key = row_key(row)

        if status == "AUTO_KEEP":
            final = "KEEP"
            decision_source = "SYSTEMATIC_TRIAGE"
            reason = clean(row["Triage_Reason"])
            comment = clean(row["Triage_Evidence_Note"])
            confidence = "SYSTEMATIC"
            evidence_pmids = ""
            needs_spot_check = "False"

        elif status == "CONTEXT_REQUIRED":
            final = "CONTEXT_REQUIRED"
            decision_source = "SYSTEMATIC_TRIAGE"
            reason = clean(row["Triage_Reason"])
            comment = clean(row["Triage_Evidence_Note"])
            confidence = "SYSTEMATIC"
            evidence_pmids = ""
            needs_spot_check = "False"

        elif status == "NOT_VALIDATED":
            final = "NOT_VALIDATED"
            decision_source = "SYSTEMATIC_TRIAGE"
            reason = clean(row["Triage_Reason"])
            comment = clean(row["Triage_Evidence_Note"])
            confidence = "SYSTEMATIC"
            evidence_pmids = ""
            needs_spot_check = "False"

        elif status == "MANUAL_REVIEW":
            m = manual_map.get(key)

            if m is None:
                errors.append(
                    f"Missing manual decision: "
                    f"{clean(row['UniProt_ID'])}/{clean(row['Alias'])}"
                )
                continue

            final = clean(m["Final_Decision"]).upper()
            decision_source = "MANUAL_REVIEW_ASSISTED"
            reason = clean(
                m["Triage_Reason"]
                if "Triage_Reason" in manual.columns
                else row["Triage_Reason"]
            )
            comment = clean(m["Curator_Comment"])
            confidence = clean(
                m["Review_Confidence"]
                if "Review_Confidence" in manual.columns
                else ""
            )
            evidence_pmids = clean(
                m["Evidence_PMIDs"]
                if "Evidence_PMIDs" in manual.columns
                else ""
            )
            needs_spot_check = clean(
                m["Needs_Spot_Check"]
                if "Needs_Spot_Check" in manual.columns
                else "False"
            )

        else:
            errors.append(
                f"Unhandled triage status: "
                f"{clean(row['UniProt_ID'])}/{clean(row['Alias'])}/{status}"
            )
            continue

        if final not in ALLOWED_FINAL:
            errors.append(
                f"Invalid final decision: "
                f"{clean(row['UniProt_ID'])}/{clean(row['Alias'])}/{final}"
            )

        decision_rows.append({
            "UniProt_ID": clean(row["UniProt_ID"]),
            "Original_Gene": clean(row["Original_Gene"]),
            "Official_Gene": clean(row["Official_Gene"]),
            "Protein_Name": clean(row["Protein_Name"]),
            "Alias": clean(row["Alias"]),
            "Source_Provenance": clean(row["Source_Provenance"]),
            "Term_Types": clean(row["Term_Types"]),
            "Triage_Status": status,
            "Triage_Reason": clean(row["Triage_Reason"]),
            "Final_Decision": final,
            "Decision_Source": decision_source,
            "Decision_Comment": comment,
            "Review_Confidence": confidence,
            "Evidence_PMIDs": evidence_pmids,
            "Needs_Spot_Check": needs_spot_check,
        })

    decisions = pd.DataFrame(decision_rows)

    if len(decisions) != len(triage):
        errors.append(
            f"Expected {len(triage)} final decision rows; "
            f"found {len(decisions)}"
        )

    final_rows = []

    for uniprot, group in decisions.groupby("UniProt_ID", sort=False):
        first = group.iloc[0]

        keep = unique_terms(
            group.loc[
                group["Final_Decision"].eq("KEEP"),
                "Alias",
            ].tolist()
        )

        context = unique_terms(
            group.loc[
                group["Final_Decision"].eq("CONTEXT_REQUIRED"),
                "Alias",
            ].tolist()
        )

        ignored = unique_terms(
            group.loc[
                group["Final_Decision"].eq("IGNORE"),
                "Alias",
            ].tolist()
        )

        not_validated = unique_terms(
            group.loc[
                group["Final_Decision"].eq("NOT_VALIDATED"),
                "Alias",
            ].tolist()
        )

        categories = {
            "KEEP": {norm(x) for x in keep},
            "CONTEXT_REQUIRED": {norm(x) for x in context},
            "IGNORE": {norm(x) for x in ignored},
            "NOT_VALIDATED": {norm(x) for x in not_validated},
        }

        labels = list(categories)
        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                overlap = categories[a] & categories[b]
                if overlap:
                    errors.append(
                        f"{clean(first['Official_Gene'])}: "
                        f"term(s) occur in both {a} and {b}: "
                        f"{sorted(overlap)}"
                    )

        active = unique_terms(keep + context)

        if not keep:
            errors.append(
                f"{clean(first['Official_Gene'])} ({uniprot}) "
                "has no final KEEP search term"
            )

        if not active:
            errors.append(
                f"{clean(first['Official_Gene'])} ({uniprot}) "
                "has no active retrieval term"
            )

        final_rows.append({
            "UniProt_ID": clean(uniprot),
            "Original_Gene": clean(first["Original_Gene"]),
            "Official_Gene": clean(first["Official_Gene"]),
            "Protein_Name": clean(first["Protein_Name"]),
            "Search_Terms": ";".join(keep),
            "Context_Required_Terms": ";".join(context),
            "Active_Retrieval_Terms": ";".join(active),
            "Ignored_Aliases": ";".join(ignored),
            "Not_Validated_Aliases": ";".join(not_validated),
            "KEEP_Term_Count": len(keep),
            "CONTEXT_REQUIRED_Term_Count": len(context),
            "Active_Retrieval_Term_Count": len(active),
            "IGNORE_Term_Count": len(ignored),
            "NOT_VALIDATED_Term_Count": len(not_validated),
        })

    final = pd.DataFrame(final_rows)

    if len(final) != EXPECTED_GPCR_COUNT:
        errors.append(
            f"Expected {EXPECTED_GPCR_COUNT} GPCR rows; found {len(final)}"
        )

    if final["UniProt_ID"].map(norm).duplicated().any():
        errors.append("Duplicate UniProt_ID rows in final alias table")

    if final["Official_Gene"].map(norm).duplicated().any():
        errors.append("Duplicate Official_Gene rows in final alias table")

    if final["Search_Terms"].map(clean).eq("").any():
        errors.append("At least one GPCR has empty Search_Terms")

    missing_official = []

    for _, row in final.iterrows():
        keep_keys = {norm(x) for x in split_terms(row["Search_Terms"])}
        if norm(row["Official_Gene"]) not in keep_keys:
            missing_official.append(clean(row["Official_Gene"]))

    if missing_official:
        errors.append(
            "Official gene symbol missing from KEEP Search_Terms for: "
            + ";".join(missing_official[:30])
        )

    collisions = active_collision_table(decisions)

    if not collisions.empty:
        keep_keep = collisions[
            collisions["Active_Decisions"].eq("KEEP")
        ]
        if not keep_keep.empty:
            errors.append(
                f"{len(keep_keep)} KEEP/KEEP cross-GPCR collisions remain"
            )

    counts = decisions["Final_Decision"].value_counts().to_dict()

    total_keep = counts.get("KEEP", 0)
    total_context = counts.get("CONTEXT_REQUIRED", 0)
    total_ignore = counts.get("IGNORE", 0)
    total_not_validated = counts.get("NOT_VALIDATED", 0)

    if (
        total_keep
        + total_context
        + total_ignore
        + total_not_validated
        != len(decisions)
    ):
        errors.append("Final decision counts do not balance")

    low_confidence = int(
        decisions["Needs_Spot_Check"]
        .map(clean)
        .str.upper()
        .isin({"TRUE", "YES", "1"})
        .sum()
    )

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 04 Integrity Report",
        "=" * 76,
        f"GPCR rows                         : {len(final)}",
        f"Candidate GPCR-term pairs         : {len(decisions)}",
        f"Final KEEP terms                  : {total_keep}",
        f"Final CONTEXT_REQUIRED terms      : {total_context}",
        f"Final IGNORE terms                : {total_ignore}",
        f"Final NOT_VALIDATED terms         : {total_not_validated}",
        f"Active retrieval terms            : {total_keep + total_context}",
        f"Manual decisions                  : {len(manual)}",
        f"Manual decisions flagged spot-check: {low_confidence}",
        f"Active cross-GPCR collision terms : {len(collisions)}",
        f"Empty final Search_Terms rows     : "
        f"{int(final['Search_Terms'].map(clean).eq('').sum())}",
        f"Integrity errors                  : {len(errors)}",
        "",
    ]

    if not collisions.empty:
        lines.append("ACTIVE CROSS-GPCR COLLISIONS")
        lines.append("-" * 76)
        lines.append(
            "Allowed only because at least one side is CONTEXT_REQUIRED."
        )
        for _, item in collisions.iterrows():
            lines.append(
                f"- {item['Term']}: {item['Genes']} "
                f"[{item['Active_Decisions']}]"
            )
        lines.append("")

    if errors:
        lines.append("ERRORS")
        lines.append("-" * 76)
        lines.extend(f"- {error}" for error in errors)

        Path(OUT_REPORT).write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8",
        )

        raise SystemExit(
            f"Step 04 failed with {len(errors)} integrity error(s). "
            f"See {OUT_REPORT}."
        )

    lines.append("ALL STEP 04 FINALIZATION INTEGRITY CHECKS PASSED")

    final.to_csv(OUT_VALIDATED, index=False)
    decisions.to_csv(OUT_DECISIONS, index=False)
    collisions.to_csv(OUT_COLLISIONS, index=False)
    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(lines))
    print("\nSaved:")
    print(f"  {OUT_VALIDATED}")
    print(f"  {OUT_DECISIONS}")
    print(f"  {OUT_COLLISIONS}")
    print(f"  {OUT_REPORT}")

if __name__ == "__main__":
    main()
