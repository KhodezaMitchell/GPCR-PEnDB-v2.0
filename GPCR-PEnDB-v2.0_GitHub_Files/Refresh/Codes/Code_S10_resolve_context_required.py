#!/usr/bin/env python3
# Purpose: Resolve publications retrieved only through context-required GPCR names.
# Input: GPCR_Cancer_PubMed_Context_Review.csv; Human_GPCR_Cancer_PubMed_Papers_Validated.csv; GPCR_Alias_Table_FINAL_VALIDATED.csv
# Output: Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv; GPCR_Cancer_Context_Manual_Review.csv

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

CONTEXT_FILE = "GPCR_Cancer_PubMed_Context_Review.csv"
AUTO_KEEP_FILE = "Human_GPCR_Cancer_PubMed_Papers_Validated.csv"
ALIAS_FILE = "GPCR_Alias_Table_FINAL_VALIDATED.csv"

OUT_ALL = "GPCR_Cancer_Context_Resolution_All.csv"
OUT_MANUAL = "GPCR_Cancer_Context_Manual_Review.csv"
OUT_REJECTED = "GPCR_Cancer_Context_Rejected.csv"
OUT_COMBINED = "Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv"
OUT_SUMMARY = "GPCR_Cancer_Context_Resolution_Summary.csv"
OUT_REPORT = "GPCR_Cancer_Context_Resolution_Integrity_Report.txt"

EXPECTED_GPCR_COUNT = 832

ALLOWED_DECISIONS = {
    "AUTO_KEEP",
    "MANUAL_REVIEW",
    "REJECT",
}

# Standardize whitespace and missing values in text.
def clean(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() == "nan":
        return ""
    return re.sub(r"\s+", " ", value)

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Choose the normalized display form of a receptor term.
def canonical_term(value):
    """
    Used only to decide whether two matched aliases are genuinely distinct.
    PAR-1 and PAR1 therefore count as one term.
    """
    return re.sub(r"[^a-z0-9]+", "", norm(value))

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

# Check whether a term mixes uppercase and lowercase.
def is_mixed_case(term):
    letters = [c for c in clean(term) if c.isalpha()]

    return (
        bool(letters)
        and any(c.islower() for c in letters)
        and any(c.isupper() for c in letters)
    )

# Check whether the term appears with the required case.
def exact_case_occurs(text, term):
    """
    Case-sensitive occurrence with alphanumeric boundaries.
    """
    term = clean(term)
    if not term:
        return False

    return (
        re.search(
            rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])",
            clean(text),
        )
        is not None
    )

# Identify an ambiguous numbered receptor name.
def generic_numbered_receptor(term):
    """
    Generic numbered descriptors are too nonspecific for auto-acceptance.
    Examples:
      G protein-coupled receptor 1
      G-protein coupled receptor 5
      orphan G protein-coupled receptor 12
    """
    t = norm(term)

    patterns = [
        r"^(?:orphan\s+)?g[\s-]*protein[\s-]*coupled\s+receptor\s+\d+[a-z]?$",
        r"^(?:orphan\s+)?gpcr\s+\d+[a-z]?$",
    ]

    return any(re.fullmatch(pattern, t) for pattern in patterns)

# Identify a descriptive receptor name.
def descriptive_receptor_term(term):
    t = norm(term)

    if generic_numbered_receptor(term):
        return False

    return (
        "receptor" in t
        or "g protein-coupled" in t
        or "g-protein-coupled" in t
    )

# Identify a receptor-specific symbol.
def receptor_nomenclature_symbol(term):
    """
    Narrow, reproducible receptor-symbol families.

    This deliberately does NOT classify all letter+number gene symbols as
    receptor-specific. It is limited to conventional historical GPCR labels
    whose prefix itself denotes a receptor family.
    """
    token = re.sub(r"[^A-Za-z0-9]", "", clean(term)).upper()

    return bool(
        re.fullmatch(
            r"(?:GPR|GPCR|CCR|CXCR|ACKR)\d+[A-Z0-9]*",
            token,
        )
    )

# Find an explicit link between alias and receptor.
def explicit_alias_receptor_relation(text, alias):
    """
    Require a local/direct alias <-> receptor naming relationship.

    Broad "same sentence contains receptor" logic is intentionally avoided
    because it would incorrectly validate examples such as an unrelated alias
    appearing in a paper that happens to mention an androgen receptor.
    """
    text = clean(text)
    alias = clean(alias)

    if not text or not alias:
        return False

    a = re.escape(alias)

    patterns = [
        rf"(?<![A-Za-z0-9]){a}(?![A-Za-z0-9])[\s-]+receptor\b",

        rf"\breceptor[\s-]+(?<![A-Za-z0-9]){a}(?![A-Za-z0-9])",

        rf"\breceptor\s*\(\s*{a}\s*\)",

        rf"\b(?:[A-Za-z0-9α-ωΑ-Ω/+\-]+[\s-]+){{0,6}}"
        rf"receptor\s*\(\s*{a}\s*\)",

        rf"(?<![A-Za-z0-9]){a}(?![A-Za-z0-9])\s*\("
        rf"\s*(?:[A-Za-z0-9α-ωΑ-Ω/+\-]+[\s-]+){{0,6}}"
        rf"receptor\s*\)",
    ]

    return any(
        re.search(pattern, text, flags=re.IGNORECASE)
        for pattern in patterns
    )

# Build active owners from the input records.
def build_active_owners(alias_df):
    """
    Determine which GPCR(s) own each currently active alias in the frozen
    vocabulary. Both KEEP and CONTEXT_REQUIRED terms count as active.
    """
    owners = defaultdict(set)

    for _, row in alias_df.iterrows():
        gene = clean(row["Official_Gene"])

        terms = (
            split_terms(row["Search_Terms"])
            + split_terms(row["Context_Required_Terms"])
        )

        for term in terms:
            owners[norm(term)].add(gene)

    return owners

# Resolve a context-required paper using receptor evidence.
def classify_context_row(row, active_owners):
    text = clean(
        f"{clean(row['Title'])}. {clean(row['Abstract'])}"
    )

    matched = split_terms(
        row["Matched_CONTEXT_REQUIRED_GPCR_Terms"]
    )

    gene = clean(row["Official_Gene"])

    usable_terms = []
    case_mismatch_terms = []

    for term in matched:
        if is_mixed_case(term) and not exact_case_occurs(text, term):
            case_mismatch_terms.append(term)
        else:
            usable_terms.append(term)

    distinct_usable = {
        canonical_term(term)
        for term in usable_terms
        if canonical_term(term)
    }

    strong_terms = []
    reasons = []

    if len(distinct_usable) >= 2:
        strong_terms.extend(usable_terms)
        reasons.append("MULTIPLE_CONTEXT_TERMS")

    for term in usable_terms:
        if generic_numbered_receptor(term):
            continue

        if descriptive_receptor_term(term):
            strong_terms.append(term)
            reasons.append("DESCRIPTIVE_RECEPTOR_TERM")
            continue

        if explicit_alias_receptor_relation(text, term):
            strong_terms.append(term)
            reasons.append("EXPLICIT_ALIAS_RECEPTOR_RELATION")
            continue

        if receptor_nomenclature_symbol(term):
            strong_terms.append(term)
            reasons.append("RECEPTOR_NOMENCLATURE_SYMBOL")
            continue

    seen = set()
    strong_terms_unique = []

    for term in strong_terms:
        key = norm(term)
        if key not in seen:
            seen.add(key)
            strong_terms_unique.append(term)

    if not strong_terms_unique:
        reason = "CONTEXT_REQUIREMENT_NOT_MET"

        if case_mismatch_terms and not usable_terms:
            reason = "CASE_MISMATCH_ONLY"

        return {
            "Context_Resolution": "REJECT",
            "Context_Resolution_Reason": reason,
            "Strong_Context_Terms": "",
            "Context_Evidence_Types": "",
            "Case_Mismatch_Terms": ";".join(case_mismatch_terms),
            "Active_Collision_Terms": "",
            "Active_Collision_Genes": "",
        }

    collision_terms = []
    collision_genes = set()
    noncollision_strong = []

    for term in strong_terms_unique:
        owners = active_owners.get(norm(term), set())

        if len(owners) > 1:
            collision_terms.append(term)
            collision_genes.update(owners)
        else:
            noncollision_strong.append(term)

    if strong_terms_unique and not noncollision_strong:
        return {
            "Context_Resolution": "MANUAL_REVIEW",
            "Context_Resolution_Reason": (
                "SUFFICIENT_CONTEXT_BUT_ACTIVE_CROSS_GPCR_COLLISION"
            ),
            "Strong_Context_Terms": ";".join(strong_terms_unique),
            "Context_Evidence_Types": ";".join(
                sorted(set(reasons))
            ),
            "Case_Mismatch_Terms": ";".join(case_mismatch_terms),
            "Active_Collision_Terms": ";".join(collision_terms),
            "Active_Collision_Genes": ";".join(
                sorted(collision_genes, key=str.casefold)
            ),
        }

    return {
        "Context_Resolution": "AUTO_KEEP",
        "Context_Resolution_Reason": (
            "CONTEXT_REQUIREMENT_SATISFIED"
        ),
        "Strong_Context_Terms": ";".join(
            noncollision_strong
        ),
        "Context_Evidence_Types": ";".join(
            sorted(set(reasons))
        ),
        "Case_Mismatch_Terms": ";".join(case_mismatch_terms),
        "Active_Collision_Terms": ";".join(collision_terms),
        "Active_Collision_Genes": ";".join(
            sorted(collision_genes, key=str.casefold)
        ),
    }

CONTEXT_REQUIRED_COLUMNS = {
    "UniProt_ID",
    "Gene",
    "Official_Gene",
    "TCGA_Code",
    "Cancer",
    "PMID",
    "Search_ID",
    "Title",
    "Abstract",
    "Cancer_Validation",
    "GPCR_Validation",
    "Matched_KEEP_GPCR_Terms",
    "Matched_CONTEXT_REQUIRED_GPCR_Terms",
    "Overall_Validation",
    "Final_KEEP_GPCR_Terms",
    "Final_CONTEXT_REQUIRED_GPCR_Terms",
}

ALIAS_REQUIRED_COLUMNS = {
    "UniProt_ID",
    "Official_Gene",
    "Search_Terms",
    "Context_Required_Terms",
    "Active_Retrieval_Terms",
}

# Load inputs from the input records.
def load_inputs():
    for filename in [
        CONTEXT_FILE,
        AUTO_KEEP_FILE,
        ALIAS_FILE,
    ]:
        if not os.path.exists(filename):
            raise FileNotFoundError(
                f"Missing required input: {filename}"
            )

    context = pd.read_csv(
        CONTEXT_FILE,
        dtype={"PMID": str},
    ).fillna("")

    auto_keep = pd.read_csv(
        AUTO_KEEP_FILE,
        dtype={"PMID": str},
    ).fillna("")

    aliases = pd.read_csv(
        ALIAS_FILE,
        dtype=str,
    ).fillna("")

    require_columns(
        context,
        CONTEXT_REQUIRED_COLUMNS,
        CONTEXT_FILE,
    )

    require_columns(
        aliases,
        ALIAS_REQUIRED_COLUMNS,
        ALIAS_FILE,
    )

    errors = []

    if len(aliases) != EXPECTED_GPCR_COUNT:
        errors.append(
            f"Expected {EXPECTED_GPCR_COUNT} GPCRs; found {len(aliases)}"
        )

    if aliases["UniProt_ID"].map(norm).duplicated().any():
        errors.append("Duplicate UniProt_ID in final alias table")

    if context.duplicated(
        subset=["UniProt_ID", "TCGA_Code", "PMID"]
    ).any():
        errors.append(
            "Duplicate GPCR-cancer-PMID rows in context review"
        )

    if not context["Overall_Validation"].eq(
        "CONTEXT_REVIEW"
    ).all():
        errors.append(
            "Context input contains rows not labeled CONTEXT_REVIEW"
        )

    if not context["Cancer_Validation"].eq("KEEP").all():
        errors.append(
            "Context input contains rows without cancer KEEP validation"
        )

    if not context["GPCR_Validation"].eq(
        "CONTEXT_REVIEW"
    ).all():
        errors.append(
            "Context input contains rows without GPCR CONTEXT_REVIEW"
        )

    if context["Matched_CONTEXT_REQUIRED_GPCR_Terms"].map(
        clean
    ).eq("").any():
        errors.append(
            "Context-review row lacks a matched CONTEXT_REQUIRED term"
        )

    if context["Matched_KEEP_GPCR_Terms"].map(clean).ne("").any():
        errors.append(
            "Context-review row unexpectedly contains a matched KEEP term"
        )

    alias_lookup = {
        clean(row["UniProt_ID"]): {
            norm(x)
            for x in split_terms(
                row["Context_Required_Terms"]
            )
        }
        for _, row in aliases.iterrows()
    }

    for _, row in context.iterrows():
        uid = clean(row["UniProt_ID"])
        allowed = alias_lookup.get(uid)

        if allowed is None:
            errors.append(
                f"Context GPCR absent from final alias table: {uid}"
            )
            break

        found = {
            norm(x)
            for x in split_terms(
                row["Matched_CONTEXT_REQUIRED_GPCR_Terms"]
            )
        }

        if not found.issubset(allowed):
            errors.append(
                f"{uid}/{clean(row['PMID'])}: matched context term "
                "not in frozen Context_Required_Terms"
            )
            break

    if errors:
        raise ValueError(
            "\n".join(f"- {x}" for x in errors)
        )

    return context, auto_keep, aliases

# Build summary from the input records.
def build_summary(resolved):
    rows = []

    for (uid, gene), group in resolved.groupby(
        ["UniProt_ID", "Official_Gene"],
        sort=False,
    ):
        counts = (
            group["Context_Resolution"]
            .value_counts()
            .to_dict()
        )

        rows.append({
            "UniProt_ID": uid,
            "Official_Gene": gene,
            "Context_Review_Rows": len(group),
            "AUTO_KEEP": counts.get("AUTO_KEEP", 0),
            "MANUAL_REVIEW": counts.get("MANUAL_REVIEW", 0),
            "REJECT": counts.get("REJECT", 0),
            "Cancer_Types_Auto_Kept": (
                group.loc[
                    group["Context_Resolution"].eq("AUTO_KEEP"),
                    "TCGA_Code",
                ]
                .nunique()
            ),
            "Unique_PMIDs_Auto_Kept": (
                group.loc[
                    group["Context_Resolution"].eq("AUTO_KEEP"),
                    "PMID",
                ]
                .nunique()
            ),
        })

    return pd.DataFrame(rows)

# Write report for the current record.
def write_report(
    context,
    auto_keep,
    resolved,
    manual,
    rejected,
    combined,
    aliases,
    summary,
):
    counts = (
        resolved["Context_Resolution"]
        .value_counts()
        .to_dict()
    )

    auto_keep_genes_before = set(
        auto_keep["Official_Gene"].map(clean)
    )
    resolved_keep_genes = set(
        resolved.loc[
            resolved["Context_Resolution"].eq("AUTO_KEEP"),
            "Official_Gene",
        ].map(clean)
    )
    combined_genes = set(
        combined["Official_Gene"].map(clean)
    )

    newly_recovered_genes = sorted(
        resolved_keep_genes - auto_keep_genes_before,
        key=str.casefold,
    )

    checks = {
        "832 GPCR alias rows": (
            len(aliases) == EXPECTED_GPCR_COUNT
        ),
        "Every context row resolved once": (
            len(resolved) == len(context)
        ),
        "Allowed resolution statuses only": (
            set(resolved["Context_Resolution"])
            .issubset(ALLOWED_DECISIONS)
        ),
        "No duplicate context associations": (
            not resolved.duplicated(
                subset=["UniProt_ID", "TCGA_Code", "PMID"]
            ).any()
        ),
        "Resolution counts balance": (
            sum(counts.get(x, 0) for x in ALLOWED_DECISIONS)
            == len(resolved)
        ),
        "Manual queue matches MANUAL_REVIEW": (
            len(manual) == counts.get("MANUAL_REVIEW", 0)
        ),
        "Rejected file matches REJECT": (
            len(rejected) == counts.get("REJECT", 0)
        ),
        "Combined file has no duplicate associations": (
            not combined.duplicated(
                subset=["UniProt_ID", "TCGA_Code", "PMID"]
            ).any()
        ),
        "Combined count is preserved automatic KEEP + Step08 AUTO_KEEP": (
            len(combined)
            == len(auto_keep) + counts.get("AUTO_KEEP", 0)
        ),
    }

    failed = [
        label
        for label, passed in checks.items()
        if not passed
    ]

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 08 Integrity Report",
        "=" * 78,
        "",
        "INPUT",
        "-" * 78,
        f"Step 07 automatic KEEP rows          : {len(auto_keep)}",
        f"Step 07 CONTEXT_REVIEW rows          : {len(context)}",
        f"GPCRs represented in context review  : {context['Official_Gene'].nunique()}",
        "",
        "CONTEXT RESOLUTION",
        "-" * 78,
        f"AUTO_KEEP                            : {counts.get('AUTO_KEEP', 0)}",
        f"MANUAL_REVIEW                        : {counts.get('MANUAL_REVIEW', 0)}",
        f"REJECT context requirement not met   : {counts.get('REJECT', 0)}",
        "",
        "COMBINED VALIDATED SET",
        "-" * 78,
        f"Rows after automatic context resolution: {len(combined)}",
        f"Unique PMIDs                          : {combined['PMID'].nunique()}",
        f"GPCRs with >=1 validated paper        : {len(combined_genes)}",
        f"New GPCRs recovered by Step 08        : {len(newly_recovered_genes)}",
        (
            "Newly recovered GPCRs                  : "
            + (
                ";".join(newly_recovered_genes)
                if newly_recovered_genes
                else "-"
            )
        ),
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
            "CONTEXT_REQUIRED aliases are not accepted solely because the "
            "alias occurs in the article. Automatic acceptance requires "
            "additional receptor-specific evidence defined before record "
            "inspection: multiple distinct context aliases, a direct "
            "alias-receptor naming relationship, a non-generic descriptive "
            "receptor term, or narrow receptor nomenclature such as GPR/CXCR/"
            "CCR/ACKR plus a number. Mixed-case aliases must preserve their "
            "stored capitalization. If sufficient evidence depends only on "
            "an active cross-GPCR collision, the record is held for manual "
            "review. Otherwise, records lacking sufficient context are "
            "rejected for that GPCR assignment."
        ),
        "",
        f"Integrity failures                     : {len(failed)}",
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
            "ALL STEP 08 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return failed

# Run this stage and write its output files.
def main():
    print("=" * 78)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 08")
    print("Systematic CONTEXT_REQUIRED Article Resolution")
    print("=" * 78)

    try:
        context, auto_keep, aliases = load_inputs()
    except Exception as exc:
        print(f"INPUT VALIDATION FAILED:\n{exc}")
        sys.exit(1)

    print(f"Automatic KEEP input : {len(auto_keep):,}")
    print(f"Context-review input : {len(context):,}")
    print(f"Context GPCRs        : {context['Official_Gene'].nunique():,}")
    print("Network requests     : 0")
    print("=" * 78)

    active_owners = build_active_owners(aliases)

    results = context.apply(
        lambda row: pd.Series(
            classify_context_row(row, active_owners)
        ),
        axis=1,
    )

    resolved = pd.concat(
        [
            context.reset_index(drop=True),
            results.reset_index(drop=True),
        ],
        axis=1,
    )

    resolved.to_csv(
        OUT_ALL,
        index=False,
        na_rep="",
    )

    manual = resolved[
        resolved["Context_Resolution"].eq("MANUAL_REVIEW")
    ].copy()

    manual.to_csv(
        OUT_MANUAL,
        index=False,
        na_rep="",
    )

    rejected = resolved[
        resolved["Context_Resolution"].eq("REJECT")
    ].copy()

    rejected.to_csv(
        OUT_REJECTED,
        index=False,
        na_rep="",
    )

    context_kept = resolved[
        resolved["Context_Resolution"].eq("AUTO_KEEP")
    ].copy()

    provenance_cols = [
        "Context_Resolution",
        "Context_Resolution_Reason",
        "Strong_Context_Terms",
        "Context_Evidence_Types",
        "Case_Mismatch_Terms",
        "Active_Collision_Terms",
        "Active_Collision_Genes",
    ]

    auto_keep_combined = auto_keep.copy()
    auto_keep_combined["Context_Resolution"] = "NOT_REQUIRED"
    auto_keep_combined["Context_Resolution_Reason"] = (
        "STEP07_FINAL_KEEP_TERM"
    )

    for col in provenance_cols[2:]:
        auto_keep_combined[col] = ""

    combined = pd.concat(
        [
            auto_keep_combined,
            context_kept,
        ],
        ignore_index=True,
        sort=False,
    )

    combined = combined.drop_duplicates(
        subset=["UniProt_ID", "TCGA_Code", "PMID"],
        keep="first",
    )

    combined.to_csv(
        OUT_COMBINED,
        index=False,
        na_rep="",
    )

    summary = build_summary(resolved)
    summary.to_csv(
        OUT_SUMMARY,
        index=False,
        na_rep="",
    )

    failed = write_report(
        context=context,
        auto_keep=auto_keep,
        resolved=resolved,
        manual=manual,
        rejected=rejected,
        combined=combined,
        aliases=aliases,
        summary=summary,
    )

    counts = (
        resolved["Context_Resolution"]
        .value_counts()
        .to_dict()
    )

    print()
    print("=" * 78)
    print(
        "STEP 08 COMPLETE"
        if not failed
        else "STEP 08 COMPLETE WITH INTEGRITY ISSUES"
    )
    print("=" * 78)
    print(f"AUTO_KEEP     : {counts.get('AUTO_KEEP', 0):,}")
    print(f"MANUAL_REVIEW : {counts.get('MANUAL_REVIEW', 0):,}")
    print(f"REJECT        : {counts.get('REJECT', 0):,}")
    print(f"Combined rows : {len(combined):,}")
    print(
        f"Combined GPCRs: "
        f"{combined['Official_Gene'].nunique():,}"
    )
    print(f"Integrity errs: {len(failed)}")
    print()
    print("Saved:")
    print(f"  {OUT_ALL}")
    print(f"  {OUT_MANUAL}")
    print(f"  {OUT_REJECTED}")
    print(f"  {OUT_COMBINED}")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_REPORT}")
    print("=" * 78)

    if failed:
        sys.exit(2)

if __name__ == "__main__":
    main()
