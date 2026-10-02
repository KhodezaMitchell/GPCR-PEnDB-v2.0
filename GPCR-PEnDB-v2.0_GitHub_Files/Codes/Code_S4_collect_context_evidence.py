#!/usr/bin/env python3
# Purpose: Collect title and abstract context for flagged GPCR names.
# Input: Alias_Audit_Recommendations.csv
# Output: Alias_Context_Review.csv

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import pandas as pd

try:
    from Bio import Entrez
except ImportError:
    Entrez = None

DEFAULT_INPUT = "Alias_Audit_Recommendations.csv"
OUT_CONTEXT = "Alias_Context_Review.csv"
OUT_INTEGRITY = "Alias_Context_Review_Integrity_Report.txt"
CACHE_FILE = "pubmed_context_cache.json"

CONTEXT_OVERLAP_PMIDS = 10
CONTEXT_UNIQUE_EARLY_PMIDS = 10
CONTEXT_UNIQUE_RECENT_PMIDS = 10

EFETCH_BATCH_SIZE = 200
MAX_RETRIES = 5
INITIAL_RETRY_WAIT = 2.0

# Standardize whitespace and missing values in text.
def clean(value):
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Split semicolon-delimited values.
def split_semicolon(value):
    return [clean(x) for x in clean(value).split(";") if clean(x)]

# Sort PubMed IDs by their numeric value.
def numeric_pmids(pmids):
    return sorted(
        {str(p).strip() for p in pmids if str(p).strip()},
        key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x),
    )

# Build a regex for an exact GPCR term.
def exact_term_pattern(term):
    return re.compile(
        rf"(?<![A-Za-z0-9]){re.escape(clean(term))}(?![A-Za-z0-9])",
        re.IGNORECASE,
    )

# Build a regex allowing controlled punctuation variants.
def flexible_term_pattern(term):
    """
    Conservative punctuation-variant matcher.

    Examples treated as equivalent:
      CXCR6   <-> CXCR-6
      GPR5    <-> GPR 5
      5-HT1A  <-> 5HT1A

    Flexibility is allowed only where the alias already contains separator
    punctuation or at letter-number / number-letter boundaries.
    """
    term = clean(term)
    if not term:
        return re.compile(r"(?!x)x")

    pieces = []
    i = 0

    while i < len(term):
        ch = term[i]

        if re.match(r"[\s\-_/\.]", ch):
            while i + 1 < len(term) and re.match(r"[\s\-_/\.]", term[i + 1]):
                i += 1
            pieces.append(r"[\s\-_/\.]*")
            i += 1
            continue

        pieces.append(re.escape(ch))

        if i + 1 < len(term):
            nxt = term[i + 1]
            if ch.isalnum() and nxt.isalnum():
                if (ch.isalpha() and nxt.isdigit()) or (ch.isdigit() and nxt.isalpha()):
                    pieces.append(r"[\s\-_/\.]*")

        i += 1

    body = "".join(pieces)

    return re.compile(
        rf"(?<![A-Za-z0-9]){body}(?![A-Za-z0-9])",
        re.IGNORECASE,
    )

# Choose the matching rule for a candidate alias.
def term_pattern(term):
    return flexible_term_pattern(term)

# Return sentence context for the current record.
def sentence_context(text, term):
    text = clean(text)
    if not text:
        return ""

    pattern = term_pattern(term)
    sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z0-9(])", text)
    matched = [clean(s) for s in sentences if pattern.search(s)]

    if matched:
        return " || ".join(matched)

    match = pattern.search(text)
    if not match:
        return ""

    start = max(0, match.start() - 180)
    end = min(len(text), match.end() + 240)
    return clean(text[start:end])

# Return parenthetical expansions for the current record.
def parenthetical_expansions(text, term):
    text = clean(text)
    if not text:
        return []

    escaped = re.escape(clean(term))
    results = []

    before = re.compile(
        rf"([A-Za-z][A-Za-z0-9α-ωΑ-Ω/,+\- ]{{2,100}}?)\s*"
        rf"\(\s*{escaped}\s*\)",
        re.IGNORECASE,
    )
    after = re.compile(
        rf"(?<![A-Za-z0-9]){escaped}(?![A-Za-z0-9])\s*"
        rf"\(\s*([^)]{{2,100}})\s*\)",
        re.IGNORECASE,
    )

    for match in before.finditer(text):
        candidate = clean(match.group(1))
        candidate = re.split(r"[.;:]", candidate)[-1].strip()
        if candidate:
            results.append(candidate)

    for match in after.finditer(text):
        candidate = clean(match.group(1))
        if candidate:
            results.append(candidate)

    seen = set()
    output = []
    for item in results:
        key = norm(item)
        if key and key not in seen:
            seen.add(key)
            output.append(item)

    return output

REQUIRED_COLUMNS = {
    "UniProt_ID",
    "Original_Gene",
    "Official_Gene",
    "Protein_Name",
    "Alias",
    "Source_Provenance",
    "Term_Types",
    "Preliminary_Status",
    "Preliminary_Reason",
    "Retrieved_PMIDs",
    "Reference_Terms",
    "Overlap_PMID_List",
    "Unique_PMID_List",
}

OPTIONAL_PASS_THROUGH = [
    "HGNC_ID",
    "HGNC_Approved_Symbol",
    "HGNC_Match_Method",
    "From_UniProt",
    "From_HGNC",
    "HGNC_Only",
    "UniProt_Only",
    "Shared_By_UniProt_HGNC",
    "Flag_X_Intrinsic_Ambiguity",
    "Flag_X_Reasons",
    "Flag_Y_Cross_GPCR",
    "Flag_Y_Reasons",
    "Flag_Z_Retrieval_Behavior",
    "Flag_Z_Reasons",
    "Other_Review_Reasons",
    "Cross_GPCR_Collision",
    "Collision_Type",
    "Collision_Genes",
    "Total_PMIDs",
    "Independent_Reference_PMIDs",
    "Overlap",
    "Unique_PMIDs",
    "Overlap_Fraction",
    "Unique_Fraction",
    "Overlap_Metrics_Exact",
]

# Load input from the input records.
def load_input(path):
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: missing input file: {path}")

    df = pd.read_csv(path, dtype=str).fillna("")

    missing = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing:
        raise SystemExit(
            "ERROR: input is missing required column(s): " + ", ".join(missing)
        )

    return df

# Sample overlapping and unique PMIDs for alias review.
def select_context_pmids(overlap_pmids, unique_pmids):
    selected = []
    seen = set()

    # Add a cleaned nonblank value to the collection.
    def add(pmids, label, limit, reverse=False):
        ordered = numeric_pmids(pmids)
        if reverse:
            ordered = list(reversed(ordered))

        for pmid in ordered[:limit]:
            if pmid not in seen:
                seen.add(pmid)
                selected.append((pmid, label))

    add(overlap_pmids, "OVERLAP", CONTEXT_OVERLAP_PMIDS)
    add(unique_pmids, "UNIQUE_EARLY", CONTEXT_UNIQUE_EARLY_PMIDS)
    add(unique_pmids, "UNIQUE_RECENT", CONTEXT_UNIQUE_RECENT_PMIDS, reverse=True)

    return selected

# Set credentials for live PubMed access.
def configure_entrez():
    if Entrez is None:
        raise SystemExit(
            "Biopython is required.\n"
            "Install with: python -m pip install biopython"
        )

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise SystemExit("ERROR: set NCBI_EMAIL before running Step 02")

    Entrez.email = email

    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if api_key:
        Entrez.api_key = api_key

    return 0.11 if api_key else 0.34

# Read previously saved PubMed results.
def load_cache():
    if not os.path.exists(CACHE_FILE):
        return {}

    try:
        with open(CACHE_FILE, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        if isinstance(data, dict):
            print(f"Loaded PubMed cache: {len(data):,} PMIDs")
            return data
    except Exception as exc:
        print(f"WARNING: could not read {CACHE_FILE}: {exc}")

    return {}

# Save a PubMed result for later reuse.
def save_cache(cache):
    temp = CACHE_FILE + ".tmp"
    with open(temp, "w", encoding="utf-8") as fh:
        json.dump(cache, fh, ensure_ascii=False)
    os.replace(temp, CACHE_FILE)

# Return abstract text for the current record.
def abstract_text(article_data):
    parts = article_data.get("Abstract", {}).get("AbstractText", [])
    output = []

    for part in parts:
        label = ""
        try:
            label = clean(part.attributes.get("Label", ""))
        except Exception:
            pass

        text = clean(str(part))
        if label and text:
            output.append(f"{label}: {text}")
        elif text:
            output.append(text)

    return clean(" ".join(output))

# Fetch titles and abstracts absent from the saved cache.
def fetch_missing_metadata(pmids, cache, sleep_seconds):
    missing = [p for p in numeric_pmids(pmids) if p not in cache]

    print(f"Unique sampled PMIDs       : {len(set(pmids)):,}")
    print(f"Already cached             : {len(set(pmids)) - len(missing):,}")
    print(f"Need PubMed metadata       : {len(missing):,}")

    if not missing:
        return cache

    total_batches = (len(missing) + EFETCH_BATCH_SIZE - 1) // EFETCH_BATCH_SIZE

    for batch_num, start in enumerate(
        range(0, len(missing), EFETCH_BATCH_SIZE), start=1
    ):
        batch = missing[start:start + EFETCH_BATCH_SIZE]
        wait = INITIAL_RETRY_WAIT

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(batch),
                    rettype="medline",
                    retmode="xml",
                )
                record = Entrez.read(handle)
                handle.close()
                break
            except Exception as exc:
                if attempt == MAX_RETRIES:
                    raise RuntimeError(
                        f"PubMed EFetch failed after {MAX_RETRIES} attempts "
                        f"for batch beginning PMID {batch[0]}: {exc}"
                    ) from exc

                print(
                    f"  EFetch failed "
                    f"(attempt {attempt}/{MAX_RETRIES}): {exc}"
                )
                print(f"  Retrying in {wait:.1f} sec...")
                time.sleep(wait)
                wait *= 2

        found = set()

        for article in record.get("PubmedArticle", []):
            citation = article.get("MedlineCitation", {})
            pmid = clean(citation.get("PMID", ""))
            data = citation.get("Article", {})

            if not pmid:
                continue

            cache[pmid] = {
                "Title": clean(data.get("ArticleTitle", "")),
                "Abstract": abstract_text(data),
            }
            found.add(pmid)

        for pmid in batch:
            if pmid not in found:
                cache[pmid] = {
                    "Title": "",
                    "Abstract": "",
                }

        save_cache(cache)

        done = min(start + len(batch), len(missing))
        print(
            f"Metadata batch {batch_num}/{total_batches} "
            f"| newly processed {done:,}/{len(missing):,}"
        )

        time.sleep(sleep_seconds)

    return cache

# Build selections from the input records.
def build_selections(flagged):
    selections = []
    all_pmids = set()
    flagged_with_pmids = 0

    for idx, row in flagged.iterrows():
        retrieved = int(float(clean(row["Retrieved_PMIDs"]) or 0))
        if retrieved <= 0:
            continue

        overlap_pmids = split_semicolon(row["Overlap_PMID_List"])
        unique_pmids = split_semicolon(row["Unique_PMID_List"])
        chosen = select_context_pmids(overlap_pmids, unique_pmids)

        if chosen:
            flagged_with_pmids += 1

        for pmid, sample_group in chosen:
            selections.append((idx, pmid, sample_group))
            all_pmids.add(pmid)

    return selections, all_pmids, flagged_with_pmids

# Build context rows from the input records.
def build_context_rows(flagged, selections, cache):
    rows = []

    for row_idx, pmid, sample_group in selections:
        row = flagged.loc[row_idx]

        alias = clean(row["Alias"])
        title = clean(cache.get(pmid, {}).get("Title", ""))
        abstract = clean(cache.get(pmid, {}).get("Abstract", ""))
        combined = clean(f"{title}. {abstract}")

        exact_pattern = exact_term_pattern(alias)
        pattern = flexible_term_pattern(alias)

        reference_terms = split_semicolon(row["Reference_Terms"])
        coanchors = [
            term for term in reference_terms
            if flexible_term_pattern(term).search(combined)
        ]

        out = {
            "UniProt_ID": clean(row["UniProt_ID"]),
            "Original_Gene": clean(row["Original_Gene"]),
            "Official_Gene": clean(row["Official_Gene"]),
            "Protein_Name": clean(row["Protein_Name"]),
            "Alias": alias,
            "Source_Provenance": clean(row["Source_Provenance"]),
            "Term_Types": clean(row["Term_Types"]),
        }

        for col in OPTIONAL_PASS_THROUGH:
            out[col] = clean(row[col]) if col in flagged.columns else ""

        out.update({
            "Preliminary_Status": clean(row["Preliminary_Status"]),
            "Preliminary_Reason": clean(row["Preliminary_Reason"]),
            "PMID": pmid,
            "Sample_Group": sample_group,
            "Title": title,
            "Matched_Context": sentence_context(combined, alias),
            "Exact_Term_Occurrence_Count": len(exact_pattern.findall(combined)),
            "Term_Occurrence_Count": len(pattern.findall(combined)),
            "Context_Match_Mode": (
                "EXACT"
                if exact_pattern.search(combined)
                else ("FLEXIBLE_VARIANT" if pattern.search(combined) else "NO_MATCH")
            ),
            "Parenthetical_Expansions": ";".join(
                parenthetical_expansions(combined, alias)
            ),
            "Cooccurring_Reference_Terms": ";".join(coanchors),
            "Independent_Support_In_Record": bool(coanchors),
            "Metadata_Retrieved": bool(title or abstract),
        })

        rows.append(out)

    return pd.DataFrame(rows)

# Build integrity report from the input records.
def build_integrity_report(
    audit,
    flagged,
    context_df,
    selections,
    flagged_with_pmids,
):
    errors = []

    allowed = {"REVIEW", "CONTEXT_REQUIRED"}

    if not set(flagged["Preliminary_Status"]).issubset(allowed):
        errors.append("Non-flagged status entered Step 02")

    expected_rows = len(selections)
    if len(context_df) != expected_rows:
        errors.append(
            f"Expected {expected_rows} context rows; found {len(context_df)}"
        )

    if not context_df.empty:
        keys = (
            context_df["UniProt_ID"].map(clean)
            + "\t"
            + context_df["Alias"].map(norm)
            + "\t"
            + context_df["PMID"].map(clean)
        )
        if keys.duplicated().any():
            errors.append("Duplicate GPCR-alias-PMID context rows detected")

    expected_pairs = set()
    for idx, pmid, _ in selections:
        row = flagged.loc[idx]
        expected_pairs.add((clean(row["UniProt_ID"]), norm(row["Alias"])))

    observed_pairs = set()
    if not context_df.empty:
        observed_pairs = set(
            zip(
                context_df["UniProt_ID"].map(clean),
                context_df["Alias"].map(norm),
            )
        )

    missing_pairs = expected_pairs - observed_pairs
    if missing_pairs:
        errors.append(
            f"{len(missing_pairs)} flagged GPCR-alias pairs are missing context rows"
        )

    metadata_missing = 0
    if not context_df.empty:
        metadata_missing = int(
            (~context_df["Metadata_Retrieved"].astype(bool)).sum()
        )

    counts = flagged["Preliminary_Status"].value_counts().to_dict()

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 02 Integrity Report",
        "=" * 72,
        f"Step 01 rows                         : {len(audit)}",
        f"Flagged aliases                     : {len(flagged)}",
        f"REVIEW                              : {counts.get('REVIEW', 0)}",
        f"CONTEXT_REQUIRED                    : {counts.get('CONTEXT_REQUIRED', 0)}",
        f"Flagged aliases with sampled PMIDs  : {flagged_with_pmids}",
        f"Context rows                        : {len(context_df)}",
        f"Unique context PMIDs                : {context_df['PMID'].nunique() if not context_df.empty else 0}",
        f"Rows missing PubMed metadata        : {metadata_missing}",
        f"Integrity errors                    : {len(errors)}",
        "",
    ]

    if errors:
        lines.append("ERRORS")
        lines.append("-" * 72)
        lines.extend(f"- {error}" for error in errors)
    else:
        lines.append("ALL STEP 02 INTEGRITY CHECKS PASSED")

    return "\n".join(lines) + "\n", errors

# Run this stage and write its output files.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT,
        help=f"Step 01 recommendation CSV (default: {DEFAULT_INPUT})",
    )
    args = parser.parse_args()

    print("=" * 72)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 02")
    print("=" * 72)

    audit = load_input(args.input)

    flagged = audit[
        audit["Preliminary_Status"].isin(["REVIEW", "CONTEXT_REQUIRED"])
    ].copy()

    print(f"Step 01 rows              : {len(audit):,}")
    print(f"Flagged aliases           : {len(flagged):,}")
    print(
        f"  REVIEW                  : "
        f"{(flagged['Preliminary_Status'] == 'REVIEW').sum():,}"
    )
    print(
        f"  CONTEXT_REQUIRED        : "
        f"{(flagged['Preliminary_Status'] == 'CONTEXT_REQUIRED').sum():,}"
    )

    selections, all_pmids, flagged_with_pmids = build_selections(flagged)

    print(f"Aliases with sampled PMIDs: {flagged_with_pmids:,}")
    print(f"Planned context rows      : {len(selections):,}")

    sleep_seconds = configure_entrez()
    cache = load_cache()
    cache = fetch_missing_metadata(all_pmids, cache, sleep_seconds)

    context_df = build_context_rows(flagged, selections, cache)
    context_df.to_csv(OUT_CONTEXT, index=False)

    report, errors = build_integrity_report(
        audit,
        flagged,
        context_df,
        selections,
        flagged_with_pmids,
    )

    with open(OUT_INTEGRITY, "w", encoding="utf-8") as fh:
        fh.write(report)

    print("\n" + "=" * 72)
    print("STEP 02 COMPLETE")
    print("=" * 72)
    print(report)
    print(f"Context output : {OUT_CONTEXT}")
    print(f"Integrity file : {OUT_INTEGRITY}")
    print(f"Metadata cache : {CACHE_FILE}")
    print("\nCode S4 collects evidence; reviewed alias decisions are applied in Code S6.")

    if errors:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
