#!/usr/bin/env python3
# Purpose: Audit candidate GPCR names using PubMed PMID overlap.
# Input: GPCR_Alias_Table.csv; HGNC_Alias_Table.csv
# Output: Alias_Candidate_Universe.csv; Alias_Audit_Recommendations.csv

# This frozen copy runs only through the staged offline runner.
from _frozen_entry import ROOT as _FROZEN_ROOT

import argparse
import csv
import os
import re
import time
from collections import defaultdict

import pandas as pd

try:
    from Bio import Entrez
except ImportError:
    Entrez = None

UNIPROT_TABLE = "GPCR_Alias_Table.csv"
HGNC_TABLE = "HGNC_Alias_Table.csv"

OUT_CANDIDATES = "Alias_Candidate_Universe.csv"
OUT_AUDIT = "Alias_Audit_Recommendations.csv"
OUT_STATS = "Alias_Audit_Stats.csv"
OUT_COLLISIONS = "Alias_Cross_GPCR_Collisions.csv"
OUT_INTEGRITY = "Alias_Audit_Integrity_Report.txt"

EXPECTED_GPCR_ROWS = 832
PUBMED_UID_LIMIT = 10000
RETMAX = PUBMED_UID_LIMIT
MAX_RETRIES = 5
INITIAL_RETRY_WAIT = 2.0

SHORT_TERM_MAX = 4
ACRONYM_RISK_MAX = 6
HIGH_REFERENCE_PRECISION = 0.80
HIGH_RETRIEVAL_COUNT = 100
HIGH_UNIQUE_FRACTION = 0.80
VERY_HIGH_RETRIEVAL_COUNT = 500
VERY_HIGH_UNIQUE_FRACTION = 0.90

VALID_PRELIMINARY = {"PASS", "REVIEW", "CONTEXT_REQUIRED"}

HGNC_TERM_COLUMNS = {
    "HGNC_Approved_Symbol": "HGNC_APPROVED_SYMBOL",
    "HGNC_Approved_Name": "HGNC_APPROVED_NAME",
    "HGNC_Alias_Symbols": "HGNC_ALIAS_SYMBOL",
    "HGNC_Previous_Symbols": "HGNC_PREVIOUS_SYMBOL",
    "HGNC_Alias_Names": "HGNC_ALIAS_NAME",
    "HGNC_Previous_Names": "HGNC_PREVIOUS_NAME",
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

# Sort unique terms reproducibly.
def unique_sorted(values):
    return sorted(set(values), key=lambda x: (str(x).casefold(), str(x)))

# Return compact term for the current record.
def compact_term(term):
    return re.sub(r"[^\w]", "", clean(term), flags=re.UNICODE)

# Check whether single token applies to the current term or record.
def is_single_token(term):
    return not bool(re.search(r"\s", clean(term)))

# Check whether short term applies to the current term or record.
def is_short_term(term):
    compact = compact_term(term)
    return bool(compact) and is_single_token(term) and len(compact) <= SHORT_TERM_MAX

# Check whether acronym risk applies to the current term or record.
def is_acronym_risk(term):
    compact = compact_term(term)
    return (
        bool(compact)
        and is_single_token(term)
        and len(compact) <= ACRONYM_RISK_MAX
        and compact.upper() == compact
        and any(ch.isalpha() for ch in compact)
    )

GENERIC_NUMBERED_RECEPTOR_PATTERNS = [
    re.compile(
        r"^(?:probable\s+|putative\s+)?"
        r"g[\s-]*protein[\s-]*coupled\s+receptor\s+\d+[a-z]?$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?:probable\s+|putative\s+)?"
        r"g[\s-]*protein[\s-]*coupled\s+receptor\s+[a-z]\d+$",
        re.IGNORECASE,
    ),
]

# Check whether generic numbered receptor applies to the current term or record.
def is_generic_numbered_receptor(term):
    value = clean(term)
    return any(p.fullmatch(value) for p in GENERIC_NUMBERED_RECEPTOR_PATTERNS)

# Sort PubMed IDs by their numeric value.
def numeric_pmids(pmids):
    return sorted(
        {str(p) for p in pmids if str(p)},
        key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x),
    )

# Join unique PMIDs in numeric order.
def join_pmids(pmids):
    return ";".join(numeric_pmids(pmids))

# Validate inputs for the current record.
def validate_inputs(uniprot_df, hgnc_df):
    req_u = {
        "UniProt_ID", "Original_Gene", "Official_Gene",
        "Protein_Name", "Search_Terms",
    }
    req_h = {
        "UniProt_ID", "Original_Gene", "Step00_Official_Gene",
        "HGNC_ID", "HGNC_Approved_Symbol", "HGNC_Approved_Name",
        "HGNC_Alias_Symbols", "HGNC_Previous_Symbols",
        "HGNC_Alias_Names", "HGNC_Previous_Names", "HGNC_Match_Method",
    }

    missing_u = sorted(req_u - set(uniprot_df.columns))
    missing_h = sorted(req_h - set(hgnc_df.columns))
    if missing_u:
        raise ValueError(f"{UNIPROT_TABLE} missing: {', '.join(missing_u)}")
    if missing_h:
        raise ValueError(f"{HGNC_TABLE} missing: {', '.join(missing_h)}")

    for name, df in [(UNIPROT_TABLE, uniprot_df), (HGNC_TABLE, hgnc_df)]:
        if df["UniProt_ID"].map(clean).duplicated().any():
            raise ValueError(f"{name} contains duplicate UniProt_ID rows")

    if set(uniprot_df["UniProt_ID"].map(clean)) != set(hgnc_df["UniProt_ID"].map(clean)):
        raise ValueError("UniProt and HGNC tables do not contain the same UniProt_ID set")

# Return add candidate for the current record.
def add_candidate(store, term, source, term_type):
    term = clean(term)
    if not term:
        return
    key = norm(term)
    if key not in store:
        store[key] = {
            "display": term,
            "sources": set(),
            "types": set(),
        }
    store[key]["sources"].add(source)
    store[key]["types"].add(term_type)

# Build candidate universe from the input records.
def build_candidate_universe(uniprot_df, hgnc_df):
    hgnc_by_acc = hgnc_df.set_index("UniProt_ID", drop=False)
    rows = []

    for _, u in uniprot_df.iterrows():
        acc = clean(u["UniProt_ID"])
        h = hgnc_by_acc.loc[acc]
        store = {}

        official = clean(u["Official_Gene"])
        original = clean(u["Original_Gene"])
        protein = clean(u["Protein_Name"])

        for term in split_terms(u["Search_Terms"]):
            if norm(term) == norm(official):
                term_type = "UNIPROT_OFFICIAL_GENE"
            elif original and norm(term) == norm(original):
                term_type = "UNIPROT_ORIGINAL_GENE"
            elif protein and norm(term) == norm(protein):
                term_type = "UNIPROT_PROTEIN_NAME"
            else:
                term_type = "UNIPROT_OTHER_TERM"
            add_candidate(store, term, "UniProt", term_type)

        for column, term_type in HGNC_TERM_COLUMNS.items():
            for term in split_terms(h[column]):
                add_candidate(store, term, "HGNC", term_type)

        symbol_disagreement = (
            bool(clean(h["HGNC_Approved_Symbol"]))
            and norm(h["HGNC_Approved_Symbol"]) != norm(official)
        )

        for item in store.values():
            sources = unique_sorted(item["sources"])
            term_types = unique_sorted(item["types"])
            rows.append({
                "UniProt_ID": acc,
                "Original_Gene": original,
                "Official_Gene": official,
                "Protein_Name": protein,
                "HGNC_ID": clean(h["HGNC_ID"]),
                "HGNC_Approved_Symbol": clean(h["HGNC_Approved_Symbol"]),
                "HGNC_Match_Method": clean(h["HGNC_Match_Method"]),
                "Alias": item["display"],
                "Source_Provenance": ";".join(sources),
                "Term_Types": ";".join(term_types),
                "From_UniProt": "UniProt" in item["sources"],
                "From_HGNC": "HGNC" in item["sources"],
                "HGNC_Only": item["sources"] == {"HGNC"},
                "UniProt_Only": item["sources"] == {"UniProt"},
                "Shared_By_UniProt_HGNC": item["sources"] == {"UniProt", "HGNC"},
                "Official_Symbol_Disagreement": symbol_disagreement,
            })

    candidates = pd.DataFrame(rows)
    candidates.sort_values(
        ["Official_Gene", "Alias"],
        key=lambda s: s.astype(str).str.casefold(),
        inplace=True,
    )
    candidates.reset_index(drop=True, inplace=True)
    return candidates

# Build ownership maps from the input records.
def build_ownership_maps(candidates):
    official_owners = defaultdict(set)
    alias_owners = defaultdict(set)

    for _, row in candidates.iterrows():
        gene = clean(row["Official_Gene"])
        official_owners[norm(gene)].add(gene)
        alias_owners[norm(row["Alias"])].add(gene)

    return official_owners, alias_owners

# Return collision info for the current record.
def collision_info(alias, current_gene, official_owners, alias_owners):
    alias_key = norm(alias)
    gene_key = norm(current_gene)

    other_official = unique_sorted(
        gene for gene in official_owners.get(alias_key, set())
        if norm(gene) != gene_key
    )
    other_alias = unique_sorted(
        gene for gene in alias_owners.get(alias_key, set())
        if norm(gene) != gene_key
    )

    types = []
    if other_official:
        types.append("ALIAS_IS_OTHER_GPCR_OFFICIAL_GENE")
    if other_alias:
        types.append("ALIAS_USED_BY_MULTIPLE_GPCRS")

    return {
        "has_collision": bool(types),
        "collision_type": ";".join(types),
        "collision_genes": ";".join(unique_sorted(set(other_official) | set(other_alias))),
    }

# Return anchor exclusion reason for the current record.
def anchor_exclusion_reason(term, gene, original, collision):
    term_key = norm(term)
    if not term_key:
        return "blank"
    if term_key in {norm(gene), norm(original)}:
        return "official_or_original_symbol"
    if is_short_term(term):
        return "short_term"
    if is_acronym_risk(term):
        return "acronym_risk"
    if is_generic_numbered_receptor(term):
        return "generic_numbered_receptor"
    if collision["has_collision"]:
        return "cross_gpcr_collision"
    if is_single_token(term) and len(compact_term(term)) < 8:
        return "short_single_token"
    return ""

# Build anchor terms from the input records.
def build_anchor_terms(uniprot_df, official_owners, alias_owners):
    """
    Use descriptive UniProt-derived terms as reference anchors; HGNC-only
    candidates do not become anchors.
    """
    anchors_by_acc = {}
    excluded_by_acc = {}

    for _, row in uniprot_df.iterrows():
        acc = clean(row["UniProt_ID"])
        gene = clean(row["Official_Gene"])
        original = clean(row["Original_Gene"])

        ordered = []
        if clean(row["Protein_Name"]):
            ordered.append(clean(row["Protein_Name"]))
        ordered.extend(split_terms(row["Search_Terms"]))

        accepted = []
        excluded = []
        seen = set()

        for term in ordered:
            key = norm(term)
            if not key or key in seen:
                continue
            seen.add(key)
            collision = collision_info(term, gene, official_owners, alias_owners)
            reason = anchor_exclusion_reason(term, gene, original, collision)
            if reason:
                excluded.append(f"{term}:{reason}")
            else:
                accepted.append(term)

        anchors_by_acc[acc] = accepted
        excluded_by_acc[acc] = excluded

    return anchors_by_acc, excluded_by_acc

# Set credentials for live PubMed access.
def configure_entrez():
    if Entrez is None:
        raise SystemExit("Biopython is required: python -m pip install biopython")

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise SystemExit(
            "Set NCBI_EMAIL before the live run, for example:\n"
            "  export NCBI_EMAIL=\"your_email@example.edu\""
        )
    Entrez.email = email

    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if api_key:
        Entrez.api_key = api_key
    return api_key

# Construct the PubMed search query for the current item.
def build_query(term):
    return (
        f'(\"{clean(term)}\"[Title/Abstract]) AND '
        f'(humans[MeSH Terms] OR human[Title/Abstract] OR patients[Title/Abstract]) '
        f'NOT Review[Publication Type]'
    )

# Search PubMed and retry transient failures.
def esearch_with_retry(query):
    wait = INITIAL_RETRY_WAIT
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=RETMAX, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            return record
        except Exception as exc:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"PubMed ESearch failed after {MAX_RETRIES} attempts.\n"
                    f"Query: {query}\nLast error: {exc}"
                ) from exc
            print(f"  PubMed request failed ({attempt}/{MAX_RETRIES}): {exc}")
            print(f"  Retrying in {wait:.1f} seconds...")
            time.sleep(wait)
            wait *= 2

# Load saved alias-search results for offline replay.
def load_search_snapshot(path):
    required = {
        "Normalized_Term", "PubMed_Query", "Total_PMIDs",
        "Retrieved_PMIDs", "PMID_List_Truncated", "Alias_PMID_List",
    }
    snapshot = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Frozen alias search snapshot has missing columns")
        for row in reader:
            query = row["PubMed_Query"]
            pmids = [x for x in row["Alias_PMID_List"].split(";") if x]
            count = int(row["Total_PMIDs"])
            retrieved = int(row["Retrieved_PMIDs"])
            if not query or query in snapshot or len(pmids) != retrieved:
                raise ValueError("Invalid or repeated frozen alias search record")
            if (count > retrieved) != (row["PMID_List_Truncated"] == "True"):
                raise ValueError("Frozen alias search truncation flag disagrees with counts")
            snapshot[query] = {"Count": count, "IdList": pmids}
    return snapshot

# Collect PMIDs for one alias from the snapshot or PubMed.
def search_pmids(term, run_cache, sleep_seconds, snapshot=None):
    key = norm(term)
    if key in run_cache:
        return run_cache[key]

    query = build_query(term)
    if snapshot is None:
        record = esearch_with_retry(query)
    else:
        if query not in snapshot:
            raise ValueError(f"Historical alias query missing from snapshot: {query}")
        record = snapshot[query]
    pmids = set(record.get("IdList", []))
    count = int(record.get("Count", len(pmids)))

    result = {
        "term": clean(term),
        "query": query,
        "count": count,
        "pmids": pmids,
        "retrieved": len(pmids),
        "truncated": count > len(pmids),
    }
    run_cache[key] = result
    time.sleep(sleep_seconds)
    return result

# Combine PMIDs from the other receptor-specific names.
def build_independent_reference(alias, anchor_terms, run_cache):
    alias_key = norm(alias)
    terms = [x for x in anchor_terms if norm(x) != alias_key]
    pmids = set()
    truncated = False
    details = []

    for term in terms:
        result = run_cache[norm(term)]
        pmids |= result["pmids"]
        truncated = truncated or result["truncated"]
        detail = f"{term}:{result['count']}"
        if result["truncated"]:
            detail += "[TRUNCATED]"
        details.append(detail)

    return {
        "terms": terms,
        "pmids": pmids,
        "truncated": truncated,
        "details": details,
    }

# Route an alias using overlap and ambiguity evidence.
def preliminary_assessment(alias, protein_name, alias_result, reference, collision, symbol_disagreement):
    """
    Return systematic risk flags and a preliminary routing recommendation.
    This is NOT the final curated KEEP/IGNORE decision.
    """
    total = alias_result["count"]
    alias_pmids = alias_result["pmids"]
    ref_pmids = reference["pmids"]
    overlap = len(alias_pmids & ref_pmids)
    unique = len(alias_pmids - ref_pmids)

    exact = not alias_result["truncated"] and not reference["truncated"]
    overlap_fraction = (overlap / total) if exact and total else 0.0 if exact else None
    unique_fraction = (unique / total) if exact and total else 0.0 if exact else None

    x_reasons = []
    y_reasons = []
    z_reasons = []
    other_reasons = []

    if len(compact_term(alias)) <= 2:
        x_reasons.append("VERY_SHORT_TERM")
    if is_short_term(alias):
        x_reasons.append("SHORT_TERM")
    if is_acronym_risk(alias):
        x_reasons.append("ACRONYM_RISK")
    if is_generic_numbered_receptor(alias):
        x_reasons.append("GENERIC_NUMBERED_RECEPTOR")

    if collision["has_collision"]:
        y_reasons.extend(collision["collision_type"].split(";"))

    if alias_result["truncated"]:
        z_reasons.append("ALIAS_PMID_LIST_TRUNCATED")
    if reference["truncated"]:
        z_reasons.append("REFERENCE_PMID_LIST_TRUNCATED")

    if total > 0 and not ref_pmids:
        z_reasons.append("NO_INDEPENDENT_REFERENCE")
    elif total > 0 and exact and overlap == 0:
        z_reasons.append("ZERO_REFERENCE_OVERLAP")

    if exact and total >= VERY_HIGH_RETRIEVAL_COUNT and unique_fraction is not None and unique_fraction >= VERY_HIGH_UNIQUE_FRACTION:
        z_reasons.append("VERY_HIGH_UNIQUE_RETRIEVAL")
    elif exact and total >= HIGH_RETRIEVAL_COUNT and unique_fraction is not None and unique_fraction >= HIGH_UNIQUE_FRACTION:
        z_reasons.append("HIGH_UNIQUE_RETRIEVAL")

    if (
        total > 0
        and not ref_pmids
        and is_single_token(alias)
        and norm(alias) == norm(protein_name)
    ):
        z_reasons.append("SINGLE_WORD_PROTEIN_NAME_WITHOUT_REFERENCE")

    if symbol_disagreement:
        other_reasons.append("UNIPROT_HGNC_APPROVED_SYMBOL_DISAGREEMENT")

    if other_reasons or y_reasons:
        status = "REVIEW"
        route_reason = "source_or_collision_review"
    elif is_generic_numbered_receptor(alias) or len(compact_term(alias)) <= 2:
        status = "REVIEW"
        route_reason = "intrinsic_term_ambiguity"
    elif alias_result["truncated"] or reference["truncated"]:
        status = "REVIEW"
        route_reason = "pubmed_truncation"
    elif total == 0:
        status = "PASS"
        route_reason = "no_pubmed_hits"
    elif is_short_term(alias) or is_acronym_risk(alias):
        if not ref_pmids:
            status = "REVIEW"
            route_reason = "short_or_acronym_without_independent_reference"
        elif overlap == 0:
            status = "REVIEW"
            route_reason = "short_or_acronym_without_reference_overlap"
        elif overlap_fraction is not None and overlap_fraction >= HIGH_REFERENCE_PRECISION:
            status = "PASS"
            route_reason = "short_or_acronym_high_reference_precision"
        else:
            status = "CONTEXT_REQUIRED"
            route_reason = "short_or_acronym_mixed_retrieval"
    elif not ref_pmids:
        status = "REVIEW"
        route_reason = "no_independent_reference"
    elif overlap == 0:
        status = "REVIEW"
        route_reason = "zero_reference_overlap"
    elif exact and total >= VERY_HIGH_RETRIEVAL_COUNT and unique_fraction >= VERY_HIGH_UNIQUE_FRACTION:
        status = "CONTEXT_REQUIRED"
        route_reason = "very_high_unique_retrieval"
    elif exact and total >= HIGH_RETRIEVAL_COUNT and unique_fraction >= HIGH_UNIQUE_FRACTION:
        status = "CONTEXT_REQUIRED"
        route_reason = "high_unique_retrieval"
    else:
        status = "PASS"
        route_reason = "acceptable_independent_reference_behavior"

    return {
        "status": status,
        "route_reason": route_reason,
        "x_reasons": unique_sorted(x_reasons),
        "y_reasons": unique_sorted(y_reasons),
        "z_reasons": unique_sorted(z_reasons),
        "other_reasons": unique_sorted(other_reasons),
        "overlap": overlap,
        "unique": unique,
        "overlap_fraction": overlap_fraction,
        "unique_fraction": unique_fraction,
        "exact": exact,
    }

# Build integrity report from the input records.
def build_integrity_report(uniprot_df, hgnc_df, candidates, audit):
    errors = []

    if len(uniprot_df) != EXPECTED_GPCR_ROWS:
        errors.append(f"Expected {EXPECTED_GPCR_ROWS} UniProt rows; found {len(uniprot_df)}")
    if len(hgnc_df) != EXPECTED_GPCR_ROWS:
        errors.append(f"Expected {EXPECTED_GPCR_ROWS} HGNC rows; found {len(hgnc_df)}")
    if candidates["UniProt_ID"].nunique() != EXPECTED_GPCR_ROWS:
        errors.append("Candidate universe does not contain 832 distinct GPCRs")
    if len(audit) != len(candidates):
        errors.append("Audit row count does not equal candidate row count")

    dup = audit["UniProt_ID"].map(norm) + "\t" + audit["Alias"].map(norm)
    if dup.duplicated().any():
        errors.append("Duplicate normalized UniProt_ID/Alias rows detected")

    invalid = sorted(set(audit["Preliminary_Status"]) - VALID_PRELIMINARY)
    if invalid:
        errors.append("Invalid preliminary statuses: " + ", ".join(invalid))

    for _, row in audit.iterrows():
        refs = {norm(x) for x in split_terms(row["Reference_Terms"])}
        if norm(row["Alias"]) in refs:
            errors.append(f"Self-reference detected: {row['Official_Gene']} / {row['Alias']}")
            break

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 01 Integrity Report",
        "=" * 72,
        f"GPCR rows                         : {len(uniprot_df)}",
        f"Candidate gene-term rows          : {len(candidates)}",
        f"HGNC-only candidate rows          : {int(candidates['HGNC_Only'].sum())}",
        f"UniProt-only candidate rows       : {int(candidates['UniProt_Only'].sum())}",
        f"Shared candidate rows             : {int(candidates['Shared_By_UniProt_HGNC'].sum())}",
        f"Cross-GPCR collision rows         : {int(audit['Cross_GPCR_Collision'].sum())}",
        f"PASS                              : {(audit['Preliminary_Status'] == 'PASS').sum()}",
        f"CONTEXT_REQUIRED                  : {(audit['Preliminary_Status'] == 'CONTEXT_REQUIRED').sum()}",
        f"REVIEW                            : {(audit['Preliminary_Status'] == 'REVIEW').sum()}",
        f"Integrity errors                  : {len(errors)}",
        "",
    ]
    if errors:
        lines.append("ERRORS")
        lines.append("-" * 72)
        lines.extend(f"- {x}" for x in errors)
    else:
        lines.append("ALL STEP 01 INTEGRITY CHECKS PASSED")
    return "\n".join(lines) + "\n", errors

# Check representative alias-routing rules.
def run_self_tests():
    assert is_short_term("ABC")
    assert is_acronym_risk("ABCDEF")
    assert not is_acronym_risk("mAbc1")
    assert is_generic_numbered_receptor("G protein-coupled receptor 5")
    assert not is_generic_numbered_receptor("C-C chemokine receptor type 7")

    no_collision = {"has_collision": False, "collision_type": "", "collision_genes": ""}
    alias_result = {"count": 10, "pmids": {str(x) for x in range(1, 11)}, "retrieved": 10, "truncated": False}
    reference = {"pmids": {str(x) for x in range(1, 10)}, "truncated": False}
    result = preliminary_assessment("ABCD", "Descriptive receptor", alias_result, reference, no_collision, False)
    assert result["status"] == "PASS"

    no_overlap = {"pmids": {"100", "101"}, "truncated": False}
    result = preliminary_assessment("Long receptor alias", "Descriptive receptor", alias_result, no_overlap, no_collision, False)
    assert result["status"] == "REVIEW"

    print("ALL STEP 01 SELF-TESTS PASSED")

# Run this stage and write its output files.
def main(snapshot_path=None):
    if not os.path.exists(UNIPROT_TABLE):
        raise SystemExit(f"Missing input: {UNIPROT_TABLE}")
    if not os.path.exists(HGNC_TABLE):
        raise SystemExit(f"Missing input: {HGNC_TABLE}")

    uniprot_df = pd.read_csv(UNIPROT_TABLE, dtype=str).fillna("")
    hgnc_df = pd.read_csv(HGNC_TABLE, dtype=str).fillna("")
    validate_inputs(uniprot_df, hgnc_df)

    candidates = build_candidate_universe(uniprot_df, hgnc_df)
    official_owners, alias_owners = build_ownership_maps(candidates)
    anchors_by_acc, excluded_anchors_by_acc = build_anchor_terms(
        uniprot_df, official_owners, alias_owners
    )

    collision_cache = {}
    for _, row in candidates.iterrows():
        key = (clean(row["UniProt_ID"]), norm(row["Alias"]))
        collision_cache[key] = collision_info(
            row["Alias"], row["Official_Gene"], official_owners, alias_owners
        )

    candidates.to_csv(OUT_CANDIDATES, index=False)

    snapshot = load_search_snapshot(snapshot_path) if snapshot_path else None
    api_key = configure_entrez() if snapshot is None else ""
    sleep_seconds = 0 if snapshot is not None else (0.11 if api_key else 0.34)

    print("=" * 76)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 01")
    print("=" * 76)
    print(f"GPCR rows                 : {len(uniprot_df)}")
    print(f"Candidate gene-term rows  : {len(candidates)}")
    print(f"Unique normalized terms   : {candidates['Alias'].map(norm).nunique()}")
    print(f"HGNC-only rows            : {int(candidates['HGNC_Only'].sum())}")
    print("S3 assigns preliminary routes; Code S6 applies the final alias decisions")
    print(f"NCBI API key              : {'yes' if api_key else 'no'}")
    print("=" * 76)

    run_cache = {}
    start = time.time()

    unique_terms = []
    seen = set()
    for term in candidates["Alias"]:
        key = norm(term)
        if key not in seen:
            seen.add(key)
            unique_terms.append(clean(term))

    for i, term in enumerate(unique_terms, start=1):
        search_pmids(term, run_cache, sleep_seconds, snapshot)
        if i == 1 or i % 250 == 0 or i == len(unique_terms):
            elapsed = (time.time() - start) / 60.0
            print(f"Searched {i}/{len(unique_terms)} unique terms | {elapsed:.1f} min")

    if snapshot is not None and len(snapshot) != len(run_cache):
        raise ValueError("Frozen alias search snapshot contains unexpected queries")

    rows = []
    total = len(candidates)

    for i, row in candidates.iterrows():
        acc = clean(row["UniProt_ID"])
        alias = clean(row["Alias"])
        gene = clean(row["Official_Gene"])
        collision = collision_cache[(acc, norm(alias))]
        alias_result = run_cache[norm(alias)]
        reference = build_independent_reference(alias, anchors_by_acc[acc], run_cache)

        assessment = preliminary_assessment(
            alias=alias,
            protein_name=clean(row["Protein_Name"]),
            alias_result=alias_result,
            reference=reference,
            collision=collision,
            symbol_disagreement=bool(row["Official_Symbol_Disagreement"]),
        )

        alias_pmids = alias_result["pmids"]
        reference_pmids = reference["pmids"]
        overlap_pmids = alias_pmids & reference_pmids
        unique_pmids = alias_pmids - reference_pmids

        out = row.to_dict()
        out.update({
            "Is_Official_Gene": norm(alias) == norm(gene),
            "Is_Original_Gene": bool(clean(row["Original_Gene"])) and norm(alias) == norm(row["Original_Gene"]),
            "Is_Protein_Name": bool(clean(row["Protein_Name"])) and norm(alias) == norm(row["Protein_Name"]),
            "Is_Short_Term": is_short_term(alias),
            "Is_Acronym_Risk": is_acronym_risk(alias),
            "Is_Generic_Numbered_Receptor": is_generic_numbered_receptor(alias),
            "Cross_GPCR_Collision": collision["has_collision"],
            "Collision_Type": collision["collision_type"],
            "Collision_Genes": collision["collision_genes"],
            "Flag_X_Intrinsic_Ambiguity": bool(assessment["x_reasons"]),
            "Flag_X_Reasons": ";".join(assessment["x_reasons"]),
            "Flag_Y_Cross_GPCR": bool(assessment["y_reasons"]),
            "Flag_Y_Reasons": ";".join(assessment["y_reasons"]),
            "Flag_Z_Retrieval_Behavior": bool(assessment["z_reasons"]),
            "Flag_Z_Reasons": ";".join(assessment["z_reasons"]),
            "Other_Review_Reasons": ";".join(assessment["other_reasons"]),
            "Total_PMIDs": alias_result["count"],
            "Retrieved_PMIDs": alias_result["retrieved"],
            "PMID_List_Truncated": alias_result["truncated"],
            "Independent_Reference_PMIDs": len(reference_pmids),
            "Reference_PMID_List_Truncated": reference["truncated"],
            "Reference_Terms": ";".join(reference["terms"]),
            "Reference_Search_Counts": ";".join(reference["details"]),
            "Excluded_Anchor_Terms": ";".join(excluded_anchors_by_acc[acc]),
            "Overlap": assessment["overlap"],
            "Unique_PMIDs": assessment["unique"],
            "Overlap_Fraction": round(assessment["overlap_fraction"], 4) if assessment["overlap_fraction"] is not None else "",
            "Unique_Fraction": round(assessment["unique_fraction"], 4) if assessment["unique_fraction"] is not None else "",
            "Overlap_Metrics_Exact": assessment["exact"],
            "Preliminary_Status": assessment["status"],
            "Preliminary_Reason": assessment["route_reason"],
            "PubMed_Query": alias_result["query"],
            "Reference_PMID_List": join_pmids(reference_pmids),
            "Alias_PMID_List": join_pmids(alias_pmids),
            "Overlap_PMID_List": join_pmids(overlap_pmids),
            "Unique_PMID_List": join_pmids(unique_pmids),
        })
        rows.append(out)

        if (i + 1) % 1000 == 0 or i + 1 == total:
            print(f"Audited {i + 1}/{total} candidate rows")

    audit = pd.DataFrame(rows)
    audit.sort_values(
        ["Official_Gene", "Alias"],
        key=lambda s: s.astype(str).str.casefold(),
        inplace=True,
    )
    audit.to_csv(OUT_AUDIT, index=False)

    stats = (
        audit.groupby(["Preliminary_Status", "Preliminary_Reason"])
        .size()
        .reset_index(name="Count")
        .sort_values(["Preliminary_Status", "Count", "Preliminary_Reason"], ascending=[True, False, True])
    )
    stats.to_csv(OUT_STATS, index=False)
    audit[audit["Cross_GPCR_Collision"]].to_csv(OUT_COLLISIONS, index=False)

    integrity_text, errors = build_integrity_report(uniprot_df, hgnc_df, candidates, audit)
    with open(OUT_INTEGRITY, "w", encoding="utf-8") as fh:
        fh.write(integrity_text)

    counts = audit["Preliminary_Status"].value_counts().to_dict()
    print("\n" + "=" * 76)
    print("STEP 01 COMPLETE")
    print("=" * 76)
    print(f"Candidate rows            : {len(audit)}")
    print(f"PASS                      : {counts.get('PASS', 0)}")
    print(f"CONTEXT_REQUIRED          : {counts.get('CONTEXT_REQUIRED', 0)}")
    print(f"REVIEW                    : {counts.get('REVIEW', 0)}")
    print(f"Cross-GPCR collisions     : {int(audit['Cross_GPCR_Collision'].sum())}")
    print(f"Integrity errors          : {len(errors)}")
    print(f"Elapsed                   : {(time.time() - start) / 60.0:.1f} min")
    print("\nCode S3 provides preliminary routes; reviewed decisions are applied in Code S6.")
    print("Next step: retrieve context for REVIEW and CONTEXT_REQUIRED terms, then curate.")

    if errors:
        raise SystemExit(f"Step 01 produced {len(errors)} integrity error(s). See {OUT_INTEGRITY}.")

# Read command-line options for this stage.
def parse_args():
    p = argparse.ArgumentParser(description="Audit combined UniProt + HGNC GPCR aliases")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--snapshot", metavar="CSV", help="Replay the saved historical PubMed search results offline")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    if args.self_test:
        run_self_tests()
    else:
        main(args.snapshot)
