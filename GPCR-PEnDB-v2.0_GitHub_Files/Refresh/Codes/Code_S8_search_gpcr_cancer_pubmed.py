#!/usr/bin/env python3
# Purpose: Search PubMed for every GPCR and cancer combination.
# Input: GPCR_Alias_Table_FINAL_VALIDATED.csv; TCGA_Cancer_Vocabulary.csv
# Output: GPCR_Cancer_PubMed_PMIDs_Long.csv; GPCR_Cancer_PubMed_Search_Summary.csv

import argparse
import hashlib
import json
import os
import re
import socket
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone

import pandas as pd
from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "").strip()

API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
if API_KEY:
    Entrez.api_key = API_KEY

ALIAS_FILE = "GPCR_Alias_Table_FINAL_VALIDATED.csv"
CANCER_FILE = "TCGA_Cancer_Vocabulary.csv"

EXPECTED_GPCR_COUNT = 832
EXPECTED_CANCER_COUNT = 33
EXPECTED_CANCER_TERM_COUNT = 122
EXPECTED_SEARCH_COUNT = 27456

MAX_GPCR_TERMS = None

CACHE_VERSION = "2026-08-21_final_fresh_aliases"
CACHE_DIR = "pubmed_cancer_search_cache_final"

OUT_SUMMARY = "GPCR_Cancer_PubMed_Search_Summary.csv"
OUT_LONG = "GPCR_Cancer_PubMed_PMIDs_Long.csv"
OUT_ERRORS = "GPCR_Cancer_PubMed_Search_Errors.csv"
OUT_REPORT = "GPCR_Cancer_PubMed_Search_Integrity_Report.txt"
OUT_PARTIAL = "GPCR_Cancer_PubMed_Search_PARTIAL.csv"

CHECKPOINT_EVERY = 100

ESEARCH_RETMAX = 10000

MIN_PUBLICATION_YEAR = 1600
MAX_PUBLICATION_YEAR = datetime.now().year + 1

MAX_RETRIES = 5
INITIAL_RETRY_WAIT = 2.0

SOCKET_TIMEOUT_SECONDS = int(
    os.environ.get("NCBI_SOCKET_TIMEOUT", "60")
)
socket.setdefaulttimeout(SOCKET_TIMEOUT_SECONDS)

REQUEST_SLEEP = 0.34 if not API_KEY else 0.11

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

# Quote a term for an exact PubMed search.
def quote_pubmed_term(term):
    """
    Escape characters that could break a quoted PubMed phrase.
    """
    term = clean(term)
    term = term.replace("\\", "\\\\")
    term = term.replace('"', '\\"')
    return f'"{term}"[Title/Abstract]'

# Return PubMed IDs in numeric order.
def numeric_pmid_sort(pmids):
    # Return a stable sort key for an item.
    def key(x):
        try:
            return int(x)
        except Exception:
            return x

    return sorted(set(pmids), key=key)

# Join unique PMIDs in numeric order.
def join_pmids(pmids):
    return ";".join(numeric_pmid_sort(pmids))

# Hash a PubMed query for its cache filename.
def query_hash(query):
    return hashlib.sha256(query.encode("utf-8")).hexdigest()

# Read command-line options for this stage.
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Search PubMed for every validated GPCR x TCGA cancer pair."
        )
    )
    parser.add_argument(
        "--preflight",
        action="store_true",
        help=(
            "validate the frozen inputs and construct every query without "
            "contacting PubMed"
        ),
    )
    return parser.parse_args()

# Load and validate inputs from the input records.
def load_and_validate_inputs():
    if not os.path.exists(ALIAS_FILE):
        raise FileNotFoundError(
            f"Missing {ALIAS_FILE}. "
            "Run Step 04 finalization before Step 06."
        )

    if not os.path.exists(CANCER_FILE):
        raise FileNotFoundError(
            f"Missing {CANCER_FILE}. "
            "Run Step 05 cancer-vocabulary construction before Step 06."
        )

    gpcr = pd.read_csv(ALIAS_FILE)
    cancer = pd.read_csv(CANCER_FILE)

    require_columns(
        gpcr,
        [
            "UniProt_ID",
            "Original_Gene",
            "Official_Gene",
            "Protein_Name",
            "Search_Terms",
            "Context_Required_Terms",
            "Active_Retrieval_Terms",
            "Ignored_Aliases",
            "Not_Validated_Aliases",
        ],
        ALIAS_FILE,
    )

    require_columns(
        cancer,
        [
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PubMed_Search_Terms",
        ],
        CANCER_FILE,
    )

    if len(gpcr) != EXPECTED_GPCR_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_GPCR_COUNT} GPCR rows; found {len(gpcr)}."
        )

    if gpcr["UniProt_ID"].astype(str).map(norm).duplicated().any():
        raise ValueError("Duplicate UniProt_ID rows found in final alias table.")

    if gpcr["Official_Gene"].astype(str).map(norm).duplicated().any():
        dup = gpcr[
            gpcr["Official_Gene"].astype(str).map(norm).duplicated(keep=False)
        ][["UniProt_ID", "Official_Gene"]]
        raise ValueError(
            "Duplicate official gene symbols found:\n"
            + dup.to_string(index=False)
        )

    empty_keep = []
    empty_active = []
    missing_official = []

    keep_term_count = 0
    context_term_count = 0
    active_term_count = 0
    ignored_term_count = 0
    not_validated_term_count = 0

    active_owners = defaultdict(list)

    for _, row in gpcr.iterrows():
        gene = clean(row["Official_Gene"])

        keep_terms = split_terms(row["Search_Terms"])
        context_terms = split_terms(row["Context_Required_Terms"])
        active_terms = split_terms(row["Active_Retrieval_Terms"])
        ignored_terms = split_terms(row["Ignored_Aliases"])
        not_validated_terms = split_terms(row["Not_Validated_Aliases"])

        keep_keys = {norm(x) for x in keep_terms}
        context_keys = {norm(x) for x in context_terms}
        active_keys = {norm(x) for x in active_terms}
        ignored_keys = {norm(x) for x in ignored_terms}
        not_validated_keys = {norm(x) for x in not_validated_terms}

        keep_term_count += len(keep_terms)
        context_term_count += len(context_terms)
        active_term_count += len(active_terms)
        ignored_term_count += len(ignored_terms)
        not_validated_term_count += len(not_validated_terms)

        if not keep_terms:
            empty_keep.append(gene)

        if not active_terms:
            empty_active.append(gene)

        expected_active = keep_keys | context_keys
        if active_keys != expected_active:
            raise ValueError(
                f"{gene}: Active_Retrieval_Terms does not equal "
                "Search_Terms + Context_Required_Terms."
            )

        categories = {
            "KEEP": keep_keys,
            "CONTEXT_REQUIRED": context_keys,
            "IGNORE": ignored_keys,
            "NOT_VALIDATED": not_validated_keys,
        }

        labels = list(categories)
        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                overlap = categories[a] & categories[b]
                if overlap:
                    raise ValueError(
                        f"{gene}: term(s) occur in both {a} and {b}: "
                        f"{sorted(overlap)}"
                    )

        if norm(gene) not in keep_keys:
            missing_official.append(gene)

        for term in keep_terms:
            active_owners[norm(term)].append((gene, "KEEP"))

        for term in context_terms:
            active_owners[norm(term)].append(
                (gene, "CONTEXT_REQUIRED")
            )

    if empty_keep:
        raise ValueError(
            "GPCRs with empty Search_Terms: " + "; ".join(empty_keep)
        )

    if empty_active:
        raise ValueError(
            "GPCRs with empty Active_Retrieval_Terms: "
            + "; ".join(empty_active)
        )

    if missing_official:
        raise ValueError(
            "Official symbol missing from final KEEP Search_Terms for: "
            + "; ".join(missing_official)
        )

    active_collisions = {}

    for term, records in active_owners.items():
        genes = {gene for gene, _ in records}
        if len(genes) > 1:
            active_collisions[term] = records

    unsafe_keep_keep = {
        term: records
        for term, records in active_collisions.items()
        if all(decision == "KEEP" for _, decision in records)
    }

    if unsafe_keep_keep:
        examples = list(unsafe_keep_keep.items())[:20]
        raise ValueError(
            "Final alias table contains cross-GPCR collisions where all "
            f"active owners are KEEP. Examples: {examples}"
        )

    if len(cancer) != EXPECTED_CANCER_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_CANCER_COUNT} cancers; found {len(cancer)}."
        )

    if cancer["TCGA_Code"].astype(str).map(norm).duplicated().any():
        raise ValueError("Duplicate TCGA codes found.")

    empty_cancer = []
    cancer_term_owners = defaultdict(set)

    for _, row in cancer.iterrows():
        code = clean(row["TCGA_Code"])
        terms = split_terms(row["PubMed_Search_Terms"])

        if not terms:
            empty_cancer.append(code)

        for term in terms:
            cancer_term_owners[norm(term)].add(code)

    if empty_cancer:
        raise ValueError(
            "Cancer categories with empty PubMed_Search_Terms: "
            + "; ".join(empty_cancer)
        )

    shared_cancer_terms = {
        term: codes
        for term, codes in cancer_term_owners.items()
        if len(codes) > 1
    }

    if shared_cancer_terms:
        examples = list(shared_cancer_terms.items())[:20]
        raise ValueError(
            "Exact PubMed cancer terms are shared across TCGA categories. "
            f"Examples: {examples}"
        )

    cancer_term_count = sum(
        len(split_terms(x))
        for x in cancer["PubMed_Search_Terms"]
    )

    if cancer_term_count != EXPECTED_CANCER_TERM_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_CANCER_TERM_COUNT} cancer PubMed terms; "
            f"found {cancer_term_count}."
        )

    expected_searches = len(gpcr) * len(cancer)
    if expected_searches != EXPECTED_SEARCH_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_SEARCH_COUNT} GPCR-cancer searches; "
            f"found {expected_searches}."
        )

    gpcr.attrs["keep_term_count"] = keep_term_count
    gpcr.attrs["context_term_count"] = context_term_count
    gpcr.attrs["active_term_count"] = active_term_count
    gpcr.attrs["ignored_term_count"] = ignored_term_count
    gpcr.attrs["not_validated_term_count"] = not_validated_term_count
    gpcr.attrs["active_collision_count"] = len(active_collisions)

    return gpcr, cancer

# Construct the PubMed search query for the current item.
def build_query(gpcr_terms, cancer_terms):
    if MAX_GPCR_TERMS is not None:
        gpcr_terms = gpcr_terms[:MAX_GPCR_TERMS]

    if not gpcr_terms:
        raise ValueError("Cannot build query with no GPCR terms.")

    if not cancer_terms:
        raise ValueError("Cannot build query with no cancer terms.")

    gpcr_clause = " OR ".join(
        quote_pubmed_term(term)
        for term in gpcr_terms
    )

    cancer_clause = " OR ".join(
        quote_pubmed_term(term)
        for term in cancer_terms
    )

    human_clause = (
        "(humans[MeSH Terms] "
        "OR human[Title/Abstract] "
        "OR humans[Title/Abstract] "
        "OR patient[Title/Abstract] "
        "OR patients[Title/Abstract])"
    )

    return (
        f"({gpcr_clause}) "
        f"AND ({cancer_clause}) "
        f"AND {human_clause} "
        f"NOT Review[Publication Type]"
    )

# Validate search inputs and queries without live requests.
def run_preflight(gpcr, cancer):
    """Construct and verify every query without making network requests."""
    pair_keys = set()
    query_hashes = set()
    query_count = 0

    for _, gpcr_row in gpcr.iterrows():
        uniprot = clean(gpcr_row["UniProt_ID"])
        gpcr_terms = split_terms(gpcr_row["Active_Retrieval_Terms"])

        if MAX_GPCR_TERMS is not None:
            gpcr_terms = gpcr_terms[:MAX_GPCR_TERMS]

        for _, cancer_row in cancer.iterrows():
            code = clean(cancer_row["TCGA_Code"])
            cancer_terms = split_terms(
                cancer_row["PubMed_Search_Terms"]
            )

            query = build_query(gpcr_terms, cancer_terms)
            pair_key = (norm(uniprot), norm(code))
            digest = query_hash(query)

            if pair_key in pair_keys:
                raise ValueError(
                    f"Duplicate GPCR-cancer pair during preflight: "
                    f"{uniprot} x {code}"
                )

            if digest in query_hashes:
                raise ValueError(
                    f"Duplicate PubMed query during preflight: "
                    f"{uniprot} x {code}"
                )

            pair_keys.add(pair_key)
            query_hashes.add(digest)
            query_count += 1

    if query_count != EXPECTED_SEARCH_COUNT:
        raise ValueError(
            f"Expected {EXPECTED_SEARCH_COUNT} constructed queries; "
            f"found {query_count}."
        )

    total_keep = sum(
        len(split_terms(x)) for x in gpcr["Search_Terms"]
    )
    total_context = sum(
        len(split_terms(x)) for x in gpcr["Context_Required_Terms"]
    )
    total_active = sum(
        len(split_terms(x)) for x in gpcr["Active_Retrieval_Terms"]
    )
    total_ignored = sum(
        len(split_terms(x)) for x in gpcr["Ignored_Aliases"]
    )
    total_not_validated = sum(
        len(split_terms(x)) for x in gpcr["Not_Validated_Aliases"]
    )
    total_cancer_terms = sum(
        len(split_terms(x)) for x in cancer["PubMed_Search_Terms"]
    )

    print("PREFLIGHT PASSED")
    print(f"GPCRs                         : {len(gpcr)}")
    print(f"KEEP GPCR terms               : {total_keep}")
    print(f"CONTEXT_REQUIRED terms        : {total_context}")
    print(f"Active raw-retrieval terms    : {total_active}")
    print(f"IGNORED terms                 : {total_ignored}")
    print(f"NOT_VALIDATED terms           : {total_not_validated}")
    print(f"Active cross-GPCR collisions  : {gpcr.attrs.get('active_collision_count', 0)}")
    print(f"TCGA cancers                  : {len(cancer)}")
    print(f"Cancer PubMed terms           : {total_cancer_terms}")
    print(f"Queries constructed           : {query_count}")
    print("Network requests               : 0")

# Locate the saved response for a query or PMID.
def cache_path(query):
    digest = query_hash(query)
    return os.path.join(CACHE_DIR, f"{digest}.json")

# Read previously saved PubMed results.
def load_cache(query):
    path = cache_path(query)

    if not os.path.exists(path):
        return None

    try:
        with open(path, "r", encoding="utf-8") as fh:
            payload = json.load(fh)
    except Exception:
        return None

    if payload.get("query") != query:
        return None

    if payload.get("cache_version") != CACHE_VERSION:
        return None

    return payload

# Save a PubMed result for later reuse.
def save_cache(query, result):
    payload = {
        "cache_version": CACHE_VERSION,
        "query": query,
        "query_sha256": query_hash(query),
        "count": int(result["count"]),
        "retrieved": int(result["retrieved"]),
        "complete": bool(result["complete"]),
        "retrieval_method": result["retrieval_method"],
        "pmids": numeric_pmid_sort(result["pmids"]),
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    with open(cache_path(query), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)

# Issue a PubMed search with retries and rate limiting.
def esearch_request(
    query,
    retmax,
    mindate=None,
    maxdate=None,
):
    """
    One NCBI ESearch request with retries and rate limiting.
    """
    wait = INITIAL_RETRY_WAIT

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            kwargs = {
                "db": "pubmed",
                "term": query,
                "retmax": retmax,
                "retmode": "xml",
            }

            if mindate is not None and maxdate is not None:
                kwargs.update({
                    "datetype": "pdat",
                    "mindate": str(mindate),
                    "maxdate": str(maxdate),
                })

            handle = Entrez.esearch(**kwargs)
            record = Entrez.read(handle)
            handle.close()

            time.sleep(REQUEST_SLEEP)
            return record

        except Exception as exc:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"PubMed ESearch failed after {MAX_RETRIES} attempts.\n"
                    f"Query: {query}\n"
                    f"mindate={mindate}, maxdate={maxdate}\n"
                    f"Last error: {exc}"
                ) from exc

            error_name = type(exc).__name__
            print(
                f"    Request failed "
                f"(attempt {attempt}/{MAX_RETRIES}) "
                f"[{error_name}]: {exc}"
            )
            print(f"    Retrying in {wait:.1f} sec...")
            time.sleep(wait)
            wait *= 2

# Retrieve PMIDs with one PubMed query.
def direct_search(query):
    record = esearch_request(
        query=query,
        retmax=ESEARCH_RETMAX,
    )

    count = int(record.get("Count", 0))
    pmids = set(record.get("IdList", []))

    return count, pmids

# Count PubMed matches within a publication-date range.
def date_range_count(query, start_year, end_year):
    record = esearch_request(
        query=query,
        retmax=0,
        mindate=start_year,
        maxdate=end_year,
    )
    return int(record.get("Count", 0))

# Retrieve PMIDs within a publication-date range.
def date_range_pmids(query, start_year, end_year):
    record = esearch_request(
        query=query,
        retmax=ESEARCH_RETMAX,
        mindate=start_year,
        maxdate=end_year,
    )

    count = int(record.get("Count", 0))
    pmids = set(record.get("IdList", []))

    return count, pmids

# Split large searches by year to recover all PMIDs.
def recover_by_year_split(query, expected_total):
    """
    Recover large PubMed result sets by recursively splitting publication
    years until each interval contains <= ESEARCH_RETMAX records.

    Returns:
        pmids
        exact_recovery (bool)
        intervals_used
    """
    recovered = set()
    intervals_used = []

    stack = [(MIN_PUBLICATION_YEAR, MAX_PUBLICATION_YEAR)]

    while stack:
        start_year, end_year = stack.pop()

        count = date_range_count(query, start_year, end_year)

        if count == 0:
            continue

        if count <= ESEARCH_RETMAX:
            returned_count, pmids = date_range_pmids(
                query,
                start_year,
                end_year,
            )

            intervals_used.append(
                {
                    "start": start_year,
                    "end": end_year,
                    "count": returned_count,
                    "retrieved": len(pmids),
                }
            )

            recovered.update(pmids)
            continue

        if start_year >= end_year:
            intervals_used.append(
                {
                    "start": start_year,
                    "end": end_year,
                    "count": count,
                    "retrieved": 0,
                    "unrecoverable": True,
                }
            )
            return recovered, False, intervals_used

        midpoint = (start_year + end_year) // 2

        stack.append((midpoint + 1, end_year))
        stack.append((start_year, midpoint))

    exact = len(recovered) == expected_total

    return recovered, exact, intervals_used

# Use the saved result or execute a PubMed search.
def search_pubmed(query):
    """
    Search one complete GPCR x cancer query.

    Uses cache when available.
    """
    cached = load_cache(query)

    if cached is not None:
        return {
            "count": int(cached["count"]),
            "pmids": set(cached.get("pmids", [])),
            "retrieved": int(cached["retrieved"]),
            "complete": bool(cached["complete"]),
            "retrieval_method": cached["retrieval_method"],
            "cache_used": True,
        }

    count, pmids = direct_search(query)

    if len(pmids) == count:
        result = {
            "count": count,
            "pmids": pmids,
            "retrieved": len(pmids),
            "complete": True,
            "retrieval_method": "DIRECT_ESEARCH",
            "cache_used": False,
        }
        save_cache(query, result)
        return result

    recovered, exact, intervals = recover_by_year_split(
        query,
        expected_total=count,
    )

    result = {
        "count": count,
        "pmids": recovered,
        "retrieved": len(recovered),
        "complete": exact,
        "retrieval_method": (
            "DATE_SPLIT_RECOVERY"
            if exact
            else "INCOMPLETE_DATE_SPLIT_RECOVERY"
        ),
        "cache_used": False,
    }

    save_cache(query, result)

    return result

# Save progress if a long search is interrupted.
def save_partial(rows):
    if rows:
        pd.DataFrame(rows).to_csv(
            OUT_PARTIAL,
            index=False,
            na_rep="",
        )

# Expand query results into GPCR, cancer, and PMID rows.
def build_long_table(summary):
    rows = []

    hits = summary[
        summary["Retrieved_PMIDs"].fillna(0).astype(int) > 0
    ]

    for _, row in hits.iterrows():
        pmids = split_terms(row["PMID_List"])

        for pmid in pmids:
            rows.append({
                "UniProt_ID": row["UniProt_ID"],
                "Gene": row["Gene"],
                "TCGA_Code": row["TCGA_Code"],
                "TCGA_Study_Name": row["TCGA_Study_Name"],
                "Cancer": row["Cancer"],
                "PMID": pmid,
                "Search_ID": row["Search_ID"],
            })

    return pd.DataFrame(
        rows,
        columns=[
            "UniProt_ID",
            "Gene",
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PMID",
            "Search_ID",
        ],
    )

# Write integrity report for the current record.
def write_integrity_report(
    gpcr,
    cancer,
    summary,
    long_df,
    error_df,
    elapsed_minutes,
):
    expected = len(gpcr) * len(cancer)

    complete = int(
        (summary["Search_Status"] == "COMPLETE").sum()
    )
    incomplete = int(
        (summary["Search_Status"] == "INCOMPLETE").sum()
    )
    errors = int(
        (summary["Search_Status"] == "ERROR").sum()
    )
    nonzero = int((summary["PubMed_Count"] > 0).sum())

    unique_pmids = (
        long_df["PMID"].nunique()
        if not long_df.empty
        else 0
    )

    recovered = int(
        (
            summary["Retrieval_Method"] == "DATE_SPLIT_RECOVERY"
        ).sum()
    )

    total_keep = sum(
        len(split_terms(x)) for x in gpcr["Search_Terms"]
    )
    total_context = sum(
        len(split_terms(x)) for x in gpcr["Context_Required_Terms"]
    )
    total_active = sum(
        len(split_terms(x)) for x in gpcr["Active_Retrieval_Terms"]
    )
    total_ignored = sum(
        len(split_terms(x)) for x in gpcr["Ignored_Aliases"]
    )
    total_not_validated = sum(
        len(split_terms(x)) for x in gpcr["Not_Validated_Aliases"]
    )

    with open(OUT_REPORT, "w", encoding="utf-8") as fh:
        fh.write("GPCR Cancer Pipeline - Step 06 Integrity Report\n")
        fh.write("=" * 68 + "\n\n")

        fh.write("INPUTS\n")
        fh.write("-" * 68 + "\n")
        fh.write(f"GPCRs                    : {len(gpcr)}\n")
        fh.write(f"TCGA cancer types        : {len(cancer)}\n")
        fh.write(f"KEEP GPCR terms          : {total_keep}\n")
        fh.write(f"CONTEXT_REQUIRED terms   : {total_context}\n")
        fh.write(f"Active retrieval terms   : {total_active}\n")
        fh.write(f"IGNORED terms            : {total_ignored}\n")
        fh.write(f"NOT_VALIDATED terms      : {total_not_validated}\n")
        fh.write(
            f"Active collision terms   : "
            f"{gpcr.attrs.get('active_collision_count', 0)}\n"
        )
        fh.write(
            f"GPCR term cap            : "
            f"{MAX_GPCR_TERMS if MAX_GPCR_TERMS is not None else 'NONE'}\n"
        )
        fh.write(
            f"Cancer PubMed terms      : "
            f"{sum(len(split_terms(x)) for x in cancer['PubMed_Search_Terms'])}\n"
        )
        fh.write(
            f"NCBI API key             : {'YES' if API_KEY else 'NO'}\n"
        )
        fh.write(
            f"Network timeout sec      : {SOCKET_TIMEOUT_SECONDS}\n\n"
        )

        fh.write("SEARCH ACCOUNTING\n")
        fh.write("-" * 68 + "\n")
        fh.write(f"Expected searches        : {expected}\n")
        fh.write(f"Observed searches        : {len(summary)}\n")
        fh.write(f"Complete searches        : {complete}\n")
        fh.write(f"Incomplete searches      : {incomplete}\n")
        fh.write(f"Error searches           : {errors}\n")
        fh.write(f"Nonzero searches         : {nonzero}\n")
        fh.write(f"Date-split recoveries    : {recovered}\n")
        fh.write(f"Unique PMIDs             : {unique_pmids}\n")
        fh.write(f"Long-table rows          : {len(long_df)}\n")
        fh.write(f"Elapsed minutes          : {elapsed_minutes:.1f}\n\n")

        fh.write("INTEGRITY CHECKS\n")
        fh.write("-" * 68 + "\n")
        fh.write(
            "832 GPCRs: "
            + ("PASS" if len(gpcr) == EXPECTED_GPCR_COUNT else "FAIL")
            + "\n"
        )
        fh.write(
            "33 TCGA cancers: "
            + ("PASS" if len(cancer) == EXPECTED_CANCER_COUNT else "FAIL")
            + "\n"
        )
        fh.write(
            "27,456 pairwise searches: "
            + ("PASS" if len(summary) == expected else "FAIL")
            + "\n"
        )
        fh.write(
            "No network/search errors: "
            + ("PASS" if errors == 0 else "FAIL")
            + "\n"
        )
        fh.write(
            "No incomplete PMID sets: "
            + ("PASS" if incomplete == 0 else "FAIL")
            + "\n"
        )
        fh.write(
            "All PubMed counts equal retrieved PMID counts: "
            + (
                "PASS"
                if (
                    summary["PubMed_Count"].fillna(-1).astype(int)
                    == summary["Retrieved_PMIDs"].fillna(-2).astype(int)
                ).all()
                else "FAIL"
            )
            + "\n"
        )

        fh.write("\nMETHOD NOTE\n")
        fh.write("-" * 68 + "\n")
        fh.write(
            "Step 06 is raw PMID retrieval only. KEEP and CONTEXT_REQUIRED "
            "terms are both allowed to retrieve candidate papers. "
            "CONTEXT_REQUIRED terms are not treated as validated standalone "
            "aliases; downstream article-level validation must enforce their "
            "context requirement before acceptance. Titles, abstracts, matched "
            "terms, and evidence classification are not performed here.\n"
        )

        if not error_df.empty:
            fh.write("\nSEARCHES REQUIRING ATTENTION\n")
            fh.write("-" * 68 + "\n")
            for _, row in error_df.head(100).iterrows():
                fh.write(
                    f"{clean(row.get('Gene'))} | "
                    f"{clean(row.get('TCGA_Code'))} | "
                    f"{clean(row.get('Cancer'))} | "
                    f"{clean(row.get('Search_Status'))} | "
                    f"{clean(row.get('Error_Message'))}\n"
                )

# Run this stage and write its output files.
def main():
    args = parse_args()

    print("=" * 76)
    print("GPCR Cancer Pipeline - Step 06 PubMed Search")
    print("=" * 76)

    try:
        gpcr, cancer = load_and_validate_inputs()
    except Exception as exc:
        print(f"INPUT VALIDATION FAILED: {exc}")
        sys.exit(1)

    total_keep = sum(
        len(split_terms(x))
        for x in gpcr["Search_Terms"]
    )
    total_context = sum(
        len(split_terms(x))
        for x in gpcr["Context_Required_Terms"]
    )
    total_active = sum(
        len(split_terms(x))
        for x in gpcr["Active_Retrieval_Terms"]
    )
    total_cancer_terms = sum(
        len(split_terms(x))
        for x in cancer["PubMed_Search_Terms"]
    )

    expected_searches = len(gpcr) * len(cancer)

    print(f"GPCRs                    : {len(gpcr)}")
    print(f"KEEP terms               : {total_keep}")
    print(f"CONTEXT_REQUIRED terms   : {total_context}")
    print(f"Active retrieval terms   : {total_active}")
    print(
        f"GPCR term cap            : "
        f"{MAX_GPCR_TERMS if MAX_GPCR_TERMS is not None else 'NONE'}"
    )
    print(f"TCGA cancers             : {len(cancer)}")
    print(f"Cancer PubMed terms      : {total_cancer_terms}")
    print(f"Pairwise searches        : {expected_searches}")
    print(f"NCBI API key             : {'YES' if API_KEY else 'NO'}")
    print(f"Network timeout          : {SOCKET_TIMEOUT_SECONDS} sec")
    print(f"Cache                    : {CACHE_DIR}")
    print("=" * 76)

    if args.preflight:
        try:
            run_preflight(gpcr, cancer)
        except Exception as exc:
            print(f"PREFLIGHT FAILED: {exc}")
            sys.exit(1)
        return

    if not clean(Entrez.email):
        raise SystemExit(
            "NCBI email is required.\n"
            'Run: export NCBI_EMAIL="your_email@example.com"'
        )

    os.makedirs(CACHE_DIR, exist_ok=True)

    rows = []
    start_time = time.time()
    search_number = 0

    for _, gpcr_row in gpcr.iterrows():
        uniprot = clean(gpcr_row["UniProt_ID"])
        gene = clean(gpcr_row["Official_Gene"])

        keep_terms = split_terms(gpcr_row["Search_Terms"])
        context_terms = split_terms(
            gpcr_row["Context_Required_Terms"]
        )
        gpcr_terms_all = split_terms(
            gpcr_row["Active_Retrieval_Terms"]
        )

        gpcr_terms_used = (
            gpcr_terms_all[:MAX_GPCR_TERMS]
            if MAX_GPCR_TERMS is not None
            else gpcr_terms_all
        )

        for _, cancer_row in cancer.iterrows():
            search_number += 1

            code = clean(cancer_row["TCGA_Code"])
            tcga_name = clean(cancer_row["TCGA_Study_Name"])
            cancer_name = clean(cancer_row["Cancer"])
            cancer_terms = split_terms(
                cancer_row["PubMed_Search_Terms"]
            )

            query = build_query(
                gpcr_terms_used,
                cancer_terms,
            )

            search_id = query_hash(query)[:20]

            elapsed_min = (time.time() - start_time) / 60.0

            if search_number > 1:
                avg_seconds = (
                    time.time() - start_time
                ) / (search_number - 1)

                remaining = expected_searches - search_number + 1
                eta_min = (avg_seconds * remaining) / 60.0
            else:
                eta_min = 0.0

            print(
                f"[{search_number:5d}/{expected_searches}] "
                f"{gene:<12} | {code:<4} | "
                f"elapsed {elapsed_min:7.1f} min | "
                f"ETA {eta_min:7.1f} min"
            )

            error_message = ""

            try:
                result = search_pubmed(query)

                status = (
                    "COMPLETE"
                    if result["complete"]
                    else "INCOMPLETE"
                )

                print(
                    f"    {cancer_name}: "
                    f"{result['count']} hits, "
                    f"{result['retrieved']} PMIDs, "
                    f"{result['retrieval_method']}, "
                    f"cache={'YES' if result['cache_used'] else 'NO'}"
                )

            except Exception as exc:
                result = {
                    "count": -1,
                    "pmids": set(),
                    "retrieved": 0,
                    "complete": False,
                    "retrieval_method": "ERROR",
                    "cache_used": False,
                }

                status = "ERROR"
                error_message = str(exc).replace("\n", " | ")

                print(f"    ERROR: {error_message}")

            rows.append({
                "Search_ID": search_id,
                "UniProt_ID": uniprot,
                "Gene": gene,
                "TCGA_Code": code,
                "TCGA_Study_Name": tcga_name,
                "Cancer": cancer_name,

                "GPCR_KEEP_Terms": ";".join(keep_terms),
                "GPCR_KEEP_Term_Count": len(keep_terms),
                "GPCR_Context_Required_Terms": ";".join(context_terms),
                "GPCR_Context_Required_Term_Count": len(context_terms),
                "GPCR_Search_Terms": ";".join(gpcr_terms_used),
                "GPCR_Term_Count": len(gpcr_terms_used),

                "Cancer_Search_Terms": ";".join(cancer_terms),
                "Cancer_Term_Count": len(cancer_terms),

                "PubMed_Query": query,
                "PubMed_Count": result["count"],
                "Retrieved_PMIDs": result["retrieved"],
                "PMID_List": join_pmids(result["pmids"]),

                "Search_Status": status,
                "Retrieval_Method": result["retrieval_method"],
                "Cache_Used": result["cache_used"],
                "Error_Message": error_message,
            })

            if (
                search_number % CHECKPOINT_EVERY == 0
                or search_number == expected_searches
            ):
                save_partial(rows)

    elapsed_minutes = (time.time() - start_time) / 60.0

    summary = pd.DataFrame(rows)

    duplicate_pairs = summary.duplicated(
        subset=["UniProt_ID", "TCGA_Code"],
        keep=False,
    )

    if duplicate_pairs.any():
        dup = summary.loc[
            duplicate_pairs,
            ["UniProt_ID", "Gene", "TCGA_Code", "Cancer"],
        ]
        print("FATAL: duplicate GPCR-cancer search pairs detected.")
        print(dup.head(50).to_string(index=False))
        sys.exit(1)

    summary.to_csv(
        OUT_SUMMARY,
        index=False,
        na_rep="",
    )

    long_df = build_long_table(summary)
    long_df.to_csv(
        OUT_LONG,
        index=False,
        na_rep="",
    )

    errors = summary[
        summary["Search_Status"].isin(["ERROR", "INCOMPLETE"])
    ].copy()

    errors.to_csv(
        OUT_ERRORS,
        index=False,
        na_rep="",
    )

    write_integrity_report(
        gpcr=gpcr,
        cancer=cancer,
        summary=summary,
        long_df=long_df,
        error_df=errors,
        elapsed_minutes=elapsed_minutes,
    )

    print()
    print("=" * 76)

    if errors.empty and len(summary) == expected_searches:
        print("STEP 06 COMPLETE - ALL INTEGRITY CHECKS PASSED")
    else:
        print("STEP 06 COMPLETE WITH UNRESOLVED SEARCH ISSUES")

    print("=" * 76)
    print(f"Expected searches     : {expected_searches}")
    print(f"Completed rows        : {len(summary)}")
    print(
        f"Nonzero searches      : "
        f"{int((summary['PubMed_Count'] > 0).sum())}"
    )
    print(
        f"Search errors         : "
        f"{int((summary['Search_Status'] == 'ERROR').sum())}"
    )
    print(
        f"Incomplete PMID sets  : "
        f"{int((summary['Search_Status'] == 'INCOMPLETE').sum())}"
    )
    print(
        f"Unique PMIDs          : "
        f"{long_df['PMID'].nunique() if not long_df.empty else 0}"
    )
    print(f"Elapsed               : {elapsed_minutes:.1f} min")
    print()
    print("Saved:")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_LONG}")
    print(f"  {OUT_ERRORS}")
    print(f"  {OUT_REPORT}")
    print(f"  {OUT_PARTIAL}")
    print("=" * 76)

    if not errors.empty:
        sys.exit(2)

if __name__ == "__main__":
    main()
