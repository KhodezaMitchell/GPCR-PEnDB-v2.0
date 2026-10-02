#!/usr/bin/env python3
# Purpose: Validate GPCR and cancer terminology in PubMed titles and abstracts.
# Input: GPCR_Cancer_PubMed_PMIDs_Long.csv; GPCR_Alias_Table_FINAL_VALIDATED.csv; TCGA_Cancer_Vocabulary.csv
# Output: Human_GPCR_Cancer_PubMed_Papers_Validated.csv; GPCR_Cancer_PubMed_Context_Review.csv

import argparse
import json
import os
import re
import socket
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "").strip()

API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
if API_KEY:
    Entrez.api_key = API_KEY

INPUT_LONG = "GPCR_Cancer_PubMed_PMIDs_Long.csv"
INPUT_SEARCH_SUMMARY = "GPCR_Cancer_PubMed_Search_Summary.csv"
ALIAS_FILE = "GPCR_Alias_Table_FINAL_VALIDATED.csv"
CANCER_FILE = "TCGA_Cancer_Vocabulary.csv"

OUT_METADATA = "PubMed_Article_Metadata.csv"
OUT_ALL = "GPCR_Cancer_PubMed_Validation_All.csv"
OUT_KEPT = "Human_GPCR_Cancer_PubMed_Papers_Validated.csv"
OUT_CONTEXT = "GPCR_Cancer_PubMed_Context_Review.csv"
OUT_SUMMARY = "GPCR_Cancer_PubMed_Validation_Summary.csv"
OUT_ERRORS = "PubMed_Article_Metadata_Errors.csv"
OUT_REPORT = "GPCR_Cancer_PubMed_Validation_Integrity_Report.txt"

CACHE_DIR = "pubmed_article_cache_final"

EXPECTED_GPCR_COUNT = 832
EXPECTED_CANCER_COUNT = 33
EXPECTED_PAIR_COUNT = EXPECTED_GPCR_COUNT * EXPECTED_CANCER_COUNT
EXPECTED_CANCER_TERM_COUNT = 122

FETCH_BATCH_SIZE = 200
MAX_RETRIES = 5
INITIAL_RETRY_WAIT = 2.0

SOCKET_TIMEOUT_SECONDS = int(
    os.environ.get("NCBI_SOCKET_TIMEOUT", "60")
)
socket.setdefaulttimeout(SOCKET_TIMEOUT_SECONDS)

REQUEST_SLEEP = 0.34 if not API_KEY else 0.11

REJECTED_CANCER_TERMS = {
    "BLCA": ["gallbladder", "gall bladder", "biliary"],
    "LGG": ["glioblastoma"],
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

# Return a numeric key for sorting PMIDs.
def pmid_sort_key(value):
    value = clean(value)
    return int(value) if value.isdigit() else value

# Read command-line options for this stage.
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Fetch PubMed metadata and validate cancer + GPCR support "
            "for every Step 06 candidate association."
        )
    )
    parser.add_argument(
        "--preflight",
        action="store_true",
        help=(
            "validate frozen inputs and cross-file relationships "
            "without contacting PubMed"
        ),
    )
    return parser.parse_args()

# Build a regex for an exact GPCR term.
def exact_term_pattern(term):
    term = clean(term)
    if not term:
        return re.compile(r"(?!x)x")

    return re.compile(
        rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])",
        re.IGNORECASE,
    )

# Build a regex allowing controlled punctuation variants.
def flexible_term_pattern(term):
    """
    Conservative GPCR-term matcher.

    Flexibility is allowed only:
      1) where the stored term already has separator punctuation, or
      2) at letter-number / number-letter boundaries.
    """
    term = clean(term)
    if not term:
        return re.compile(r"(?!x)x")

    pieces = []
    i = 0

    while i < len(term):
        ch = term[i]

        if re.match(r"[\s\-_/\.]", ch):
            while (
                i + 1 < len(term)
                and re.match(r"[\s\-_/\.]", term[i + 1])
            ):
                i += 1
            pieces.append(r"[\s\-_/\.]*")
            i += 1
            continue

        pieces.append(re.escape(ch))

        if i + 1 < len(term):
            nxt = term[i + 1]
            if ch.isalnum() and nxt.isalnum():
                if (
                    (ch.isalpha() and nxt.isdigit())
                    or (ch.isdigit() and nxt.isalpha())
                ):
                    pieces.append(r"[\s\-_/\.]*")

        i += 1

    body = "".join(pieces)

    return re.compile(
        rf"(?<![A-Za-z0-9]){body}(?![A-Za-z0-9])",
        re.IGNORECASE,
    )

# Find matched GPCR terms in a title or abstract.
def gpcr_term_matches(text, terms):
    """
    Return:
        all_matches,
        exact_matches,
        flexible_only_matches
    """
    text = clean(text)

    all_matches = []
    exact_matches = []
    flexible_only_matches = []

    for term in terms:
        exact = exact_term_pattern(term)
        flexible = flexible_term_pattern(term)

        if exact.search(text):
            all_matches.append(term)
            exact_matches.append(term)
        elif flexible.search(text):
            all_matches.append(term)
            flexible_only_matches.append(term)

    return (
        all_matches,
        exact_matches,
        flexible_only_matches,
    )

# Check whether text contains a complete phrase.
def contains_phrase(text, phrase):
    """
    Preserve previous cancer-term validation behavior.

    Short alphanumeric terms (<=5 characters) are whole-word matched.
    Longer terms are case-insensitive substring matched.
    """
    text = clean(text).lower()
    phrase = clean(phrase).lower()

    if not phrase:
        return False

    if len(phrase) <= 5 and phrase.replace("-", "").isalnum():
        return (
            re.search(
                r"\b" + re.escape(phrase) + r"\b",
                text,
            )
            is not None
        )

    return phrase in text

# Find accepted and excluded cancer terms in a record.
def matched_cancer_terms(text, terms):
    return [
        term
        for term in terms
        if contains_phrase(text, term)
    ]

# Validate cancer for the current record.
def validate_cancer(code, title, abstract, accepted_terms):
    text = f"{clean(title)} {clean(abstract)}"

    rejected_terms = REJECTED_CANCER_TERMS.get(code, [])

    rejected_matches = matched_cancer_terms(
        text,
        rejected_terms,
    )
    accepted_matches = matched_cancer_terms(
        text,
        accepted_terms,
    )

    if rejected_matches:
        return (
            "REJECT",
            "rejected_cancer_term",
            accepted_matches,
            rejected_matches,
        )

    if accepted_matches:
        return (
            "KEEP",
            "accepted_cancer_term",
            accepted_matches,
            [],
        )

    return (
        "REJECT",
        "no_accepted_cancer_term_in_title_abstract",
        [],
        [],
    )

# Validate gpcr for the current record.
def validate_gpcr(title, abstract, keep_terms, context_terms):
    text = f"{clean(title)} {clean(abstract)}"

    (
        keep_matches,
        keep_exact,
        keep_flexible,
    ) = gpcr_term_matches(text, keep_terms)

    (
        context_matches,
        context_exact,
        context_flexible,
    ) = gpcr_term_matches(text, context_terms)

    if keep_matches:
        return {
            "Decision": "KEEP",
            "Reason": "matched_final_keep_gpcr_term",
            "KEEP_Matches": keep_matches,
            "KEEP_Exact": keep_exact,
            "KEEP_Flexible": keep_flexible,
            "Context_Matches": context_matches,
            "Context_Exact": context_exact,
            "Context_Flexible": context_flexible,
        }

    if context_matches:
        return {
            "Decision": "CONTEXT_REVIEW",
            "Reason": "context_required_gpcr_term_only",
            "KEEP_Matches": [],
            "KEEP_Exact": [],
            "KEEP_Flexible": [],
            "Context_Matches": context_matches,
            "Context_Exact": context_exact,
            "Context_Flexible": context_flexible,
        }

    return {
        "Decision": "REJECT",
        "Reason": "no_active_gpcr_term_in_title_abstract",
        "KEEP_Matches": [],
        "KEEP_Exact": [],
        "KEEP_Flexible": [],
        "Context_Matches": [],
        "Context_Exact": [],
        "Context_Flexible": [],
    }

# Combine receptor and cancer validation decisions.
def overall_validation(cancer_decision, gpcr_decision):
    if cancer_decision == "ERROR" or gpcr_decision == "ERROR":
        return "ERROR"

    if cancer_decision != "KEEP":
        return "REJECT"

    if gpcr_decision == "KEEP":
        return "KEEP"

    if gpcr_decision == "CONTEXT_REVIEW":
        return "CONTEXT_REVIEW"

    return "REJECT"

# Load inputs from the input records.
def load_inputs():
    required_files = [
        INPUT_LONG,
        INPUT_SEARCH_SUMMARY,
        ALIAS_FILE,
        CANCER_FILE,
    ]

    missing = [f for f in required_files if not os.path.exists(f)]
    if missing:
        raise FileNotFoundError(
            "Missing required file(s): " + ", ".join(missing)
        )

    long_df = pd.read_csv(
        INPUT_LONG,
        dtype={"PMID": str},
    ).fillna("")

    search_summary = pd.read_csv(
        INPUT_SEARCH_SUMMARY,
        dtype={"PMID_List": str},
    ).fillna("")

    aliases = pd.read_csv(ALIAS_FILE).fillna("")
    cancers = pd.read_csv(CANCER_FILE).fillna("")

    require_columns(
        long_df,
        [
            "UniProt_ID",
            "Gene",
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PMID",
            "Search_ID",
        ],
        INPUT_LONG,
    )

    require_columns(
        search_summary,
        [
            "Search_ID",
            "UniProt_ID",
            "Gene",
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "GPCR_KEEP_Terms",
            "GPCR_Context_Required_Terms",
            "GPCR_Search_Terms",
            "Cancer_Search_Terms",
            "PubMed_Query",
            "PubMed_Count",
            "Retrieved_PMIDs",
            "PMID_List",
            "Search_Status",
        ],
        INPUT_SEARCH_SUMMARY,
    )

    require_columns(
        aliases,
        [
            "UniProt_ID",
            "Official_Gene",
            "Original_Gene",
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
        cancers,
        [
            "TCGA_Code",
            "TCGA_Study_Name",
            "Cancer",
            "PubMed_Search_Terms",
        ],
        CANCER_FILE,
    )

    errors = []

    if len(aliases) != EXPECTED_GPCR_COUNT:
        errors.append(
            f"Expected {EXPECTED_GPCR_COUNT} GPCR rows; found {len(aliases)}"
        )

    if len(cancers) != EXPECTED_CANCER_COUNT:
        errors.append(
            f"Expected {EXPECTED_CANCER_COUNT} cancers; found {len(cancers)}"
        )

    if len(search_summary) != EXPECTED_PAIR_COUNT:
        errors.append(
            f"Expected {EXPECTED_PAIR_COUNT} search rows; "
            f"found {len(search_summary)}"
        )

    if aliases["UniProt_ID"].map(norm).duplicated().any():
        errors.append("Duplicate UniProt_ID rows in alias table")

    if cancers["TCGA_Code"].map(norm).duplicated().any():
        errors.append("Duplicate TCGA_Code rows in cancer vocabulary")

    if search_summary["Search_ID"].map(clean).duplicated().any():
        errors.append("Duplicate Search_ID values in search summary")

    if search_summary.duplicated(
        subset=["UniProt_ID", "TCGA_Code"]
    ).any():
        errors.append("Duplicate GPCR-cancer pairs in search summary")

    long_df["PMID"] = long_df["PMID"].astype(str).map(clean)

    if long_df["PMID"].eq("").any():
        errors.append("Blank PMID found in long table")

    if not long_df["PMID"].str.fullmatch(r"\d+").all():
        errors.append("Nonnumeric PMID found in long table")

    if long_df.duplicated(
        subset=["UniProt_ID", "TCGA_Code", "PMID"]
    ).any():
        errors.append("Duplicate GPCR-cancer-PMID rows in long table")

    if long_df.duplicated(
        subset=["Search_ID", "PMID"]
    ).any():
        errors.append("Duplicate Search_ID-PMID rows in long table")

    bad_status = search_summary[
        search_summary["Search_Status"].map(clean) != "COMPLETE"
    ]
    if not bad_status.empty:
        errors.append(
            f"{len(bad_status)} Step 06 searches are not COMPLETE"
        )

    pubmed_counts = pd.to_numeric(
        search_summary["PubMed_Count"],
        errors="coerce",
    )
    retrieved_counts = pd.to_numeric(
        search_summary["Retrieved_PMIDs"],
        errors="coerce",
    )

    if pubmed_counts.isna().any() or retrieved_counts.isna().any():
        errors.append(
            "Nonnumeric PubMed_Count/Retrieved_PMIDs in search summary"
        )
    else:
        if not (
            pubmed_counts.astype(int)
            == retrieved_counts.astype(int)
        ).all():
            errors.append(
                "PubMed_Count and Retrieved_PMIDs disagree"
            )

        if int(retrieved_counts.sum()) != len(long_df):
            errors.append(
                "Sum of Retrieved_PMIDs does not equal long-table row count"
            )

    alias_ids = set(aliases["UniProt_ID"].map(clean))
    summary_ids = set(search_summary["UniProt_ID"].map(clean))
    cancer_codes = set(cancers["TCGA_Code"].map(clean))
    summary_codes = set(search_summary["TCGA_Code"].map(clean))

    if alias_ids != summary_ids:
        errors.append(
            "Search-summary GPCR set does not match final alias table"
        )

    if cancer_codes != summary_codes:
        errors.append(
            "Search-summary cancer set does not match cancer vocabulary"
        )

    alias_lookup = {}

    for _, row in aliases.iterrows():
        uid = clean(row["UniProt_ID"])
        keep = split_terms(row["Search_Terms"])
        context = split_terms(row["Context_Required_Terms"])
        active = split_terms(row["Active_Retrieval_Terms"])

        keep_keys = {norm(x) for x in keep}
        context_keys = {norm(x) for x in context}
        active_keys = {norm(x) for x in active}

        if not keep:
            errors.append(
                f"{uid}: empty final KEEP Search_Terms"
            )

        if norm(row["Official_Gene"]) not in keep_keys:
            errors.append(
                f"{uid}: official gene symbol missing from KEEP terms"
            )

        if active_keys != (keep_keys | context_keys):
            errors.append(
                f"{uid}: active terms != KEEP + CONTEXT_REQUIRED"
            )

        if keep_keys & context_keys:
            errors.append(
                f"{uid}: KEEP and CONTEXT_REQUIRED overlap"
            )

        alias_lookup[uid] = {
            "KEEP": keep,
            "CONTEXT": context,
            "ACTIVE": active,
        }

    for _, row in search_summary.iterrows():
        uid = clean(row["UniProt_ID"])
        frozen = alias_lookup.get(uid)

        if frozen is None:
            continue

        used_keep = {norm(x) for x in split_terms(row["GPCR_KEEP_Terms"])}
        used_context = {
            norm(x)
            for x in split_terms(row["GPCR_Context_Required_Terms"])
        }
        used_active = {norm(x) for x in split_terms(row["GPCR_Search_Terms"])}

        if used_keep != {norm(x) for x in frozen["KEEP"]}:
            errors.append(
                f"{uid}: Step 06 KEEP terms differ from frozen alias table"
            )
            break

        if used_context != {norm(x) for x in frozen["CONTEXT"]}:
            errors.append(
                f"{uid}: Step 06 CONTEXT_REQUIRED terms differ from frozen alias table"
            )
            break

        if used_active != {norm(x) for x in frozen["ACTIVE"]}:
            errors.append(
                f"{uid}: Step 06 active terms differ from frozen alias table"
            )
            break

    cancer_term_count = sum(
        len(split_terms(x))
        for x in cancers["PubMed_Search_Terms"]
    )
    if cancer_term_count != EXPECTED_CANCER_TERM_COUNT:
        errors.append(
            f"Expected {EXPECTED_CANCER_TERM_COUNT} cancer terms; "
            f"found {cancer_term_count}"
        )

    expected_pmids_by_search = {}
    identity_by_search = {}

    for _, row in search_summary.iterrows():
        sid = clean(row["Search_ID"])
        pmids = set(split_terms(row["PMID_List"]))

        if len(pmids) != int(float(clean(row["Retrieved_PMIDs"]) or 0)):
            errors.append(
                f"{sid}: PMID_List count != Retrieved_PMIDs"
            )
            break

        expected_pmids_by_search[sid] = pmids
        identity_by_search[sid] = (
            clean(row["UniProt_ID"]),
            clean(row["Gene"]),
            clean(row["TCGA_Code"]),
            clean(row["TCGA_Study_Name"]),
            clean(row["Cancer"]),
        )

    observed_pmids_by_search = defaultdict(set)

    for _, row in long_df.iterrows():
        sid = clean(row["Search_ID"])

        if sid not in identity_by_search:
            errors.append(
                f"Long-table Search_ID absent from summary: {sid}"
            )
            break

        observed_identity = (
            clean(row["UniProt_ID"]),
            clean(row["Gene"]),
            clean(row["TCGA_Code"]),
            clean(row["TCGA_Study_Name"]),
            clean(row["Cancer"]),
        )

        if observed_identity != identity_by_search[sid]:
            errors.append(
                f"Identity mismatch for Search_ID {sid}"
            )
            break

        observed_pmids_by_search[sid].add(clean(row["PMID"]))

    if not errors:
        for sid, expected_pmids in expected_pmids_by_search.items():
            if observed_pmids_by_search.get(sid, set()) != expected_pmids:
                errors.append(
                    f"Summary/long-table PMID mismatch for Search_ID {sid}"
                )
                break

    if errors:
        raise ValueError(
            "\n".join(f"- {x}" for x in errors)
        )

    return long_df, search_summary, aliases, cancers

# Validate search inputs and queries without live requests.
def run_preflight(long_df, search_summary, aliases, cancers):
    unique_pmids = long_df["PMID"].nunique()

    keep_terms = sum(
        len(split_terms(x))
        for x in aliases["Search_Terms"]
    )
    context_terms = sum(
        len(split_terms(x))
        for x in aliases["Context_Required_Terms"]
    )
    active_terms = sum(
        len(split_terms(x))
        for x in aliases["Active_Retrieval_Terms"]
    )
    cancer_terms = sum(
        len(split_terms(x))
        for x in cancers["PubMed_Search_Terms"]
    )

    expected_batches = (
        unique_pmids + FETCH_BATCH_SIZE - 1
    ) // FETCH_BATCH_SIZE

    print("PREFLIGHT PASSED", flush=True)
    print(f"GPCRs                         : {len(aliases)}", flush=True)
    print(f"KEEP GPCR terms               : {keep_terms}", flush=True)
    print(f"CONTEXT_REQUIRED GPCR terms   : {context_terms}", flush=True)
    print(f"Active GPCR terms             : {active_terms}", flush=True)
    print(f"TCGA cancers                  : {len(cancers)}", flush=True)
    print(f"Cancer PubMed terms           : {cancer_terms}", flush=True)
    print(f"GPCR-cancer search rows       : {len(search_summary)}", flush=True)
    print(f"Candidate associations        : {len(long_df)}", flush=True)
    print(f"Unique PMIDs                  : {unique_pmids}", flush=True)
    print(f"Maximum metadata batches      : {expected_batches}", flush=True)
    print("Cross-file PMID agreement      : PASS", flush=True)
    print("Network requests               : 0", flush=True)

# Locate the saved response for a query or PMID.
def cache_path(pmid):
    return os.path.join(
        CACHE_DIR,
        f"{clean(pmid)}.json",
    )

# Read a PubMed article from the local cache.
def load_cached_article(pmid):
    path = cache_path(pmid)

    if not os.path.exists(path):
        return None

    try:
        with open(path, "r", encoding="utf-8") as fh:
            article = json.load(fh)
    except Exception:
        return None

    if clean(article.get("PMID")) != clean(pmid):
        return None

    return article

# Save parsed article metadata to the local cache.
def save_cached_article(article):
    path = cache_path(article["PMID"])

    payload = dict(article)
    payload["Cached_At_UTC"] = (
        datetime.now(timezone.utc).isoformat()
    )

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(
            payload,
            fh,
            indent=2,
            ensure_ascii=False,
        )

# Combine structured abstract sections into text.
def _abstract_from_container(container):
    try:
        parts = container.get("Abstract", {}).get("AbstractText", [])
        return clean(" ".join(str(x) for x in parts))
    except Exception:
        return ""

# Collect authors from PubMed metadata.
def _collect_author_names(obj):
    names = []

    # Traverse nested PubMed content to collect text.
    def walk(x):
        if isinstance(x, dict):
            last = clean(x.get("LastName", ""))
            initials = clean(x.get("Initials", ""))
            collective = clean(x.get("CollectiveName", ""))

            if last:
                names.append(f"{last} {initials}".strip())
                return

            if collective:
                names.append(collective)
                return

            for value in x.values():
                walk(value)

        elif isinstance(x, (list, tuple)):
            for value in x:
                walk(value)

    walk(obj)

    seen = set()
    out = []
    for name in names:
        key = norm(name)
        if key and key not in seen:
            seen.add(key)
            out.append(name)

    return "; ".join(out)

# Extract the publication title.
def extract_title(article):
    try:
        return clean(
            article["MedlineCitation"]["Article"]["ArticleTitle"]
        )
    except Exception:
        return ""

# Extract the publication abstract.
def extract_abstract(article):
    try:
        return _abstract_from_container(
            article["MedlineCitation"]["Article"]
        )
    except Exception:
        return ""

# Extract the journal name.
def extract_journal(article):
    try:
        return clean(
            article["MedlineCitation"]["Article"]["Journal"]["Title"]
        )
    except Exception:
        return ""

# Extract or validate a publication year.
def extract_year(article):
    try:
        pub_date = (
            article["MedlineCitation"]["Article"]["Journal"]
            ["JournalIssue"]["PubDate"]
        )
        return clean(
            pub_date.get("Year", pub_date.get("MedlineDate", ""))
        )
    except Exception:
        return ""

# Extract the publication author list.
def extract_authors(article):
    try:
        author_list = (
            article["MedlineCitation"]["Article"]
            .get("AuthorList", [])
        )
        return _collect_author_names(author_list)
    except Exception:
        return ""

# Extract fields from a PubMed article response.
def parse_pubmed_article(article):
    pmid = clean(article["MedlineCitation"]["PMID"])

    return {
        "PMID": pmid,
        "Title": extract_title(article),
        "Abstract": extract_abstract(article),
        "Journal": extract_journal(article),
        "Year": extract_year(article),
        "Authors": extract_authors(article),
        "Record_Type": "PUBMED_ARTICLE",
    }

# Extract fields from a PubMed book-article response.
def parse_pubmed_book_article(article):
    """
    PubmedBookArticle uses BookDocument instead of MedlineCitation.
    Treat it as a valid PubMed record and validate its title/abstract normally.
    """
    doc = article.get("BookDocument", {})
    book = doc.get("Book", {}) or {}

    pmid = clean(doc.get("PMID", ""))
    title = clean(doc.get("ArticleTitle", ""))
    book_title = clean(book.get("BookTitle", ""))

    if not title:
        title = book_title

    abstract = _abstract_from_container(doc)

    pub_date = book.get("PubDate", {}) or {}
    year = clean(
        pub_date.get("Year", pub_date.get("MedlineDate", ""))
    )

    authors = _collect_author_names(doc.get("AuthorList", []))

    return {
        "PMID": pmid,
        "Title": title,
        "Abstract": abstract,
        "Journal": book_title,
        "Year": year,
        "Authors": authors,
        "Record_Type": "PUBMED_BOOK_ARTICLE",
    }

# Fetch metadata for a batch of PMIDs.
def efetch_batch(pmids):
    pmids = [clean(x) for x in pmids if clean(x)]
    wait = INITIAL_RETRY_WAIT

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(pmids),
                retmode="xml",
            )
            record = Entrez.read(handle)
            handle.close()

            time.sleep(REQUEST_SLEEP)

            articles = []

            for raw in record.get("PubmedArticle", []):
                article = parse_pubmed_article(raw)
                if article["PMID"]:
                    save_cached_article(article)
                    articles.append(article)

            for raw in record.get("PubmedBookArticle", []):
                article = parse_pubmed_book_article(raw)
                if article["PMID"]:
                    save_cached_article(article)
                    articles.append(article)

            return articles

        except Exception as exc:
            if attempt == MAX_RETRIES:
                raise RuntimeError(
                    f"PubMed EFetch failed after {MAX_RETRIES} attempts "
                    f"for {len(pmids)} PMID(s). "
                    f"First PMIDs: {pmids[:10]}. "
                    f"Last error: {exc}"
                ) from exc

            print(
                f"    EFetch failed "
                f"(attempt {attempt}/{MAX_RETRIES}) "
                f"[{type(exc).__name__}]: {exc}",
                flush=True,
            )
            print(
                f"    Retrying in {wait:.1f} sec...",
                flush=True,
            )
            time.sleep(wait)
            wait *= 2

# Divide PMIDs into batches for metadata retrieval.
def chunks(items, size):
    for start in range(0, len(items), size):
        yield items[start:start + size]

# Read cached PubMed metadata and fetch missing records.
def fetch_all_metadata(unique_pmids):
    metadata = {}
    missing = []

    for pmid in unique_pmids:
        cached = load_cached_article(pmid)

        if cached:
            metadata[pmid] = {
                "PMID": clean(cached.get("PMID")),
                "Title": clean(cached.get("Title")),
                "Abstract": clean(cached.get("Abstract")),
                "Journal": clean(cached.get("Journal")),
                "Year": clean(cached.get("Year")),
                "Authors": clean(cached.get("Authors")),
                "Record_Type": clean(
                    cached.get("Record_Type", "PUBMED_ARTICLE")
                ),
            }
        else:
            missing.append(pmid)

    print(f"Metadata cache hits   : {len(metadata)}", flush=True)
    print(f"PMIDs to fetch        : {len(missing)}", flush=True)

    error_rows = []
    total_batches = (
        len(missing) + FETCH_BATCH_SIZE - 1
    ) // FETCH_BATCH_SIZE

    for batch_number, batch in enumerate(
        chunks(missing, FETCH_BATCH_SIZE),
        start=1,
    ):
        print(
            f"[Metadata batch {batch_number}/{total_batches}] "
            f"{len(batch)} PMIDs",
            flush=True,
        )

        try:
            returned = efetch_batch(batch)
            returned_ids = set()

            for article in returned:
                pmid = clean(article["PMID"])
                returned_ids.add(pmid)
                metadata[pmid] = article

            for pmid in batch:
                if pmid not in returned_ids:
                    error_rows.append({
                        "PMID": pmid,
                        "Error_Type": "NOT_RETURNED_BY_EFETCH",
                        "Error_Message": (
                            "PMID requested but no PubmedArticle "
                            "record was returned."
                        ),
                    })

        except Exception:
            print(
                "    Batch failed after retries. "
                "Retrying individual PMIDs.",
                flush=True,
            )

            for pmid in batch:
                try:
                    returned = efetch_batch([pmid])

                    if returned:
                        metadata[pmid] = returned[0]
                    else:
                        error_rows.append({
                            "PMID": pmid,
                            "Error_Type": "NOT_RETURNED_BY_EFETCH",
                            "Error_Message": (
                                "No PubmedArticle record returned."
                            ),
                        })

                except Exception as single_exc:
                    error_rows.append({
                        "PMID": pmid,
                        "Error_Type": "EFETCH_ERROR",
                        "Error_Message": str(single_exc),
                    })

    return metadata, error_rows

# Build lookups from the input records.
def build_lookups(aliases, cancers, search_summary):
    gpcr_lookup = {}

    for _, row in aliases.iterrows():
        uid = clean(row["UniProt_ID"])

        gpcr_lookup[uid] = {
            "Official_Gene": clean(row["Official_Gene"]),
            "Original_Gene": clean(row["Original_Gene"]),
            "Protein_Name": clean(row["Protein_Name"]),
            "KEEP_Terms": split_terms(row["Search_Terms"]),
            "CONTEXT_Terms": split_terms(
                row["Context_Required_Terms"]
            ),
        }

    cancer_lookup = {}

    for _, row in cancers.iterrows():
        code = clean(row["TCGA_Code"])

        cancer_lookup[code] = {
            "TCGA_Study_Name": clean(row["TCGA_Study_Name"]),
            "Cancer": clean(row["Cancer"]),
            "Accepted_Terms": split_terms(
                row["PubMed_Search_Terms"]
            ),
            "Rejected_Terms": REJECTED_CANCER_TERMS.get(
                code,
                [],
            ),
        }

    search_lookup = {}

    for _, row in search_summary.iterrows():
        sid = clean(row["Search_ID"])

        search_lookup[sid] = {
            "GPCR_KEEP_Terms": clean(
                row["GPCR_KEEP_Terms"]
            ),
            "GPCR_Context_Required_Terms": clean(
                row["GPCR_Context_Required_Terms"]
            ),
            "GPCR_Search_Terms": clean(
                row["GPCR_Search_Terms"]
            ),
            "Cancer_Search_Terms": clean(
                row["Cancer_Search_Terms"]
            ),
            "PubMed_Query": clean(row["PubMed_Query"]),
        }

    return gpcr_lookup, cancer_lookup, search_lookup

# Check candidate papers for the intended GPCR and cancer.
def validate_candidates(
    long_df,
    metadata,
    gpcr_lookup,
    cancer_lookup,
    search_lookup,
):
    rows = []
    total = len(long_df)

    for position, (_, candidate) in enumerate(
        long_df.iterrows(),
        start=1,
    ):
        if (
            position == 1
            or position % 5000 == 0
            or position == total
        ):
            print(
                f"Validating candidate {position}/{total}",
                flush=True,
            )

        pmid = clean(candidate["PMID"])
        uid = clean(candidate["UniProt_ID"])
        code = clean(candidate["TCGA_Code"])
        sid = clean(candidate["Search_ID"])

        gpcr = gpcr_lookup[uid]
        cancer = cancer_lookup[code]
        search = search_lookup[sid]

        article = metadata.get(pmid)

        base = {
            "UniProt_ID": uid,
            "Gene": clean(candidate["Gene"]),
            "Official_Gene": gpcr["Official_Gene"],
            "Original_Gene": gpcr["Original_Gene"],
            "Protein_Name": gpcr["Protein_Name"],
            "TCGA_Code": code,
            "TCGA_Study_Name": clean(candidate["TCGA_Study_Name"]),
            "Cancer": clean(candidate["Cancer"]),
            "PMID": pmid,
            "Search_ID": sid,
        }

        if article is None:
            rows.append({
                **base,
                "Title": "",
                "Abstract": "",
                "Journal": "",
                "Year": "",
                "Authors": "",
                "PubMed_Record_Type": "",
                "Cancer_Validation": "ERROR",
                "Cancer_Validation_Reason": "missing_pubmed_metadata",
                "Matched_Accepted_Cancer_Terms": "",
                "Matched_Rejected_Cancer_Terms": "",
                "GPCR_Validation": "ERROR",
                "GPCR_Validation_Reason": "missing_pubmed_metadata",
                "Matched_KEEP_GPCR_Terms": "",
                "Matched_KEEP_GPCR_Terms_Exact": "",
                "Matched_KEEP_GPCR_Terms_Flexible": "",
                "Matched_CONTEXT_REQUIRED_GPCR_Terms": "",
                "Matched_CONTEXT_REQUIRED_GPCR_Terms_Exact": "",
                "Matched_CONTEXT_REQUIRED_GPCR_Terms_Flexible": "",
                "Overall_Validation": "ERROR",
                "Final_KEEP_GPCR_Terms": ";".join(gpcr["KEEP_Terms"]),
                "Final_CONTEXT_REQUIRED_GPCR_Terms": ";".join(
                    gpcr["CONTEXT_Terms"]
                ),
                "Accepted_Cancer_Terms": ";".join(
                    cancer["Accepted_Terms"]
                ),
                "Rejected_Cancer_Terms": ";".join(
                    cancer["Rejected_Terms"]
                ),
                "PubMed_Query": search["PubMed_Query"],
            })
            continue

        cancer_decision, cancer_reason, accepted, rejected = (
            validate_cancer(
                code=code,
                title=article["Title"],
                abstract=article["Abstract"],
                accepted_terms=cancer["Accepted_Terms"],
            )
        )

        gpcr_result = validate_gpcr(
            title=article["Title"],
            abstract=article["Abstract"],
            keep_terms=gpcr["KEEP_Terms"],
            context_terms=gpcr["CONTEXT_Terms"],
        )

        overall = overall_validation(
            cancer_decision,
            gpcr_result["Decision"],
        )

        rows.append({
            **base,
            "Title": clean(article["Title"]),
            "Abstract": clean(article["Abstract"]),
            "Journal": clean(article["Journal"]),
            "Year": clean(article["Year"]),
            "Authors": clean(article["Authors"]),
            "PubMed_Record_Type": clean(
                article.get("Record_Type", "PUBMED_ARTICLE")
            ),

            "Cancer_Validation": cancer_decision,
            "Cancer_Validation_Reason": cancer_reason,
            "Matched_Accepted_Cancer_Terms": ";".join(accepted),
            "Matched_Rejected_Cancer_Terms": ";".join(rejected),

            "GPCR_Validation": gpcr_result["Decision"],
            "GPCR_Validation_Reason": gpcr_result["Reason"],

            "Matched_KEEP_GPCR_Terms": ";".join(
                gpcr_result["KEEP_Matches"]
            ),
            "Matched_KEEP_GPCR_Terms_Exact": ";".join(
                gpcr_result["KEEP_Exact"]
            ),
            "Matched_KEEP_GPCR_Terms_Flexible": ";".join(
                gpcr_result["KEEP_Flexible"]
            ),

            "Matched_CONTEXT_REQUIRED_GPCR_Terms": ";".join(
                gpcr_result["Context_Matches"]
            ),
            "Matched_CONTEXT_REQUIRED_GPCR_Terms_Exact": ";".join(
                gpcr_result["Context_Exact"]
            ),
            "Matched_CONTEXT_REQUIRED_GPCR_Terms_Flexible": ";".join(
                gpcr_result["Context_Flexible"]
            ),

            "Overall_Validation": overall,

            "Final_KEEP_GPCR_Terms": ";".join(
                gpcr["KEEP_Terms"]
            ),
            "Final_CONTEXT_REQUIRED_GPCR_Terms": ";".join(
                gpcr["CONTEXT_Terms"]
            ),
            "Accepted_Cancer_Terms": ";".join(
                cancer["Accepted_Terms"]
            ),
            "Rejected_Cancer_Terms": ";".join(
                cancer["Rejected_Terms"]
            ),
            "PubMed_Query": search["PubMed_Query"],
        })

    return pd.DataFrame(rows)

# Summarize validation for each GPCR and cancer pair.
def build_pair_summary(search_summary, validated):
    grouped = (
        validated
        .groupby("Search_ID", as_index=False)
        .agg(
            Candidate_PMIDs=("PMID", "size"),
            Cancer_KEEP=(
                "Cancer_Validation",
                lambda s: int((s == "KEEP").sum()),
            ),
            GPCR_KEEP=(
                "GPCR_Validation",
                lambda s: int((s == "KEEP").sum()),
            ),
            GPCR_CONTEXT_REVIEW=(
                "GPCR_Validation",
                lambda s: int((s == "CONTEXT_REVIEW").sum()),
            ),
            Overall_KEEP=(
                "Overall_Validation",
                lambda s: int((s == "KEEP").sum()),
            ),
            Overall_CONTEXT_REVIEW=(
                "Overall_Validation",
                lambda s: int((s == "CONTEXT_REVIEW").sum()),
            ),
            Overall_REJECT=(
                "Overall_Validation",
                lambda s: int((s == "REJECT").sum()),
            ),
            Validation_Errors=(
                "Overall_Validation",
                lambda s: int((s == "ERROR").sum()),
            ),
        )
    )

    summary = search_summary.merge(
        grouped,
        on="Search_ID",
        how="left",
        validate="one_to_one",
    )

    fill_cols = [
        "Candidate_PMIDs",
        "Cancer_KEEP",
        "GPCR_KEEP",
        "GPCR_CONTEXT_REVIEW",
        "Overall_KEEP",
        "Overall_CONTEXT_REVIEW",
        "Overall_REJECT",
        "Validation_Errors",
    ]

    for col in fill_cols:
        summary[col] = (
            summary[col]
            .fillna(0)
            .astype(int)
        )

    return summary

# Write report for the current record.
def write_report(
    long_df,
    metadata_df,
    validated,
    kept,
    context_review,
    pair_summary,
    error_df,
    aliases,
    cancers,
    elapsed_minutes,
):
    candidate_rows = len(long_df)
    unique_pmids = long_df["PMID"].nunique()

    overall_counts = (
        validated["Overall_Validation"]
        .value_counts()
        .to_dict()
    )

    cancer_counts = (
        validated["Cancer_Validation"]
        .value_counts()
        .to_dict()
    )

    gpcr_counts = (
        validated["GPCR_Validation"]
        .value_counts()
        .to_dict()
    )

    candidate_pmid_set = set(long_df["PMID"].map(clean))
    metadata_pmid_set = set(metadata_df["PMID"].map(clean))

    book_record_count = int(
        metadata_df["Record_Type"]
        .map(clean)
        .eq("PUBMED_BOOK_ARTICLE")
        .sum()
    )

    checks = {
        "832 GPCRs": len(aliases) == EXPECTED_GPCR_COUNT,
        "33 TCGA cancers": len(cancers) == EXPECTED_CANCER_COUNT,
        "27,456 GPCR-cancer summary rows": (
            len(pair_summary) == EXPECTED_PAIR_COUNT
        ),
        "All unique PMIDs have metadata": (
            metadata_pmid_set == candidate_pmid_set
        ),
        "No duplicate metadata rows": (
            not metadata_df["PMID"].map(clean).duplicated().any()
        ),
        "Candidate row count preserved": (
            len(validated) == candidate_rows
        ),
        "No duplicate validated associations": (
            not validated.duplicated(
                subset=["UniProt_ID", "TCGA_Code", "PMID"]
            ).any()
        ),
        "All candidate rows resolved": (
            sum(overall_counts.get(x, 0) for x in [
                "KEEP",
                "CONTEXT_REVIEW",
                "REJECT",
                "ERROR",
            ])
            == candidate_rows
        ),
        "No validation errors": (
            overall_counts.get("ERROR", 0) == 0
        ),
        "No metadata fetch errors": error_df.empty,
        "Kept rows equal overall KEEP": (
            len(kept) == overall_counts.get("KEEP", 0)
        ),
        "Context-review rows equal overall CONTEXT_REVIEW": (
            len(context_review)
            == overall_counts.get("CONTEXT_REVIEW", 0)
        ),
        "Pair-summary candidate count preserved": (
            int(pair_summary["Candidate_PMIDs"].sum())
            == candidate_rows
        ),
    }

    failed = [label for label, passed in checks.items() if not passed]

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 07 Integrity Report",
        "=" * 76,
        "",
        "INPUT ACCOUNTING",
        "-" * 76,
        f"Candidate associations           : {candidate_rows}",
        f"Unique candidate PMIDs           : {unique_pmids}",
        f"Metadata rows recovered          : {len(metadata_df)}",
        f"PubMed book/chapter records      : {book_record_count}",
        f"GPCR x cancer pair rows          : {len(pair_summary)}",
        "",
        "CANCER VALIDATION",
        "-" * 76,
        f"KEEP                             : {cancer_counts.get('KEEP', 0)}",
        f"REJECT                           : {cancer_counts.get('REJECT', 0)}",
        f"ERROR                            : {cancer_counts.get('ERROR', 0)}",
        "",
        "GPCR VALIDATION",
        "-" * 76,
        f"KEEP via final KEEP term         : {gpcr_counts.get('KEEP', 0)}",
        f"CONTEXT_REVIEW                   : {gpcr_counts.get('CONTEXT_REVIEW', 0)}",
        f"REJECT no active term in text    : {gpcr_counts.get('REJECT', 0)}",
        f"ERROR                            : {gpcr_counts.get('ERROR', 0)}",
        "",
        "OVERALL VALIDATION",
        "-" * 76,
        f"KEEP                             : {overall_counts.get('KEEP', 0)}",
        f"CONTEXT_REVIEW                   : {overall_counts.get('CONTEXT_REVIEW', 0)}",
        f"REJECT                           : {overall_counts.get('REJECT', 0)}",
        f"ERROR                            : {overall_counts.get('ERROR', 0)}",
        "",
        "INTEGRITY CHECKS",
        "-" * 76,
    ]

    for label, passed in checks.items():
        lines.append(
            f"{label}: {'PASS' if passed else 'FAIL'}"
        )

    lines += [
        "",
        "METHOD NOTE",
        "-" * 76,
        (
            "Cancer validation is performed on Title + Abstract using the "
            "preserved accepted/rejected cancer-term logic. GPCR validation "
            "is performed independently on the same text. A paper is accepted "
            "automatically only when at least one final KEEP GPCR term is "
            "present. Papers supported only by CONTEXT_REQUIRED GPCR terms are "
            "not accepted automatically and are exported to a separate "
            "context-review file. PubMed book/chapter citations returned as "
            "PubmedBookArticle are parsed through BookDocument and validated "
            "under the same Title + Abstract rules rather than being counted "
            "as metadata-fetch errors."
        ),
        "",
        "TECHNICAL SETTINGS",
        "-" * 76,
        f"Metadata fetch batch size        : {FETCH_BATCH_SIZE}",
        f"Network timeout seconds          : {SOCKET_TIMEOUT_SECONDS}",
        f"NCBI API key                     : {'YES' if API_KEY else 'NO'}",
        f"Elapsed minutes                  : {elapsed_minutes:.1f}",
        f"Integrity failures               : {len(failed)}",
    ]

    if failed:
        lines += [
            "",
            "FAILED CHECKS",
            "-" * 76,
        ]
        lines.extend(f"- {x}" for x in failed)
    else:
        lines += [
            "",
            "ALL STEP 07 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return failed

# Run this stage and write its output files.
def main():
    args = parse_args()
    start_time = time.time()

    print("=" * 76, flush=True)
    print(
        "GPCR Cancer Pipeline - Final Fresh Run Step 07",
        flush=True,
    )
    print(
        "Fetch Metadata + Validate Cancer and GPCR Support",
        flush=True,
    )
    print("=" * 76, flush=True)

    try:
        long_df, search_summary, aliases, cancers = load_inputs()
    except Exception as exc:
        print(
            f"INPUT VALIDATION FAILED:\n{exc}",
            flush=True,
        )
        sys.exit(1)

    long_df["PMID"] = (
        long_df["PMID"]
        .astype(str)
        .map(clean)
    )

    unique_pmids = sorted(
        set(long_df["PMID"]),
        key=pmid_sort_key,
    )

    print(f"Candidate rows         : {len(long_df):,}", flush=True)
    print(f"Unique PMIDs           : {len(unique_pmids):,}", flush=True)
    print(f"GPCR x cancer pairs    : {len(search_summary):,}", flush=True)
    print(f"TCGA cancers           : {len(cancers)}", flush=True)
    print(f"Metadata batch size    : {FETCH_BATCH_SIZE}", flush=True)
    print(
        f"Network timeout        : {SOCKET_TIMEOUT_SECONDS} sec",
        flush=True,
    )
    print(f"Metadata cache         : {CACHE_DIR}", flush=True)
    print("=" * 76, flush=True)

    if args.preflight:
        run_preflight(
            long_df,
            search_summary,
            aliases,
            cancers,
        )
        return

    if not Entrez.email:
        raise SystemExit(
            "NCBI email is required.\n"
            'Run: export NCBI_EMAIL="your_email@example.com"'
        )

    os.makedirs(CACHE_DIR, exist_ok=True)

    metadata, fetch_errors = fetch_all_metadata(unique_pmids)

    metadata_df = pd.DataFrame(
        [
            metadata[pmid]
            for pmid in sorted(
                metadata,
                key=pmid_sort_key,
            )
        ],
        columns=[
            "PMID",
            "Title",
            "Abstract",
            "Journal",
            "Year",
            "Authors",
            "Record_Type",
        ],
    )

    metadata_df.to_csv(
        OUT_METADATA,
        index=False,
        na_rep="",
    )

    error_df = pd.DataFrame(
        fetch_errors,
        columns=[
            "PMID",
            "Error_Type",
            "Error_Message",
        ],
    )

    error_df.to_csv(
        OUT_ERRORS,
        index=False,
        na_rep="",
    )

    gpcr_lookup, cancer_lookup, search_lookup = build_lookups(
        aliases,
        cancers,
        search_summary,
    )

    validated = validate_candidates(
        long_df,
        metadata,
        gpcr_lookup,
        cancer_lookup,
        search_lookup,
    )

    validated.to_csv(
        OUT_ALL,
        index=False,
        na_rep="",
    )

    kept = validated[
        validated["Overall_Validation"].eq("KEEP")
    ].copy()

    kept.to_csv(
        OUT_KEPT,
        index=False,
        na_rep="",
    )

    context_review = validated[
        validated["Overall_Validation"].eq("CONTEXT_REVIEW")
    ].copy()

    context_review.to_csv(
        OUT_CONTEXT,
        index=False,
        na_rep="",
    )

    pair_summary = build_pair_summary(
        search_summary,
        validated,
    )

    pair_summary.to_csv(
        OUT_SUMMARY,
        index=False,
        na_rep="",
    )

    elapsed_minutes = (
        time.time() - start_time
    ) / 60.0

    failed = write_report(
        long_df,
        metadata_df,
        validated,
        kept,
        context_review,
        pair_summary,
        error_df,
        aliases,
        cancers,
        elapsed_minutes,
    )

    overall_counts = (
        validated["Overall_Validation"]
        .value_counts()
        .to_dict()
    )

    print()
    print("=" * 76, flush=True)

    if not failed:
        print(
            "STEP 07 COMPLETE - ALL INTEGRITY CHECKS PASSED",
            flush=True,
        )
        exit_code = 0
    else:
        print(
            "STEP 07 COMPLETE WITH UNRESOLVED ISSUES",
            flush=True,
        )
        exit_code = 2

    print("=" * 76, flush=True)
    print(
        f"Metadata recovered    : "
        f"{len(metadata_df):,}/{len(unique_pmids):,}",
        flush=True,
    )
    print(
        f"Candidate associations: {len(validated):,}",
        flush=True,
    )
    print(
        f"Overall KEEP          : "
        f"{overall_counts.get('KEEP', 0):,}",
        flush=True,
    )
    print(
        f"CONTEXT_REVIEW        : "
        f"{overall_counts.get('CONTEXT_REVIEW', 0):,}",
        flush=True,
    )
    print(
        f"REJECT                : "
        f"{overall_counts.get('REJECT', 0):,}",
        flush=True,
    )
    print(
        f"ERROR                 : "
        f"{overall_counts.get('ERROR', 0):,}",
        flush=True,
    )
    print(
        f"Elapsed               : {elapsed_minutes:.1f} min",
        flush=True,
    )
    print()
    print("Saved:", flush=True)
    print(f"  {OUT_METADATA}", flush=True)
    print(f"  {OUT_ALL}", flush=True)
    print(f"  {OUT_KEPT}", flush=True)
    print(f"  {OUT_CONTEXT}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_ERRORS}", flush=True)
    print(f"  {OUT_REPORT}", flush=True)
    print("=" * 76, flush=True)

    sys.exit(exit_code)

if __name__ == "__main__":
    main()
