#!/usr/bin/env python3
# Purpose: Collect additional GPCR terminology from HGNC.
# Input: GPCR_Alias_Table.csv
# Output: HGNC_Alias_Table.csv; GPCR_Alias_Source_Comparison.csv

import os
import re
import sys
from collections import OrderedDict

import pandas as pd
import requests

INPUT_FILE = "GPCR_Alias_Table.csv"

HGNC_FILE = "hgnc_complete_set.txt"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/"
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
)

HGNC_OUTPUT = "HGNC_Alias_Table.csv"
COMPARISON_OUTPUT = "GPCR_Alias_Source_Comparison.csv"

REQUEST_TIMEOUT = 120

# Standardize whitespace in a text value.
def clean_text(value):
    if value is None:
        return ""
    value = str(value).strip()
    if not value or value.lower() == "nan":
        return ""
    return re.sub(r"\s+", " ", value)

# Split an HGNC pipe-delimited name field.
def split_pipe(value):
    """HGNC multi-valued TSV fields are pipe-delimited."""
    value = clean_text(value)
    if not value:
        return []
    return [clean_text(x) for x in value.split("|") if clean_text(x)]

# Split semicolon-delimited values.
def split_semicolon(value):
    value = clean_text(value)
    if not value:
        return []
    return [clean_text(x) for x in value.split(";") if clean_text(x)]

# Deduplicate terms while preserving display spelling.
def unique_preserve_case(terms):
    """
    Case-insensitive deduplication while preserving the first observed spelling.
    """
    seen = OrderedDict()
    for term in terms:
        term = clean_text(term)
        if not term:
            continue
        key = term.casefold()
        if key not in seen:
            seen[key] = term
    return list(seen.values())

# Join unique terms in their selected order.
def joined(terms):
    return ";".join(unique_preserve_case(terms))

# Download an HGNC source table when requested.
def download_hgnc():
    if os.path.exists(HGNC_FILE) and os.path.getsize(HGNC_FILE) > 1000:
        print(f"Using existing HGNC file: {HGNC_FILE}")
        return

    print("Downloading current HGNC complete set once...")
    r = requests.get(HGNC_URL, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()

    with open(HGNC_FILE, "wb") as fh:
        fh.write(r.content)

    print(f"Saved: {HGNC_FILE}")

print("=" * 72)
print("GPCR Cancer Pipeline - Fresh Final Run")
print("Step 00_1 - HGNC Alias Collection and UniProt Comparison")
print("=" * 72)

if not os.path.exists(INPUT_FILE):
    sys.exit(f"ERROR: Missing input file: {INPUT_FILE}")

download_hgnc()

gpcr = pd.read_csv(INPUT_FILE, dtype=str).fillna("")
hgnc = pd.read_csv(HGNC_FILE, sep="\t", dtype=str, low_memory=False).fillna("")

required_gpcr = {"UniProt_ID", "Original_Gene", "Official_Gene", "Protein_Name", "Search_Terms"}
missing = required_gpcr - set(gpcr.columns)
if missing:
    sys.exit(f"ERROR: {INPUT_FILE} is missing columns: {sorted(missing)}")

required_hgnc = {
    "hgnc_id", "symbol", "name", "status",
    "alias_symbol", "prev_symbol", "alias_name", "prev_name", "uniprot_ids"
}
missing_hgnc = required_hgnc - set(hgnc.columns)
if missing_hgnc:
    sys.exit(f"ERROR: HGNC file is missing columns: {sorted(missing_hgnc)}")

hgnc = hgnc[hgnc["status"].astype(str).str.strip().eq("Approved")].copy()

by_uniprot = {}
for idx, row in hgnc.iterrows():
    for acc in split_pipe(row["uniprot_ids"]):
        by_uniprot.setdefault(acc, []).append(idx)

by_symbol = {}
for idx, row in hgnc.iterrows():
    symbol = clean_text(row["symbol"])
    if symbol:
        by_symbol.setdefault(symbol.casefold(), []).append(idx)

hgnc_rows = []
compare_rows = []

n_uniprot = 0
n_symbol = 0
n_unmatched = 0
n_ambiguous = 0

for i, row in gpcr.iterrows():
    acc = clean_text(row["UniProt_ID"])
    original_gene = clean_text(row["Original_Gene"])
    official_gene = clean_text(row["Official_Gene"])

    match_indices = by_uniprot.get(acc, [])
    match_method = ""

    if len(match_indices) == 1:
        match_method = "UNIPROT_ID"
        n_uniprot += 1

    elif len(match_indices) > 1:
        symbol_matches = [
            idx for idx in match_indices
            if clean_text(hgnc.loc[idx, "symbol"]).casefold() == official_gene.casefold()
        ]
        if len(symbol_matches) == 1:
            match_indices = symbol_matches
            match_method = "UNIPROT_ID+SYMBOL"
            n_uniprot += 1
        else:
            match_method = "AMBIGUOUS_UNIPROT"
            n_ambiguous += 1

    else:
        symbol_matches = by_symbol.get(official_gene.casefold(), []) if official_gene else []
        if len(symbol_matches) == 1:
            match_indices = symbol_matches
            match_method = "OFFICIAL_SYMBOL_FALLBACK"
            n_symbol += 1
        elif len(symbol_matches) > 1:
            match_indices = symbol_matches
            match_method = "AMBIGUOUS_SYMBOL"
            n_ambiguous += 1
        else:
            match_indices = []
            match_method = "UNMATCHED"
            n_unmatched += 1

    if len(match_indices) == 1:
        h = hgnc.loc[match_indices[0]]

        approved_symbol = clean_text(h["symbol"])
        approved_name = clean_text(h["name"])
        alias_symbols = split_pipe(h["alias_symbol"])
        previous_symbols = split_pipe(h["prev_symbol"])
        alias_names = split_pipe(h["alias_name"])
        previous_names = split_pipe(h["prev_name"])

        hgnc_terms = unique_preserve_case(
            [approved_symbol, approved_name]
            + alias_symbols
            + previous_symbols
            + alias_names
            + previous_names
        )

        hgnc_id = clean_text(h["hgnc_id"])
    else:
        approved_symbol = ""
        approved_name = ""
        alias_symbols = []
        previous_symbols = []
        alias_names = []
        previous_names = []
        hgnc_terms = []
        hgnc_id = ""

    uniprot_terms = unique_preserve_case(split_semicolon(row["Search_Terms"]))

    up_map = {x.casefold(): x for x in uniprot_terms}
    hg_map = {x.casefold(): x for x in hgnc_terms}

    shared_keys = sorted(set(up_map) & set(hg_map))
    up_only_keys = sorted(set(up_map) - set(hg_map))
    hg_only_keys = sorted(set(hg_map) - set(up_map))

    shared = [up_map[k] for k in shared_keys]
    uniprot_only = [up_map[k] for k in up_only_keys]
    hgnc_only = [hg_map[k] for k in hg_only_keys]

    combined = unique_preserve_case(uniprot_terms + hgnc_terms)

    hgnc_rows.append({
        "UniProt_ID": acc,
        "Original_Gene": original_gene,
        "Step00_Official_Gene": official_gene,
        "HGNC_ID": hgnc_id,
        "HGNC_Approved_Symbol": approved_symbol,
        "HGNC_Approved_Name": approved_name,
        "HGNC_Alias_Symbols": joined(alias_symbols),
        "HGNC_Previous_Symbols": joined(previous_symbols),
        "HGNC_Alias_Names": joined(alias_names),
        "HGNC_Previous_Names": joined(previous_names),
        "HGNC_All_Terms": joined(hgnc_terms),
        "HGNC_Match_Method": match_method,
    })

    compare_rows.append({
        "UniProt_ID": acc,
        "Original_Gene": original_gene,
        "UniProt_Official_Gene": official_gene,
        "HGNC_Approved_Symbol": approved_symbol,
        "HGNC_Match_Method": match_method,
        "UniProt_Term_Count": len(uniprot_terms),
        "HGNC_Term_Count": len(hgnc_terms),
        "Shared_Term_Count": len(shared),
        "UniProt_Only_Count": len(uniprot_only),
        "HGNC_Only_Count": len(hgnc_only),
        "Shared_Terms": joined(shared),
        "UniProt_Only_Terms": joined(uniprot_only),
        "HGNC_Only_Terms": joined(hgnc_only),
        "Combined_Candidate_Terms": joined(combined),
    })

    print(
        f"[{i+1:3d}/{len(gpcr)}] {acc:10s} "
        f"{official_gene:15s} HGNC={match_method:24s} "
        f"UP={len(uniprot_terms):2d} HGNC={len(hgnc_terms):2d} "
        f"new={len(hgnc_only):2d}"
    )

pd.DataFrame(hgnc_rows).to_csv(HGNC_OUTPUT, index=False)
pd.DataFrame(compare_rows).to_csv(COMPARISON_OUTPUT, index=False)

print("\n" + "=" * 72)
print("STEP 00_1 COMPLETE")
print("=" * 72)
print(f"GPCR rows                 : {len(gpcr)}")
print(f"Matched by UniProt        : {n_uniprot}")
print(f"Matched by symbol fallback: {n_symbol}")
print(f"Ambiguous matches         : {n_ambiguous}")
print(f"Unmatched                 : {n_unmatched}")
print(f"HGNC table                : {HGNC_OUTPUT}")
print(f"Comparison table          : {COMPARISON_OUTPUT}")
print("\nCode S2 collects names; reviewed alias decisions are applied in Code S6.")
print("Next: inspect HGNC-only and UniProt-only terms before defining alias rules.")
