#!/usr/bin/env python3
# Purpose: Collect GPCR names and aliases from UniProt.
# Input: gpcrpendb_results_1782027334.2077546.tsv
# Output: GPCR_Alias_Table.csv

import json
import os
import re
import time
from collections import Counter

import pandas as pd
import requests

INPUT_FILE = "gpcrpendb_results_1782027334.2077546.tsv"
OUTPUT_FILE = "GPCR_Alias_Table.csv"

CACHE_DIR = "uniprot_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

REQUEST_TIMEOUT = 30
MAX_RETRIES = 5
INITIAL_WAIT = 1.0

# Normalize an alias while keeping its display spelling.
def normalize_alias(term):
    if term is None:
        return None
    term = re.sub(r"\s+", " ", str(term)).strip()
    if not term or term.lower() == "nan" or term == "-":
        return None
    return term

# Locate the saved response for a query or PMID.
def cache_path(acc):
    return os.path.join(CACHE_DIR, f"{acc}.json")

# Read the saved UniProt record or fetch it when absent.
def get_uniprot(acc):
    p = cache_path(acc)
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as fh:
            return json.load(fh), True

    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"

    wait = INITIAL_WAIT
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                data = r.json()
                with open(p, "w", encoding="utf-8") as fh:
                    json.dump(data, fh, indent=2, ensure_ascii=False)
                return data, False
            raise RuntimeError(f"HTTP {r.status_code}")
        except Exception:
            if attempt == MAX_RETRIES:
                return None, False
            time.sleep(wait)
            wait *= 2

# Add a cleaned nonblank value to the collection.
def add(term, bag, stats):
    stats["raw"] += 1
    t = normalize_alias(term)
    if t is None:
        stats["blank"] += 1
        return
    if t in bag:
        stats["dup"] += 1
    bag.add(t)

# Extract gene and protein names from a UniProt record.
def parse_uniprot(data, stats):
    terms = set()
    official = ""
    protein_name = ""

    for gene in data.get("genes", []):
        if "geneName" in gene:
            official = gene["geneName"]["value"]
            add(official, terms, stats)
        for key in ("synonyms","orderedLocusNames","orfNames"):
            for x in gene.get(key, []):
                add(x["value"], terms, stats)

    pdesc = data.get("proteinDescription", {})

    rn = pdesc.get("recommendedName", {})
    if "fullName" in rn:
        protein_name = rn["fullName"]["value"]
        add(protein_name, terms, stats)

    for s in rn.get("shortNames", []):
        add(s["value"], terms, stats)

    for alt in pdesc.get("alternativeNames", []):
        if "fullName" in alt:
            add(alt["fullName"]["value"], terms, stats)
        for s in alt.get("shortNames", []):
            add(s["value"], terms, stats)

    return official, protein_name, sorted(terms, key=str.casefold)

start = time.time()

print("="*60)
print("GPCR Cancer Pipeline v1.0")
print("Step 00 - Build UniProt Alias Table")
print("="*60)

df = pd.read_csv(INPUT_FILE, sep="\t")
human = df[df["Species"].astype(str).str.strip() == "sapiens"].copy()

gpcrs = (
    human[["ID","Gene"]]
    .dropna()
    .drop_duplicates("ID")
    .sort_values("Gene")
    .reset_index(drop=True)
)

print(f"Human GPCRs: {len(gpcrs)}")

stats = Counter()
rows = []
failed = 0
cached = 0

for i, row in gpcrs.iterrows():
    acc = str(row["ID"]).strip()
    original = normalize_alias(row["Gene"]) or ""

    print(f"[{i+1:3d}/{len(gpcrs)}] {acc} {original}")

    data, from_cache = get_uniprot(acc)
    if from_cache:
        cached += 1

    if data is None:
        failed += 1
        rows.append({
            "UniProt_ID": acc,
            "Original_Gene": original,
            "Official_Gene": "",
            "Protein_Name": "",
            "Search_Terms": original
        })
        continue

    official, protein, terms = parse_uniprot(data, stats)

    add(original, set(terms), stats)
    if original and original not in terms:
        terms.append(original)

    terms = sorted(set(terms), key=str.casefold)

    rows.append({
        "UniProt_ID": acc,
        "Original_Gene": original,
        "Official_Gene": official,
        "Protein_Name": protein,
        "Search_Terms": ";".join(terms)
    })

out = pd.DataFrame(rows)
out.to_csv(OUTPUT_FILE, index=False)

elapsed = time.time() - start

print("\nFinished.")
print(f"Output file               : {OUTPUT_FILE}")
print(f"Cached UniProt records    : {cached}")
print(f"Failed downloads          : {failed}")
print(f"Raw aliases seen          : {stats['raw']}")
print(f"Duplicate aliases removed : {stats['dup']}")
print(f"Blank aliases removed     : {stats['blank']}")
print(f"Elapsed time              : {elapsed:0.1f} sec")

if failed:
    raise SystemExit(
        f"Incomplete UniProt collection: {failed} download(s) failed. "
        "Repeat the prepare command to reuse successful cached records and retry."
    )
