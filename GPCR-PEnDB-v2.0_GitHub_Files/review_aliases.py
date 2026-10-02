#!/usr/bin/env python3
"""Record optional keyboard decisions for a pending refresh alias-review queue.

Place beside run_refresh.py. Uses only the Python standard library and saved
evidence in --work; the curator supplies all decisions and comments.
"""

import argparse
from collections import defaultdict
import csv
import json
import os
from pathlib import Path
import tempfile
import textwrap


KEYS = ("UniProt_ID", "Alias")
QUEUE = "Alias_Manual_Review_Queue.csv"
SHEET = "REVIEW_DECISIONS_TO_COMPLETE.csv"
CHOICES = {
    "K": "KEEP", "KEEP": "KEEP",
    "I": "IGNORE", "IGNORE": "IGNORE",
    "C": "CONTEXT_REQUIRED", "CONTEXT_REQUIRED": "CONTEXT_REQUIRED",
}
ALLOWED = set(CHOICES.values())


def clean(value):
    return " ".join(str(value or "").split())


def key(row):
    return tuple(clean(row[column]) for column in KEYS)


def read_csv(path, required):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate CSV column names: {path}")
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f"{path.name} lacks columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Malformed CSV row: {path}")
    return fields, rows


def index_rows(rows, name):
    indexed = {}
    for row in rows:
        identifier = key(row)
        if not all(identifier) or identifier in indexed:
            raise ValueError(f"Incomplete or duplicate UniProt_ID + Alias in {name}")
        indexed[identifier] = row
    return indexed


def save_csv(path, fields, rows):
    """Replace the sheet only after a complete CSV has been written."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def completed(row):
    return clean(row["Final_Decision"]).upper() in ALLOWED and bool(clean(row["Curator_Comment"]))


def display(label, value):
    value = clean(value)
    if value:
        print(textwrap.fill(f"{label}: {value}", width=100, subsequent_indent="  "))


def show_alias(row, number, total):
    print(f"\n[{number}/{total}] {row.get('Official_Gene', '')} | {row['UniProt_ID']} | {row['Alias']}")
    fields = (
        ("Receptor", "Protein_Name"), ("Source", "Source_Provenance"),
        ("Term types", "Term_Types"), ("Triage reason", "Triage_Reason"),
        ("Evidence note", "Triage_Evidence_Note"),
        ("Intrinsic ambiguity", "Flag_X_Reasons"),
        ("Cross-GPCR flags", "Flag_Y_Reasons"),
        ("Retrieval flags", "Flag_Z_Reasons"),
        ("Other review reasons", "Other_Review_Reasons"),
        ("Other receptor names in collision", "Collision_Genes"),
        ("Total alias PMIDs", "Total_PMIDs"),
        ("Independent reference PMIDs", "Independent_Reference_PMIDs"),
        ("Overlapping PMIDs", "Overlap"), ("Alias-only PMIDs", "Unique_PMIDs"),
        ("Alias PMID list truncated", "PMID_List_Truncated"),
        ("Reference PMID list truncated", "Reference_PMID_List_Truncated"),
        ("Reference terms", "Reference_Terms"),
    )
    for label, field in fields:
        display(label, row.get(field, ""))


def show_samples(samples, cache, full=False):
    if not samples:
        print("No sampled evidence is available for this alias in Alias_Context_Review.csv.")
        return
    print(f"\nSaved sampled records: {len(samples)}")
    for row in samples:
        pmid = clean(row.get("PMID", ""))
        record = cache.get(pmid, {})
        if not isinstance(record, dict):
            record = {}
        print(f"\nPMID {pmid} | {row.get('Sample_Group', '')} | https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
        display("Title", row.get("Title", "") or record.get("Title", ""))
        if full:
            display("Abstract", record.get("Abstract", "") or "Unavailable in the saved cache.")
        else:
            display("Matched context", row.get("Matched_Context", ""))
            display("Parenthetical expansions", row.get("Parenthetical_Expansions", ""))
            display("Cooccurring reference terms", row.get("Cooccurring_Reference_Terms", ""))
            display("Independent support flag", row.get("Independent_Support_In_Record", ""))


def review(work):
    state_path = work / ".refresh_run_state.json"
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if "6" in state.get("completed", {}) or "alias" in state.get("reviews", {}):
            raise ValueError("Alias decisions were already applied. Start a new run to change them.")
        if "5" not in state.get("completed", {}):
            raise ValueError("Complete prepare through Code S5 before reviewing aliases.")

    _, queue = read_csv(work / QUEUE, KEYS)
    expected = index_rows(queue, QUEUE)
    target = work / SHEET
    if target.exists():
        fields, rows = read_csv(target, (*KEYS, "Final_Decision", "Curator_Comment"))
    else:
        fields = [*KEYS]
        if queue and "Official_Gene" in queue[0]:
            fields.append("Official_Gene")
        fields.extend(["Final_Decision", "Curator_Comment"])
        rows = [{column: row.get(column, "") if column in (*KEYS, "Official_Gene") else ""
                 for column in fields} for row in queue]
    decisions = index_rows(rows, SHEET)
    if set(decisions) != set(expected):
        raise ValueError("Decision identifiers must match the pending review queue exactly.")
    for row in rows:
        decision = clean(row["Final_Decision"]).upper()
        if decision and decision not in ALLOWED:
            raise ValueError(f"Invalid existing decision for {key(row)}: {decision}")

    samples = defaultdict(list)
    context_path = work / "Alias_Context_Review.csv"
    if context_path.exists():
        _, evidence = read_csv(context_path, (*KEYS, "PMID"))
        for row in evidence:
            if key(row) in expected:
                samples[key(row)].append(row)
    cache_path = work / "pubmed_context_cache.json"
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
    if not isinstance(cache, dict):
        raise ValueError("pubmed_context_cache.json must contain a PMID-to-record object.")
    if not target.exists():
        save_csv(target, fields, rows)

    print(f"Decision CSV: {target}")
    print(f"Already complete: {sum(completed(row) for row in rows)}/{len(rows)}")
    print("K = KEEP; I = IGNORE; C = CONTEXT_REQUIRED; A = saved full abstracts; S = skip; Q = quit.")
    print("Every decision requires a curator comment. Each completed decision is saved immediately.")
    try:
        for number, alias in enumerate(queue, 1):
            row = decisions[key(alias)]
            if completed(row):
                continue
            show_alias(alias, number, len(queue))
            show_samples(samples[key(alias)], cache)
            display("Existing decision", row["Final_Decision"])
            display("Existing comment", row["Curator_Comment"])
            while True:
                choice = clean(input("Decision [K/I/C, A, S, Q]: ")).upper()
                if choice == "Q":
                    print("Paused. Completed decisions are saved; run the same helper to continue.")
                    return 0
                if choice == "A":
                    show_samples(samples[key(alias)], cache, full=True)
                    continue
                if choice == "S":
                    break
                decision = CHOICES.get(choice) or (clean(row["Final_Decision"]).upper() if not choice else "")
                if decision not in ALLOWED:
                    print("Enter K, I, or C for a decision, or A, S, or Q.")
                    continue
                while True:
                    comment = clean(input("Curator comment (Enter keeps an existing comment): ")) or clean(row["Curator_Comment"])
                    if comment:
                        break
                    print("A curator comment is required.")
                row["Final_Decision"] = decision
                row["Curator_Comment"] = comment
                save_csv(target, fields, rows)
                print(f"Saved: {decision}")
                break
    except (EOFError, KeyboardInterrupt):
        print("\nPaused. Completed decisions are saved; run the same helper to continue.")
        return 0

    remaining = sum(not completed(row) for row in rows)
    if remaining:
        print(f"\n{remaining} aliases remain incomplete. Run the helper again or complete the CSV before resume.")
        return 3
    print("\nAll alias decisions are complete. Run run_refresh.py resume with this decision CSV.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True, help="Refresh output folder used for prepare")
    args = parser.parse_args()
    try:
        return review(args.work.expanduser().resolve())
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Alias review stopped: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
