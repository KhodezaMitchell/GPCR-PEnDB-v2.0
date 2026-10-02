# GPCR-PEnDB manuscript and refresh workflows

This package reproduces the manuscript results from saved evidence or starts a new analysis of the same 832 human GPCRs and 33 TCGA cancer types.

Use Linux and install dependencies:

```bash
python3 -m pip install -r requirements.txt pandas==2.2.3
```

Tested with Python 3.12.14, pandas 2.2.3, requests 2.34.2, and Biopython 1.88.

## Reproduce the manuscript results

From this folder, run:

```bash
python3 run_frozen.py --work ../reviewer_run
```

Use a new, empty output folder. The calculation runs offline without NCBI credentials. Success reports PASS for all 15 stages and compares generated tables with the supplied data: 458 GPCRs, 33 cancer types, 3,680 associations, 18,892 records, 15,010 unique PMIDs, 16 categories, and zero integrity failures. Data S7 differs only in its `Cache_Used` flag.

S1 checks rebuilt term membership, then uses the archived display table for the exact downstream handoff. The historical S3 replay uses a recorded-search snapshot reconstructed from the audit output; its original raw PubMed responses were not preserved.

## Start a fresh analysis

Supply a current [HGNC complete-set TSV](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt), set your email, then collect evidence:

```bash
export NCBI_EMAIL="your.email@example.edu"
python3 run_refresh.py prepare --work ../refresh_run --hgnc /path/to/current_hgnc_complete_set.txt
```

The runner pauses after S5. Complete `../refresh_run/REVIEW_DECISIONS_TO_COMPLETE.csv` with KEEP, IGNORE, or CONTEXT_REQUIRED and a curator comment for every queued alias. Then continue:

```bash
python3 run_refresh.py resume --work ../refresh_run --decisions ../refresh_run/REVIEW_DECISIONS_TO_COMPLETE.csv
```

New publication-review cases cause additional pauses. Complete the indicated sheet with KEEP or REJECT and comments. KEEP decisions also require `Supporting_GPCR_Terms` from the current receptor terminology that occur in the title or abstract. Continue with the corresponding command:

```bash
python3 run_refresh.py resume --work ../refresh_run --context-decisions ../refresh_run/CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv
python3 run_refresh.py resume --work ../refresh_run --ambiguity-decisions ../refresh_run/AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv
```

Previously accepted alias decisions remain available; their flag need not be repeated. Repeat the same prepare or resume command after an interruption. Completed stages are verified and skipped; saved responses are reused. Use a new folder if inputs, code, or already-applied decisions change. Current sources and new reviews can produce different counts. Logs, code/input/output hashes, and reviewed decisions record each refresh. `NCBI_API_KEY` is optional.

## Files and mapping

`Codes/` and `run_frozen.py` preserve the manuscript workflow. `Refresh/Codes/` and `run_refresh.py` support new runs; changes are listed in `Refresh/CHANGES.md`. `Inputs/`, `Caches/`, and `Intermediates/` contain the historical evidence and decisions. `Data/` contains Data S1–S17. S5 and S9 are supplied reference tables, not generated results.

| Code | Purpose | Data relationship |
| --- | --- | --- |
| S1 | UniProt terminology | Intermediate tables |
| S2 | HGNC terminology | Intermediate tables |
| S3 | Candidate names and alias audit | S1–S2 |
| S4 | Alias context evidence | Intermediate tables |
| S5 | Alias triage and review queue | Intermediate tables |
| S6 | Final terminology decisions | S3–S4 |
| S7 | Cancer vocabulary | Reference S5 |
| S8 | PubMed retrieval | S6–S7 |
| S9 | Publication validation | Intermediate tables; exclusions in reference S5 |
| S10 | Context resolution | Intermediate tables |
| S11 | Final publication review | Intermediate tables |
| S12 | Ambiguity audit | S8 |
| S13 | Category assignment | S10–S11; rules in reference S9 |
| S14 | Database tables | S12–S15 |
| S15 | Metrics and integrity | S16–S17 |

Literature-search date: 22 August 2026. Classification revision: 26 September 2026. Verification: `WORKFLOW_VERIFICATION.txt`. Code hashes: `CODE_VERSIONS.json`.
