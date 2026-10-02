# GPCR-PEnDB 2.0 Cancer Literature Workflow

This repository provides the supplementary Python workflow and supporting inputs for the cancer module of GPCR-PEnDB 2.0. The workflow assesses GPCR terminology, retrieves PubMed publications, screens receptor and cancer matches, resolves ambiguous assignments, assigns association categories, and constructs database tables.

## Code and archived study files

The repository includes the original workflow code, runners, the optional `review_aliases.py` helper, smaller input files, dependencies, and verification documentation. The Python files are unchanged from the complete supplementary workflow archive.

The full supplement is available through the companion [Zenodo record](https://doi.org/10.5281/zenodo.23093282) as `GPCRPEnDB_2.0_Supplementary_files.zip`. It contains:

- `GPCR_PEnDB_Manuscript_Supplement.zip`: the complete workflow, including the `Data/`, `Caches/`, and `Intermediates/` folders.
- `File_S1_GPCR_PEnDB_v2_Cancer_Networks.zip`: the supplementary cancer network pages.

To reproduce the manuscript results, download and extract the outer supplement, then extract `GPCR_PEnDB_Manuscript_Supplement.zip`. Run the frozen workflow from its extracted `GPCR_PEnDB_Manuscript_Supplement` folder. The archived data and evidence required for this run are supplied in Zenodo.

To start an independent analysis, run the Refresh workflow from this GitHub repository with a current HGNC complete-set file and your NCBI email. Refresh pauses at designated stages for documented human review decisions.

The file mapping and verification reports below describe the complete supplementary package.

## Workflow instructions

This package reproduces the manuscript results from saved evidence or starts a new analysis of the same 832 human GPCRs and 33 TCGA cancer types.

## Setup and command location

Use Linux and Python 3. Run setup and Refresh commands from the repository folder containing `run_refresh.py` and `review_aliases.py`. For frozen reproduction, use the complete workflow folder extracted from Zenodo.

```bash
python3 -m pip install -r requirements.txt pandas==2.2.3
```

The supplied workflow was tested with Python 3.12.14, pandas 2.2.3, requests 2.34.2, and Biopython 1.88. Install the dependencies before attempting offline reproduction.

The optional keyboard review helper, `review_aliases.py`, is supplied beside `run_refresh.py` and uses only the Python standard library.

## Reproduce the manuscript results

From the complete `GPCR_PEnDB_Manuscript_Supplement` folder extracted from Zenodo, run:

```bash
python3 run_frozen.py --work ../reviewer_run
```

Use a new, empty output folder. `../reviewer_run` is a folder named `reviewer_run` beside the package folder. The calculation runs offline without NCBI credentials and uses the recorded study decisions. Success reports PASS for all 15 stages and compares generated tables with the supplied data: 458 GPCRs, 33 cancer types, 3,680 associations, 18,892 records, 15,010 unique PMIDs, 16 categories, and zero integrity failures. Data S7 differs only in its `Cache_Used` flag.

S1 checks rebuilt term membership, then uses the archived display table for the exact downstream handoff. The historical S3 replay uses a recorded-search snapshot reconstructed from the audit output; its original raw PubMed responses were not preserved.

## Start a fresh analysis

A refresh requires internet access and your email address for NCBI requests. It retrieves current source information and requires curator decisions for flagged cases. Current sources and new reviews can produce different counts from the manuscript.

### 1. Download the current HGNC complete set and prepare

Download the complete HGNC TSV using `wget`:

```bash
wget -O current_hgnc_complete_set.txt \
  "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
```

Alternatively, download the [HGNC complete-set TSV](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt) manually and save it in the package folder as `current_hgnc_complete_set.txt`.

Replace the example email with your own, then run:

```bash
export NCBI_EMAIL="your.email@example.com"

python3 run_refresh.py prepare \
  --work ../refresh_run \
  --hgnc current_hgnc_complete_set.txt
```

`NCBI_API_KEY` is optional. If you open a new terminal for a later command, set `NCBI_EMAIL` again in that terminal.

### 2. Find the output folder

`--work` selects the output folder. In these examples, `../refresh_run` is a folder named `refresh_run` beside the package folder. The runner creates it for a new analysis. Use a new, empty folder outside the package; output folders inside the package are rejected.

Generated CSV tables, copied source inputs, caches, review sheets, and integrity reports are saved in the output folder. Stage attempt logs are saved in its `Logs/` subfolder. The archived manuscript tables are supplied in the full Zenodo package's `Data/` folder.

You may choose another output folder, such as `../refresh_run_02`. Use that same path in every prepare, review, and resume command for that run. Relative paths are interpreted from the directory where you run the command, so keep running the examples from the package folder.

### 3. Complete alias review after S5

The runner pauses after S5 and creates:

- `../refresh_run/Alias_Manual_Review_Queue.csv`: queued aliases and their assessment information.
- `../refresh_run/REVIEW_DECISIONS_TO_COMPLETE.csv`: the decision sheet to complete.

Review the queued aliases using their receptor information, flags, and sampled publications. Saved evidence is available in `Alias_Context_Review.csv` and `pubmed_context_cache.json` in the output folder. Record one decision and a curator comment for every queued alias.

| Keyboard choice | Value saved in `Final_Decision` | Meaning |
| --- | --- | --- |
| K | KEEP | Evidence supports use of the name for the intended receptor without additional receptor-specific terminology. |
| I | IGNORE | Retrieved records do not support safe assignment of the name to the intended receptor. |
| C | CONTEXT_REQUIRED | The name has supporting evidence but requires additional receptor-specific terminology in retrieved publications. |

`NOT_VALIDATED` is handled by the triage rules and is not a choice in this manual review sheet.

Choose either of the following review methods. Both use the same decision CSV, so you can switch between them before submitting the decisions to `resume`.

**Option A: edit the CSV directly.** Open `REVIEW_DECISIONS_TO_COMPLETE.csv`, fill `Final_Decision` with `KEEP`, `IGNORE`, or `CONTEXT_REQUIRED`, and complete `Curator_Comment` for every row. Preserve the column names, `UniProt_ID`, `Alias`, and all rows. Save as UTF-8 CSV, retaining the `.csv` extension, in the output folder.

**Option B: use the keyboard review helper.** From the package folder, run:

```bash
python3 review_aliases.py --work ../refresh_run
```

The helper displays alias information and saved sampled titles and context excerpts. Enter `K`, `I`, or `C`, then a curator comment. Each completed decision is saved immediately to `../refresh_run/REVIEW_DECISIONS_TO_COMPLETE.csv`. Existing completed decisions are preserved and skipped when the helper restarts.

- `A` displays full abstracts available in the saved cache.
- `S` skips the current alias, leaving it pending.
- `Q` exits; rerun the same command to continue pending decisions.

The helper uses saved evidence and requires no internet access or NCBI credentials. The curator supplies the decisions; Code S6 applies them. This optional interface is provided for refresh review and does not alter the recorded manuscript decisions or the numbered scripts. Its location beside `run_refresh.py` keeps it outside the runner's stage-code signature.

### 4. Resume with the completed alias CSV

Save the completed CSV as `REVIEW_DECISIONS_TO_COMPLETE.csv` in the refresh output folder. From the same package folder where you ran `prepare`, run:

```bash
python3 run_refresh.py resume \
  --work ../refresh_run \
  --decisions ../refresh_run/REVIEW_DECISIONS_TO_COMPLETE.csv
```

The runner checks that all queued identifiers are present exactly once and that every row has an allowed decision and a curator comment. Once alias decisions have been accepted, their flag need not be repeated in later resume commands. Use a new run to change already-applied decisions.

### 5. Complete any later publication reviews

New publication-review cases can cause additional pauses at S11 and S12. Complete the indicated sheet with `KEEP` or `REJECT` and a curator comment. `KEEP` decisions also require `Supporting_GPCR_Terms` from the current receptor terminology that occur in the title or abstract. Separate multiple supporting terms with semicolons. Preserve all identifier columns and rows.

For an S11 context-review pause, complete `CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv`, then run:

```bash
python3 run_refresh.py resume \
  --work ../refresh_run \
  --context-decisions ../refresh_run/CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv
```

For an S12 ambiguity-review pause, complete `AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv`, then run:

```bash
python3 run_refresh.py resume \
  --work ../refresh_run \
  --ambiguity-decisions ../refresh_run/AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv
```

These publication review sheets are completed directly as CSV files; `review_aliases.py` handles the alias review after S5.

### Interruptions and final results

After an interruption, repeat the prepare or resume command you were using, with the same inputs and output folder. Reuse the HGNC file downloaded for that run rather than downloading a different version. Completed stages are verified and skipped; saved responses are reused. Use a new output folder if source inputs, stage code, or already-applied decisions change. Keep the output folder and its run state together.

Success reports `COMPLETE: all 15 refresh stages passed` and gives the results folder. Final classifications include `CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv`; database tables include `CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv` and `CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv`. Final metrics and checks are in `CORRECTED_FINAL12_Manuscript_Metrics.txt` and `CORRECTED_FINAL12_Integrity_Report.txt`. These files are all in the refresh output folder.

Logs, code/input/output hashes, and reviewed decisions record each refresh. Preserve these files with the generated results to document the run.

## Files and mapping

`Codes/` and `run_frozen.py` preserve the manuscript workflow. `Refresh/Codes/` and `run_refresh.py` support new runs; changes are listed in `Refresh/CHANGES.md`. The optional `review_aliases.py` helper is supplied beside the runners. In the complete Zenodo workflow package, `Inputs/`, `Caches/`, and `Intermediates/` contain the historical evidence and decisions, and `Data/` contains Data S1–S17. This repository includes the smaller `Inputs/` files; download the larger evidence and data folders from Zenodo. Data S5 and S9 are supplied reference tables, not generated results.

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

Literature-search date for the manuscript: 22 August 2026. Classification revision: 26 September 2026. Archived workflow verification: `WORKFLOW_VERIFICATION.txt`. Packaged workflow code hashes: `CODE_VERSIONS.json`. These records describe the supplied workflow; each refresh also records its own run information.
