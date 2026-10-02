#!/usr/bin/env python3
"""Replay the frozen GPCR-PEnDB calculation with bundled historical evidence."""

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
INPUTS = ROOT / "Inputs"
INTERMEDIATES = ROOT / "Intermediates"
STATE_FILE = ".frozen_run_state.json"
RUN_VERSION = 2

DATA_OUTPUTS = {
    "Data_S1_Alias_Candidate_Universe.csv": "Alias_Candidate_Universe.csv",
    "Data_S2_Alias_Audit_Recommendations.csv": "Alias_Audit_Recommendations.csv",
    "Data_S3_Alias_Final_Decisions.csv": "Alias_Final_Decisions.csv",
    "Data_S4_GPCR_Alias_Table_FINAL_VALIDATED.csv": "GPCR_Alias_Table_FINAL_VALIDATED.csv",
    "Data_S6_GPCR_Cancer_PubMed_PMIDs_Long.csv": "GPCR_Cancer_PubMed_PMIDs_Long.csv",
    "Data_S7_GPCR_Cancer_PubMed_Search_Summary.csv": "GPCR_Cancer_PubMed_Search_Summary.csv",
    "Data_S8_Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv": "Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv",
    "Data_S10_CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv": "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv",
    "Data_S11_CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Associations_Retained.csv": "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Associations_Retained.csv",
    "Data_S12_CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv": "CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv",
    "Data_S13_CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv": "CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv",
    "Data_S14_CORRECTED_FINAL16_DB_Page_GPCR_to_Cancers.csv": "CORRECTED_FINAL16_DB_Page_GPCR_to_Cancers.csv",
    "Data_S15_CORRECTED_FINAL16_DB_Page_Cancer_to_GPCRs.csv": "CORRECTED_FINAL16_DB_Page_Cancer_to_GPCRs.csv",
    "Data_S16_CORRECTED_FINAL12_Manuscript_Metrics.txt": "CORRECTED_FINAL12_Manuscript_Metrics.txt",
    "Data_S17_CORRECTED_FINAL12_Integrity_Report.txt": "CORRECTED_FINAL12_Integrity_Report.txt",
}

DATA_STAGE = {
    1: 3, 2: 3, 3: 6, 4: 6, 6: 8, 7: 8, 8: 12,
    10: 13, 11: 13, 12: 14, 13: 14, 14: 14, 15: 14,
    16: 15, 17: 15,
}

INTERMEDIATE_STAGE = {
    "Alias_Audit_Integrity_Report.txt": 3,
    "Alias_Audit_Stats.csv": 3,
    "Alias_Cross_GPCR_Collisions.csv": 3,
    "Alias_Context_Review.csv": 4,
    "Alias_Context_Review_Integrity_Report.txt": 4,
    "Alias_Manual_Review_Queue.csv": 5,
    "Alias_Systematic_Triage.csv": 5,
    "TCGA_Cancer_PubMed_Terms_Long.csv": 7,
    "TCGA_Cancer_UI_Aliases_Long.csv": 7,
    "TCGA_Cancer_Vocabulary.csv": 7,
    "PubMed_Article_Metadata.csv": 9,
    "GPCR_Cancer_PubMed_Context_Review.csv": 9,
    "Human_GPCR_Cancer_PubMed_Papers_Validated.csv": 9,
    "GPCR_Cancer_Context_Manual_Review.csv": 10,
    "GPCR_Cancer_Context_Resolution_All.csv": 10,
    "Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv": 10,
    "Human_GPCR_Cancer_PubMed_Papers_FINAL.csv": 11,
    "FINAL13_Ambiguity_Rejected.csv": 12,
    "CORRECTED_FINAL16_DB_Summary_Counts.csv": 14,
    "CORRECTED_FINAL12_Category_Summary.csv": 15,
    "CORRECTED_FINAL12_Overall_Summary.csv": 15,
}

SAVED_INPUTS = {
    "gpcrpendb_results_1782027334.2077546.tsv": INPUTS / "gpcrpendb_results_1782027334.2077546.tsv",
    "hgnc_complete_set.txt": INPUTS / "hgnc_complete_set.txt",
    "Alias_Manual_Review_Assistant_Decisions.csv": INPUTS / "Alias_Manual_Review_Assistant_Decisions.csv",
    "uniprot_cache": ROOT / "Caches" / "uniprot_cache",
    "pubmed_context_cache.json": ROOT / "Caches" / "pubmed_context_cache.json",
    "pubmed_cancer_search_cache_final": ROOT / "Caches" / "pubmed_cancer_search_cache_final",
    "pubmed_article_cache_final": ROOT / "Caches" / "pubmed_article_cache_final",
}


def same_csv_except_cache_flag(expected, observed):
    """Ignore the provenance flag that changes when saved searches are replayed."""
    with expected.open(newline="", encoding="utf-8") as a, observed.open(newline="", encoding="utf-8") as b:
        ea, ob = csv.DictReader(a), csv.DictReader(b)
        if ea.fieldnames != ob.fieldnames:
            return False
        for left, right in zip(ea, ob):
            left.pop("Cache_Used", None)
            right.pop("Cache_Used", None)
            if left != right:
                return False
        return next(ea, None) is None and next(ob, None) is None


def same_uniprot_terms(expected, observed):
    """Check receptor identity and term membership despite tied display spellings."""
    with expected.open(newline="", encoding="utf-8") as a, observed.open(newline="", encoding="utf-8") as b:
        ea, ob = list(csv.DictReader(a)), list(csv.DictReader(b))
    if len(ea) != 832 or len(ob) != len(ea):
        return False
    for left, right in zip(ea, ob):
        if any(left[k] != right[k] for k in ("UniProt_ID", "Original_Gene", "Official_Gene", "Protein_Name")):
            return False
        if {x.casefold() for x in left["Search_Terms"].split(";")} != {x.casefold() for x in right["Search_Terms"].split(";")}:
            return False
    return True


def stage_inputs(work):
    """Link bundled input evidence into an otherwise empty work directory."""
    for name, source in SAVED_INPUTS.items():
        if not source.exists():
            raise FileNotFoundError(source)
        (work / name).symlink_to(source.resolve(), target_is_directory=source.is_dir())


def check_staged_inputs(work):
    """Reject a resumed work directory with changed archived inputs."""
    for name, source in SAVED_INPUTS.items():
        link = work / name
        if not link.is_symlink() or link.resolve() != source.resolve():
            raise RuntimeError(f"Frozen input link has changed: {name}")


def run_step(work, number):
    """Run one code with saved evidence and an enforced network guard."""
    env = dict(os.environ)
    env.update(NCBI_EMAIL="frozen-replay@example.invalid", PYTHONHASHSEED="0")
    env["GPCR_PENDB_FROZEN_RUNNER"] = str(ROOT)
    env["PYTHONPATH"] = str(ROOT / "Offline_Guard") + os.pathsep + env.get("PYTHONPATH", "")
    matches = list((ROOT / "Codes").glob(f"Code_S{number}_*.py"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one Code S{number}, found {len(matches)}")
    command = [sys.executable, str(matches[0])]
    if number == 3:
        command += ["--snapshot", str(INPUTS / "Alias_PubMed_Search_Snapshot.csv")]
    print(f"Running Code S{number}: {matches[0].name}", flush=True)
    with (work / f"Code_S{number}_run.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    if number == 1:
        rebuilt = work / "GPCR_Alias_Table.csv"
        canonical = INPUTS / "Frozen_Terminology" / "GPCR_Alias_Table.csv"
        if not same_uniprot_terms(canonical, rebuilt):
            raise RuntimeError("S1 term membership or receptor identity differs from frozen input")
        rebuilt.rename(work / "GPCR_Alias_Table_Rebuilt_From_Cache.csv")
        shutil.copyfile(canonical, rebuilt)
        print("S1: rebuilt terms match; archived display table selected for exact downstream replay.", flush=True)


def verify_stage(work, number):
    """Compare this stage's recreated tables to the packaged frozen tables."""
    checked = 0
    if number == 1:
        canonical = INPUTS / "Frozen_Terminology" / "GPCR_Alias_Table.csv"
        rebuilt = work / "GPCR_Alias_Table_Rebuilt_From_Cache.csv"
        if not rebuilt.is_file() or not same_uniprot_terms(canonical, rebuilt):
            raise RuntimeError("S1 receptor identities or term membership differ")
        if (work / "GPCR_Alias_Table.csv").read_bytes() != canonical.read_bytes():
            raise RuntimeError("S1 downstream alias table differs from frozen input")
        checked += 2
    if number == 2:
        canonical = INPUTS / "Frozen_Terminology" / "HGNC_Alias_Table.csv"
        if (work / "HGNC_Alias_Table.csv").read_bytes() != canonical.read_bytes():
            raise RuntimeError("S2 HGNC alias table differs from frozen input")
        checked += 1
    for published, actual in DATA_OUTPUTS.items():
        index = int(published.split("_")[1][1:])
        if DATA_STAGE[index] != number:
            continue
        expected_path, actual_path = DATA / published, work / actual
        if not actual_path.is_file():
            raise RuntimeError(f"Missing reproduced output: {actual}")
        if index == 7:
            okay = same_csv_except_cache_flag(expected_path, actual_path)
        else:
            okay = expected_path.read_bytes() == actual_path.read_bytes()
        if not okay:
            raise RuntimeError(f"Stage S{number} differs from packaged Data: {published}")
        checked += 1
    for name, stage in INTERMEDIATE_STAGE.items():
        if stage != number:
            continue
        expected_path, actual_path = INTERMEDIATES / name, work / name
        if not actual_path.is_file():
            raise RuntimeError(f"Missing reproduced intermediate: {name}")
        if name.endswith(".csv") and expected_path.read_bytes() != actual_path.read_bytes():
            raise RuntimeError(f"Stage S{number} differs from packaged intermediate: {name}")
        checked += 1
    print(f"PASS Code S{number}: {checked} frozen stage checks", flush=True)


def save_state(work, completed):
    """Record which verified stage can be resumed."""
    target = work / STATE_FILE
    temporary = work / (STATE_FILE + ".tmp")
    temporary.write_text(
        json.dumps({"version": RUN_VERSION, "root": str(ROOT), "completed": completed}),
        encoding="utf-8",
    )
    temporary.replace(target)


def load_state(work):
    """Accept only an existing work area made by this frozen runner."""
    path = work / STATE_FILE
    if not path.is_file():
        raise RuntimeError("Work directory is not empty and has no frozen runner state")
    state = json.loads(path.read_text(encoding="utf-8"))
    if state.get("version") != RUN_VERSION or state.get("root") != str(ROOT):
        raise RuntimeError("Work directory belongs to a different frozen package or runner version")
    completed = state.get("completed")
    if type(completed) is not int or not 0 <= completed <= 15:
        raise RuntimeError("Frozen runner state is invalid")
    check_staged_inputs(work)
    for number in range(1, completed + 1):
        verify_stage(work, number)
    if completed < 15 and (work / f"Code_S{completed + 1}_run.log").exists():
        raise RuntimeError("Next stage previously failed; start a new empty work directory")
    return completed


def verify(work):
    """Compare recalculated tables to the packaged manuscript-facing files."""
    checked = 0
    for published, actual in DATA_OUTPUTS.items():
        expected_path, actual_path = DATA / published, work / actual
        if not actual_path.is_file():
            raise RuntimeError(f"Missing reproduced output: {actual}")
        if published.startswith("Data_S7_"):
            okay = same_csv_except_cache_flag(expected_path, actual_path)
        else:
            okay = expected_path.read_bytes() == actual_path.read_bytes()
        if not okay:
            raise RuntimeError(f"Output differs from packaged Data: {published}")
        checked += 1
    for expected_path in sorted(INTERMEDIATES.iterdir()):
        observed = work / expected_path.name
        if not observed.is_file():
            raise RuntimeError(f"Missing reproduced intermediate: {expected_path.name}")
        if expected_path.suffix == ".csv" and expected_path.read_bytes() != observed.read_bytes():
            raise RuntimeError(f"Intermediate differs: {expected_path.name}")
    print(f"PASS: {checked} Data outputs, S1 term membership, and bundled CSV intermediates match.")
    print("Data S7 differs only in Cache_Used, which records replay provenance.")
    print("Data S5 and S9 are static reference tables; no code regenerates them.")


def main():
    """Run all stages, or resume one verified stage at a time."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True, help="New or previously verified frozen work directory")
    parser.add_argument("--until", default="15", metavar="S<N>", help="Last stage to run, S1 through S15 (default: S15)")
    args = parser.parse_args()
    try:
        until = int(args.until.upper().removeprefix("S"))
    except ValueError:
        parser.error("--until must be S1 through S15")
    if not 1 <= until <= 15:
        parser.error("--until must be S1 through S15")
    work = args.work.resolve()
    if work.exists():
        completed = load_state(work) if any(work.iterdir()) else 0
    else:
        work.mkdir(parents=True)
        completed = 0
    if completed == 0 and not any(work.iterdir()):
        stage_inputs(work)
        save_state(work, 0)
    if completed > until:
        parser.error(f"Work directory already completed S{completed}; requested S{until}")
    for number in range(completed + 1, until + 1):
        run_step(work, number)
        verify_stage(work, number)
        save_state(work, number)
    if until == 15:
        verify(work)
    else:
        print(f"PAUSED after verified Code S{until}. Inspect the stage log and outputs in {work}.")
        print(f"Continue with: python {ROOT / 'run_frozen.py'} --work {work} --until S{until + 1}")


if __name__ == "__main__":
    main()
