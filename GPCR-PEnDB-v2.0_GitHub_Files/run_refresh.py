#!/usr/bin/env python3
"""Collect current sources, pause for documented reviews, and resume safely."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CODES = ROOT / "Refresh" / "Codes"
sys.path.insert(0, str(CODES))
from _refresh_review_support import ensure_template, read_csv, validate_sheet

STATE_FILE = ".refresh_run_state.json"
RUN_VERSION = 1
EXPORT_NAME = "gpcrpendb_results_1782027334.2077546.tsv"
HGNC_FIELDS = {
    "hgnc_id", "symbol", "name", "status", "alias_symbol", "prev_symbol",
    "alias_name", "prev_name", "uniprot_ids",
}
ALIAS_QUEUE = "Alias_Manual_Review_Queue.csv"
ALIAS_TEMPLATE = "REVIEW_DECISIONS_TO_COMPLETE.csv"
REVIEW_FILES = {
    "alias": (6, "Alias_Manual_Review_Assistant_Decisions.csv"),
    "context": (11, "Context_Reviewed_Decisions.csv"),
    "ambiguity": (12, "Ambiguity_Reviewed_Decisions.csv"),
}
STAGE_OUTPUTS = {
    1: ["GPCR_Alias_Table.csv"],
    2: ["HGNC_Alias_Table.csv", "GPCR_Alias_Source_Comparison.csv"],
    3: ["Alias_Candidate_Universe.csv", "Alias_Audit_Recommendations.csv", "Alias_Audit_Stats.csv",
        "Alias_Cross_GPCR_Collisions.csv", "Alias_Audit_Integrity_Report.txt", "Alias_PubMed_Search_Snapshot.csv"],
    4: ["Alias_Context_Review.csv", "Alias_Context_Review_Integrity_Report.txt"],
    5: ["Alias_Systematic_Triage.csv", ALIAS_QUEUE, "Alias_Systematic_Triage_Summary.txt",
        "Alias_Systematic_Triage_Integrity_Report.txt"],
    6: ["GPCR_Alias_Table_FINAL_VALIDATED.csv", "Alias_Final_Decisions.csv",
        "Alias_Final_Active_Collisions.csv", "Alias_Finalization_Integrity_Report.txt"],
    7: ["TCGA_Cancer_Vocabulary.csv", "TCGA_Cancer_PubMed_Terms_Long.csv",
        "TCGA_Cancer_UI_Aliases_Long.csv", "TCGA_Cancer_Vocabulary_Integrity_Report.txt"],
    8: ["GPCR_Cancer_PubMed_Search_Summary.csv", "GPCR_Cancer_PubMed_PMIDs_Long.csv",
        "GPCR_Cancer_PubMed_Search_Errors.csv", "GPCR_Cancer_PubMed_Search_Integrity_Report.txt"],
    9: ["PubMed_Article_Metadata.csv", "GPCR_Cancer_PubMed_Validation_All.csv",
        "Human_GPCR_Cancer_PubMed_Papers_Validated.csv", "GPCR_Cancer_PubMed_Context_Review.csv",
        "GPCR_Cancer_PubMed_Validation_Summary.csv", "PubMed_Article_Metadata_Errors.csv",
        "GPCR_Cancer_PubMed_Validation_Integrity_Report.txt"],
    10: ["GPCR_Cancer_Context_Resolution_All.csv", "GPCR_Cancer_Context_Manual_Review.csv",
         "GPCR_Cancer_Context_Rejected.csv", "Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv",
         "GPCR_Cancer_Context_Resolution_Summary.csv", "GPCR_Cancer_Context_Resolution_Integrity_Report.txt"],
    11: ["Human_GPCR_Cancer_PubMed_Papers_FINAL.csv", "GPCR_Cancer_Context_Manual_Resolution.csv",
         "GPCR_Cancer_Context_Finalization_Integrity_Report.txt"],
    12: ["FINAL13_Ambiguity_Audit_All.csv", "FINAL13_Ambiguity_Rejected.csv", "FINAL13_Ambiguity_Manual_Review.csv",
         "Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv", "FINAL13_Ambiguous_Official_Symbol_Summary.csv",
         "FINAL13_Integrity_Report.txt"],
    13: ["CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv",
         "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Associations_Retained.csv",
         "CORRECTED_FINAL_GPCR_Cancer_Association_Classification_Summary.csv",
         "CORRECTED_FINAL_GPCR_Cancer_Association_Classification_Integrity_Report.txt"],
    14: ["CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv", "CORRECTED_FINAL16_DB_Associations_GPCR_Cancer.csv",
         "CORRECTED_FINAL16_DB_Page_GPCR_to_Cancers.csv", "CORRECTED_FINAL16_DB_Page_Cancer_to_GPCRs.csv",
         "CORRECTED_FINAL16_DB_Cancer_Search_Aliases.csv", "CORRECTED_FINAL16_DB_Summary_Counts.csv",
         "CORRECTED_FINAL16_GPCR_Cancer_Database_Integrity_Report.txt"],
    15: ["CORRECTED_FINAL12_Overall_Summary.csv", "CORRECTED_FINAL12_Category_Summary.csv",
         "CORRECTED_FINAL12_GPCR_Summary.csv", "CORRECTED_FINAL12_Cancer_Summary.csv",
         "CORRECTED_FINAL12_GPCR_Cancer_Count_Distribution.csv", "CORRECTED_FINAL12_Association_Paper_Count_Distribution.csv",
         "CORRECTED_FINAL12_Year_Summary.csv", "CORRECTED_FINAL12_Manuscript_Metrics.txt", "CORRECTED_FINAL12_Integrity_Report.txt"],
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code(number):
    matches = list(CODES.glob(f"Code_S{number}_*.py"))
    if len(matches) != 1:
        raise ValueError(f"Expected one Code S{number}; found {len(matches)}")
    return matches[0]


def code_signature():
    files = [Path(__file__).resolve(), *sorted(CODES.glob("*.py")), ROOT / "Refresh" / "_run_stage.py"]
    return {str(path.relative_to(ROOT)): sha256(path) for path in files}


def receptor_ids(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"ID", "Gene", "Species"}.issubset(reader.fieldnames or []):
            raise ValueError("GPCR export must have ID, Gene, and Species columns")
        rows = [row for row in reader if row["Species"].strip() == "sapiens"]
    if any(not row["ID"].strip() or not row["Gene"].strip() for row in rows):
        raise ValueError("Human GPCR identifiers and original names must be complete")
    identifiers = [row["ID"].strip() for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate human receptor IDs in the GPCR export")
    return set(identifiers)


def check_inputs(hgnc, export):
    expected = receptor_ids(ROOT / "Inputs" / EXPORT_NAME)
    if len(expected) != 832 or receptor_ids(export) != expected:
        raise ValueError("The refresh uses the same 832 human receptor IDs as the manuscript")
    with Path(hgnc).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = HGNC_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError("HGNC input lacks: " + ", ".join(sorted(missing)))
        if not any(row.get("status", "").strip() == "Approved" for row in reader):
            raise ValueError("HGNC complete set has no approved records")
    if Path(hgnc).stat().st_size <= 1000:
        raise ValueError("Supply the complete HGNC TSV, not a small header or subset")


def save_state(work, state):
    target = work / STATE_FILE
    temporary = work / (STATE_FILE + ".tmp")
    temporary.write_text(json.dumps(state, indent=2), encoding="utf-8")
    temporary.replace(target)


def load_state(work):
    path = work / STATE_FILE
    if not path.is_file():
        raise ValueError("This folder has no refresh state; start prepare in a new empty folder")
    state = json.loads(path.read_text(encoding="utf-8"))
    if state.get("version") != RUN_VERSION or state.get("scripts") != code_signature():
        raise ValueError("Refresh code differs from the version used to start this run")
    for name, digest in state["inputs"].items():
        if not (work / name).is_file() or sha256(work / name) != digest:
            raise ValueError(f"Staged source input changed: {name}")
    for stage, record in state["completed"].items():
        for name, digest in record["outputs"].items():
            if not (work / name).is_file() or sha256(work / name) != digest:
                raise ValueError(f"Completed Code S{stage} output changed: {name}")
    for kind, record in state.get("reviews", {}).items():
        stage, name = REVIEW_FILES[kind]
        if str(stage) in state["completed"] and sha256(work / name) != record["sha256"]:
            raise ValueError(f"Already-applied {kind} decisions changed; use a new run")
    return state


def start_or_continue(work, hgnc, export):
    if work == ROOT or ROOT in work.parents:
        raise ValueError("Keep refresh outputs outside the supplementary package")
    check_inputs(hgnc, export)
    hashes = {"hgnc_complete_set.txt": sha256(hgnc), EXPORT_NAME: sha256(export)}
    if work.exists() and any(work.iterdir()):
        state = load_state(work)
        if state["inputs"] != hashes:
            raise ValueError("Prepare inputs differ from the inputs already used in this folder")
        return state
    work.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(hgnc, work / "hgnc_complete_set.txt")
    shutil.copyfile(export, work / EXPORT_NAME)
    state = {
        "version": RUN_VERSION, "scripts": code_signature(), "inputs": hashes,
        "started_at_utc": utc_now(), "completed": {}, "attempts": {}, "reviews": {},
    }
    save_state(work, state)
    return state


def register_review(work, state, kind, source):
    stage, filename = REVIEW_FILES[kind]
    target = work / filename
    if source is not None:
        source = Path(source).resolve()
        if not source.is_file():
            raise ValueError(f"Decision file does not exist: {source}")
        digest = sha256(source)
        prior = state["reviews"].get(kind)
        if str(stage) in state["completed"] and (not prior or prior["sha256"] != digest):
            raise ValueError(f"These {kind} decisions differ from decisions already applied")
        if source != target.resolve():
            shutil.copyfile(source, target)
        state["reviews"][kind] = {"file": filename, "sha256": digest, "accepted_at_utc": utc_now()}
        save_state(work, state)
    return target if kind in state["reviews"] else None


def run_step(work, state, number):
    if str(number) in state["completed"]:
        print(f"SKIP Code S{number}: completed outputs verified", flush=True)
        return 0
    env = dict(os.environ)
    env.pop("GPCR_PENDB_FROZEN_RUNNER", None)
    paths = [entry for entry in env.get("PYTHONPATH", "").split(os.pathsep)
             if entry and "Offline_Guard" not in Path(entry).parts]
    env["PYTHONPATH"] = os.pathsep.join(paths)
    env["GPCR_PENDB_REFRESH_XML_CACHE"] = str(work / "Xml_Parser_Cache")
    env.update(PYTHONHASHSEED="0", PYTHONUNBUFFERED="1")
    command = [sys.executable, str(ROOT / "Refresh" / "_run_stage.py"), str(code(number))]
    if number in (11, 12):
        kind = "context" if number == 11 else "ambiguity"
        if kind in state["reviews"]:
            command += ["--decisions", str(work / REVIEW_FILES[kind][1])]
    attempt = state["attempts"].get(str(number), 0) + 1
    state["attempts"][str(number)] = attempt
    state["running"] = {"stage": number, "started_at_utc": utc_now()}
    save_state(work, state)
    logs = work / "Logs"
    logs.mkdir(exist_ok=True)
    log_path = logs / f"Code_S{number}_attempt_{attempt:03d}.log"
    print(f"Running Code S{number}: {code(number).name}", flush=True)
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
    shutil.copyfile(log_path, work / f"Code_S{number}_run.log")
    state.pop("running", None)
    if result.returncode:
        state["blocked"] = {"stage": number, "exit_code": result.returncode, "log": str(log_path.relative_to(work))}
        save_state(work, state)
        if result.returncode == 3 and number in (11, 12):
            kind = "context" if number == 11 else "ambiguity"
            template = "CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv" if number == 11 else "AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv"
            print(f"PAUSED Code S{number}: complete {work / template}")
            print(f"Then resume with --{kind}-decisions pointing to that completed file.")
            return 3
        raise RuntimeError(f"Code S{number} stopped. Inspect {log_path}; repeat the same command after resolving it.")
    outputs = {}
    names = list(STAGE_OUTPUTS[number])
    if number == 12 and (work / "FINAL13_Ambiguity_Manual_Resolution.csv").is_file():
        names.append("FINAL13_Ambiguity_Manual_Resolution.csv")
    for name in names:
        if not (work / name).is_file():
            raise RuntimeError(f"Code S{number} did not produce expected output: {name}")
        outputs[name] = sha256(work / name)
    state["completed"][str(number)] = {
        "outputs": outputs, "completed_at_utc": utc_now(), "log": str(log_path.relative_to(work)),
    }
    state.pop("blocked", None)
    save_state(work, state)
    print(f"PASS Code S{number}: outputs recorded", flush=True)
    return 0


def prepare(work, state):
    for number in range(1, 6):
        run_step(work, state, number)
    queue = read_csv(work / ALIAS_QUEUE)
    ensure_template(queue, work / ALIAS_TEMPLATE, ("UniProt_ID", "Alias"), support=False)
    print(f"PAUSED: {len(queue)} aliases need current reviewed decisions.")
    print(f"Complete {work / ALIAS_TEMPLATE}, then run resume with --decisions.")


def resume(work, state, args):
    if any(str(number) not in state["completed"] for number in range(1, 6)):
        raise ValueError("Complete prepare through Code S5 before resume")
    queue = read_csv(work / ALIAS_QUEUE)
    alias_path = args.decisions or (
        work / REVIEW_FILES["alias"][1] if "alias" in state["reviews"] else None
    )
    if alias_path is None and queue.empty:
        alias_path = work / ALIAS_TEMPLATE
    if alias_path is None:
        raise ValueError("Supply --decisions with the completed current alias review sheet")
    validate_sheet(queue, alias_path, {"KEEP", "IGNORE", "CONTEXT_REQUIRED"},
                   ("UniProt_ID", "Alias"), support=False)
    register_review(work, state, "alias", alias_path)
    register_review(work, state, "context", args.context_decisions)
    register_review(work, state, "ambiguity", args.ambiguity_decisions)
    for number in range(6, 16):
        status = run_step(work, state, number)
        if status == 3:
            return status
    state["finished_at_utc"] = utc_now()
    save_state(work, state)
    print(f"COMPLETE: all 15 refresh stages passed. Results and logs: {work}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    first = modes.add_parser("prepare", help="Collect current sources through S5; repeat to continue interruptions")
    first.add_argument("--work", type=Path, required=True)
    first.add_argument("--hgnc", type=Path, required=True, help="Current complete-set HGNC TSV")
    first.add_argument("--gpcr-export", type=Path, default=ROOT / "Inputs" / EXPORT_NAME)
    second = modes.add_parser("resume", help="Apply current reviews and continue; repeat after interruptions")
    second.add_argument("--work", type=Path, required=True)
    second.add_argument("--decisions", type=Path, help="Completed current alias decisions")
    second.add_argument("--context-decisions", type=Path, help="Completed new context-review decisions")
    second.add_argument("--ambiguity-decisions", type=Path, help="Completed new ambiguity-review decisions")
    args = parser.parse_args(argv)
    if not os.environ.get("NCBI_EMAIL", "").strip():
        parser.error("Set NCBI_EMAIL to your email address for live PubMed requests")
    work = args.work.resolve()
    try:
        if args.mode == "prepare":
            state = start_or_continue(work, args.hgnc, args.gpcr_export)
            prepare(work, state)
            return 0
        return resume(work, load_state(work), args)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Refresh stopped: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
