"""Focused tests for fresh-run recovery, raw responses, and curator decisions.

These use synthetic records and simulated services; they make no live requests.
"""
import ast
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
CODES = ROOT / 'Refresh' / 'Codes'
sys.path.insert(0, str(CODES))
from _refresh_review_support import validate_sheet


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load('refresh_test_runner', ROOT / 'run_refresh.py')
s3 = load('refresh_test_s3', runner.code(3))
s12 = load('refresh_test_s12', runner.code(12))


@contextlib.contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def stage(number, work, decisions=None, frozen=False):
    script = runner.code(number)
    env = dict(os.environ)
    if frozen:
        script = next((ROOT / 'Codes').glob(f'Code_S{number}_*.py'))
        env['GPCR_PENDB_FROZEN_RUNNER'] = str(ROOT)
        env['PYTHONPATH'] = str(ROOT / 'Offline_Guard')
    command = [sys.executable, str(script)]
    if decisions is not None:
        command += ['--decisions', str(decisions)]
    return subprocess.run(command, cwd=work, env=env, capture_output=True, text=True)


def aliases(work, uid, gene, keep, context=''):
    pd.DataFrame([{
        'UniProt_ID': uid, 'Official_Gene': gene, 'Protein_Name': gene + ' receptor',
        'Search_Terms': keep, 'Context_Required_Terms': context,
    }]).to_csv(work / 'GPCR_Alias_Table_FINAL_VALIDATED.csv', index=False)


def context_inputs(work, gene='ACKR4'):
    base = {
        'UniProt_ID': 'TEST_LPAR1', 'Official_Gene': 'LPAR1', 'TCGA_Code': 'OV', 'PMID': '1',
        'Title': 'LPAR1 in ovarian cancer', 'Abstract': 'LPAR1 expression was measured.',
        'Matched_KEEP_GPCR_Terms': 'LPAR1', 'Matched_CONTEXT_REQUIRED_GPCR_Terms': '',
        'Context_Resolution': 'NOT_REQUIRED', 'Context_Resolution_Reason': '',
        'Strong_Context_Terms': '', 'Active_Collision_Terms': '', 'Active_Collision_Genes': '',
    }
    pd.DataFrame([base]).to_csv(work / 'Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv', index=False)
    row = dict(base, UniProt_ID='TEST_' + gene, Official_Gene=gene, PMID='2',
               Title='CCBP2 in ovarian cancer', Abstract='CCBP2 was examined in tumor samples.',
               Matched_KEEP_GPCR_Terms='', Matched_CONTEXT_REQUIRED_GPCR_Terms='CCBP2',
               Context_Resolution='MANUAL_REVIEW', Active_Collision_Terms='CCBP2',
               Active_Collision_Genes='ACKR2;ACKR4')
    if gene == 'ACKR2':
        row.update(Matched_CONTEXT_REQUIRED_GPCR_Terms='CCR9', Active_Collision_Terms='CCR9',
                   Active_Collision_Genes='ACKR2;CCR9')
    pd.DataFrame([row]).to_csv(work / 'GPCR_Cancer_Context_Manual_Review.csv', index=False)
    aliases(work, row['UniProt_ID'], gene, gene, 'CCBP2')


def ambiguity_inputs(work):
    rows = []
    for gene, uid, pmid, abstract in [
        ('SCTR', 'TEST_SCTR', '2', 'SCTR was observed in these samples.'),
        ('LPAR1', 'TEST_LPAR1', '1', 'LPAR1 expression was measured.'),
    ]:
        row = {field: '' for field in s12.REQUIRED_COLUMNS}
        row.update(UniProt_ID=uid, Official_Gene=gene, TCGA_Code='OV', PMID=pmid,
                   Title=gene + ' observations in ovarian cancer', Abstract=abstract,
                   Validated_Search_Terms=gene, Final_KEEP_GPCR_Terms=gene,
                   Final_CONTEXT_REQUIRED_GPCR_Terms='', Final_Literature_Status='KEEP')
        rows.append(row)
    pd.DataFrame(rows).to_csv(work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv', index=False)
    aliases(work, 'TEST_SCTR', 'SCTR', 'SCTR;secretin receptor')


def fill_sheet(path, decision='KEEP', supporting=''):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    frame['Final_Decision'] = decision
    frame['Curator_Comment'] = 'Synthetic test: curator assessed the intended receptor.'
    if 'Supporting_GPCR_Terms' in frame.columns:
        frame['Supporting_GPCR_Terms'] = supporting
    frame.to_csv(path, index=False)


class RefreshTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='gpcr_refresh_test_')
        self.work = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_scientific_functions_and_classifier_preserved(self):
        checks = {1: ['normalize_alias', 'parse_uniprot'],
                  3: ['build_query', 'preliminary_assessment', 'build_independent_reference'],
                  12: ['classify_ambiguity_row', 'rho_strict_decision']}
        for number, names in checks.items():
            original = next((ROOT / 'Codes').glob(f'Code_S{number}_*.py'))
            before = {node.name: ast.dump(node) for node in ast.parse(original.read_text()).body
                      if isinstance(node, ast.FunctionDef)}
            after = {node.name: ast.dump(node) for node in ast.parse(runner.code(number).read_text()).body
                     if isinstance(node, ast.FunctionDef)}
            for name in names:
                self.assertEqual(before[name], after[name])
        original = ast.parse(next((ROOT / 'Codes').glob('Code_S13_*.py')).read_text())
        original.body = [node for node in original.body
                         if not (isinstance(node, ast.ImportFrom) and node.module == '_frozen_entry')]
        self.assertEqual(ast.dump(original), ast.dump(ast.parse(runner.code(13).read_text())))

    def test_raw_search_response_persists_and_detects_tampering(self):
        xml = b'''<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart><IdList><Id>123</Id><Id>456</Id></IdList><TranslationSet/><QueryTranslation>LPAR1</QueryTranslation></eSearchResult>'''
        with working_directory(self.work):
            query = s3.build_query('LPAR1')
            s3.Entrez.local_cache = str(self.work / 'xml_parser_cache')
            with patch.object(s3.Entrez, 'esearch', side_effect=lambda **kwargs: io.BytesIO(xml)) as request:
                record = s3.esearch_with_retry(query)
                self.assertEqual(int(record['Count']), 2)
                self.assertEqual(request.call_count, 1)
            again = load('fresh_s3_cache_test', runner.code(3))
            with patch.object(again.Entrez, 'esearch', side_effect=AssertionError('unexpected request')):
                record = again.esearch_with_retry(query)
                self.assertEqual(list(record['IdList']), ['123', '456'])
            raw_path, metadata_path = again.esearch_cache_paths(query)
            self.assertEqual(raw_path.read_bytes(), xml)
            self.assertTrue(metadata_path.is_file())
            cache = {}
            again.search_pmids('LPAR1', cache, 0)
            again.export_search_snapshot(cache)
            snapshot = again.load_search_snapshot('Alias_PubMed_Search_Snapshot.csv')
            self.assertEqual(snapshot[query]['IdList'], ['123', '456'])
            raw_path.write_bytes(xml + b'changed')
            with self.assertRaises(ValueError):
                again.esearch_with_retry(query)

    def test_failed_uniprot_download_stops_stage(self):
        (self.work / runner.EXPORT_NAME).write_text('ID\tGene\tSpecies\nQ92633\tLPAR1\tsapiens\n')
        with working_directory(self.work), patch.object(requests, 'get', return_value=SimpleNamespace(status_code=503)), \
             patch('time.sleep'), contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit):
                runpy.run_path(str(runner.code(1)), run_name='__main__')

    def test_input_validation_before_collection(self):
        hgnc, export = ROOT / 'Inputs/hgnc_complete_set.txt', ROOT / 'Inputs' / runner.EXPORT_NAME
        runner.check_inputs(hgnc, export)
        invalid = self.work / 'invalid_hgnc.txt'
        invalid.write_text('symbol\nSCTR\n')
        with self.assertRaises(ValueError):
            runner.check_inputs(invalid, export)
        changed = self.work / 'different_roster.tsv'
        changed.write_text('ID\tGene\tSpecies\nP12345\tDIFFERENT\tsapiens\n')
        with self.assertRaises(ValueError):
            runner.check_inputs(hgnc, changed)
        self.assertFalse((self.work / runner.STATE_FILE).exists())

    def test_restart_skips_completed_stages_and_keeps_attempt_logs(self):
        job = self.work / 'job'
        hgnc, export = ROOT / 'Inputs/hgnc_complete_set.txt', ROOT / 'Inputs' / runner.EXPORT_NAME
        state = runner.start_or_continue(job, hgnc, export)
        calls = []
        fail_second = [True]
        def simulated_stage(command, cwd, **kwargs):
            number = int(Path(command[2]).name.split('_')[1][1:])
            calls.append(number)
            if number == 2 and fail_second[0]:
                return SimpleNamespace(returncode=1)
            for name in runner.STAGE_OUTPUTS[number]:
                (Path(cwd) / name).write_text('synthetic checkpoint output\n')
            return SimpleNamespace(returncode=0)
        with patch.object(runner.subprocess, 'run', side_effect=simulated_stage):
            runner.run_step(job, state, 1)
            with self.assertRaises(RuntimeError):
                runner.run_step(job, state, 2)
            state = runner.start_or_continue(job, hgnc, export)
            fail_second[0] = False
            runner.run_step(job, state, 1)
            runner.run_step(job, state, 2)
        self.assertEqual(calls, [1, 2, 2])
        self.assertTrue((job / 'Logs/Code_S2_attempt_001.log').exists())
        self.assertTrue((job / 'Logs/Code_S2_attempt_002.log').exists())
        runner.load_state(job)
        (job / 'GPCR_Alias_Table.csv').write_text('tampered')
        with self.assertRaises(ValueError):
            runner.load_state(job)

    def test_original_collision_rule_has_same_final_records(self):
        old, current = self.work / 'old', self.work / 'current'
        old.mkdir(); current.mkdir()
        context_inputs(old, 'ACKR2'); context_inputs(current, 'ACKR2')
        baseline, refresh = stage(11, old, frozen=True), stage(11, current)
        self.assertEqual(baseline.returncode, 0, baseline.stdout + baseline.stderr)
        self.assertEqual(refresh.returncode, 0, refresh.stdout + refresh.stderr)
        self.assertEqual((old / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv').read_bytes(),
                         (current / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv').read_bytes())

    def test_runner_worker_executes_the_real_context_stage(self):
        context_inputs(self.work, 'ACKR2')
        state = {'completed': {}, 'attempts': {}, 'reviews': {}}
        self.assertEqual(runner.run_step(self.work, state, 11), 0)
        self.assertIn('11', state['completed'])

    def test_new_context_case_requires_review_then_accepts_documented_keep(self):
        context_inputs(self.work)
        result = stage(11, self.work)
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertFalse((self.work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv').exists())
        sheet = self.work / 'CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv'
        fill_sheet(sheet, supporting='CCBP2')
        result = stage(11, self.work, sheet)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        final = pd.read_csv(self.work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv')
        self.assertEqual(len(final), 2)
        self.assertIn('REFRESH_MANUAL_CONTEXT_REVIEW', set(final['Final_Validation_Source']))

    def test_new_context_reject_and_invalid_support(self):
        context_inputs(self.work)
        stage(11, self.work)
        sheet = self.work / 'CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv'
        fill_sheet(sheet, supporting='imaginary receptor')
        result = stage(11, self.work, sheet)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv').exists())
        fill_sheet(sheet, decision='REJECT')
        result = stage(11, self.work, sheet)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(len(pd.read_csv(self.work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv')), 1)

    def test_ambiguity_pause_then_keep_and_reject(self):
        for decision, expected in [('KEEP', 2), ('REJECT', 1)]:
            work = self.work / decision
            work.mkdir()
            ambiguity_inputs(work)
            result = stage(12, work)
            self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
            self.assertFalse((work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv').exists())
            sheet = work / 'AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv'
            fill_sheet(sheet, decision=decision, supporting='SCTR' if decision == 'KEEP' else '')
            result = stage(12, work, sheet)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            final = pd.read_csv(work / 'Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv')
            self.assertEqual(len(final), expected)
            self.assertTrue((work / 'FINAL13_Ambiguity_Manual_Resolution.csv').is_file())

    def test_incomplete_or_changed_review_identifiers_rejected(self):
        queue = pd.DataFrame([{'UniProt_ID': 'TEST', 'Alias': 'term'}])
        path = self.work / 'decisions.csv'
        pd.DataFrame([{'UniProt_ID': 'TEST', 'Alias': 'term', 'Final_Decision': 'KEEP',
                       'Curator_Comment': ''}]).to_csv(path, index=False)
        with self.assertRaises(ValueError):
            validate_sheet(queue, path, {'KEEP'}, ('UniProt_ID', 'Alias'), support=False)
        frame = pd.read_csv(path).fillna('')
        frame['Alias'], frame['Curator_Comment'] = 'other', 'A comment'
        frame.to_csv(path, index=False)
        with self.assertRaises(ValueError):
            validate_sheet(queue, path, {'KEEP'}, ('UniProt_ID', 'Alias'), support=False)


if __name__ == '__main__':
    unittest.main(verbosity=2)
