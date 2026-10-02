# Refresh changes

The original manuscript scripts, offline runner, Data S1–S17, archived inputs, evidence, and decisions are preserved. The refresh scripts start from the supplied online versions of Code S1–S15. Their file hashes are recorded in `../CODE_VERSIONS.json`.

| File | Refresh change |
| --- | --- |
| Code S1 | Stop after failed UniProt downloads; successful cached responses remain reusable. |
| Code S3 | Save raw ESearch XML and request metadata for every query, check cached-response hashes, and export the completed search snapshot. Queries, PMID limits, reference construction, and alias assessment rules are unchanged. |
| Code S11 | Preserve the ACKR2/CCR9/CCR10 rejection rule. Queue other cases for current KEEP/REJECT decisions and carry documented KEEP records forward with supporting receptor terms. |
| Code S12 | Preserve automatic ambiguity rules. Accept current decisions for unresolved rows and retain both automatic and curator decisions in the audit trail. |
| `_refresh_review_support.py` | Validate review identifiers, comments, decisions, and supporting terms in current terminology and article text. |
| `_run_stage.py` | Execute each stage with a writable Biopython XML parser cache in the work folder. |
| `../run_refresh.py` | Check inputs before requests; record code/input/output hashes, dated attempts, and reviewed files; verify and skip completed stages; retain failed-attempt logs; pause at required reviews. |

Codes S2, S4–S10, and S13–S15 are unchanged from the supplied online scripts. The 832-receptor roster, 33 cancer types, human-related search criteria, review exclusion, alias assessment thresholds, and 16-category classification rules remain unchanged.

The preserved manuscript component reproduces the archived calculation. Refresh-specific review handling is a documented extension for new evidence, not a claim that these additions were part of the historical run.

Run focused checks, which use synthetic records and simulated services without live requests:

```bash
python3 Refresh/tests/test_refresh.py
```
