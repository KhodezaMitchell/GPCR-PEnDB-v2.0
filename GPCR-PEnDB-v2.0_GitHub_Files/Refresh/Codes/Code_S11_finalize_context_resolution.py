#!/usr/bin/env python3
# Purpose: Preserve the collision rule and apply new documented context reviews.
# Input: Context-resolved records; manual queue; optional current review decisions
# Output: Human_GPCR_Cancer_PubMed_Papers_FINAL.csv and review/integrity records

import argparse
from pathlib import Path

import pandas as pd

from _refresh_review_support import (
    PAPER_KEYS, REVIEW_EXIT, clean, ensure_template, key, read_csv, terms,
    validate_sheet, validate_supporting_terms,
)

INPUT_VALIDATED = 'Human_GPCR_Cancer_PubMed_Papers_Context_Resolved.csv'
INPUT_MANUAL = 'GPCR_Cancer_Context_Manual_Review.csv'
OUT_FINAL = 'Human_GPCR_Cancer_PubMed_Papers_FINAL.csv'
OUT_MANUAL_DECISIONS = 'GPCR_Cancer_Context_Manual_Resolution.csv'
OUT_REPORT = 'GPCR_Cancer_Context_Finalization_Integrity_Report.txt'
TEMPLATE = 'CONTEXT_REVIEW_DECISIONS_TO_COMPLETE.csv'
ALLOWED_COLLISION_TERMS = {'ccr9', 'ccr10', 'chemokine receptor ccr-9'}
ALLOWED_COLLISION_GENES = {'ackr2', 'ccr9', 'ccr10'}
REQUIRED = {
    *PAPER_KEYS, 'Official_Gene', 'Matched_KEEP_GPCR_Terms',
    'Matched_CONTEXT_REQUIRED_GPCR_Terms', 'Context_Resolution',
    'Context_Resolution_Reason', 'Strong_Context_Terms',
    'Active_Collision_Terms', 'Active_Collision_Genes',
}


def original_collision_rule(row):
    """Use precisely the existing ACKR2/CCR9/CCR10 rejection conditions."""
    context = {term.casefold() for term in terms(row['Matched_CONTEXT_REQUIRED_GPCR_Terms'])}
    collisions = {term.casefold() for term in terms(row['Active_Collision_Terms'])}
    owners = {term.casefold() for term in terms(row['Active_Collision_Genes'])}
    return (
        clean(row['Official_Gene']).casefold() == 'ackr2'
        and not terms(row['Matched_KEEP_GPCR_Terms'])
        and bool(context) and context.issubset(ALLOWED_COLLISION_TERMS)
        and bool(collisions) and collisions.issubset(ALLOWED_COLLISION_TERMS)
        and bool(owners) and owners.issubset(ALLOWED_COLLISION_GENES)
        and bool({'ccr9', 'ccr10'} & owners)
        and clean(row['Context_Resolution']) == 'MANUAL_REVIEW'
    )


def initial_final(validated):
    """Retain the established evidence-source and supporting-term assignments."""
    final = validated.copy()
    final['Final_Literature_Status'] = 'KEEP'
    sources, supporting = [], []
    for _, row in final.iterrows():
        resolution = clean(row.get('Context_Resolution', ''))
        if resolution == 'AUTO_KEEP':
            sources.append('STEP08_CONTEXT_RESOLVED')
            supporting.append(clean(row.get('Strong_Context_Terms', '')))
        elif resolution == 'NOT_REQUIRED':
            sources.append('STEP07_FINAL_KEEP_TERM')
            supporting.append(clean(row.get('Matched_KEEP_GPCR_Terms', '')))
        else:
            sources.append('VALIDATED_PRIOR_STEP')
            supporting.append(clean(row.get('Matched_KEEP_GPCR_Terms', '')))
    final['Final_Validation_Source'] = sources
    final['Validated_Search_Terms'] = supporting
    return final


def main(decisions_path=None):
    Path(OUT_FINAL).unlink(missing_ok=True)
    validated, manual = read_csv(INPUT_VALIDATED), read_csv(INPUT_MANUAL)
    if not REQUIRED.issubset(manual.columns):
        raise ValueError('Context review input lacks: ' + ', '.join(sorted(REQUIRED - set(manual.columns))))
    if manual.duplicated(subset=list(PAPER_KEYS)).any():
        raise ValueError('Duplicate context-review association rows')
    if len(manual) and not manual['Context_Resolution'].map(clean).eq('MANUAL_REVIEW').all():
        raise ValueError('Context queue contains rows not held for manual review')
    original = manual.apply(original_collision_rule, axis=1) if len(manual) else pd.Series(False, index=manual.index, dtype=bool)
    pending = manual.loc[~original].copy()
    decisions = {}
    if len(pending):
        ensure_template(pending, TEMPLATE)
    if decisions_path:
        decisions = validate_sheet(pending, decisions_path, {'KEEP', 'REJECT'})
    if len(pending) and not decisions_path:
        Path(OUT_FINAL).unlink(missing_ok=True)
        Path(OUT_REPORT).write_text(
            f'PAUSED: {len(pending)} new context-review rows require documented decisions.\n'
            f'Complete {TEMPLATE}; KEEP rows also need Supporting_GPCR_Terms.\n',
            encoding='utf-8',
        )
        print(Path(OUT_REPORT).read_text())
        raise SystemExit(REVIEW_EXIT)

    final = initial_final(validated)
    accepted, reviewed = [], []
    for index, row in manual.iterrows():
        if original.loc[index]:
            decision, supporting = 'REJECT', []
            comment = (
                'ACKR2 assignment is supported only by active CCR9/CCR10 collision '
                'terminology, with no final KEEP ACKR2 term in the article.'
            )
            reason = 'ACKR2_CCR9_CCR10_COLLISION'
        else:
            item = decisions[key(row)]
            decision, supporting, comment = item['decision'], item['supporting'], item['comment']
            reason = 'CURATOR_' + decision
            if decision == 'KEEP':
                validate_supporting_terms(row, supporting)
                record = row.to_dict()
                record.update({
                    'Final_Literature_Status': 'KEEP',
                    'Final_Validation_Source': 'REFRESH_MANUAL_CONTEXT_REVIEW',
                    'Validated_Search_Terms': ';'.join(supporting),
                    'Manual_Curator_Comment': comment,
                })
                accepted.append(record)
        record = row.to_dict()
        record.update({
            'Manual_Final_Decision': decision, 'Manual_Final_Reason': reason,
            'Manual_Final_Comment': comment, 'Manual_Supporting_GPCR_Terms': ';'.join(supporting),
        })
        reviewed.append(record)
    if accepted:
        final = pd.concat([final, pd.DataFrame(accepted)], ignore_index=True)
    if final['Validated_Search_Terms'].map(clean).eq('').any():
        raise ValueError('A final record has no supporting receptor term')
    if final.duplicated(subset=list(PAPER_KEYS)).any():
        raise ValueError('Duplicate final association records')
    reviewed_frame = pd.DataFrame(reviewed, columns=list(manual.columns) + [
        'Manual_Final_Decision', 'Manual_Final_Reason', 'Manual_Final_Comment',
        'Manual_Supporting_GPCR_Terms',
    ])
    reviewed_frame.to_csv(OUT_MANUAL_DECISIONS, index=False, na_rep='')
    final.to_csv(OUT_FINAL, index=False, na_rep='')
    lines = [
        'GPCR Cancer Refresh - Code S11 Integrity Report',
        f'Previously context-resolved rows: {len(validated)}',
        f'Original-rule rejections: {int(original.sum())}',
        f'New reviewed rows: {len(pending)}',
        f'New reviewed KEEP rows: {len(accepted)}',
        f'Final records: {len(final)}',
        'Unresolved review rows: 0',
        'No duplicate final records: PASS',
        'Every final record has supporting receptor terminology: PASS',
        'ALL CODE S11 INTEGRITY CHECKS PASSED',
    ]
    Path(OUT_REPORT).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--decisions', type=Path, help='Completed decisions for new context-review cases')
    args = parser.parse_args()
    main(args.decisions)
