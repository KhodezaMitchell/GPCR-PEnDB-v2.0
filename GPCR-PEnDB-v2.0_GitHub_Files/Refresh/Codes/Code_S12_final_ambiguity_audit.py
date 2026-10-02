#!/usr/bin/env python3
# Purpose: Audit short official symbols and repair ambiguous assignments.
# Input: Human_GPCR_Cancer_PubMed_Papers_FINAL.csv; GPCR_Alias_Table_FINAL_VALIDATED.csv
# Output: Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from functools import lru_cache
from pathlib import Path

import pandas as pd

from _refresh_review_support import apply_ambiguity_reviews, REVIEW_EXIT

INPUT_FILE = "Human_GPCR_Cancer_PubMed_Papers_FINAL.csv"

OUT_AUDIT = "FINAL13_Ambiguity_Audit_All.csv"
OUT_REJECTED = "FINAL13_Ambiguity_Rejected.csv"
OUT_MANUAL = "FINAL13_Ambiguity_Manual_Review.csv"
OUT_CORRECTED = "Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv"
OUT_SYMBOL_SUMMARY = "FINAL13_Ambiguous_Official_Symbol_Summary.csv"
OUT_REPORT = "FINAL13_Integrity_Report.txt"

# Standardize whitespace and missing values in text.
def clean(value):
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()

# Normalize a term for case-insensitive comparison.
def norm(value):
    return clean(value).casefold()

# Split a semicolon-delimited term field into names.
def split_terms(value):
    seen = set()
    out = []

    for item in clean(value).split(";"):
        term = clean(item)
        key = norm(term)

        if term and key not in seen:
            seen.add(key)
            out.append(term)

    return out

# Normalize a phrase for contextual comparison.
def phrase_norm(value):
    value = re.sub(
        r"[^a-z0-9]+",
        " ",
        clean(value).lower(),
    )
    return " ".join(value.split())

# Build a regex allowing controlled punctuation variants.
@lru_cache(maxsize=50000)
def flexible_term_pattern(term):
    """
    Conservative punctuation normalization used throughout the fresh run.
    """
    term = clean(term)

    if not term:
        return re.compile(r"(?!x)x")

    pieces = []
    i = 0

    while i < len(term):
        ch = term[i]

        if re.match(r"[\s\-_/\.]", ch):
            while (
                i + 1 < len(term)
                and re.match(r"[\s\-_/\.]", term[i + 1])
            ):
                i += 1

            pieces.append(r"[\s\-_/\.]*")
            i += 1
            continue

        pieces.append(re.escape(ch))

        if i + 1 < len(term):
            nxt = term[i + 1]

            if (
                ch.isalnum()
                and nxt.isalnum()
                and (
                    (ch.isalpha() and nxt.isdigit())
                    or (ch.isdigit() and nxt.isalpha())
                )
            ):
                pieces.append(r"[\s\-_/\.]*")

        i += 1

    return re.compile(
        rf"(?<![A-Za-z0-9])"
        rf"{''.join(pieces)}"
        rf"(?![A-Za-z0-9])",
        re.IGNORECASE,
    )

# Separate a title or abstract into sentences.
@lru_cache(maxsize=50000)
def split_sentences(title, abstract):
    title = clean(title)
    abstract = clean(abstract)

    text = title
    if abstract:
        text = (
            f"{title}. {abstract}"
            if title
            else abstract
        )

    return tuple(
        sentence.strip()
        for sentence in re.split(
            r"(?<=[.!?])\s+",
            text,
        )
        if sentence.strip()
    )

# Find sentences containing the target gene symbol.
def symbol_sentences(row, symbol):
    pattern = flexible_term_pattern(symbol)

    return [
        sentence
        for sentence in split_sentences(
            clean(row["Title"]),
            clean(row["Abstract"]),
        )
        if pattern.search(sentence)
    ]

# Identify short official symbols needing extra evidence.
def is_high_risk_short_official_symbol(symbol):
    """
    Reproducible audit-universe rule:
      current official symbol = letters only, length 2-4.
    """
    return bool(
        re.fullmatch(
            r"[A-Za-z]{2,4}",
            clean(symbol),
        )
    )

# Check whether only a short symbol supports the assignment.
def supported_only_by_official_symbol(row):
    official = norm(
        row["Official_Gene"]
    )

    terms = split_terms(
        row["Validated_Search_Terms"]
    )

    return (
        bool(terms)
        and all(
            norm(term) == official
            for term in terms
        )
    )

ANCHOR_STOPWORDS = {
    "protein",
    "receptor",
    "receptors",
    "gene",
    "genes",
    "human",
    "homolog",
    "family",
    "member",
    "type",
    "subtype",
    "coupled",
    "associated",
    "related",
    "factor",
    "antigen",
    "activity",
    "isoform",
    "precursor",
    "putative",
    "growth",
    "stimulating",
    "extracellular",
    "dependent",
    "inhibitory",
    "glucose",
    "gastric",
    "polypeptide",
    "preferring",
    "pigment",
    "visual",
}

# Tokenize the target receptor name for matching.
def tokenize_target_term(term):
    return [
        token.lower()
        for token in re.findall(
            r"[A-Za-z][A-Za-z0-9]{4,}",
            clean(term),
        )
        if token.lower()
        not in ANCHOR_STOPWORDS
    ]

# Collect active terms for the ambiguity audit.
def build_term_catalog(df):
    catalog = {}

    for gene, group in df.groupby(
        "Official_Gene",
        sort=False,
    ):
        row = group.iloc[0]

        terms = (
            split_terms(
                row["Final_KEEP_GPCR_Terms"]
            )
            + split_terms(
                row["Final_CONTEXT_REQUIRED_GPCR_Terms"]
            )
        )

        seen = set()
        ordered = []

        for term in terms:
            key = norm(term)

            if key and key not in seen:
                seen.add(key)
                ordered.append(term)

        catalog[clean(gene)] = ordered

    return catalog

# Count receptor anchor terms across publications.
def build_anchor_document_frequency(catalog):
    dfreq = Counter()

    for _, terms in catalog.items():
        tokens = set()

        for term in terms:
            tokens.update(
                tokenize_target_term(term)
            )

        for token in tokens:
            dfreq[token] += 1

    return dfreq

# Find receptor-specific anchor tokens.
def target_anchor_tokens(
    gene,
    catalog,
    anchor_dfreq,
):
    """
    Rare terminology-derived target anchors.

    A token is considered a strong anchor if it occurs in final terminology
    for <=5 GPCRs.
    """
    out = set()

    for term in catalog.get(
        gene,
        [],
    ):
        if norm(term) == norm(gene):
            continue

        for token in tokenize_target_term(
            term
        ):
            if anchor_dfreq[token] <= 5:
                out.add(token)

    return sorted(out)

# Find another unambiguous target receptor name.
def find_alternate_target_term(
    row,
    gene,
    catalog,
):
    text = (
        clean(row["Title"])
        + " "
        + clean(row["Abstract"])
    )

    matches = []

    for term in catalog.get(
        gene,
        [],
    ):
        if norm(term) == norm(gene):
            continue

        compact = re.sub(
            r"[^A-Za-z0-9]",
            "",
            term,
        )

        if (
            " " not in term
            and len(compact) <= 4
        ):
            continue

        if flexible_term_pattern(
            term
        ).search(text):
            matches.append(term)

    return matches

# Find anchor tokens in the publication text.
def find_target_anchor_tokens(
    row,
    gene,
    catalog,
    anchor_dfreq,
):
    text = (
        clean(row["Title"])
        + " "
        + clean(row["Abstract"])
    )

    matches = []

    for token in target_anchor_tokens(
        gene,
        catalog,
        anchor_dfreq,
    ):
        if re.search(
            rf"\b{re.escape(token)}\b",
            text,
            flags=re.IGNORECASE,
        ):
            matches.append(token)

    return sorted(
        set(matches),
        key=str.casefold,
    )

# Compare an acronym with its preceding expansion.
def infer_preceding_initialism(
    text,
    acronym,
):
    """
    Find expansions of the form:
        metastatic lymph node ratio (MLNR)

    The shortest suffix whose word initials reproduce the acronym is returned.
    """
    text = re.sub(
        r"<[^>]+>",
        " ",
        clean(text),
    )

    expansions = []

    for match in re.finditer(
        rf"\(\s*{re.escape(acronym)}\s*\)",
        text,
        flags=re.IGNORECASE,
    ):
        prefix = text[:match.start()]
        chunk = prefix[-250:]

        raw_tokens = [
            item.group(0)
            for item in re.finditer(
                r"[A-Za-z][A-Za-z0-9-]*",
                chunk,
            )
        ]

        for size in range(
            1,
            min(8, len(raw_tokens)) + 1,
        ):
            suffix = raw_tokens[-size:]

            initials = ""

            for token in suffix:
                parts = re.findall(
                    r"[A-Za-z]+",
                    token,
                )

                initials += "".join(
                    part[0]
                    for part in parts
                    if part
                )

            if norm(initials) == norm(acronym):
                expansions.append(
                    " ".join(suffix)
                )
                break

    return expansions

# Measure whether an acronym expansion names the receptor.
def expansion_target_similarity(
    expansion,
    gene,
    catalog,
):
    expansion_n = phrase_norm(
        expansion
    )

    best_score = 0.0
    best_term = ""

    for term in catalog.get(
        gene,
        [],
    ):
        if norm(term) == norm(gene):
            continue

        term_n = phrase_norm(
            term
        )

        if not term_n:
            continue

        sequence_score = (
            SequenceMatcher(
                None,
                expansion_n,
                term_n,
            ).ratio()
        )

        expansion_tokens = set(
            expansion_n.split()
        )
        term_tokens = set(
            term_n.split()
        )

        union = (
            expansion_tokens
            | term_tokens
        )

        jaccard = (
            len(
                expansion_tokens
                & term_tokens
            )
            / len(union)
            if union
            else 0.0
        )

        score = max(
            sequence_score,
            jaccard,
        )

        if score > best_score:
            best_score = score
            best_term = term

    return best_score, best_term

COLLISION_PATTERNS = {
    "SMO": re.compile(
        r"(?:"
        r"\bspermine oxidase\b|"
        r"\bpolyamine catabolic\b|"
        r"\bpolyamine metabolism\b|"
        r"\bPAOh1\b|"
        r"\bsequential[- ]minimal[- ]optimization\b|"
        r"\bSVM-SMO\b|"
        r"\bsemiconductor metal oxide\b|"
        r"\bSMOTE\b|"
        r"classifiers?\s*\([^)]*\bSMO\b|"
        r"synchronous metastases? to the ovary"
        r")",
        re.IGNORECASE,
    ),

    "RRH": re.compile(
        r"(?:"
        r"\brobotic radical hysterectomy\b|"
        r"\bradical hysterectomy\b|"
        r"\bregional referral hospital\b|"
        r"\bGulu RRH\b|"
        r"\bRRH\+PLH\b|"
        r"\bright hemicolectomy\b|"
        r"\brubber workers\b|"
        r"\bRRH\s*[=:]?\s*\d+(?:\.\d+)?\b|"
        r"medium/high level of exposure\s*\(RRH\)"
        r")",
        re.IGNORECASE,
    ),

    "MLNR": re.compile(
        r"(?:"
        r"\bmetastatic lymph nodes? ratio\b|"
        r"\bmodified lymph node ratio\b|"
        r"\blymph nodes? ratio\b|"
        r"\bMCF7-MLNr\b|"
        r"\bmelphalan-resistant\b|"
        r"melphalan\s*\(MLNR\)|"
        r"\bL-PAM\b"
        r")",
        re.IGNORECASE,
    ),

    "RGR": re.compile(
        r"(?:"
        r"\brelative growth rate\b|"
        r"\brobot-assisted gastric resection\b|"
        r"\bRGR peptides?\b|"
        r"\bpeptides? RGR\b"
        r")",
        re.IGNORECASE,
    ),

    "CASR": re.compile(
        r"(?:"
        r"\bage[- ]standardized (?:death )?rate\b|"
        r"\bage[- ]standardised (?:death )?rate\b|"
        r"\bstandard population\s*\(CASR\)|"
        r"\bChinese population\s*\(CASR\)|"
        r"\bCASR\s+(?:was|of|and)\s+\d"
        r")",
        re.IGNORECASE,
    ),

    "SCTR": re.compile(
        r"(?:"
        r"\bshort-course radiotherapy\b|"
        r"\bSCTR-delay\b|"
        r"\bSaudi Food and Drug Authority\b|"
        r"\bSCTR\s+\d{5,}\b"
        r")",
        re.IGNORECASE,
    ),

    "NMBR": re.compile(
        r"(?:"
        r"\bNottingham Modification of Scarff Bloom-Richardson\b|"
        r"\bNMBR grade\b"
        r")",
        re.IGNORECASE,
    ),

    "GCGR": re.compile(
        r"\bGCGR-E55\b",
        re.IGNORECASE,
    ),

    "TSHR": re.compile(
        r"(?:"
        r"\bhost range \(?tshr\)? phenotype\b|"
        r"\btshr phenotype\b"
        r")",
        re.IGNORECASE,
    ),
}

SMO_TARGET_CONTEXT = re.compile(
    r"\b(?:"
    r"hedgehog|"
    r"sonic hedgehog|"
    r"SHH|"
    r"GLI-?1|"
    r"GLI-?2|"
    r"PTCH1?|"
    r"patched|"
    r"cyclopamine|"
    r"vismodegib|"
    r"sonidegib|"
    r"glasdegib|"
    r"GDC-0449|"
    r"LDE225|"
    r"SAG"
    r")\b",
    re.IGNORECASE,
)

GENE_CONTEXT = re.compile(
    r"\b(?:"
    r"gene|genes|"
    r"mRNA|RNA-seq|RNA sequencing|"
    r"transcript|transcripts|transcriptome|"
    r"expression|express|expressed|"
    r"methylation|methylated|promoter|"
    r"mutation|mutations|variant|variants|"
    r"polymorphism|polymorphisms|haplotype|"
    r"SNP|locus|loci|"
    r"DEG|DEGs|IRG|IRGs|"
    r"differentially expressed|"
    r"hub|hub gene|"
    r"signature|biomarker|marker|markers|"
    r"prognostic|predictor|"
    r"candidate|candidates|"
    r"knockdown|silencing|"
    r"overexpression|upregulated|downregulated|"
    r"receptor|GPCR|G protein-coupled|"
    r"protein|proteins|"
    r"agonist|agonists|"
    r"antagonist|antagonists|"
    r"inhibitor|inhibitors|inhibit|inhibition|"
    r"ligand|ligands|radioligand|radioligands|"
    r"target|targets|targeted|targeting|"
    r"theragnostic|"
    r"immunostaining|immunohistochemistry|"
    r"transcription factor|"
    r"oncogene|oncogenes|"
    r"deletion|deletions|"
    r"autoantibody|autoantibodies|"
    r"risk score|activity|activation|"
    r"diagnostic performance|marker panel"
    r")\b",
    re.IGNORECASE,
)

# Check for explicit gene or receptor context.
def has_gene_context(
    row,
    symbol,
):
    return any(
        GENE_CONTEXT.search(sentence)
        for sentence in symbol_sentences(
            row,
            symbol,
        )
    )

# Check for a structured list of gene symbols.
def has_structured_gene_list_context(
    row,
    symbol,
):
    """
    Detect contexts such as:
        GDF3, CASR, CLDN10, and COL9A2
    where the short symbol is functioning as a gene identifier.
    """
    for sentence in symbol_sentences(
        row,
        symbol,
    ):
        plain = re.sub(
            r"<[^>]+>",
            " ",
            sentence,
        )

        tokens = set(
            re.findall(
                r"(?<![A-Za-z0-9])"
                r"([A-Z][A-Z0-9-]{1,9})"
                r"(?![A-Za-z0-9])",
                plain,
            )
        )

        tokens -= {
            "CI",
            "OR",
            "AND",
            "THE",
            "AUC",
            "PCR",
            "OS",
            "HR",
            "RR",
            "P",
        }

        if (
            symbol in tokens
            and len(tokens) >= 3
        ):
            return True

    return False

# Apply the strict rule for ambiguous RHO mentions.
def rho_strict_decision(row):
    """
    Fresh-run article-level RHO audit.

    The final RHO candidate set was inspected after QC. Candidate records using
    "RHO" were dominated by Rho GTPase / RHOA-RHOC terminology. Candidate
    records using "Rhodopsin" were generic rhodopsin-family classifications,
    structural comparators, retinal toxicity/retinopathy contexts, or assays
    using rhodopsin rather than evidence that RHO itself was associated with
    the queried cancer.

    No RHO candidate in this fresh final set demonstrated receptor-specific
    cancer evidence; therefore RHO rows are removed from the cancer-association
    resource.
    """
    return (
        "REJECT",
        "RHO_NO_RECEPTOR_SPECIFIC_CANCER_EVIDENCE",
        (
            "Fresh-run article-level audit found RHO/Rhodopsin usage to be "
            "Rho-GTPase terminology or generic rhodopsin-family/comparator/"
            "retinal context rather than receptor-specific cancer evidence."
        ),
    )

# Check whether a short-symbol assignment has support.
def classify_ambiguity_row(
    row,
    catalog,
    anchor_dfreq,
):
    gene = clean(
        row["Official_Gene"]
    )

    if gene == "RHO":
        decision, reason, comment = (
            rho_strict_decision(row)
        )

        return {
            "Ambiguity_Decision": decision,
            "Ambiguity_Reason": reason,
            "Ambiguity_Comment": comment,
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": "",
            "Best_Expansion_Target_Similarity": "",
        }

    if not is_high_risk_short_official_symbol(
        gene
    ):
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "NOT_SHORT_ALPHABETIC_OFFICIAL_SYMBOL",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": "",
            "Best_Expansion_Target_Similarity": "",
        }

    if not supported_only_by_official_symbol(
        row
    ):
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "INDEPENDENT_VALIDATED_GPCR_TERM_PRESENT",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": "",
            "Best_Expansion_Target_Similarity": "",
        }

    text = (
        clean(row["Title"])
        + " "
        + clean(row["Abstract"])
    )

    target_terms = find_alternate_target_term(
        row,
        gene,
        catalog,
    )

    if target_terms:
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "INDEPENDENT_TARGET_TERM",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": ";".join(
                target_terms
            ),
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": "",
            "Best_Expansion_Target_Similarity": "",
        }

    anchors = find_target_anchor_tokens(
        row,
        gene,
        catalog,
        anchor_dfreq,
    )

    if anchors:
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "TARGET_SPECIFIC_ANCHOR",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": ";".join(
                anchors
            ),
            "Acronym_Expansions_Found": "",
            "Best_Expansion_Target_Similarity": "",
        }

    expansions = infer_preceding_initialism(
        text,
        gene,
    )

    incompatible = []
    best_similarity = 0.0

    for expansion in expansions:
        score, best_term = (
            expansion_target_similarity(
                expansion,
                gene,
                catalog,
            )
        )

        best_similarity = max(
            best_similarity,
            score,
        )

        if score < 0.65:
            incompatible.append(
                expansion
            )

    if incompatible:
        return {
            "Ambiguity_Decision": "REJECT",
            "Ambiguity_Reason": "NON_GPCR_ACRONYM_EXPANSION",
            "Ambiguity_Comment": (
                "Short official symbol is explicitly expanded to a "
                "non-target meaning."
            ),
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": ";".join(
                incompatible
            ),
            "Best_Expansion_Target_Similarity": (
                f"{best_similarity:.4f}"
            ),
        }

    pattern = COLLISION_PATTERNS.get(
        gene
    )

    if (
        pattern is not None
        and pattern.search(text)
    ):
        return {
            "Ambiguity_Decision": "REJECT",
            "Ambiguity_Reason": "RECURRENT_NON_GPCR_COLLISION_CONTEXT",
            "Ambiguity_Comment": (
                "Fresh-run corpus contains a recurrent non-GPCR use of "
                "this short official symbol."
            ),
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": ";".join(
                expansions
            ),
            "Best_Expansion_Target_Similarity": (
                f"{best_similarity:.4f}"
                if expansions
                else ""
            ),
        }

    if (
        gene == "SMO"
        and SMO_TARGET_CONTEXT.search(text)
    ):
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "SMO_HEDGEHOG_TARGET_CONTEXT",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": ";".join(
                expansions
            ),
            "Best_Expansion_Target_Similarity": (
                f"{best_similarity:.4f}"
                if expansions
                else ""
            ),
        }

    if has_gene_context(
        row,
        gene,
    ):
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "EXPLICIT_GENE_RECEPTOR_CONTEXT",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": ";".join(
                expansions
            ),
            "Best_Expansion_Target_Similarity": (
                f"{best_similarity:.4f}"
                if expansions
                else ""
            ),
        }

    if has_structured_gene_list_context(
        row,
        gene,
    ):
        return {
            "Ambiguity_Decision": "KEEP",
            "Ambiguity_Reason": "STRUCTURED_GENE_LIST_CONTEXT",
            "Ambiguity_Comment": "",
            "Alternate_Target_Terms_Found": "",
            "Target_Anchor_Tokens_Found": "",
            "Acronym_Expansions_Found": ";".join(
                expansions
            ),
            "Best_Expansion_Target_Similarity": (
                f"{best_similarity:.4f}"
                if expansions
                else ""
            ),
        }

    return {
        "Ambiguity_Decision": "MANUAL_REVIEW",
        "Ambiguity_Reason": "SHORT_SYMBOL_CONTEXT_UNRESOLVED",
        "Ambiguity_Comment": (
            "No independent target evidence, explicit non-target expansion, "
            "known collision context, or explicit gene/receptor context."
        ),
        "Alternate_Target_Terms_Found": "",
        "Target_Anchor_Tokens_Found": "",
        "Acronym_Expansions_Found": ";".join(
            expansions
        ),
        "Best_Expansion_Target_Similarity": (
            f"{best_similarity:.4f}"
            if expansions
            else ""
        ),
    }

REQUIRED_COLUMNS = {
    "UniProt_ID",
    "Official_Gene",
    "TCGA_Code",
    "Cancer",
    "PMID",
    "Title",
    "Abstract",
    "Validated_Search_Terms",
    "Final_KEEP_GPCR_Terms",
    "Final_CONTEXT_REQUIRED_GPCR_Terms",
    "Final_Literature_Status",
}

# Load input from the input records.
def load_input():
    if not os.path.exists(
        INPUT_FILE
    ):
        raise FileNotFoundError(
            f"Missing required input: {INPUT_FILE}"
        )

    df = pd.read_csv(
        INPUT_FILE,
        dtype={"PMID": str},
        low_memory=False,
    ).fillna("")

    missing = sorted(
        REQUIRED_COLUMNS
        - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Input is missing required column(s): "
            + ", ".join(missing)
        )

    if df.empty:
        raise ValueError(
            "Input final literature table is empty."
        )

    if not df[
        "Final_Literature_Status"
    ].map(clean).eq("KEEP").all():
        raise ValueError(
            "Input contains rows not marked "
            "Final_Literature_Status=KEEP."
        )

    if df.duplicated(
        [
            "UniProt_ID",
            "TCGA_Code",
            "PMID",
        ]
    ).any():
        raise ValueError(
            "Duplicate GPCR-cancer-PMID rows in input."
        )

    if df[
        "Validated_Search_Terms"
    ].map(clean).eq("").any():
        raise ValueError(
            "Input contains blank Validated_Search_Terms."
        )

    return df

# Build symbol summary from the input records.
def build_symbol_summary(
    audited,
    input_df,
    corrected,
):
    symbols = sorted(
        {
            clean(gene)
            for gene in input_df[
                "Official_Gene"
            ].unique()
            if is_high_risk_short_official_symbol(
                gene
            )
        }
        | {"RHO"},
        key=str.casefold,
    )

    rows = []

    for gene in symbols:
        before = input_df[
            input_df[
                "Official_Gene"
            ].eq(gene)
        ]

        audit_gene = audited[
            audited[
                "Official_Gene"
            ].eq(gene)
        ]

        after = corrected[
            corrected[
                "Official_Gene"
            ].eq(gene)
        ]

        rows.append({
            "Official_Gene": gene,
            "Rows_Before": len(before),
            "Cancer_Types_Before":
                before["TCGA_Code"].nunique(),
            "Audited_Rows": len(audit_gene),
            "Kept_Audited_Rows": int(
                audit_gene[
                    "Ambiguity_Decision"
                ].eq("KEEP").sum()
            ),
            "Rejected_Rows": int(
                audit_gene[
                    "Ambiguity_Decision"
                ].eq("REJECT").sum()
            ),
            "Manual_Review_Rows": int(
                audit_gene[
                    "Ambiguity_Decision"
                ].eq(
                    "MANUAL_REVIEW"
                ).sum()
            ),
            "Rows_After": len(after),
            "Cancer_Types_After":
                after["TCGA_Code"].nunique(),
        })

    return pd.DataFrame(rows)

# Write report for the current record.
def write_report(
    input_df,
    audited,
    rejected,
    manual,
    corrected,
    symbol_summary,
):
    input_keys = set(
        map(
            tuple,
            input_df[
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ].values,
        )
    )

    corrected_keys = set(
        map(
            tuple,
            corrected[
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ].values,
        )
    )

    rejected_keys = set(
        map(
            tuple,
            rejected[
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ].values,
        )
    )

    manual_keys = set(
        map(
            tuple,
            manual[
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ].values,
        )
    )

    resolved_partition = (
        corrected_keys
        | rejected_keys
        | manual_keys
    )

    checks = {
        "No duplicate input association rows":
            len(input_keys)
            == len(input_df),

        "Every input row resolved exactly once":
            (
                len(
                    corrected_keys
                    & rejected_keys
                ) == 0
                and len(
                    corrected_keys
                    & manual_keys
                ) == 0
                and len(
                    rejected_keys
                    & manual_keys
                ) == 0
                and resolved_partition
                == input_keys
            ),

        "Corrected output has no duplicates":
            not corrected.duplicated(
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ).any(),

        "All corrected rows remain KEEP literature":
            corrected[
                "Final_Literature_Status"
            ].map(clean).eq("KEEP").all(),

        "All rejected rows have REJECT decision":
            rejected[
                "Ambiguity_Decision"
            ].eq("REJECT").all(),

        "Manual queue matches MANUAL_REVIEW":
            manual[
                "Ambiguity_Decision"
            ].eq("MANUAL_REVIEW").all(),

        "No unresolved ambiguity rows":
            len(manual) == 0,
    }

    failed = [
        label
        for label, passed in checks.items()
        if not passed
    ]

    lines = [
        "GPCR Cancer Pipeline - Final Fresh Run Step 13 Integrity Report",
        "=" * 78,
        "",
        "INPUT",
        "-" * 78,
        f"Step 09 final literature rows         : {len(input_df)}",
        f"GPCRs before ambiguity repair         : "
        f"{input_df['Official_Gene'].nunique()}",
        f"GPCR-cancer pairs before repair       : "
        f"{input_df[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]}",
        f"Unique PMIDs before repair            : "
        f"{input_df['PMID'].nunique()}",
        "",
        "AMBIGUITY AUDIT",
        "-" * 78,
        f"Rows explicitly audited               : {len(audited)}",
        f"Rows rejected                          : {len(rejected)}",
        f"Rows unresolved/manual review          : {len(manual)}",
        "",
        "CORRECTED VALIDATED LITERATURE",
        "-" * 78,
        f"Rows retained                          : {len(corrected)}",
        f"GPCRs retained                         : "
        f"{corrected['Official_Gene'].nunique()}",
        f"GPCR-cancer pairs retained             : "
        f"{corrected[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]}",
        f"Unique PMIDs retained                  : "
        f"{corrected['PMID'].nunique()}",
        "",
        "INTEGRITY CHECKS",
        "-" * 78,
    ]

    for label, passed in checks.items():
        lines.append(
            f"{label}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    lines += [
        "",
        "METHOD NOTE",
        "-" * 78,
        (
            "Step 13 re-audits short alphabetic official symbols because "
            "official-symbol ownership does not guarantee lexical specificity. "
            "The audit is performed at the article level using final target "
            "terminology, acronym-expansion checks, recurrent non-GPCR "
            "collision contexts identified in the fresh-run corpus, and "
            "explicit gene/receptor usage. RHO is handled by a strict "
            "fresh-run article audit because both RHO and Rhodopsin candidate "
            "records lacked receptor-specific cancer evidence."
        ),
        "",
        f"Integrity failures                    : {len(failed)}",
    ]

    if failed:
        lines += [
            "",
            "FAILED CHECKS",
            "-" * 78,
        ]
        lines.extend(
            f"- {item}"
            for item in failed
        )
    else:
        lines += [
            "",
            "ALL STEP 13 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return failed

# Run this stage and write its output files.
def main(decisions_path=None):
    Path(OUT_CORRECTED).unlink(missing_ok=True)
    print("=" * 78)
    print("GPCR Cancer Pipeline - Final Fresh Run Step 13")
    print("Final Ambiguous Official-Symbol Audit and Repair")
    print("=" * 78)

    try:
        df = load_input()

        catalog = build_term_catalog(
            df
        )

        anchor_dfreq = (
            build_anchor_document_frequency(
                catalog
            )
        )

    except Exception as exc:
        print(
            f"STEP 13 INPUT VALIDATION FAILED:\n{exc}"
        )
        sys.exit(1)

    short_symbols = sorted(
        {
            clean(gene)
            for gene in df[
                "Official_Gene"
            ].unique()
            if is_high_risk_short_official_symbol(
                gene
            )
        },
        key=str.casefold,
    )

    print(
        f"Input rows                    : "
        f"{len(df):,}"
    )
    print(
        f"Input GPCRs                   : "
        f"{df['Official_Gene'].nunique():,}"
    )
    print(
        f"Short alphabetic symbols      : "
        f"{len(short_symbols):,}"
    )
    print(
        "Symbols                       : "
        + ";".join(short_symbols)
    )
    print("Network requests              : 0")
    print("=" * 78)

    resolution_rows = []
    audit_indices = []

    for number, (index, row) in enumerate(
        df.iterrows(),
        start=1,
    ):
        gene = clean(
            row["Official_Gene"]
        )

        explicitly_audited = (
            gene == "RHO"
            or (
                is_high_risk_short_official_symbol(
                    gene
                )
                and supported_only_by_official_symbol(
                    row
                )
            )
        )

        if explicitly_audited:
            audit_indices.append(index)

            result = classify_ambiguity_row(
                row,
                catalog,
                anchor_dfreq,
            )
        else:
            result = {
                "Ambiguity_Decision": "KEEP",
                "Ambiguity_Reason": (
                    "OUTSIDE_EXPLICIT_AMBIGUITY_AUDIT"
                ),
                "Ambiguity_Comment": "",
                "Alternate_Target_Terms_Found": "",
                "Target_Anchor_Tokens_Found": "",
                "Acronym_Expansions_Found": "",
                "Best_Expansion_Target_Similarity": "",
            }

        resolution_rows.append(
            result
        )

        if (
            number == 1
            or number % 5000 == 0
            or number == len(df)
        ):
            print(
                f"Auditing "
                f"{number:,}/{len(df):,}"
            )

    resolution = pd.DataFrame(
        resolution_rows,
        index=df.index,
    )

    resolved = pd.concat(
        [
            df,
            resolution,
        ],
        axis=1,
    )

    resolved = apply_ambiguity_reviews(
        resolved, decisions_path, "AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv"
    )

    audited = resolved.loc[
        audit_indices
    ].copy()

    rejected = resolved[
        resolved[
            "Ambiguity_Decision"
        ].eq("REJECT")
    ].copy()

    manual = resolved[
        resolved[
            "Ambiguity_Decision"
        ].eq("MANUAL_REVIEW")
    ].copy()

    corrected = resolved[
        resolved[
            "Ambiguity_Decision"
        ].eq("KEEP")
    ].copy()

    corrected[
        "Final_Ambiguity_QC_Status"
    ] = "KEEP"

    corrected[
        "Final_Ambiguity_QC_Reason"
    ] = corrected[
        "Ambiguity_Reason"
    ]

    symbol_summary = (
        build_symbol_summary(
            audited,
            df,
            corrected,
        )
    )

    audited.to_csv(
        OUT_AUDIT,
        index=False,
        na_rep="",
    )

    rejected.to_csv(
        OUT_REJECTED,
        index=False,
        na_rep="",
    )

    manual.to_csv(
        OUT_MANUAL,
        index=False,
        na_rep="",
    )

    symbol_summary.to_csv(
        OUT_SYMBOL_SUMMARY,
        index=False,
        na_rep="",
    )

    failed = write_report(
        input_df=df,
        audited=audited,
        rejected=rejected,
        manual=manual,
        corrected=corrected,
        symbol_summary=symbol_summary,
    )

    if len(manual) == 0:
        corrected.to_csv(
            OUT_CORRECTED,
            index=False,
            na_rep="",
        )
    else:
        Path(OUT_CORRECTED).unlink(missing_ok=True)
        print(
            "\nCorrected final literature file was NOT written because "
            "manual-review rows remain."
        )

    counts = (
        audited[
            "Ambiguity_Decision"
        ]
        .value_counts()
        .to_dict()
    )

    print()
    print("=" * 78)
    print(
        "STEP 13 COMPLETE"
        if not failed
        else "STEP 13 COMPLETE WITH INTEGRITY ISSUES"
    )
    print("=" * 78)
    print(
        f"Explicitly audited rows       : "
        f"{len(audited):,}"
    )
    print(
        f"Audited KEEP                 : "
        f"{counts.get('KEEP', 0):,}"
    )
    print(
        f"Audited REJECT               : "
        f"{counts.get('REJECT', 0):,}"
    )
    print(
        f"Manual review                : "
        f"{counts.get('MANUAL_REVIEW', 0):,}"
    )
    print(
        f"Corrected literature rows    : "
        f"{len(corrected):,}"
    )
    print(
        f"Corrected GPCRs              : "
        f"{corrected['Official_Gene'].nunique():,}"
    )
    print(
        f"Corrected GPCR-cancer pairs  : "
        f"{corrected[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]:,}"
    )
    print(
        f"Corrected unique PMIDs       : "
        f"{corrected['PMID'].nunique():,}"
    )
    print(
        f"Integrity errors             : "
        f"{len(failed)}"
    )
    print()
    print("Saved:")
    print(f"  {OUT_AUDIT}")
    print(f"  {OUT_REJECTED}")
    print(f"  {OUT_MANUAL}")
    print(f"  {OUT_SYMBOL_SUMMARY}")
    print(f"  {OUT_REPORT}")

    if len(manual) == 0:
        print(f"  {OUT_CORRECTED}")

    if len(manual):
        print("PAUSED: complete AMBIGUITY_REVIEW_DECISIONS_TO_COMPLETE.csv and resume with --ambiguity-decisions.")
        sys.exit(REVIEW_EXIT)

    if failed:
        sys.exit(2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Audit ambiguous symbols and apply current manual decisions")
    parser.add_argument("--decisions", type=Path, help="Completed ambiguity review decisions")
    args = parser.parse_args()
    main(args.decisions)
