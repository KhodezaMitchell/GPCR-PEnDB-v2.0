#!/usr/bin/env python3
# Purpose: Assign biological and clinical categories using GPCR-centered text context.
# Input: Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv
# Output: CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv; CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Associations_Retained.csv

import os
import re
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

INPUT_FILE = "Human_GPCR_Cancer_PubMed_Papers_FINAL_CORRECTED.csv"

OUT_ALL = "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Classified.csv"
OUT_RETAINED = "CORRECTED_FINAL_Human_GPCR_Cancer_PubMed_Associations_Retained.csv"
OUT_SUMMARY = "CORRECTED_FINAL_GPCR_Cancer_Association_Classification_Summary.csv"
OUT_REPORT = "CORRECTED_FINAL_GPCR_Cancer_Association_Classification_Integrity_Report.txt"

CATEGORY_KEYWORDS = {
    "Proliferation/Apoptosis": [
        "proliferation", "proliferative", "cell proliferation",
        "apoptosis", "apoptotic", "cell cycle", "cell viability",
        "senescence", "cell death",
    ],
    "Migration/Invasion/Metastasis": [
        "migration", "migratory", "invasion", "invasive",
        "metastasis", "metastatic", "motility", "emt",
        "epithelial-mesenchymal", "epithelial mesenchymal",
    ],
    "Tumor Growth/Progression": [
        "tumor growth", "tumour growth", "cancer progression",
        "tumor progression", "tumour progression", "disease progression",
        "progression gene", "tumorigenesis", "tumourigenesis",
        "tumorigenic", "tumourigenic", "oncogenic", "oncogenesis",
        "malignant transformation", "cancer development",
        "tumor development", "tumour development", "stemness",
        "angiogenesis", "colony formation", "antitumor", "anti-tumor",
        "antitumour", "anti-tumour", "tumor suppressor",
        "tumour suppressor", "tumor-suppressive", "tumour-suppressive",
        "tumor-promoting", "tumour-promoting", "malignant progression",
    ],
    "Biomarker/Diagnosis": [
        "biomarker", "marker", "markers", "diagnostic", "diagnosis",
        "signature", "risk score", "risk model", "predictor",
        "predictive", "prediction", "classifier", "roc", "auc",
    ],
    "Prognosis": [
        "prognosis", "prognostic", "poor outcome",
        "clinical outcome", "recurrence", "recurrent",
    ],
    "Survival": [
        "survival", "overall survival", "disease-free survival",
        "progression-free survival", "relapse-free survival",
        "hazard ratio", "mortality",
    ],
    "Therapeutic Target": [
        "therapeutic target", "drug target", "potential target",
        "treatment target", "agonist", "antagonist", "inhibitor",
        "blockade", "therapy", "treatment", "immunotherapeutic",
        "targeting", "targeted therapy", "knockdown", "silencing",
        "therapeutic potential", "candidate therapeutic",
        "monoclonal antibody", "antibody therapy",
        "drug discovery", "drug development",
    ],
    "Expression": [
        "expression", "expressed", "expressing", "overexpression",
        "overexpressed", "upregulated", "downregulated", "up-regulated",
        "down-regulated", "increased expression", "decreased expression",
        "elevated expression", "reduced expression", "mrna", "transcript",
        "transcripts", "transcriptional", "protein level", "express",
        "upregulation", "downregulation", "immunohistochemistry",
        "immunohistochemical", "immunostain", "staining",
        "immunoreactivity", "positivity", "receptor-positive",
        "receptor positive",
    ],
    "Signaling": [
        "signaling", "signalling", "signal transduction",
        "signaling protein", "signalling protein", "downstream signaling",
        "downstream signalling", "upstream signaling",
        "upstream signalling", "phosphorylation", "phosphorylated",
        "g protein", "g-protein", "beta-arrestin", "β-arrestin",
        "camp", "pi3k", "akt", "mapk", "erk", "rho gtpase",
        "second messenger", "activation", "activated", "activates",
        "inhibition", "inhibited", "inhibits", "pathway", "pathways",
        "axis", "crosstalk", "cross-talk", "constitutive activity",
        "receptor activity", "receptor activities", "internalization",
        "internalisation", "endocytosis", "ubiquitination",
        "ubiquitin-mediated", "desensitization", "desensitisation",
        "trafficking",
    ],
    "Genetic Association": [
        "mutation", "mutant", "variant", "variants", "polymorphism",
        "snp", "snps", "snv", "snvs", "copy number", "cnv",
        "methylation", "hypermethylation", "genotype", "genotypes",
        "allele", "alleles", "gwas", "genome-wide association",
        "amplification", "deletion", "fusion", "rearrangement",
        "somatic alteration", "genomic alteration", "mutated", "locus",
        "loci", "genetic alteration", "genetic association",
        "predisposition gene", "predisposition genes",
    ],
    "Cancer Risk/Susceptibility": [
        "cancer risk", "risk of cancer", "susceptibility", "risk locus",
        "risk loci", "risk allele", "risk of pancreatic cancer",
        "risk of breast cancer", "risk of prostate cancer",
        "risk of colorectal cancer", "risk of colon cancer",
        "risk of rectal cancer", "risk of esophageal cancer",
        "risk of oesophageal cancer", "risk of gastric cancer",
        "risk of ovarian cancer", "risk of lung cancer",
        "risk of bladder cancer", "risk of thyroid cancer",
        "risk of endometrial cancer",
    ],
    "Drug Response/Resistance": [
        "drug resistance", "chemoresistance", "radioresistance",
        "resistant", "resistance", "drug sensitivity", "sensitivity",
        "sensitized", "sensitised", "treatment response",
        "therapy response", "cisplatin", "temozolomide", "5-fu",
        "fluorouracil", "gemcitabine", "paclitaxel", "tamoxifen",
        "trastuzumab",
    ],
    "Immune Microenvironment": [
        "immune", "macrophage", "t cell", "b cell", "nk cell",
        "myeloid", "dendritic", "tumor microenvironment",
        "tumour microenvironment", "stromal", "immune infiltration",
        "immunosuppressive", "immune checkpoint", "til", "tils",
        "tumor-infiltrating lymphocyte", "tumour-infiltrating lymphocyte",
        "lymphocyte infiltration", "immune cell infiltration",
    ],
    "Imaging/Theranostics": [
        "imaging", "targeted imaging", "molecular imaging", "pet",
        "pet/ct", "positron emission tomography", "spect",
        "single photon emission computed tomography", "radiotracer",
        "radiotracers", "radioligand", "radioligands",
        "radiopharmaceutical", "radiopharmaceuticals", "radionuclide",
        "radiolabeled", "radiolabelled", "theranostic", "theranostics",
        "scintigraphy",
    ],
    "Computational/Network Evidence": [
        "hub gene", "hub genes", "ppi", "ppi network",
        "protein-protein interaction network",
        "gene-gene interaction network", "interaction network",
        "network analysis", "wgcna", "cytohubba", "molecular docking",
        "virtual screening", "core gene", "core genes", "central deg",
        "central degs", "bioinformatics analysis",
        "bioinformatic analysis", "computational analysis",
    ],
}

OTHER_CATEGORY = "Other/Unspecified"

ALL_CATEGORIES = list(CATEGORY_KEYWORDS) + [OTHER_CATEGORY]

IMMUNE_CATEGORY = "Immune Microenvironment"

IMMUNE_DIRECT_PHRASES = [
    "tumor microenvironment",
    "tumour microenvironment",
    "cancer microenvironment",
    "immune microenvironment",
    "immune infiltration",
    "immune cell infiltration",
    "lymphocyte infiltration",
    "tumor-infiltrating lymphocyte",
    "tumour-infiltrating lymphocyte",
    "immune checkpoint",
    "immunosuppressive",
    "tumor immunity",
    "tumour immunity",
    "antitumor immunity",
    "anti-tumor immunity",
    "anti-tumour immunity",
    "immune evasion",
    "immune surveillance",
    "tumor-associated macrophage",
    "tumour-associated macrophage",
    "cancer-associated fibroblast",
    "tumor stroma",
    "tumour stroma",
]

IMMUNE_CONTEXT_REQUIRED_TERMS = [
    "immune",
    "macrophage",
    "t cell",
    "b cell",
    "nk cell",
    "myeloid",
    "dendritic",
    "stromal",
    "til",
    "tils",
]

IMMUNE_CONTEXT_QUALIFIERS = [
    "infiltrat",
    "recruit",
    "polariz",
    "exhaust",
    "suppress",
    "checkpoint",
    "immunotherap",
    "cytotoxic",
    "regulatory",
    "antitumor",
    "anti-tumor",
    "anti-tumour",
    "tumor-associated",
    "tumour-associated",
    "cancer-associated",
    "microenvironment",
    "immune response",
    "immune activation",
    "immune regulation",
    "immune escape",
    "immune evasion",
    "immune surveillance",
    "immunity",
    "antigen",
    "lymph node",
    "t-cell therapy",
    "t cell therapy",
    "tcr-t",
    "car t",
]

IMMUNE_MODEL_MARKERS = [
    "jurkat",
    "cem t cell",
    "cem t-cell",
    "molt-4",
    "molt 4",
    "hut 78",
    "sup-t1",
    "sup t1",
]

IMMUNE_FALSE_LEXICAL_PHRASES = [
    "macrophage migration inhibitory factor",
    "stromal cell-derived factor",
    "stromal cell derived factor",
    "stromal-derived factor",
    "stromal derived factor",
    "acute myeloid leukemia",
    "acute myeloid leukaemia",
    "chronic myeloid leukemia",
    "chronic myeloid leukaemia",
    "myeloid leukemia",
    "myeloid leukaemia",
    "b-cell lymphoma",
    "b cell lymphoma",
    "b-cell non-hodgkin",
    "b cell non-hodgkin",
    "t-cell lymphoma",
    "t cell lymphoma",
    "t-cell leukemia",
    "t cell leukemia",
    "t-cell leukaemia",
    "t cell leukaemia",
    "t-cell factor",
    "t cell factor",
    "t-cell death-associated gene",
    "t cell death-associated gene",
    "t-cell death associated gene",
    "t cell death associated gene",
]

IMMUNE_BROAD_PATTERNS = {
    "immune": re.compile(r"\bimmune\b", re.IGNORECASE),
    "macrophage": re.compile(r"\bmacrophages?\b", re.IGNORECASE),
    "t cell": re.compile(r"\bt[\s-]+cells?\b", re.IGNORECASE),
    "b cell": re.compile(r"\bb[\s-]+cells?\b", re.IGNORECASE),
    "nk cell": re.compile(r"\bnk[\s-]+cells?\b", re.IGNORECASE),
    "myeloid": re.compile(r"\bmyeloid\b", re.IGNORECASE),
    "dendritic": re.compile(r"\bdendritic\b", re.IGNORECASE),
    "stromal": re.compile(r"\bstromal\b", re.IGNORECASE),
    "til": re.compile(r"\btils?\b", re.IGNORECASE),
}

BIOMARKER_CATEGORY = "Biomarker/Diagnosis"

BIOMARKER_STRONG_TERMS = [
    "biomarker",
    "diagnostic",
    "diagnosis",
    "signature",
    "risk score",
    "risk model",
    "predictor",
    "predictive",
    "prediction",
    "classifier",
    "roc",
    "auc",
]

BIOMARKER_MARKER_QUALIFIERS = [
    "diagnos",
    "prognos",
    "patient",
    "clinical",
    "survival",
    "outcome",
    "risk",
    "detection",
    "screening",
    "sensitivity",
    "specificity",
    "circulating",
    "serum",
    "plasma",
    "blood",
    "urine",
    "tumor aggressiveness",
    "tumour aggressiveness",
    "molecular marker",
    "cancer marker",
    "tumor marker",
    "tumour marker",
]

BIOMARKER_FALSE_MARKER_PHRASES = [
    "endosomal marker",
    "cell surface marker",
    "stem cell marker",
    "stem-cell marker",
    "differentiation marker",
    "proliferation marker",
    "lineage marker",
    "memory marker",
    "emt marker",
    "organelle marker",
]

SURVIVAL_CATEGORY = "Survival"

SURVIVAL_STRONG_PHRASES = [
    "overall survival",
    "disease-free survival",
    "progression-free survival",
    "relapse-free survival",
    "hazard ratio",
    "mortality",
]

SURVIVAL_FALSE_PATTERNS = [
    re.compile(r"\b(?:cancer\s+|tumou?r\s+)?cells?\s+survival\b", re.I),
    re.compile(r"\bcellular\s+survival\b", re.I),
    re.compile(
        r"\b(?:cancer\s+|tumou?r\s+)?cells?\s+"
        r"(?:proliferation|growth|viability)\s+and\s+survival\b",
        re.I,
    ),
    re.compile(
        r"\bcellular\s+(?:proliferation|growth|viability)\s+"
        r"and\s+survival\b",
        re.I,
    ),
    re.compile(r"\bsurvival\s+of\s+(?:the\s+)?(?:cancer\s+|tumou?r\s+)?cells?\b", re.I),
    re.compile(r"\bpro[\s-]?survival\b", re.I),
    re.compile(r"\bsurvival\s+(?:signaling|signalling|pathways?|factor)\b", re.I),
]

IMAGING_CATEGORY = "Imaging/Theranostics"

IMAGING_STRONG_TERMS = [
    "targeted imaging",
    "molecular imaging",
    "diagnostic imaging",
    "receptor imaging",
    "nuclear imaging",
    "noninvasive imaging",
    "non-invasive imaging",
    "imaging agent",
    "imaging probe",
    "bioimaging",
    "fluorescent probe",
    "luminescent probe",
    "contrast agent",
    "pet/ct",
    "positron emission tomography",
    "spect",
    "single photon emission computed tomography",
    "radiotracer",
    "radiotracers",
    "radioligand",
    "radioligands",
    "radiopharmaceutical",
    "radiopharmaceuticals",
    "radionuclide",
    "radiolabeled",
    "radiolabelled",
    "theranostic",
    "theranostics",
    "scintigraphy",
]

IMAGING_CONTEXT_QUALIFIERS = [
    "tumor",
    "tumour",
    "cancer",
    "xenograft",
    "metast",
    "diagnos",
    "screen",
    "tracer",
    "probe",
    "target",
    "in vivo",
    "patient",
    "clinical",
    "therapeutic",
    "biodistribution",
]

PRIMARY_PRIORITY = [
    "Therapeutic Target",
    "Drug Response/Resistance",
    "Imaging/Theranostics",
    "Biomarker/Diagnosis",
    "Prognosis",
    "Survival",
    "Tumor Growth/Progression",
    "Migration/Invasion/Metastasis",
    "Proliferation/Apoptosis",
    "Cancer Risk/Susceptibility",
    "Genetic Association",
    "Signaling",
    "Expression",
    "Immune Microenvironment",
    "Computational/Network Evidence",
    OTHER_CATEGORY,
]

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

# Split article sentences for the current record.
@lru_cache(maxsize=50000)
def split_article_sentences(title, abstract):
    title = clean(title)
    abstract = clean(abstract)

    text = title
    if abstract:
        text = f"{title}. {abstract}" if title else abstract

    return tuple(
        sentence.strip()
        for sentence in re.split(
            r"(?<=[.!?])\s+",
            text,
        )
        if sentence.strip()
    )

# Build a regex allowing controlled punctuation variants.
@lru_cache(maxsize=20000)
def flexible_term_pattern(term):
    """
    Preserve the conservative punctuation/case matching used in the validated
    literature pipeline. This classification step cannot redefine which
    terms validated the paper.
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

# Find active GPCR terms in the selected text.
def matching_gpcr_terms(sentence, terms):
    return [
        term
        for term in terms
        if flexible_term_pattern(term).search(sentence)
    ]

ANAPHORIC_PATTERN = re.compile(
    r"^(?:"
    r"this|these|such|it|its|they|their|"
    r"these\s+(?:findings|results|data|observations)|"
    r"this\s+(?:finding|result|observation|expression|association|effect|"
    r"receptor|gene|protein|pathway|axis|signal|signaling|signalling|"
    r"variant|mutation)|"
    r"the\s+(?:receptor|gene|protein|expression|finding|result|results|"
    r"effect|effects|association|relationship|pathway|axis|signal|"
    r"signaling|signalling|variant|mutation|marker|biomarker)"
    r")\b",
    re.IGNORECASE,
)

# Check whether a sentence refers back to the GPCR.
def begins_with_anaphoric_reference(sentence):
    return bool(
        ANAPHORIC_PATTERN.search(clean(sentence))
    )

SHORT_KEYWORD_PATTERNS = {}
LONG_KEYWORDS = {}

for category, words in CATEGORY_KEYWORDS.items():
    SHORT_KEYWORD_PATTERNS[category] = []
    LONG_KEYWORDS[category] = []

    for word in words:
        word_l = word.lower()
        compact = word_l.replace("-", "")

        if len(word_l) <= 5 and compact.isalnum():
            SHORT_KEYWORD_PATTERNS[category].append(
                (
                    word,
                    re.compile(
                        r"\b"
                        + re.escape(word_l)
                        + r"\b",
                        re.IGNORECASE,
                    ),
                )
            )
        else:
            LONG_KEYWORDS[category].append(
                (word, word_l)
            )

# Find biological and clinical categories in context.
def classify_context_text(text, document_title=""):
    text_clean = clean(text)
    text_l = text_clean.lower()

    categories = []
    keywords = []

    for category in CATEGORY_KEYWORDS:
        hits = []

        if category == IMMUNE_CATEGORY:
            hits = immune_microenvironment_hits(
                text_clean,
                document_title=document_title,
            )

            if hits:
                categories.append(category)
                keywords.extend(hits)

            continue

        if category == BIOMARKER_CATEGORY:
            hits = biomarker_diagnosis_hits(
                text_clean,
                document_title=document_title,
            )

            if hits:
                categories.append(category)
                keywords.extend(hits)

            continue

        if category == SURVIVAL_CATEGORY:
            hits = survival_hits(text_clean)

            if hits:
                categories.append(category)
                keywords.extend(hits)

            continue

        if category == IMAGING_CATEGORY:
            hits = imaging_theranostics_hits(
                text_clean,
                document_title=document_title,
            )

            if hits:
                categories.append(category)
                keywords.extend(hits)

            continue

        for word, pattern in SHORT_KEYWORD_PATTERNS[
            category
        ]:
            if pattern.search(text_clean):
                hits.append(word)

        for word, word_l in LONG_KEYWORDS[
            category
        ]:
            if word_l in text_l:
                hits.append(word)

        if hits:
            categories.append(category)
            keywords.extend(hits)

    if not categories:
        categories = [OTHER_CATEGORY]

    return (
        categories,
        sorted(set(keywords), key=str.casefold),
    )

# Find evidence of a cancer immune microenvironment.
def immune_microenvironment_hits(text, document_title=""):
    """Return immune-category triggers supported by sentence-level context."""

    sentences = [
        clean(sentence)
        for sentence in re.split(r"\s*\|\|\s*", text)
        if clean(sentence)
    ]

    accepted_hits = set()
    broad_hits_across_context = set()

    for sentence in sentences:
        sentence_l = sentence.lower()

        direct_hits = [
            phrase
            for phrase in IMMUNE_DIRECT_PHRASES
            if phrase in sentence_l
        ]

        if direct_hits:
            accepted_hits.update(direct_hits)
            continue

        lexical_screened_sentence = sentence_l

        for phrase in IMMUNE_FALSE_LEXICAL_PHRASES:
            lexical_screened_sentence = lexical_screened_sentence.replace(
                phrase,
                " ",
            )

        broad_hits = [
            term
            for term, pattern in IMMUNE_BROAD_PATTERNS.items()
            if pattern.search(lexical_screened_sentence)
        ]

        broad_hits_across_context.update(broad_hits)

        qualifier_hits = [
            qualifier
            for qualifier in IMMUNE_CONTEXT_QUALIFIERS
            if qualifier in sentence_l
        ]

        model_usage = any(
            marker in sentence_l
            for marker in IMMUNE_MODEL_MARKERS
        )

        if broad_hits and (not model_usage or qualifier_hits):
            accepted_hits.update(broad_hits)

            if model_usage:
                accepted_hits.update(qualifier_hits)

    title_l = clean(document_title).lower()
    title_direct_hits = [
        phrase
        for phrase in IMMUNE_DIRECT_PHRASES
        if phrase in title_l
    ]

    if broad_hits_across_context and title_direct_hits:
        accepted_hits.update(broad_hits_across_context)
        accepted_hits.update(title_direct_hits)

    return sorted(accepted_hits, key=str.casefold)

# Match a keyword at word boundaries.
def boundary_keyword_present(text, keyword):
    """Match a keyword as a complete token or phrase, allowing plurals."""

    if keyword in {"roc", "auc", "pet", "spect"}:
        return bool(
            re.search(
                r"\b" + re.escape(keyword) + r"\b",
                text,
                re.IGNORECASE,
            )
        )

    return keyword in text.lower()

# Find clinical biomarker or diagnostic evidence.
def biomarker_diagnosis_hits(text, document_title=""):
    """Classify biomarker contexts without treating every marker as clinical."""

    sentences = [
        clean(sentence)
        for sentence in re.split(r"\s*\|\|\s*", text)
        if clean(sentence)
    ]

    accepted_hits = set()
    title_l = clean(document_title).lower()
    title_context = any(
        qualifier in title_l
        for qualifier in BIOMARKER_MARKER_QUALIFIERS
    ) or any(
        boundary_keyword_present(title_l, term)
        for term in BIOMARKER_STRONG_TERMS
    )

    for sentence in sentences:
        sentence_l = sentence.lower()

        strong_hits = [
            term
            for term in BIOMARKER_STRONG_TERMS
            if boundary_keyword_present(sentence_l, term)
        ]

        if strong_hits:
            accepted_hits.update(strong_hits)

        screened_sentence = sentence_l

        for phrase in BIOMARKER_FALSE_MARKER_PHRASES:
            screened_sentence = screened_sentence.replace(phrase, " ")

        marker_present = bool(
            re.search(r"\bmarkers?\b", screened_sentence, re.I)
        )

        marker_context = title_context or any(
            qualifier in sentence_l
            for qualifier in BIOMARKER_MARKER_QUALIFIERS
        )

        if marker_present and marker_context:
            accepted_hits.add("marker")

    return sorted(accepted_hits, key=str.casefold)

# Distinguish patient survival from cellular survival.
def survival_hits(text):
    """Distinguish subject survival from cellular survival mechanisms."""

    sentences = [
        clean(sentence)
        for sentence in re.split(r"\s*\|\|\s*", text)
        if clean(sentence)
    ]

    accepted_hits = set()

    for sentence in sentences:
        sentence_l = sentence.lower()

        strong_hits = [
            phrase
            for phrase in SURVIVAL_STRONG_PHRASES
            if phrase in sentence_l
        ]

        if strong_hits:
            accepted_hits.update(strong_hits)
            continue

        screened_sentence = sentence_l

        for pattern in SURVIVAL_FALSE_PATTERNS:
            screened_sentence = pattern.sub(" ", screened_sentence)

        if re.search(r"\bsurvivals?\b", screened_sentence, re.I):
            accepted_hits.add("survival")

    return sorted(accepted_hits, key=str.casefold)

# Find cancer imaging or theranostic evidence.
def imaging_theranostics_hits(text, document_title=""):
    """Require a cancer/targeting context for standalone imaging language."""

    sentences = [
        clean(sentence)
        for sentence in re.split(r"\s*\|\|\s*", text)
        if clean(sentence)
    ]

    accepted_hits = set()
    title_l = clean(document_title).lower()
    title_context = any(
        qualifier in title_l
        for qualifier in IMAGING_CONTEXT_QUALIFIERS
    )

    for sentence in sentences:
        sentence_l = sentence.lower()

        strong_hits = [
            term
            for term in IMAGING_STRONG_TERMS
            if boundary_keyword_present(sentence_l, term)
        ]

        if re.search(r"\bpet\b", sentence, re.I):
            strong_hits.append("pet")

        if strong_hits:
            accepted_hits.update(strong_hits)
            continue

        imaging_present = bool(
            re.search(r"\bimaging\b", sentence_l, re.I)
        )
        imaging_context = title_context or any(
            qualifier in sentence_l
            for qualifier in IMAGING_CONTEXT_QUALIFIERS
        )

        if imaging_present and imaging_context:
            accepted_hits.add("imaging")

    return sorted(accepted_hits, key=str.casefold)

# Select the primary category from the assigned labels.
def primary_category(categories):
    for category in PRIMARY_PRIORITY:
        if category in categories:
            return category

    return OTHER_CATEGORY

# Assign categories using GPCR-centered publication text.
def classify_row(row):
    aliases = split_terms(
        row["Validated_Search_Terms"]
    )

    sentences = split_article_sentences(
        clean(row["Title"]),
        clean(row["Abstract"]),
    )

    gpcr_sentence_indices = []
    matched_aliases = set()

    for index, sentence in enumerate(sentences):
        hits = matching_gpcr_terms(
            sentence,
            aliases,
        )

        if hits:
            gpcr_sentence_indices.append(index)
            matched_aliases.update(hits)

    if not gpcr_sentence_indices:
        raise ValueError(
            f"No validated GPCR term found for "
            f"{clean(row['Official_Gene'])}/"
            f"{clean(row['TCGA_Code'])}/"
            f"{clean(row['PMID'])}"
        )

    context_indices = set(
        gpcr_sentence_indices
    )
    extension_indices = set()

    for index in gpcr_sentence_indices:
        if (
            index + 1 < len(sentences)
            and begins_with_anaphoric_reference(
                sentences[index + 1]
            )
        ):
            context_indices.add(index + 1)
            extension_indices.add(index + 1)

        if (
            index > 0
            and begins_with_anaphoric_reference(
                sentences[index]
            )
        ):
            context_indices.add(index - 1)
            extension_indices.add(index - 1)

    context_indices = sorted(context_indices)
    extension_indices = sorted(extension_indices)

    context_sentences = [
        sentences[index]
        for index in context_indices
    ]

    categories, keywords = classify_context_text(
        " || ".join(context_sentences),
        document_title=clean(row["Title"]),
    )

    out = dict(row)

    out["Classification_Version"] = (
        "GPCR_PEnDB_16_categories_final"
    )

    out["Matched_Aliases"] = ";".join(
        sorted(matched_aliases, key=str.casefold)
    )
    out["Matched_Alias_Count"] = len(
        matched_aliases
    )

    out["GPCR_Sentence_Indices"] = ";".join(
        str(index)
        for index in gpcr_sentence_indices
    )

    out["Context_Extension_Indices"] = ";".join(
        str(index)
        for index in extension_indices
    )

    out["Context_Extension_Used"] = (
        "YES"
        if extension_indices
        else "NO"
    )

    out["Supporting_Sentence_Count"] = len(
        context_sentences
    )

    out["Supporting_Sentences"] = " || ".join(
        context_sentences
    )

    out["Association_Categories"] = "; ".join(
        categories
    )

    out["Association_Category_Count"] = len(
        categories
    )

    out["Primary_Association_Category"] = (
        primary_category(categories)
    )

    out["Matched_Keywords"] = "; ".join(
        keywords
    )

    out["Association_Status"] = "Retain"
    out["Review_Reason"] = ""

    return out

# Check that required columns are present.
def require_columns(df):
    required = {
        "UniProt_ID",
        "Gene",
        "Official_Gene",
        "Original_Gene",
        "Protein_Name",
        "TCGA_Code",
        "TCGA_Study_Name",
        "Cancer",
        "PMID",
        "Title",
        "Abstract",
        "Journal",
        "Year",
        "Authors",
        "Search_ID",
        "Validated_Search_Terms",
        "Final_Literature_Status",
        "Final_Validation_Source",
        "PubMed_Query",
    }

    missing = sorted(
        required - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Input is missing required column(s): "
            + ", ".join(missing)
        )

# Validate input for the current record.
def validate_input(df):
    require_columns(df)

    if df.empty:
        raise ValueError("Input is empty.")

    if not df[
        "Final_Literature_Status"
    ].map(clean).eq("KEEP").all():
        raise ValueError(
            "Input contains a row not marked "
            "Final_Literature_Status=KEEP."
        )

    if df[
        "Validated_Search_Terms"
    ].map(clean).eq("").any():
        raise ValueError(
            "Input contains blank Validated_Search_Terms."
        )

    if df.duplicated(
        ["UniProt_ID", "TCGA_Code", "PMID"]
    ).any():
        raise ValueError(
            "Input contains duplicate GPCR-cancer-PMID rows."
        )

# Build summary from the input records.
def build_summary(classified):
    rows = []

    for category in ALL_CATEGORIES:
        count = int(
            classified[
                "Association_Categories"
            ]
            .fillna("")
            .apply(
                lambda value: category in [
                    item.strip()
                    for item in value.split(";")
                    if item.strip()
                ]
            )
            .sum()
        )

        rows.append({
            "Summary_Type":
                "Association_Category_All_Matches",
            "Value": category,
            "Count": count,
        })

    for category in PRIMARY_PRIORITY:
        rows.append({
            "Summary_Type":
                "Primary_Association_Category",
            "Value": category,
            "Count": int(
                classified[
                    "Primary_Association_Category"
                ].eq(category).sum()
            ),
        })

    rows.extend([
        {
            "Summary_Type":
                "Context_Extension",
            "Value": "Used",
            "Count": int(
                classified[
                    "Context_Extension_Used"
                ].eq("YES").sum()
            ),
        },
        {
            "Summary_Type":
                "Context_Extension",
            "Value": "Not used",
            "Count": int(
                classified[
                    "Context_Extension_Used"
                ].eq("NO").sum()
            ),
        },
        {
            "Summary_Type":
                "Association_Status",
            "Value": "Retain",
            "Count": len(classified),
        },
    ])

    return pd.DataFrame(rows)

# Write report for the current record.
def write_report(input_df, classified):
    observed_categories = {
        item
        for value in classified[
            "Association_Categories"
        ]
        for item in split_terms(value)
    }

    lpar1_10632375 = classified[
        classified["Official_Gene"].eq("LPAR1")
        & classified["Cancer"].eq("Ovarian cancer")
        & classified["PMID"].astype(str).eq("10632375")
    ]

    lpar1_30353652 = classified[
        classified["Official_Gene"].eq("LPAR1")
        & classified["Cancer"].eq("Ovarian cancer")
        & classified["PMID"].astype(str).eq("30353652")
    ]

    calcr_8143880 = classified[
        classified["Official_Gene"].eq("CALCR")
        & classified["TCGA_Code"].eq("THCA")
        & classified["PMID"].astype(str).eq("8143880")
    ]

    ackr3_39018005 = classified[
        classified["Official_Gene"].eq("ACKR3")
        & classified["TCGA_Code"].eq("COAD")
        & classified["PMID"].astype(str).eq("39018005")
    ]

    gpr68_32948783 = classified[
        classified["Official_Gene"].eq("GPR68")
        & classified["TCGA_Code"].eq("OV")
        & classified["PMID"].astype(str).eq("32948783")
    ]

    cxcr1_37541772 = classified[
        classified["Official_Gene"].eq("CXCR1")
        & classified["TCGA_Code"].eq("PAAD")
        & classified["PMID"].astype(str).eq("37541772")
    ]

    cxcr4_20563261 = classified[
        classified["Official_Gene"].eq("CXCR4")
        & classified["TCGA_Code"].eq("PRAD")
        & classified["PMID"].astype(str).eq("20563261")
    ]

    lpar1_24657741 = classified[
        classified["Official_Gene"].eq("LPAR1")
        & classified["TCGA_Code"].eq("OV")
        & classified["PMID"].astype(str).eq("24657741")
    ]

    tshr_25786552 = classified[
        classified["Official_Gene"].eq("TSHR")
        & classified["TCGA_Code"].eq("THCA")
        & classified["PMID"].astype(str).eq("25786552")
    ]

    lgr5_26003047 = classified[
        classified["Official_Gene"].eq("LGR5")
        & classified["TCGA_Code"].eq("COAD")
        & classified["PMID"].astype(str).eq("26003047")
    ]

    fzd10_19137009 = classified[
        classified["Official_Gene"].eq("FZD10")
        & classified["TCGA_Code"].eq("SARC")
        & classified["PMID"].astype(str).eq("19137009")
    ]

    lpar1_29029500 = classified[
        classified["Official_Gene"].eq("LPAR1")
        & classified["TCGA_Code"].eq("LIHC")
        & classified["PMID"].astype(str).eq("29029500")
    ]

    chrm3_40541171 = classified[
        classified["Official_Gene"].eq("CHRM3")
        & classified["TCGA_Code"].eq("GBM")
        & classified["PMID"].astype(str).eq("40541171")
    ]

    sstr2_29747035 = classified[
        classified["Official_Gene"].eq("SSTR2")
        & classified["TCGA_Code"].eq("PRAD")
        & classified["PMID"].astype(str).eq("29747035")
    ]

    # Check whether a record has the specified category.
    def row_has_category(frame, category):
        return (
            len(frame) == 1
            and category in split_terms(
                frame.iloc[0]["Association_Categories"]
            )
        )

    checks = {
        "Input/output row count preserved": (
            len(input_df) == len(classified)
        ),
        "All validated GPCRs preserved": (
            input_df["Official_Gene"].nunique()
            == classified[
                "Official_Gene"
            ].nunique()
        ),
        "All validated GPCR-cancer pairs preserved": (
            input_df[
                ["UniProt_ID", "TCGA_Code"]
            ].drop_duplicates().shape[0]
            == classified[
                ["UniProt_ID", "TCGA_Code"]
            ].drop_duplicates().shape[0]
        ),
        "No duplicate output association keys": (
            not classified.duplicated(
                [
                    "UniProt_ID",
                    "TCGA_Code",
                    "PMID",
                ]
            ).any()
        ),
        "Every row retained": (
            classified[
                "Association_Status"
            ].eq("Retain").all()
        ),
        "Every row has GPCR-supporting context": (
            classified[
                "Supporting_Sentence_Count"
            ].astype(int).ge(1).all()
        ),
        "Every row has >=1 category": (
            classified[
                "Association_Category_Count"
            ].astype(int).ge(1).all()
        ),
        "Only frozen 16 categories used": (
            observed_categories.issubset(
                set(ALL_CATEGORIES)
            )
        ),
        "Regression 10632375 excludes immune category": (
            len(lpar1_10632375) == 1
            and not row_has_category(
                lpar1_10632375,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 30353652 retains immune category": (
            row_has_category(
                lpar1_30353652,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 8143880 excludes TT-cell false match": (
            len(calcr_8143880) == 1
            and not row_has_category(
                calcr_8143880,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 39018005 excludes MIF-name false match": (
            len(ackr3_39018005) == 1
            and not row_has_category(
                ackr3_39018005,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 32948783 excludes T-cell gene-name match": (
            len(gpr68_32948783) == 1
            and not row_has_category(
                gpr68_32948783,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 37541772 retains immune recruitment": (
            row_has_category(
                cxcr1_37541772,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 20563261 retains stromal-cell context": (
            row_has_category(
                cxcr4_20563261,
                IMMUNE_CATEGORY,
            )
        ),
        "Regression 24657741 excludes endosomal-marker match": (
            len(lpar1_24657741) == 1
            and not row_has_category(
                lpar1_24657741,
                BIOMARKER_CATEGORY,
            )
        ),
        "Regression 25786552 retains clinical marker": (
            row_has_category(
                tshr_25786552,
                BIOMARKER_CATEGORY,
            )
        ),
        "Regression 26003047 excludes stem-cell marker": (
            len(lgr5_26003047) == 1
            and not row_has_category(
                lgr5_26003047,
                BIOMARKER_CATEGORY,
            )
        ),
        "Regression 19137009 excludes cell-survival match": (
            len(fzd10_19137009) == 1
            and not row_has_category(
                fzd10_19137009,
                SURVIVAL_CATEGORY,
            )
        ),
        "Regression 29029500 retains patient survival": (
            row_has_category(
                lpar1_29029500,
                SURVIVAL_CATEGORY,
            )
        ),
        "Regression 40541171 excludes calcium imaging": (
            len(chrm3_40541171) == 1
            and not row_has_category(
                chrm3_40541171,
                IMAGING_CATEGORY,
            )
        ),
        "Regression 29747035 retains PET imaging": (
            row_has_category(
                sstr2_29747035,
                IMAGING_CATEGORY,
            )
        ),
    }

    failed = [
        label
        for label, passed in checks.items()
        if not passed
    ]

    other = classified[
        classified[
            "Association_Categories"
        ].eq(OTHER_CATEGORY)
    ]

    lines = [
        "GPCR Cancer Pipeline - Code S13 Classification Integrity Report",
        "=" * 78,
        "",
        "INPUT / PRESERVATION",
        "-" * 78,
        f"Corrected input rows              : {len(input_df)}",
        f"Corrected input GPCRs             : {input_df['Official_Gene'].nunique()}",
        f"Corrected input GPCR-cancer pairs : "
        f"{input_df[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]}",
        f"Corrected input unique PMIDs      : {input_df['PMID'].nunique()}",
        "",
        "CATEGORY CLASSIFICATION",
        "-" * 78,
        f"Classified/retained rows            : {len(classified)}",
        f"Classified GPCRs                    : {classified['Official_Gene'].nunique()}",
        f"Classified GPCR-cancer pairs        : "
        f"{classified[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]}",
        f"Unique PMIDs                        : {classified['PMID'].nunique()}",
        f"Predefined categories               : {len(CATEGORY_KEYWORDS)}",
        f"Fallback category                   : {OTHER_CATEGORY}",
        f"Total categories                    : {len(ALL_CATEGORIES)}",
        f"Other/Unspecified paper rows        : {len(other)}",
        f"GPCRs with >=1 Other paper          : {other['Official_Gene'].nunique()}",
        f"Rows using controlled adjacent context: "
        f"{int(classified['Context_Extension_Used'].eq('YES').sum())}",
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
            "Code S13 is an organizational classification layer and does not "
            "revalidate the corrected literature input from Code S12. Categories are assigned from "
            "sentences containing the validated GPCR-supporting terms. "
            "One adjacent sentence may be included only through a direct "
            "anaphoric/reference relationship. For Immune Microenvironment, "
            "specific immune-microenvironment phrases are accepted directly. "
            "Broad immune and stromal cell terms use boundary-aware matching; "
            "misleading lexical phrases are masked, and named experimental "
            "cell lines require an additional immune-context qualifier. "
            "Generic marker terms require biomarker/clinical context, cellular "
            "survival phrases are excluded from the Survival category, and "
            "standalone imaging language requires a cancer, targeting, or "
            "diagnostic imaging context. "
            "Records with "
            "no match to the "
            "15 predefined biological/clinical categories are retained as "
            "Other/Unspecified rather than being discarded."
        ),
        "",
        f"Integrity failures                  : {len(failed)}",
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
            "ALL CODE S13 INTEGRITY CHECKS PASSED",
        ]

    Path(OUT_REPORT).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    return failed

# Run this stage and write its output files.
def main():
    print("=" * 78)
    print("GPCR Cancer Pipeline - Code S13")
    print("GPCR-Centered 16-Category Classification")
    print("=" * 78)

    if not os.path.exists(INPUT_FILE):
        raise SystemExit(
            f"ERROR: missing {INPUT_FILE}"
        )

    df = pd.read_csv(
        INPUT_FILE,
        dtype={"PMID": str},
        low_memory=False,
    ).fillna("")

    try:
        validate_input(df)
    except Exception as exc:
        print(
            f"INPUT VALIDATION FAILED:\n{exc}"
        )
        sys.exit(1)

    print(f"Input rows       : {len(df):,}")
    print(
        f"Input GPCRs      : "
        f"{df['Official_Gene'].nunique():,}"
    )
    print(
        f"Input pairs      : "
        f"{df[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]:,}"
    )
    print("Network requests : 0")
    print("=" * 78)

    rows = []

    for number, (_, row) in enumerate(
        df.iterrows(),
        start=1,
    ):
        if (
            number == 1
            or number % 2500 == 0
            or number == len(df)
        ):
            print(
                f"Classifying "
                f"{number:,}/{len(df):,}"
            )

        try:
            rows.append(
                classify_row(row)
            )
        except Exception as exc:
            raise RuntimeError(
                f"Classification failed at row "
                f"{number}: {exc}"
            ) from exc

    classified = pd.DataFrame(rows)

    summary = build_summary(
        classified
    )

    failed = write_report(
        df,
        classified,
    )

    classified.to_csv(
        OUT_ALL,
        index=False,
        na_rep="",
    )

    classified.to_csv(
        OUT_RETAINED,
        index=False,
        na_rep="",
    )

    summary.to_csv(
        OUT_SUMMARY,
        index=False,
        na_rep="",
    )

    print()
    print("=" * 78)
    print(
        "CODE S13 COMPLETE"
        if not failed
        else "CODE S13 COMPLETE WITH INTEGRITY ISSUES"
    )
    print("=" * 78)
    print(
        f"Classified rows : "
        f"{len(classified):,}"
    )
    print(
        f"GPCRs           : "
        f"{classified['Official_Gene'].nunique():,}"
    )
    print(
        f"Pairs           : "
        f"{classified[['UniProt_ID','TCGA_Code']].drop_duplicates().shape[0]:,}"
    )
    print(
        f"Unique PMIDs    : "
        f"{classified['PMID'].nunique():,}"
    )
    print(
        f"Other rows      : "
        f"{int(classified['Association_Categories'].eq(OTHER_CATEGORY).sum()):,}"
    )
    print(
        f"Integrity errors: "
        f"{len(failed)}"
    )

    if failed:
        sys.exit(2)

if __name__ == "__main__":
    main()
