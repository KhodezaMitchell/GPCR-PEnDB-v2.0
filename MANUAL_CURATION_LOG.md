# GPCR-PEnDB 2.0 Manual Curation Log

Documentation version: 1.0  
Documentation date: 26 October 2026

This log documents manual corrections applied to the cancer-network release after the saved computational workflow outputs were generated. It explains the difference between the preserved supplementary CSV files and the final website release. The documentation date identifies this log version, rather than the date of each individual curation action.

## Release counts and preserved outputs

The supplementary CSV files preserve the saved workflow outputs and have not been overwritten to incorporate the subsequent manual corrections.

| Metric | Saved workflow outputs | Final website release |
|---|---:|---:|
| GPCR–cancer–PMID records | 18,892 | 18,890 |
| Distinct PMIDs | 15,010 | 15,009 |
| GPCRs with retained literature | 458 | 458 |
| Cancer types | 33 | 33 |
| GPCR–cancer associations | 3,680 | 3,680 |
| Classification categories | 16 | 16 |

The automated ambiguity audit excluded 3,253 of the 22,145 records retained after initial validation and contextual review, producing 18,892 records. The two additional manual exclusions below produced 18,890 records. Thus, the ambiguity audit and final quality-control checks together excluded 3,255 records:

22,145 − 3,253 − 2 = 18,890.

A record is identified by UniProt accession, TCGA cancer code, and PMID. Removing one assignment does not necessarily remove that publication from every other association.

## Manual corrections

All three corrections concern LPAR1, UniProt accession Q92633.

| ID | PMID | Cancer code | Action | Reason and effect |
|---|---|---|---|---|
| MC-01 | 38396744 | COAD | Remove the LPAR1–colon assignment; retain the LPAR1–PRAD assignment. | The study concerns PC-3 human prostate cancer cells. An incidental MACC1 reference led to an incorrect colon assignment. This removes one record, without reducing the number of distinct PMIDs in the resource. |
| MC-02 | 31257480 | PRAD | Remove the LPAR1–prostate assignment. | The publication concerns acupuncture, moxibustion, endometrial receptivity, and infertility treatment. Its circ-LPAR1 mention does not establish a prostate-cancer connection. This removes one record and one distinct PMID from the resource. |
| MC-03 | 39267662 | PAAD | Repair the displayed abstract and supporting evidence text. | Restore text damaged during formatting/HTML handling, including the malformed “PENPP2” fragment and missing text around statistical comparisons. This is a text-restoration correction; the record remains included and no record or PMID count changes. |

Publication identifiers:

- MC-01: [Global Proteomics Analysis of Lysophosphatidic Acid Signaling in PC-3 Human Prostate Cancer Cells: Role of CCN1](https://pubmed.ncbi.nlm.nih.gov/38396744/).
- MC-02: [Circular RNA sequencing reveals the molecular mechanism of the effects of acupuncture and moxibustion on endometrial receptivity in patients undergoing infertility treatment](https://pubmed.ncbi.nlm.nih.gov/31257480/).
- MC-03: [Insights into autotaxin- and lysophosphatidate-mediated signaling in the pancreatic ductal adenocarcinoma tumor microenvironment: a survey of pathway gene expression](https://pubmed.ncbi.nlm.nih.gov/39267662/).

For MC-03, the restored supporting evidence retains the comparisons P≤0.003, P<0.02, and P<0.01 in their corresponding sentences. The complete corrected abstract and evidence text are stored in `cancer_networks/data/pairs/Q92633__PAAD.json` in the accompanying network archive.

## Website category alignment

The website category assignments were synchronized with the saved submitted classifications in `Data_S12_CORRECTED_FINAL16_DB_Papers_GPCR_Cancer.csv` for all retained records. The original CSV was preserved. This alignment copied existing assignments rather than performing a new literature search or a new classification run.

The accompanying archive's `CATEGORY_ALIGNMENT_REPORT.json` records 1,667 retained records whose category sets changed relative to the earlier website build. Primary-category and matched-keyword fields were synchronized as well. Final category counts include Expression 11,819; Signaling 8,465; Therapeutic Target 7,281; and Other/Unspecified 932. The alignment does not change the final record or distinct-PMID totals.

## Using the saved tables and corrected release

To reconcile the saved publication-level outputs with the final website release:

1. Exclude only the keys `(Q92633, COAD, 38396744)` and `(Q92633, PRAD, 31257480)`.
2. Preserve `(Q92633, PRAD, 38396744)`.
3. Use the restored abstract and evidence text for `(Q92633, PAAD, 39267662)` from the accompanying corrected JSON file.
4. Use the saved submitted category assignments for the remaining records and recompute the relevant website summaries.

After these corrections, LPAR1 has 115 records representing 113 distinct PMIDs across 15 cancer types. Its ovarian-cancer association retains 17 PMIDs.

These changes do not alter the resource's scope: an association is a retained textual connection between a GPCR and a cancer type, which may include a background mention. It does not by itself establish a direct experimental, functional, or causal relationship. Review reduction relies on PubMed publication-type indexing, subject to the limitations described in the manuscript.
