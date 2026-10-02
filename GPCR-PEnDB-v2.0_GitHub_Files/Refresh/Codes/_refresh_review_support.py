"""Validate fresh-run review sheets without changing automatic scientific rules."""

from pathlib import Path
import re

import pandas as pd


PAPER_KEYS = ("UniProt_ID", "TCGA_Code", "PMID")
REVIEW_EXIT = 3


def clean(value):
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def terms(value):
    result, seen = [], set()
    for part in clean(value).split(";"):
        term = clean(part)
        if term and term.casefold() not in seen:
            result.append(term)
            seen.add(term.casefold())
    return result


def key(row, columns=PAPER_KEYS):
    return tuple(clean(row[column]) for column in columns)


def read_csv(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def check_keys(frame, columns):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError("Missing identifier columns: " + ", ".join(sorted(missing)))
    keys = [key(row, columns) for _, row in frame.iterrows()]
    if any(not all(item) for item in keys) or len(keys) != len(set(keys)):
        raise ValueError("Review identifiers must be complete and unique")
    return set(keys)


def ensure_template(queue, path, columns=PAPER_KEYS, support=True):
    """Create a review sheet once; preserve an existing curator's entries."""
    path = Path(path)
    expected = check_keys(queue, columns)
    if path.exists():
        if check_keys(read_csv(path), columns) != expected:
            raise ValueError(f"Existing review sheet has changed identifiers: {path}")
        return
    identifiers = list(columns)
    if "Official_Gene" in queue.columns and "Official_Gene" not in identifiers:
        identifiers.append("Official_Gene")
    sheet = queue[identifiers].copy()
    sheet["Final_Decision"] = ""
    if support:
        sheet["Supporting_GPCR_Terms"] = ""
    sheet["Curator_Comment"] = ""
    sheet.to_csv(path, index=False)


def validate_sheet(queue, path, allowed, columns=PAPER_KEYS, support=True):
    """Require exactly one documented decision per pending row."""
    frame = read_csv(path)
    fields = {"Final_Decision", "Curator_Comment"}
    if support:
        fields.add("Supporting_GPCR_Terms")
    if not fields.issubset(frame.columns):
        raise ValueError("Decision file lacks: " + ", ".join(sorted(fields - set(frame.columns))))
    if check_keys(queue, columns) != check_keys(frame, columns):
        raise ValueError("Decision rows must match the pending review queue exactly")
    decisions = {}
    for _, row in frame.iterrows():
        decision = clean(row["Final_Decision"]).upper()
        comment = clean(row["Curator_Comment"])
        if decision not in allowed or not comment:
            raise ValueError("Each row needs an allowed decision and a curator comment")
        supporting = terms(row.get("Supporting_GPCR_Terms", ""))
        if support and decision == "KEEP" and not supporting:
            raise ValueError("Each KEEP decision needs Supporting_GPCR_Terms")
        decisions[key(row, columns)] = {
            "decision": decision, "comment": comment, "supporting": supporting,
        }
    return decisions


def flexible_pattern(term):
    """Use the established punctuation and letter/number boundary treatment."""
    pieces, index = [], 0
    while index < len(term):
        ch = term[index]
        if re.match(r"[\s\-_/\.]", ch):
            while index + 1 < len(term) and re.match(r"[\s\-_/\.]", term[index + 1]):
                index += 1
            pieces.append(r"[\s\-_/\.]*")
        else:
            pieces.append(re.escape(ch))
            if index + 1 < len(term):
                nxt = term[index + 1]
                if ch.isalnum() and nxt.isalnum() and (
                    (ch.isalpha() and nxt.isdigit()) or (ch.isdigit() and nxt.isalpha())
                ):
                    pieces.append(r"[\s\-_/\.]*")
        index += 1
    return re.compile(r"(?<![A-Za-z0-9])" + "".join(pieces) + r"(?![A-Za-z0-9])", re.I)


def validate_supporting_terms(row, supporting, alias_file="GPCR_Alias_Table_FINAL_VALIDATED.csv"):
    """Check that curator-named support is current receptor terminology in the article."""
    aliases = read_csv(alias_file)
    matched = aliases[aliases["UniProt_ID"].map(clean).eq(clean(row["UniProt_ID"]))]
    if len(matched) != 1:
        raise ValueError("Reviewed receptor does not have one current terminology record")
    receptor = matched.iloc[0]
    allowed = set()
    for field in ("Official_Gene", "Protein_Name", "Search_Terms", "Context_Required_Terms"):
        allowed.update(term.casefold() for term in terms(receptor.get(field, "")))
    text = clean(row.get("Title", "")) + " " + clean(row.get("Abstract", ""))
    for term in supporting:
        if term.casefold() not in allowed:
            raise ValueError(f"Supporting term is outside this receptor's current terminology: {term}")
        if not flexible_pattern(term).search(text):
            raise ValueError(f"Supporting term is absent from the supplied title/abstract: {term}")


def apply_ambiguity_reviews(resolved, decisions_path, template_path):
    """Resolve only the rows held for manual review by the existing S12 rules."""
    pending = resolved[resolved["Ambiguity_Decision"].eq("MANUAL_REVIEW")].copy()
    if pending.empty:
        if decisions_path:
            validate_sheet(pending, decisions_path, {"KEEP", "REJECT"})
        return resolved
    ensure_template(pending, template_path)
    if not decisions_path:
        return resolved
    decisions = validate_sheet(pending, decisions_path, {"KEEP", "REJECT"})
    reviewed = []
    result = resolved.copy()
    for index, row in pending.iterrows():
        item = decisions[key(row)]
        if item["decision"] == "KEEP":
            validate_supporting_terms(row, item["supporting"])
        result.loc[index, "Automatic_Ambiguity_Decision"] = row["Ambiguity_Decision"]
        result.loc[index, "Automatic_Ambiguity_Reason"] = row["Ambiguity_Reason"]
        result.loc[index, "Ambiguity_Decision"] = item["decision"]
        result.loc[index, "Ambiguity_Reason"] = "CURATOR_" + item["decision"]
        result.loc[index, "Ambiguity_Comment"] = item["comment"]
        result.loc[index, "Reviewed_Supporting_GPCR_Terms"] = ";".join(item["supporting"])
        if item["decision"] == "KEEP":
            result.loc[index, "Validated_Search_Terms"] = ";".join(item["supporting"])
        reviewed.append(result.loc[index].to_dict())
    pd.DataFrame(reviewed).to_csv("FINAL13_Ambiguity_Manual_Resolution.csv", index=False, na_rep="")
    return result
