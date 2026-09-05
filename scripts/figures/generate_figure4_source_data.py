#!/usr/bin/env python3

import csv
import re
from pathlib import Path
from collections import Counter

H5 = Path("results/H5/H5_62_variant_conservation_evidence.tsv")
H4 = Path("results/H4/H4_complete_62_variant_evidence.tsv")
INTEGRATED = Path(
    "results/H_integration/H2_H6_H3_62_candidate_complete_evidence_matrix.tsv"
)

OUT = Path(
    "results/figures/source_data/Figure4_candidate_conservation.tsv"
)

OUT.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Read authoritative H5 conservation evidence
# ---------------------------------------------------------------------
with H5.open(newline="") as fh:
    h5_rows = list(csv.DictReader(fh, delimiter="\t"))

# ---------------------------------------------------------------------
# Read H4 SIFT evidence
# ---------------------------------------------------------------------
with H4.open(newline="") as fh:
    h4_rows = list(csv.DictReader(fh, delimiter="\t"))

h4_index = {
    (r["SYMBOL"], r["CHROM"], r["POS"], r["REF"], r["ALT"]): r
    for r in h4_rows
}

# ---------------------------------------------------------------------
# Read integrated evidence matrix
# Used only to recover publication-facing protein changes.
# ---------------------------------------------------------------------
with INTEGRATED.open(newline="") as fh:
    integrated_rows = list(csv.DictReader(fh, delimiter="\t"))

integrated_index = {
    (r["SYMBOL"], r["CHROM"], r["POS"], r["REF"], r["ALT"]): r
    for r in integrated_rows
}

# ---------------------------------------------------------------------
# Build compact protein labels
#
# Examples:
#   HGVSp = I/L      + position 73  -> I73L
#   HGVSp = R/Q      + position 129 -> R129Q
#   HGVSp = E/*      + position 128 -> E128*
#   HGVSp = */L      + position 396 -> *396L
#
# RefSeq-resolved records may already contain compact forms such as I388M.
# ---------------------------------------------------------------------
def compact_change(h5_row, integrated_row):
    ref = h5_row["reference_aa"]
    pos = h5_row["protein_position"]
    hgvsp = (integrated_row.get("HGVSp") or "").strip()

    if not pos:
        raise ValueError(
            f"Missing protein position for "
            f'{h5_row["SYMBOL"]} chr{h5_row["CHROM"]}:{h5_row["POS"]}'
        )

    # Slash-style form: I/L, R/Q, E/*, */L
    if "/" in hgvsp:
        left, right = hgvsp.split("/", 1)
        left = left.strip()
        right = right.strip()

        if left and right:
            return f"{left}{pos}{right}"

    # Already compact form, e.g. I388M, P424L, N386S
    compact_match = re.fullmatch(
        r"([A-Z*])(\d+)([A-Z*])",
        hgvsp
    )
    if compact_match:
        return hgvsp

    # Forms such as p.I388M
    compact_match = re.fullmatch(
        r"p\.([A-Z*])(\d+)([A-Z*])",
        hgvsp
    )
    if compact_match:
        return "".join(compact_match.groups())

    raise ValueError(
        f"Cannot derive compact protein change for "
        f'{h5_row["SYMBOL"]} chr{h5_row["CHROM"]}:{h5_row["POS"]}; '
        f"HGVSp={hgvsp!r}"
    )


# ---------------------------------------------------------------------
# Merge H5 + H4 + publication protein labels
# ---------------------------------------------------------------------
out_rows = []

for r in h5_rows:
    key = (
        r["SYMBOL"],
        r["CHROM"],
        r["POS"],
        r["REF"],
        r["ALT"],
    )

    h4 = h4_index.get(key)
    integrated = integrated_index.get(key)

    if h4 is None:
        raise RuntimeError(f"No matching H4 record for {key}")

    if integrated is None:
        raise RuntimeError(f"No matching integrated record for {key}")

    protein_change = compact_change(r, integrated)

    sift_prediction = (h4["SIFT_prediction"] or "").strip()
    sift_status = (h4["SIFT_status"] or "").strip()

    if sift_prediction:
        sift_category = sift_prediction
    elif sift_status == "SIFT_UNAVAILABLE":
        sift_category = "SIFT unavailable"
    elif sift_status == "NOT_APPLICABLE_HIGH_IMPACT":
        sift_category = "Not applicable"
    else:
        sift_category = sift_status

    out_rows.append(
        {
            "Gene": r["SYMBOL"],
            "Protein_change": protein_change,
            "Chromosome": r["CHROM"],
            "Position": r["POS"],
            "REF": r["REF"],
            "ALT": r["ALT"],
            "Consequence": r["Consequence"],
            "Conservation_status": r["rate4site_status"],
            "Rate4Site_score": r["rate4site_score"],
            "Comparative_non_gap": r["comparative_non_gap"],
            "Comparative_total": r["comparative_total"],
            "SIFT_category": sift_category,
        }
    )

# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------
assert len(out_rows) == 62, (
    f"Expected 62 candidates, found {len(out_rows)}"
)

status_counts = Counter(
    r["Conservation_status"] for r in out_rows
)

expected_status = {
    "ASSESSABLE": 57,
    "NOT_ASSESSABLE": 4,
    "STOP_LOST_SPECIAL_CASE": 1,
}

assert status_counts == expected_status, (
    f"Unexpected conservation status counts: {dict(status_counts)}"
)

assessable = [
    r for r in out_rows
    if r["Conservation_status"] == "ASSESSABLE"
]

assert len(assessable) == 57
assert all(r["Rate4Site_score"] != "" for r in assessable)

# ---------------------------------------------------------------------
# Validate independently confirmed strongest conservation signals
# ---------------------------------------------------------------------
score_lookup = {
    (r["Gene"], r["Protein_change"]): float(r["Rate4Site_score"])
    for r in assessable
}

expected_scores = {
    ("ST8SIA4", "I73L"): -0.9582,
    ("ST8SIA2", "R129Q"): -0.9269,
    ("ST8SIA4", "R203Q"): -0.8719,
    ("CMAS", "E128*"): -0.8396,
    ("NEU3", "E94K"): -0.6812,
    ("CMAH", "I388M"): -0.6392,
    ("CMAH", "P424L"): -0.5806,
    ("NANS", "P302A"): -0.5498,
    ("ST8SIA2", "S63L"): -0.5438,
    ("NEU3", "F159C"): -0.5026,
}

for key, expected in expected_scores.items():
    observed = score_lookup.get(key)

    if observed is None:
        raise AssertionError(
            f"Missing expected candidate {key}"
        )

    assert abs(observed - expected) < 1e-6, (
        f"{key}: expected {expected}, observed {observed}"
    )

# ---------------------------------------------------------------------
# Validate non-assessable / special-case identities
# ---------------------------------------------------------------------
status_lookup = {
    (r["Gene"], r["Protein_change"]): r["Conservation_status"]
    for r in out_rows
}

expected_nonassessable = {
    ("NEU1", "L275P"),
    ("ST3GAL4", "P4L"),
    ("ST3GAL4", "R16Q"),
    ("ST3GAL4", "C38R"),
}

observed_nonassessable = {
    key
    for key, status in status_lookup.items()
    if status == "NOT_ASSESSABLE"
}

assert observed_nonassessable == expected_nonassessable, (
    f"Unexpected non-assessable set: {observed_nonassessable}"
)

assert status_lookup.get(
    ("ST3GAL4", "*396L")
) == "STOP_LOST_SPECIAL_CASE"

# ---------------------------------------------------------------------
# Write frozen Figure 4 source data
# ---------------------------------------------------------------------
fields = [
    "Gene",
    "Protein_change",
    "Chromosome",
    "Position",
    "REF",
    "ALT",
    "Consequence",
    "Conservation_status",
    "Rate4Site_score",
    "Comparative_non_gap",
    "Comparative_total",
    "SIFT_category",
]

with OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(out_rows)

print(f"Saved: {OUT}")
print()

print("CONSERVATION STATUS:")
for status in [
    "ASSESSABLE",
    "NOT_ASSESSABLE",
    "STOP_LOST_SPECIAL_CASE",
]:
    print(f"{status}: {status_counts[status]}")

print()
print("FIGURE 4 SOURCE DATA GENERATION: PASS")
