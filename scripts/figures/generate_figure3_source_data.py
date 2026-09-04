#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path

INPUT = Path("results/H4/H4_complete_62_variant_evidence.tsv")
OUTDIR = Path("results/figures/source_data")

PANEL_A = OUTDIR / "Figure3_panelA_sift_categories.tsv"
PANEL_B = OUTDIR / "Figure3_panelB_gene_deleterious.tsv"

OUTDIR.mkdir(parents=True, exist_ok=True)

with INPUT.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

assert len(rows) == 62, f"Expected 62 candidates, found {len(rows)}"

# ------------------------------------------------------------------
# Panel A: overall SIFT classification
# ------------------------------------------------------------------

prediction_counts = Counter(r["SIFT_prediction"] for r in rows)
status_counts = Counter(r["SIFT_status"] for r in rows)

expected_predictions = {
    "deleterious": 6,
    "deleterious_low_confidence": 6,
    "tolerated": 38,
    "tolerated_low_confidence": 6,
    "": 6,
}

for key, expected in expected_predictions.items():
    observed = prediction_counts[key]
    assert observed == expected, (
        f"SIFT prediction count mismatch for {key!r}: "
        f"expected {expected}, observed {observed}"
    )

assert status_counts["SIFT_AVAILABLE"] == 56
assert status_counts["SIFT_UNAVAILABLE"] == 4
assert status_counts["NOT_APPLICABLE_HIGH_IMPACT"] == 2

panel_a_rows = [
    ("deleterious", 6),
    ("deleterious_low_confidence", 6),
    ("tolerated", 38),
    ("tolerated_low_confidence", 6),
]

with PANEL_A.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["SIFT_category", "Count"])
    w.writerows(panel_a_rows)

# ------------------------------------------------------------------
# Panel B: gene distribution of deleterious calls
# ------------------------------------------------------------------

gene_counts = defaultdict(lambda: {
    "deleterious": 0,
    "deleterious_low_confidence": 0,
})

for r in rows:
    pred = r["SIFT_prediction"]
    if pred in gene_counts[r["SYMBOL"]]:
        gene_counts[r["SYMBOL"]][pred] += 1

panel_b_rows = []
for gene, counts in gene_counts.items():
    total = counts["deleterious"] + counts["deleterious_low_confidence"]
    if total > 0:
        panel_b_rows.append((
            gene,
            counts["deleterious"],
            counts["deleterious_low_confidence"],
            total,
        ))

panel_b_rows.sort(key=lambda x: (-x[3], x[0]))

expected_gene_totals = {
    "ST3GAL4": 5,
    "NEU3": 2,
    "NEU1": 2,
    "ST8SIA2": 1,
    "CMAH": 1,
    "GNE": 1,
}

observed_gene_totals = {r[0]: r[3] for r in panel_b_rows}
assert observed_gene_totals == expected_gene_totals, (
    f"Gene-level deleterious distribution mismatch:\n"
    f"Expected: {expected_gene_totals}\n"
    f"Observed: {observed_gene_totals}"
)

with PANEL_B.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([
        "Gene",
        "Deleterious",
        "Deleterious_low_confidence",
        "Total_flagged",
    ])
    w.writerows(panel_b_rows)

print(f"Saved: {PANEL_A}")
print(f"Saved: {PANEL_B}")
print()
print("FIGURE 3 SOURCE DATA GENERATION: PASS")
