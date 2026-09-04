#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SRC = Path("results/figures/source_data")
OUT = Path("results/figures")
OUT.mkdir(parents=True, exist_ok=True)

panel_a_file = SRC / "Figure3_panelA_sift_categories.tsv"
panel_b_file = SRC / "Figure3_panelB_gene_deleterious.tsv"

png_out = OUT / "Figure3_predicted_functional_impact.png"
svg_out = OUT / "Figure3_predicted_functional_impact.svg"


# ---------------------------------------------------------------------
# Load frozen Figure 3 source data
# ---------------------------------------------------------------------

with panel_a_file.open(newline="") as fh:
    panel_a = list(csv.DictReader(fh, delimiter="\t"))

with panel_b_file.open(newline="") as fh:
    panel_b = list(csv.DictReader(fh, delimiter="\t"))

assert sum(int(r["Count"]) for r in panel_a) == 56
assert sum(int(r["Total_flagged"]) for r in panel_b) == 12


# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------

fig, (axA, axB) = plt.subplots(
    1, 2,
    figsize=(11.5, 5.2),
    gridspec_kw={"width_ratios": [1.05, 1.0]},
)

# =====================================================================
# Panel A — Overall SIFT classification
# =====================================================================

category_labels = {
    "deleterious": "Deleterious",
    "deleterious_low_confidence": "Deleterious\nlow confidence",
    "tolerated": "Tolerated",
    "tolerated_low_confidence": "Tolerated\nlow confidence",
}

cats = [category_labels[r["SIFT_category"]] for r in panel_a]
counts = [int(r["Count"]) for r in panel_a]

bars = axA.bar(
    np.arange(len(cats)),
    counts,
    width=0.68,
)

for bar, n in zip(bars, counts):
    axA.text(
        bar.get_x() + bar.get_width() / 2,
        n + 0.8,
        str(n),
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

axA.set_xticks(np.arange(len(cats)))
axA.set_xticklabels(cats, fontsize=8.5)
axA.set_ylabel("Number of candidates")
axA.set_ylim(0, 43)
axA.set_title("A  Overall SIFT classification", loc="left", fontweight="bold")

axA.text(
    0.5,
    -0.20,
    "56 missense candidates with SIFT predictions; "
    "4 RefSeq-resolved CMAH variants unavailable; "
    "2 stop variants not applicable",
    transform=axA.transAxes,
    ha="center",
    va="top",
    fontsize=8,
    color="dimgrey",
    clip_on=False,
)

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)


# =====================================================================
# Panel B — Gene distribution of deleterious predictions
# =====================================================================

genes = [r["Gene"] for r in panel_b]
deleterious = [int(r["Deleterious"]) for r in panel_b]
low_conf = [int(r["Deleterious_low_confidence"]) for r in panel_b]

y = np.arange(len(genes))

b1 = axB.barh(
    y,
    deleterious,
    height=0.62,
    label="Deleterious",
)

b2 = axB.barh(
    y,
    low_conf,
    left=deleterious,
    height=0.62,
    label="Deleterious, low confidence",
)

for i, total in enumerate(
    [d + lc for d, lc in zip(deleterious, low_conf)]
):
    axB.text(
        total + 0.08,
        i,
        str(total),
        va="center",
        ha="left",
        fontsize=8.5,
        fontweight="bold",
    )

axB.set_yticks(y)
axB.set_yticklabels(genes, fontstyle="italic")
axB.invert_yaxis()
axB.set_xlabel("Number of candidates")
axB.set_xlim(0, 5.8)
axB.set_title(
    "B  Genes containing deleterious predictions",
    loc="left",
    fontweight="bold",
)

axB.legend(
    frameon=False,
    fontsize=8,
    loc="lower right",
)

axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)


# ---------------------------------------------------------------------
# Final layout
# ---------------------------------------------------------------------

fig.subplots_adjust(
    left=0.09,
    right=0.98,
    bottom=0.17,
    top=0.90,
    wspace=0.34,
)

fig.savefig(png_out, dpi=300, bbox_inches="tight")
fig.savefig(svg_out, bbox_inches="tight")

print(f"Saved: {png_out}")
print(f"Saved: {svg_out}")
print()
print("FIGURE 3 PLOT GENERATION: PASS")
