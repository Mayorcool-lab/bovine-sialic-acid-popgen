#!/usr/bin/env python3
"""
Figure 1: Variant landscape across the 13-gene bovine sialic-acid network.

Reads exclusively from the frozen, versioned source-data tables produced by
scripts/figures/generate_figure1_source_data.py. This script performs no
independent recalculation from the raw Phase G outputs; it is a pure
plotting step, kept separate from data derivation for reproducibility.

Inputs:
  results/figures/source_data/Figure1_panelB_gene_landscape.tsv
  results/figures/source_data/Figure1_panelC_consequences.tsv
  results/figures/source_data/Figure1_panelD_candidates.tsv

Outputs:
  results/figures/Figure1_variant_landscape.png
  results/figures/Figure1_variant_landscape.svg
"""

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    "font.size": 9,
    "font.family": "sans-serif",
    "axes.linewidth": 0.8,
})

SOURCE_DIR = "results/figures/source_data"
OUT_DIR = "results/figures"

PANEL_B_FILE = f"{SOURCE_DIR}/Figure1_panelB_gene_landscape.tsv"
PANEL_C_FILE = f"{SOURCE_DIR}/Figure1_panelC_consequences.tsv"
PANEL_D_FILE = f"{SOURCE_DIR}/Figure1_panelD_candidates.tsv"

panel_b = {}
with open(PANEL_B_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        panel_b[r["gene"]] = (
            float(r["analysed_kb"]),
            int(r["target_gene_variants"]),
            float(r["variants_per_kb"]),
            int(r["protein_altering_candidates"]),
        )

if len(panel_b) != 13:
    raise RuntimeError(f"Expected 13 genes in Panel B source data, got {len(panel_b)}")

panel_d = {}
with open(PANEL_D_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        panel_d[r["gene"]] = {
            "missense": int(r["missense"]),
            "missense_splice": int(r["missense_splice"]),
            "stop_gained": int(r["stop_gained"]),
            "stop_lost": int(r["stop_lost"]),
        }

total_candidates = sum(sum(v.values()) for v in panel_d.values())
if total_candidates != 62:
    raise RuntimeError(f"Expected 62 total candidates in Panel D source data, got {total_candidates}")

panel_c = {}
with open(PANEL_C_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        panel_c[r["consequence_group"]] = int(r["count"])

if sum(panel_c.values()) != 9724:
    raise RuntimeError(f"Expected 9724 total variants in Panel C source data, got {sum(panel_c.values())}")

print("Data loaded and validated OK")

tiers = {
    "Core pathway": ["GNE", "NANS", "NANP", "CMAS", "SLC35A1", "CMAH"],
    "Polysialylation": ["ST8SIA2", "ST8SIA4"],
    "Terminal sialylation": ["ST6GAL1", "ST3GAL1", "ST3GAL4"],
    "Sialic acid removal": ["NEU1", "NEU3"],
}
tier_colors = {
    "Core pathway": "#4C72B0",
    "Polysialylation": "#DD8452",
    "Terminal sialylation": "#55A868",
    "Sialic acid removal": "#C44E52",
}

gene_order = []
gene_tier_color = {}
for tier, genes in tiers.items():
    for g in genes:
        if g not in panel_b:
            raise RuntimeError(f"Gene {g} in fixed tier structure not found in Panel B source data")
        gene_order.append(g)
        gene_tier_color[g] = tier_colors[tier]

consequence_types = ["missense", "missense_splice", "stop_gained", "stop_lost"]
consequence_colors = {
    "missense": "#4C72B0",
    "missense_splice": "#8172B2",
    "stop_gained": "#C44E52",
    "stop_lost": "#CCB974",
}
consequence_labels = {
    "missense": "Missense",
    "missense_splice": "Missense + splice-region",
    "stop_gained": "Stop-gained",
    "stop_lost": "Stop-lost",
}

fig = plt.figure(figsize=(12.5, 9.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.38, wspace=0.28)

axA = fig.add_subplot(gs[0, 0])
axA.set_xlim(0, 10)
axA.set_ylim(0, 10)
axA.axis("off")
axA.set_title("A   Pathway overview", loc="left", fontweight="bold", fontsize=12, pad=18)

tier_y = {"Core pathway": 8.3, "Polysialylation": 5.9, "Terminal sialylation": 4.0, "Sialic acid removal": 1.8}
tier_step = {"Core pathway": 1.5, "Polysialylation": 2.2, "Terminal sialylation": 2.6, "Sialic acid removal": 2.6}

for tier, genes in tiers.items():
    color = tier_colors[tier]
    y = tier_y[tier]
    axA.text(-0.3, y, tier, fontsize=8, va="center", ha="right", fontweight="bold", color=color)
    x0 = 0.3
    step = tier_step[tier]
    for i, g in enumerate(genes):
        x = x0 + i * step
        box = FancyBboxPatch((x, y - 0.4), step - 0.3, 0.8,
                              boxstyle="round,pad=0.02,rounding_size=0.08",
                              linewidth=1, edgecolor=color, facecolor=color, alpha=0.75)
        axA.add_patch(box)
        axA.text(x + (step - 0.3) / 2, y, g, ha="center", va="center",
                  fontsize=7.5, color="white", fontweight="bold")
    if tier != "Sialic acid removal":
        axA.annotate("", xy=(x0, y - 0.9), xytext=(x0, y - 0.5),
                      arrowprops=dict(arrowstyle="-|>", color="grey", lw=1.2))

axA.text(0.3, 9.5, "Sialic acid biosynthesis, activation, transport,\nsialylation and removal",
          fontsize=7.5, style="italic", color="dimgrey")

axB = fig.add_subplot(gs[0, 1])
axB.set_title("B   Gene-level variant landscape", loc="left", fontweight="bold", fontsize=12)

x = np.arange(len(gene_order))
counts = [panel_b[g][1] for g in gene_order]
density = [panel_b[g][2] for g in gene_order]
colors = [gene_tier_color[g] for g in gene_order]

axB.bar(x, counts, color=colors, edgecolor="black", linewidth=0.4, width=0.65)
axB.set_ylabel("Target-gene variants (n)")
axB.set_xticks(x)
axB.set_xticklabels(gene_order, rotation=55, ha="right", fontsize=7.5)
axB.spines[["top", "right"]].set_visible(False)

ax2 = axB.twinx()
ax2.plot(x, density, marker="o", color="black", linewidth=1.2, markersize=4)
ax2.set_ylabel("Variants / kb")
ax2.spines[["top"]].set_visible(False)

tier_patches = [mpatches.Patch(color=c, label=t) for t, c in tier_colors.items()]
line_patch = plt.Line2D([0], [0], color="black", marker="o", markersize=4, linewidth=1.2, label="Variants/kb")
axB.legend(handles=tier_patches + [line_patch], fontsize=6.5, loc="upper left", frameon=False, ncol=1)

axC = fig.add_subplot(gs[1, 0])
axC.set_title("C   Consequence composition of target-gene variants", loc="left", fontweight="bold", fontsize=12)

labels_sorted = sorted(panel_c.items(), key=lambda kv: kv[1])
labels_c = [l for l, v in labels_sorted]
vals_c = [v for l, v in labels_sorted]
palette_map = {
    "Intronic": "#8C8C8C",
    "Upstream/downstream": "#B0B0B0",
    "UTR": "#CFCFCF",
    "Synonymous": "#55A868",
    "Protein-altering": "#4C72B0",
    "Splice-associated": "#DD8452",
}
colors_c = [palette_map.get(l, "#999999") for l in labels_c]

ypos = range(len(labels_c))
bars_c = axC.barh(ypos, vals_c, color=colors_c, edgecolor="black", linewidth=0.5, height=0.65)
axC.set_yticks(ypos)
axC.set_yticklabels([l.replace("/", "/\n") if len(l) > 12 else l for l in labels_c], fontsize=8)
axC.set_xscale("log")
axC.set_xlabel("Variants (n), log scale")
axC.spines[["top", "right"]].set_visible(False)

for bar, v in zip(bars_c, vals_c):
    pct = 100 * v / sum(vals_c)
    axC.text(bar.get_width() * 1.15, bar.get_y() + bar.get_height() / 2,
              f"n={v:,} ({pct:.1f}%)", va="center", fontsize=7.5)

axC.set_xlim(right=axC.get_xlim()[1] * 4)

axD = fig.add_subplot(gs[1, 1])
axD.set_title("D   Protein-altering candidates across genes", loc="left", fontweight="bold", fontsize=12)

bottoms = np.zeros(len(gene_order))
for ctype in consequence_types:
    vals = [panel_d.get(g, {}).get(ctype, 0) for g in gene_order]
    axD.bar(x, vals, bottom=bottoms, color=consequence_colors[ctype],
            edgecolor="black", linewidth=0.4, width=0.65, label=consequence_labels[ctype])
    bottoms += np.array(vals)

axD.set_ylabel("Protein-altering candidates (n)")
axD.set_xticks(x)
axD.set_xticklabels(gene_order, rotation=55, ha="right", fontsize=7.5)
axD.spines[["top", "right"]].set_visible(False)
axD.legend(fontsize=6.5, loc="upper left", frameon=False)
axD.text(0.5, -0.42, "62 candidates total, including 4 RefSeq-resolved\nCMAH missense variants (Section 2.4.1)",
          transform=axD.transAxes, ha="center", fontsize=7, style="italic", color="dimgrey")

fig.text(0.25, 0.015, "9,724 VEP-annotated target-gene variants (58 VEP-native protein-altering)",
         ha="center", fontsize=7.5, style="italic", color="dimgrey")

plt.savefig(f"{OUT_DIR}/Figure1_variant_landscape.png", dpi=350, bbox_inches="tight", facecolor="white")
plt.savefig(f"{OUT_DIR}/Figure1_variant_landscape.svg", bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_DIR}/Figure1_variant_landscape.png")
print(f"Saved: {OUT_DIR}/Figure1_variant_landscape.svg")
print("\nFIGURE 1 PLOT GENERATION: PASS")
