#!/usr/bin/env python3
"""
Figure 2: Population- and breed-level distribution of protein-altering variants.

Reads exclusively from the frozen source-data tables produced by:
  scripts/figures/generate_figure2_source_data.py

Inputs:
  results/figures/source_data/Figure2_panelA_population.tsv
  results/figures/source_data/Figure2_panelB_breed_heatmap.tsv

Outputs:
  results/figures/Figure2_population_breed_distribution.png
  results/figures/Figure2_population_breed_distribution.svg
"""

import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


plt.rcParams.update({
    "font.size": 9,
    "font.family": "sans-serif",
    "axes.linewidth": 0.8,
})

SOURCE_DIR = "results/figures/source_data"
OUT_DIR = "results/figures"

PANEL_A_FILE = f"{SOURCE_DIR}/Figure2_panelA_population.tsv"
PANEL_B_FILE = f"{SOURCE_DIR}/Figure2_panelB_breed_heatmap.tsv"


# ---------------------------------------------------------------------
# Fixed biological/display ordering
# ---------------------------------------------------------------------

gene_order = [
    "GNE", "NANS", "NANP", "CMAS", "SLC35A1", "CMAH",
    "ST8SIA2", "ST8SIA4",
    "ST6GAL1", "ST3GAL1", "ST3GAL4",
    "NEU1", "NEU3",
]

sa_breeds = [
    "Afrikaner",
    "Bonsmara",
    "Nguni",
    "Tuli",
]

eu_breeds = [
    "Barrosa",
    "Deep_Red",
    "Dutch_Belted",
    "Dutch_Friesian",
    "Eastern_Finncattle",
    "Groningen_White_Headed",
    "Holstein_Friesian",
    "Mertolenga",
    "Meuse_Rhine_Yssel",
    "Mirandesa",
    "Northern_Finncattle",
    "Western_Finncattle",
]

breed_order = sa_breeds + eu_breeds


# ---------------------------------------------------------------------
# Load Panel A
# ---------------------------------------------------------------------

with open(PANEL_A_FILE, newline="") as f:
    pop_rows = list(csv.DictReader(f, delimiter="\t"))

if len(pop_rows) != 62:
    raise RuntimeError(f"Expected 62 Panel A rows, got {len(pop_rows)}")

shared = [
    r for r in pop_rows
    if r["SA_status"] == "PRESENT" and r["EU_status"] == "PRESENT"
]

if len(shared) != 31:
    raise RuntimeError(f"Expected 31 shared callable candidates, got {len(shared)}")

sa_af = np.array([float(r["SA_AF"]) for r in shared])
eu_af = np.array([float(r["EU_AF"]) for r in shared])

for r in shared:
    r["_delta"] = abs(float(r["SA_AF"]) - float(r["EU_AF"]))

top5 = sorted(shared, key=lambda r: r["_delta"], reverse=True)[:5]


# ---------------------------------------------------------------------
# Load Panel B
# ---------------------------------------------------------------------

with open(PANEL_B_FILE, newline="") as f:
    breed_rows = list(csv.DictReader(f, delimiter="\t"))

if len(breed_rows) != 992:
    raise RuntimeError(f"Expected 992 Panel B rows, got {len(breed_rows)}")

observed_breeds = set(r["breed"] for r in breed_rows)
if observed_breeds != set(breed_order):
    raise RuntimeError(
        "Breed set mismatch.\n"
        f"Expected: {breed_order}\n"
        f"Observed: {sorted(observed_breeds)}"
    )

candidate_meta = {}
candidate_rows = defaultdict(dict)

for r in breed_rows:
    cid = r["candidate_id"]

    candidate_meta[cid] = {
        "gene": r["gene"],
        "protein_change": r["protein_change"],
        "variant_label": r["variant_label"],
    }

    if r["breed"] in candidate_rows[cid]:
        raise RuntimeError(
            f"Duplicate candidate-breed combination: {cid} / {r['breed']}"
        )

    candidate_rows[cid][r["breed"]] = r

if len(candidate_meta) != 62:
    raise RuntimeError(
        f"Expected 62 unique candidates in Panel B, got {len(candidate_meta)}"
    )

for cid in candidate_meta:
    if len(candidate_rows[cid]) != 16:
        raise RuntimeError(
            f"Candidate {cid} has {len(candidate_rows[cid])} breed rows, expected 16"
        )


def protein_position(protein_change):
    """
    Extract numeric protein position from compact notation:
      Y136H -> 136
      E128* -> 128
      *396L -> 396
    """
    digits = "".join(ch for ch in protein_change if ch.isdigit())
    if not digits:
        raise RuntimeError(
            f"Cannot extract protein position from {protein_change}"
        )
    return int(digits)


gene_rank = {g: i for i, g in enumerate(gene_order)}

candidate_order = sorted(
    candidate_meta,
    key=lambda cid: (
        gene_rank[candidate_meta[cid]["gene"]],
        protein_position(candidate_meta[cid]["protein_change"]),
        cid,
    )
)

row_labels = [
    candidate_meta[cid]["variant_label"]
    for cid in candidate_order
]


# ---------------------------------------------------------------------
# Construct AF matrix and status masks
# ---------------------------------------------------------------------

n_rows = len(candidate_order)
n_cols = len(breed_order)

af_matrix = np.full((n_rows, n_cols), np.nan)
source_missing = np.zeros((n_rows, n_cols), dtype=bool)
subthreshold = np.zeros((n_rows, n_cols), dtype=bool)
very_small_n = np.zeros((n_rows, n_cols), dtype=bool)

for i, cid in enumerate(candidate_order):
    for j, breed in enumerate(breed_order):
        r = candidate_rows[cid][breed]
        status = r["status"]

        if status == "NOT_DETECTED_IN_SOURCE":
            source_missing[i, j] = True
            continue

        if r["allele_frequency"] == "":
            raise RuntimeError(
                f"Estimable status lacks AF: {cid} / {breed} / {status}"
            )

        af_matrix[i, j] = float(r["allele_frequency"])

        if status == "ABSENT_AFTER_QC":
            subthreshold[i, j] = True
        elif status != "PRESENT":
            raise RuntimeError(
                f"Unexpected breed status: {status}"
            )

        if r["sample_size_flag"] == "VERY_SMALL_N":
            very_small_n[i, j] = True


if np.sum(subthreshold) != 16:
    raise RuntimeError(
        f"Expected 16 sub-threshold cells, got {np.sum(subthreshold)}"
    )

if np.sum(source_missing) != 276:
    raise RuntimeError(
        f"Expected 276 source-unavailable cells, got {np.sum(source_missing)}"
    )

if np.sum(very_small_n) != 82:
    raise RuntimeError(
        f"Expected 82 very-small-N cells, got {np.sum(very_small_n)}"
    )


# ---------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------

fig = plt.figure(figsize=(14.5, 19))
gs = fig.add_gridspec(
    2, 1,
    height_ratios=[1.0, 4.2],
    hspace=0.12,
)

# Panel A uses two independent axes: scatter left, status summary right.
gsA = gs[0, 0].subgridspec(
    1, 2,
    width_ratios=[1.05, 1.0],
    wspace=0.38,
)

# =====================================================================
# Panel A: population-level allele frequencies
# =====================================================================

axA = fig.add_subplot(gsA[0, 0])
axA.set_title(
    "A   Population-level allele-frequency distribution",
    loc="left",
    fontweight="bold",
    fontsize=12,
)

axA.scatter(
    eu_af,
    sa_af,
    s=38,
    edgecolor="black",
    linewidth=0.45,
    alpha=0.85,
)

axA.plot(
    [0, 1],
    [0, 1],
    linestyle="--",
    linewidth=1,
    color="grey",
)

axA.set_xlim(-0.03, 1.03)
axA.set_ylim(-0.03, 1.03)
axA.set_aspect("equal", adjustable="box")

axA.set_xlabel("European cohort allele frequency")
axA.set_ylabel("South African cohort allele frequency")

axA.spines[["top", "right"]].set_visible(False)

# Label the five largest absolute population AF differences.
label_offsets = [
    (28, -12),   # NEU3 Y136H
    (28, 38),    # ST3GAL4 L59F
    (28, 20),    # ST3GAL4 R87C
    (28, 2),     # ST3GAL4 P43L
    (28, -18),   # NANS H16Q
]

for r, offset in zip(top5, label_offsets):
    x = float(r["EU_AF"])
    y = float(r["SA_AF"])

    axA.annotate(
        r["variant_label"],
        xy=(x, y),
        xytext=offset,
        textcoords="offset points",
        fontsize=7.5,
        arrowprops=dict(
            arrowstyle="-",
            linewidth=0.6,
            color="grey",
        ),
    )

axA.text(
    0.04,
    0.96,
    "31 candidates callable in both cohorts\n"
    "Dashed line: equal allele frequency\n"
    "Labels: five largest absolute ΔAF values",
    transform=axA.transAxes,
    va="top",
    ha="left",
    fontsize=7.5,
    color="dimgrey",
    bbox=dict(
        facecolor="white",
        edgecolor="none",
        alpha=0.9,
        pad=1.5,
    ),
)


# Population-status summary axis.
inset = fig.add_subplot(gsA[0, 1])

status_labels = [
    "Callable in both",
    "SA present;\nEU source unavailable",
    "EU present;\nSA source unavailable",
    "EU present;\nSA sub-threshold",
]

status_counts = [31, 21, 6, 4]

y = np.arange(len(status_labels))

bars = inset.barh(
    y,
    status_counts,
    edgecolor="black",
    linewidth=0.45,
)

inset.set_yticks(y)
inset.set_yticklabels(status_labels, fontsize=7.5)
inset.invert_yaxis()
inset.set_xlabel("Candidates (n)", fontsize=8)
inset.set_title(
    "Representation status across all 62 candidates",
    fontsize=9,
    fontweight="bold",
    pad=6,
)
inset.spines[["top", "right"]].set_visible(False)

for bar, n in zip(bars, status_counts):
    inset.text(
        n + 0.5,
        bar.get_y() + bar.get_height() / 2,
        str(n),
        va="center",
        fontsize=8,
    )

inset.set_xlim(0, 35)


# =====================================================================
# Panel B: full 62 × 16 breed-level heatmap
# =====================================================================

axB = fig.add_subplot(gs[1, 0])
axB.set_title(
    "B   Breed-level allele-frequency distribution",
    loc="left",
    fontweight="bold",
    fontsize=12,
    pad=12,
)

# White -> blue AF gradient; source-unavailable cells grey.
cmap = LinearSegmentedColormap.from_list(
    "af_gradient",
    ["#FFFFFF", "#D7E6F5", "#7FA9D1", "#315F8C"],
)
cmap.set_bad("#D9D9D9")

masked_af = np.ma.masked_invalid(af_matrix)

im = axB.imshow(
    masked_af,
    aspect="auto",
    interpolation="none",
    cmap=cmap,
    vmin=0,
    vmax=1,
)

axB.set_xticks(np.arange(n_cols))
axB.set_xticklabels(
    [b.replace('_', ' ') for b in breed_order],
    rotation=60,
    ha="right",
    rotation_mode="anchor",
    fontsize=7.5,
)

axB.set_yticks(np.arange(n_rows))
axB.set_yticklabels(
    row_labels,
    fontsize=6.2,
)

axB.tick_params(axis="both", length=0)

# Thin cell boundaries.
axB.set_xticks(
    np.arange(-0.5, n_cols, 1),
    minor=True,
)
axB.set_yticks(
    np.arange(-0.5, n_rows, 1),
    minor=True,
)
axB.grid(
    which="minor",
    color="white",
    linewidth=0.35,
)
axB.tick_params(which="minor", bottom=False, left=False)

# Strong divider between South African and European breeds.
axB.axvline(
    len(sa_breeds) - 0.5,
    color="black",
    linewidth=1.6,
)

# Cohort labels and brackets below breed names.
sa_mid = (len(sa_breeds) - 1) / 2
eu_mid = len(sa_breeds) + (len(eu_breeds) - 1) / 2

axB.text(
    sa_mid,
    -0.105,
    "South Africa",
    transform=axB.get_xaxis_transform(),
    ha="center",
    va="top",
    fontweight="bold",
    fontsize=8.5,
    clip_on=False,
)

axB.text(
    eu_mid,
    -0.105,
    "Europe",
    transform=axB.get_xaxis_transform(),
    ha="center",
    va="top",
    fontweight="bold",
    fontsize=8.5,
    clip_on=False,
)

# Simple cohort brackets.
bracket_y = -0.095

axB.plot(
    [-0.35, len(sa_breeds) - 0.65],
    [bracket_y, bracket_y],
    transform=axB.get_xaxis_transform(),
    color="black",
    linewidth=0.8,
    clip_on=False,
)
axB.plot(
    [len(sa_breeds) - 0.35, n_cols - 0.65],
    [bracket_y, bracket_y],
    transform=axB.get_xaxis_transform(),
    color="black",
    linewidth=0.8,
    clip_on=False,
)


# Overlay sub-threshold estimates with ×.
sub_y, sub_x = np.where(subthreshold)
axB.scatter(
    sub_x,
    sub_y,
    marker="x",
    s=22,
    linewidth=0.9,
    color="black",
    zorder=4,
)

# Overlay VERY_SMALL_N estimates with open circles.
small_y, small_x = np.where(very_small_n)
axB.scatter(
    small_x,
    small_y,
    marker="o",
    s=28,
    facecolors="none",
    edgecolors="black",
    linewidth=0.75,
    zorder=5,
)


# Horizontal separators between genes.
previous_gene = None
for i, cid in enumerate(candidate_order):
    gene = candidate_meta[cid]["gene"]

    if previous_gene is not None and gene != previous_gene:
        axB.axhline(
            i - 0.5,
            color="black",
            linewidth=0.9,
        )

    previous_gene = gene


# Gene labels on far left.
gene_start_end = defaultdict(list)
for i, cid in enumerate(candidate_order):
    gene_start_end[candidate_meta[cid]["gene"]].append(i)

for gene in gene_order:
    indices = gene_start_end[gene]
    midpoint = (min(indices) + max(indices)) / 2

    axB.text(
        -2.15,
        midpoint,
        gene,
        ha="right",
        va="center",
        fontsize=7,
        fontweight="bold",
        clip_on=False,
    )


# AF colour bar.
cbar = fig.colorbar(
    im,
    ax=axB,
    fraction=0.018,
    pad=0.025,
)
cbar.set_label("Allele frequency", rotation=90)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])


# Heatmap-state legend.
source_patch = mpatches.Patch(
    facecolor="#D9D9D9",
    edgecolor="black",
    linewidth=0.4,
    label="Not represented in source VCF",
)

sub_handle = plt.Line2D(
    [0],
    [0],
    marker="x",
    color="black",
    linestyle="None",
    markersize=6,
    label="Sub-threshold site; AF estimated",
)

small_handle = plt.Line2D(
    [0],
    [0],
    marker="o",
    markerfacecolor="none",
    markeredgecolor="black",
    linestyle="None",
    markersize=6,
    label="Very small sample (n < 5)",
)

axB.legend(
    handles=[
        source_patch,
        sub_handle,
        small_handle,
    ],
    loc="upper center",
    bbox_to_anchor=(0.5, -0.145),
    frameon=False,
    fontsize=7.5,
    ncol=3,
)


fig.subplots_adjust(
    left=0.20,
    right=0.94,
    top=0.97,
    bottom=0.13,
)

plt.savefig(
    f"{OUT_DIR}/Figure2_population_breed_distribution.png",
    dpi=350,
    bbox_inches="tight",
    facecolor="white",
)

plt.savefig(
    f"{OUT_DIR}/Figure2_population_breed_distribution.svg",
    bbox_inches="tight",
    facecolor="white",
)

print(f"Saved: {OUT_DIR}/Figure2_population_breed_distribution.png")
print(f"Saved: {OUT_DIR}/Figure2_population_breed_distribution.svg")
print("\nFIGURE 2 PLOT GENERATION: PASS")
