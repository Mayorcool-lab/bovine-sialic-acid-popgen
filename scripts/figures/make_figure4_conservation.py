#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


SOURCE = Path(
    "results/figures/source_data/Figure4_candidate_conservation.tsv"
)

OUT_PNG = Path(
    "results/figures/Figure4_evolutionary_conservation.png"
)

OUT_SVG = Path(
    "results/figures/Figure4_evolutionary_conservation.svg"
)


# ---------------------------------------------------------------------
# Read frozen Figure 4 source data
# ---------------------------------------------------------------------
with SOURCE.open(newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

assert len(rows) == 62

assessable = [
    r for r in rows
    if r["Conservation_status"] == "ASSESSABLE"
]

nonassessable = [
    r for r in rows
    if r["Conservation_status"] == "NOT_ASSESSABLE"
]

special = [
    r for r in rows
    if r["Conservation_status"] == "STOP_LOST_SPECIAL_CASE"
]

assert len(assessable) == 57
assert len(nonassessable) == 4
assert len(special) == 1


# ---------------------------------------------------------------------
# Gene order
# ---------------------------------------------------------------------
gene_order = [
    "GNE",
    "NANS",
    "NANP",
    "CMAS",
    "SLC35A1",
    "CMAH",
    "ST8SIA2",
    "ST8SIA4",
    "ST6GAL1",
    "ST3GAL1",
    "ST3GAL4",
    "NEU1",
    "NEU3",
]

gene_rank = {
    gene: i
    for i, gene in enumerate(gene_order)
}


# ---------------------------------------------------------------------
# Sort candidates
# Within each gene: stronger conservation first
# ---------------------------------------------------------------------
assessable.sort(
    key=lambda r: (
        gene_rank[r["Gene"]],
        float(r["Rate4Site_score"]),
    )
)


# ---------------------------------------------------------------------
# SIFT display categories
# ---------------------------------------------------------------------
category_style = {
    "deleterious": {
        "label": "Deleterious",
        "marker": "o",
        "color": "#B22222",
    },
    "deleterious_low_confidence": {
        "label": "Deleterious, low confidence",
        "marker": "s",
        "color": "#E67E22",
    },
    "tolerated": {
        "label": "Tolerated",
        "marker": "o",
        "color": "#4C78A8",
    },
    "tolerated_low_confidence": {
        "label": "Tolerated, low confidence",
        "marker": "s",
        "color": "#72A0C1",
    },
    "SIFT unavailable": {
        "label": "SIFT unavailable",
        "marker": "D",
        "color": "#6C6C6C",
    },
    "Not applicable": {
        "label": "SIFT not applicable",
        "marker": "X",
        "color": "#111111",
    },
}

observed_categories = {
    r["SIFT_category"]
    for r in assessable
}

unknown_categories = (
    observed_categories - set(category_style)
)

if unknown_categories:
    raise RuntimeError(
        f"Unexpected SIFT categories: "
        f"{sorted(unknown_categories)}"
    )


# ---------------------------------------------------------------------
# Create y positions with gaps between genes
# ---------------------------------------------------------------------
plot_rows = []
gene_positions = defaultdict(list)

y = 0
previous_gene = None

for r in assessable:
    gene = r["Gene"]

    if (
        previous_gene is not None
        and gene != previous_gene
    ):
        y += 0.7

    plot_rows.append((y, r))
    gene_positions[gene].append(y)

    y += 1
    previous_gene = gene

max_y = max(
    ypos
    for ypos, _ in plot_rows
)


# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
fig, ax = plt.subplots(
    figsize=(11.5, 14.5)
)


# ---------------------------------------------------------------------
# Zero reference
# ---------------------------------------------------------------------
ax.axvline(
    0,
    linewidth=0.9,
    linestyle="--",
    color="0.55",
    zorder=1,
)


# ---------------------------------------------------------------------
# Horizontal stems
# ---------------------------------------------------------------------
for ypos, r in plot_rows:
    score = float(
        r["Rate4Site_score"]
    )

    ax.plot(
        [0, score],
        [ypos, ypos],
        linewidth=0.8,
        color="0.78",
        zorder=1,
    )


# ---------------------------------------------------------------------
# Candidate points
# ---------------------------------------------------------------------
for category, style in category_style.items():

    subset = [
        (ypos, r)
        for ypos, r in plot_rows
        if r["SIFT_category"] == category
    ]

    if not subset:
        continue

    xs = [
        float(r["Rate4Site_score"])
        for ypos, r in subset
    ]

    ys = [
        ypos
        for ypos, r in subset
    ]

    ax.scatter(
        xs,
        ys,
        s=45,
        marker=style["marker"],
        color=style["color"],
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
    )


# ---------------------------------------------------------------------
# Protein-change labels
# ---------------------------------------------------------------------
yticks = [
    ypos
    for ypos, r in plot_rows
]

yticklabels = [
    r["Protein_change"]
    for ypos, r in plot_rows
]

ax.set_yticks(yticks)

ax.set_yticklabels(
    yticklabels,
    fontsize=8,
)

ax.tick_params(
    axis="y",
    length=0,
    pad=5,
)


# ---------------------------------------------------------------------
# X-axis limits based only on conservation scores
# ---------------------------------------------------------------------
xmin = min(
    float(r["Rate4Site_score"])
    for r in assessable
)

xmax = max(
    float(r["Rate4Site_score"])
    for r in assessable
)

xrange_ = xmax - xmin

ax.set_xlim(
    xmin - (0.12 * xrange_),
    xmax + (0.15 * xrange_),
)


# ---------------------------------------------------------------------
# Gene labels
# Dedicated left-hand column separate from protein-change labels
# ---------------------------------------------------------------------
for gene in gene_order:

    ys = gene_positions.get(gene)

    if not ys:
        continue

    centre = (
        sum(ys) / len(ys)
    )

    ax.text(
        -0.16,
        centre,
        gene,
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
        clip_on=False,
    )


# ---------------------------------------------------------------------
# Gene separators
# ---------------------------------------------------------------------
ordered_gene_ranges = []

for gene in gene_order:

    ys = gene_positions.get(gene)

    if ys:
        ordered_gene_ranges.append(
            (
                gene,
                min(ys),
                max(ys),
            )
        )

for i in range(
    len(ordered_gene_ranges) - 1
):

    _, _, current_max = (
        ordered_gene_ranges[i]
    )

    _, next_min, _ = (
        ordered_gene_ranges[i + 1]
    )

    separator = (
        current_max + next_min
    ) / 2

    ax.axhline(
        separator,
        linewidth=0.5,
        color="0.88",
        zorder=0,
    )


# ---------------------------------------------------------------------
# Axes
# ---------------------------------------------------------------------
ax.set_xlabel(
    "Rate4Site conservation score\n"
    "(lower values indicate stronger evolutionary conservation)",
    fontsize=10,
)

ax.set_ylabel("")

ax.set_title(
    "Evolutionary conservation of protein-altering candidate residues",
    fontsize=12,
    fontweight="bold",
    pad=12,
)

ax.set_ylim(
    max_y + 1,
    -1,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)


# ---------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------
legend_handles = []

for category, style in category_style.items():

    if category not in observed_categories:
        continue

    legend_handles.append(
        Line2D(
            [0],
            [0],
            marker=style["marker"],
            linestyle="",
            markerfacecolor=style["color"],
            markeredgecolor="white",
            markersize=7,
            label=style["label"],
        )
    )

ax.legend(
    handles=legend_handles,
    title="SIFT classification",
    loc="lower right",
    frameon=False,
    fontsize=8,
    title_fontsize=9,
)


# ---------------------------------------------------------------------
# Explicit accounting for candidates lacking Rate4Site scores
#
# H5 distinguishes NEU1 L275P from the three ST3GAL4 candidates:
# - NEU1 L275P: candidate position falls within an alignment gap
# - ST3GAL4 P4L/R16Q/C38R: insufficient comparative sequence coverage
# ---------------------------------------------------------------------
alignment_gap = [
    r for r in nonassessable
    if (
        r["Gene"] == "NEU1"
        and r["Protein_change"] == "L275P"
    )
]

insufficient_coverage = [
    r for r in nonassessable
    if not (
        r["Gene"] == "NEU1"
        and r["Protein_change"] == "L275P"
    )
]

assert len(alignment_gap) == 1
assert len(insufficient_coverage) == 3

assert {
    (r["Gene"], r["Protein_change"])
    for r in insufficient_coverage
} == {
    ("ST3GAL4", "P4L"),
    ("ST3GAL4", "R16Q"),
    ("ST3GAL4", "C38R"),
}

alignment_gap_text = ", ".join(
    f'{r["Gene"]} {r["Protein_change"]}'
    for r in alignment_gap
)

insufficient_coverage_text = ", ".join(
    f'{r["Gene"]} {r["Protein_change"]}'
    for r in insufficient_coverage
)

special_text = ", ".join(
    f'{r["Gene"]} {r["Protein_change"]}'
    for r in special
)

fig.text(
    0.50,
    0.012,
    "Rate4Site assessable: 57/62. "
    f"Alignment gap: {alignment_gap_text}. "
    f"Insufficient comparative coverage: "
    f"{insufficient_coverage_text}. "
    f"Stop-lost special case: {special_text}.",
    ha="center",
    va="bottom",
    fontsize=8,
    color="dimgrey",
)


# ---------------------------------------------------------------------
# Final validation
# ---------------------------------------------------------------------
assert len(yticks) == 57

scores = [
    float(r["Rate4Site_score"])
    for r in assessable
]

assert (
    abs(
        min(scores) - (-0.9582)
    )
    < 1e-6
)


# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------
fig.subplots_adjust(
    left=0.27,
    right=0.97,
    top=0.95,
    bottom=0.085,
)


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------
fig.savefig(
    OUT_PNG,
    dpi=300,
    bbox_inches="tight",
)

fig.savefig(
    OUT_SVG,
    bbox_inches="tight",
)

plt.close(fig)


print(
    f"Saved: {OUT_PNG}"
)

print(
    f"Saved: {OUT_SVG}"
)

print()

print(
    "FIGURE 4 PLOT GENERATION: PASS"
)
