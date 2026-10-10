#!/usr/bin/env python3
"""Figure 9: Breed-, country-, and region-level FST Manhattan plots."""

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from H13_supplementary_fst import (
    load_genes,
    load_fst,
    chromosome_layout,
    annotate_genes,
    CHROMOSOMES,
    LABEL_ROWS,
)

OUT = ROOT / "results/H_integration/H13_supplementary_fst"
OUT.mkdir(parents=True, exist_ok=True)

LEVELS = ("breed", "country", "region")
PANELS = {"breed": "A", "country": "B", "region": "C"}

genes = load_genes()
datasets = {}

# Load and validate each population-grouping analysis.
for level in LEVELS:
    fst, n_input, threshold = load_fst(level)
    offsets, centres, running = chromosome_layout(fst)
    gene_summary = annotate_genes(fst, genes, offsets)

    datasets[level] = {
        "fst": fst,
        "n_input": n_input,
        "threshold": threshold,
        "offsets": offsets,
        "centres": centres,
        "running": running,
        "gene_summary": gene_summary,
    }

# Confirm shared genomic positions and chromosome layout.
reference = pd.read_csv(
    ROOT / "data/external/aiyesa_popgen/breed_all.weir.fst",
    sep=r"\s+",
    usecols=["CHROM", "POS"],
)

for level in ("country", "region"):
    comparison = pd.read_csv(
        ROOT / f"data/external/aiyesa_popgen/{level}_all.weir.fst",
        sep=r"\s+",
        usecols=["CHROM", "POS"],
    )

    if not reference.equals(comparison):
        raise ValueError(
            f"{level} genomic positions differ from breed positions"
        )

    if datasets[level]["offsets"] != datasets["breed"]["offsets"]:
        raise ValueError(
            f"{level} chromosome offsets differ from breed offsets"
        )

print("All three datasets have matching genomic positions and layout.")

# Figure layout: one gene-label strip and three aligned FST panels.
fig = plt.figure(figsize=(18, 15))

gs = fig.add_gridspec(
    4, 1,
    height_ratios=[1.0, 3.5, 3.5, 3.5],
    left=0.065,
    right=0.985,
    bottom=0.075,
    top=0.865,
    hspace=0.13,
)

ax_genes = fig.add_subplot(gs[0])
axes = {}

for i, level in enumerate(LEVELS):
    axes[level] = fig.add_subplot(
        gs[i + 1],
        sharex=ax_genes,
    )

for level in LEVELS:
    ax = axes[level]
    d = datasets[level]
    fst = d["fst"]
    threshold = d["threshold"]

    outliers = fst.loc[fst["OUTLIER"]]
    in_genes = fst.loc[fst["GENE"] != ""]
    gene_outliers = in_genes.loc[in_genes["OUTLIER"]]

    for chrom in CHROMOSOMES:
        subset = fst.loc[fst["CHROM"] == chrom]

        colour = (
            "#444444" if chrom % 2 else "#A0A0A0"
        )

        ax.scatter(
            subset["X"],
            subset["FST_PLOT"],
            s=1.0,
            alpha=0.45,
            color=colour,
            linewidths=0,
            rasterized=True,
        )

    ax.scatter(
        outliers["X"],
        outliers["FST_PLOT"],
        s=4,
        color="#D62728",
        alpha=0.85,
        linewidths=0,
        rasterized=True,
        zorder=3,
    )

    ax.scatter(
        in_genes["X"],
        in_genes["FST_PLOT"],
        s=12,
        color="#087E8B",
        alpha=0.95,
        edgecolors="white",
        linewidths=0.2,
        rasterized=True,
        zorder=4,
    )

    if not gene_outliers.empty:
        ax.scatter(
            gene_outliers["X"],
            gene_outliers["FST_PLOT"],
            s=30,
            facecolors="#FFBF00",
            edgecolors="black",
            linewidths=0.5,
            rasterized=True,
            zorder=5,
        )

    ax.axhline(
        threshold,
        color="#2855B8",
        linestyle="--",
        linewidth=1.2,
    )

    ax.set_xlim(0, d["running"])
    ax.set_ylim(0, 1.0)
    ax.set_ylabel(r"Per-site $F_{ST}$", fontsize=10)

    ax.grid(axis="y", alpha=0.12)
    ax.spines[["top", "right"]].set_visible(False)

    ax.text(
        0.008,
        0.96,
        f"({PANELS[level]}) {level.capitalize()}-level FST",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        fontweight="bold",
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.85,
        },
    )

    if level != "region":
        ax.tick_params(
            axis="x",
            which="both",
            bottom=False,
            labelbottom=False,
        )

    print(
        f"{level}: "
        f"finite SNPs={len(fst):,}, "
        f"threshold={threshold:.6f}, "
        f"pathway SNPs={len(in_genes)}, "
        f"pathway outliers={len(gene_outliers)}"
    )

# Show gene annotations once above all panels.
breed = datasets["breed"]

ax_genes.set_xlim(0, breed["running"])
ax_genes.set_ylim(0, 1)
ax_genes.axis("off")

for g in breed["gene_summary"].itertuples(index=False):
    x = g.gene_mid_x
    label_y = LABEL_ROWS.get(g.gene, 0.88)

    ax_genes.plot(
        [x, x],
        [0.04, 0.10],
        color="#087E8B",
        linewidth=1.0,
        alpha=0.8,
    )

    ax_genes.text(
        x,
        label_y,
        g.gene,
        fontsize=8,
        ha="center",
        va="bottom",
        color="#125E68",
    )

axes["region"].set_xticks(
    [breed["centres"][c] for c in CHROMOSOMES]
)

axes["region"].set_xticklabels(
    [str(c) for c in CHROMOSOMES],
    fontsize=8,
)

axes["region"].set_xlabel(
    "Bovine autosome",
    fontsize=11,
)

fig.suptitle(
    "Genetic differentiation across bovine autosomes "
    "at three population-grouping levels",
    fontsize=14,
    y=0.975,
)

legend_handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markerfacecolor="#D62728",
        markeredgecolor="none",
        markersize=5,
        label="Top 0.1% of per-site FST",
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markerfacecolor="#087E8B",
        markeredgecolor="none",
        markersize=6,
        label="SNP within pathway gene",
    ),
    Line2D(
        [0], [0],
        color="#2855B8",
        linestyle="--",
        linewidth=1.2,
        label="Group-specific empirical 99.9th percentile",
    ),
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.935),
    ncol=3,
    frameon=False,
    fontsize=9,
)

fig.text(
    0.5,
    0.028,
    "Negative FST estimates are displayed as zero. "
    "Empirical thresholds: breed = 0.452; "
    "country = 0.475; region = 0.580.",
    ha="center",
    fontsize=9,
    color="#444444",
)

prefix = "figure9_breed_country_region_fst"

for extension in ("png", "pdf", "svg"):
    kwargs = {"facecolor": "white"}

    if extension == "png":
        kwargs["dpi"] = 300

    fig.savefig(
        OUT / f"{prefix}.{extension}",
        **kwargs,
    )

plt.close(fig)

print(
    "\nThree-panel Figure 9 saved to: "
    f"{OUT.relative_to(ROOT)}"
)
