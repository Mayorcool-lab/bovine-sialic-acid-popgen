#!/usr/bin/env python3
"""Breed-level FST Manhattan plot with readable pathway-gene annotations."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/external/aiyesa_popgen/breed_all.weir.fst"
GENES = ROOT / "data/reference/candidate_gene_coordinates_provisional.tsv"
OUT = ROOT / "results/H_integration/H12_figure9"
OUT.mkdir(parents=True, exist_ok=True)

# Read original FST data without modifying source files.
fst = pd.read_csv(INPUT, sep=r"\s+")
genes = pd.read_csv(GENES, sep="\t")

required_fst = {"CHROM", "POS", "WEIR_AND_COCKERHAM_FST"}
required_genes = {"gene_symbol", "chromosome", "start", "end"}

if not required_fst.issubset(fst.columns):
    raise ValueError("Missing required FST columns")
if not required_genes.issubset(genes.columns):
    raise ValueError("Missing required gene-coordinate columns")
if len(genes) != 13 or genes["gene_symbol"].nunique() != 13:
    raise ValueError("Expected exactly 13 distinct pathway genes")

fst["CHROM"] = pd.to_numeric(fst["CHROM"], errors="coerce")
fst["POS"] = pd.to_numeric(fst["POS"], errors="coerce")
fst["FST_RAW"] = pd.to_numeric(
    fst["WEIR_AND_COCKERHAM_FST"], errors="coerce"
)

n_input = len(fst)

valid = (
    np.isfinite(fst["CHROM"])
    & np.isfinite(fst["POS"])
    & np.isfinite(fst["FST_RAW"])
)

fst = fst.loc[valid].copy()

if not (
    fst["CHROM"].between(1, 29).all()
    and (fst["CHROM"] % 1 == 0).all()
    and (fst["POS"] > 0).all()
    and (fst["POS"] % 1 == 0).all()
):
    raise ValueError("Invalid chromosome or position values")

fst["CHROM"] = fst["CHROM"].astype(int)
fst["POS"] = fst["POS"].astype(int)

# Preserve raw negative estimates; clip only for plotting.
fst["FST_PLOT"] = fst["FST_RAW"].clip(lower=0)

threshold = float(fst["FST_PLOT"].quantile(0.999))
fst["OUTLIER"] = fst["FST_PLOT"] >= threshold

# Cumulative coordinates follow Aiyesa's original plotting convention.
chromosomes = list(range(1, 30))
max_positions = fst.groupby("CHROM")["POS"].max()

missing = set(chromosomes) - set(max_positions.index)
if missing:
    raise ValueError(f"Missing chromosomes: {sorted(missing)}")

offsets = {}
centres = {}
running = 0

for chrom in chromosomes:
    length = int(max_positions.loc[chrom])
    offsets[chrom] = running
    centres[chrom] = running + length / 2
    running += length

fst["X"] = fst["POS"] + fst["CHROM"].map(offsets)
fst["GENE"] = ""

gene_rows = []

for g in genes.itertuples(index=False):
    chrom = int(g.chromosome)
    start = int(g.start)
    end = int(g.end)

    if chrom not in offsets or start < 1 or end < start:
        raise ValueError(f"Invalid coordinates for {g.gene_symbol}")

    mask = (
        (fst["CHROM"] == chrom)
        & fst["POS"].between(start, end)
    )

    if (fst.loc[mask, "GENE"] != "").any():
        raise ValueError(f"Overlapping pathway annotations: {g.gene_symbol}")

    fst.loc[mask, "GENE"] = g.gene_symbol

    gene_rows.append({
        "gene": g.gene_symbol,
        "chromosome": chrom,
        "start": start,
        "end": end,
        "fst_snp_count": int(mask.sum()),
        "fst_outlier_count": int((mask & fst["OUTLIER"]).sum()),
        "gene_mid_x": offsets[chrom] + (start + end) / 2,
    })

gene_summary = pd.DataFrame(gene_rows)
gene_summary.to_csv(
    OUT / "figure9_gene_fst_summary.tsv",
    sep="\t",
    index=False
)

outliers = fst.loc[fst["OUTLIER"]]
in_genes = fst.loc[fst["GENE"] != ""]
gene_outliers = in_genes.loc[in_genes["OUTLIER"]]

# Separate gene annotations from the Manhattan plot.
fig = plt.figure(figsize=(18, 8.5))

gs = fig.add_gridspec(
    2, 1,
    height_ratios=[1.5, 6],
    left=0.065,
    right=0.985,
    bottom=0.10,
    top=0.80,
    hspace=0.03
)

ax_genes = fig.add_subplot(gs[0])
ax = fig.add_subplot(gs[1], sharex=ax_genes)

# Plot all valid FST positions.
for chrom in chromosomes:
    subset = fst.loc[fst["CHROM"] == chrom]
    colour = "#444444" if chrom % 2 else "#A0A0A0"

    ax.scatter(
        subset["X"],
        subset["FST_PLOT"],
        s=1.2,
        alpha=0.5,
        color=colour,
        linewidths=0,
        rasterized=True
    )

ax.scatter(
    outliers["X"],
    outliers["FST_PLOT"],
    s=5,
    color="#D62728",
    alpha=0.8,
    linewidths=0,
    rasterized=True,
    zorder=3
)

ax.scatter(
    in_genes["X"],
    in_genes["FST_PLOT"],
    s=15,
    color="#087E8B",
    alpha=0.95,
    edgecolors="white",
    linewidths=0.2,
    rasterized=True,
    zorder=4
)

if not gene_outliers.empty:
    ax.scatter(
        gene_outliers["X"],
        gene_outliers["FST_PLOT"],
        s=34,
        facecolors="#FFBF00",
        edgecolors="black",
        linewidths=0.5,
        rasterized=True,
        zorder=5
    )

ax.axhline(
    threshold,
    color="#2855B8",
    linestyle="--",
    linewidth=1.2
)

ax.set_xlim(0, running)
ax.set_ylim(0, max(1.0, float(fst["FST_PLOT"].max()) * 1.04))
ax.set_xticks([centres[c] for c in chromosomes])
ax.set_xticklabels([str(c) for c in chromosomes], fontsize=8)
ax.set_xlabel("Bovine autosome", fontsize=11)
ax.set_ylabel(
    r"Per-site $F_{ST}$ (negative estimates shown as zero)",
    fontsize=10
)
ax.grid(axis="y", alpha=0.12)
ax.spines[["top", "right"]].set_visible(False)

# Horizontal gene labels in staggered rows.
# Use chromosome centres for crowded chromosome-specific groups.
ax_genes.set_xlim(0, running)
ax_genes.set_ylim(0, 1)
ax_genes.axis("off")

gene_summary = gene_summary.sort_values(
    ["chromosome", "gene_mid_x"]
).copy()

# Manual rows for nearby labels; other genes use the default top row.
label_rows = {
    "GNE": 0.88,
    "NANS": 0.58,
    "CMAH": 0.88,
    "NEU1": 0.58,
}

for g in gene_summary.itertuples(index=False):
    x = g.gene_mid_x
    label_y = label_rows.get(g.gene, 0.88)

    # Short position tick only; no line connecting or crossing gene labels.
    ax_genes.plot(
        [x, x],
        [0.04, 0.10],
        color="#087E8B",
        linewidth=1.0,
        alpha=0.8
    )

    ax_genes.text(
        x,
        label_y,
        g.gene,
        rotation=0,
        fontsize=8,
        fontweight="normal",
        ha="center",
        va="bottom",
        color="#125E68"
    )

fig.suptitle(
    "Breed-level genetic differentiation across bovine autosomes",
    fontsize=14,
    y=0.965
)

legend_handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markerfacecolor="#D62728",
        markeredgecolor="none",
        markersize=5,
        label="Top 0.1% of per-site FST"
    ),
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markerfacecolor="#087E8B",
        markeredgecolor="none",
        markersize=6,
        label="SNP within pathway gene"
    ),
    Line2D(
        [0], [0],
        color="#2855B8",
        linestyle="--",
        linewidth=1.2,
        label=f"Empirical 99.9th percentile ({threshold:.3f})"
    ),
]

if not gene_outliers.empty:
    legend_handles.append(
        Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            markerfacecolor="#FFBF00",
            markeredgecolor="black",
            markersize=6,
            label="Pathway-gene SNP in top 0.1%"
        )
    )

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.925),
    ncol=3,
    frameon=False,
    fontsize=9
)

fig.savefig(
    OUT / "figure9_breed_fst_v4.png",
    dpi=300,
    facecolor="white"
)
fig.savefig(
    OUT / "figure9_breed_fst_v4.pdf",
    facecolor="white"
)
fig.savefig(
    OUT / "figure9_breed_fst_v4.svg",
    facecolor="white"
)
plt.close(fig)

print("Figure 9 version 4")
print(f"Input FST rows: {n_input:,}")
print(f"Finite FST rows: {len(fst):,}")
print(f"Negative FST estimates: {(fst['FST_RAW'] < 0).sum():,}")
print(f"Empirical 99.9th percentile: {threshold:.6f}")
print(f"Genome-wide empirical outlier SNPs: {len(outliers):,}")
print(f"Pathway-gene SNPs: {len(in_genes):,}")
print(f"Pathway-gene empirical outliers: {len(gene_outliers):,}")
print(
    "Genes with at least one FST SNP: "
    f"{(gene_summary['fst_snp_count'] > 0).sum()}/13"
)
print(f"Output directory: {OUT.relative_to(ROOT)}")
