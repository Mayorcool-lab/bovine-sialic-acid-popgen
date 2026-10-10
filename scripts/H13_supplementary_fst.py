#!/usr/bin/env python3
"""Country- and region-level supplementary FST Manhattan plots."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/external/aiyesa_popgen"
GENES = ROOT / "data/reference/candidate_gene_coordinates_provisional.tsv"
OUT = ROOT / "results/H_integration/H13_supplementary_fst"

CHROMOSOMES = list(range(1, 30))
LEVELS = ("country", "region")

REQUIRED_FST = {"CHROM", "POS", "WEIR_AND_COCKERHAM_FST"}
REQUIRED_GENES = {"gene_symbol", "chromosome", "start", "end"}

LABEL_ROWS = {
    "GNE": 0.88,
    "NANS": 0.58,
    "CMAH": 0.88,
    "NEU1": 0.58,
}


def load_genes():
    genes = pd.read_csv(GENES, sep="\t")

    if not REQUIRED_GENES.issubset(genes.columns):
        raise ValueError("Missing required gene-coordinate columns")

    if len(genes) != 13 or genes["gene_symbol"].nunique() != 13:
        raise ValueError("Expected exactly 13 distinct pathway genes")

    for col in ("chromosome", "start", "end"):
        genes[col] = pd.to_numeric(genes[col], errors="raise")

    if not (
        genes["chromosome"].between(1, 29).all()
        and (genes["chromosome"] % 1 == 0).all()
        and (genes["start"] >= 1).all()
        and (genes["start"] % 1 == 0).all()
        and (genes["end"] >= genes["start"]).all()
        and (genes["end"] % 1 == 0).all()
    ):
        raise ValueError("Invalid gene coordinates")

    return genes


def load_fst(level):
    path = DATA / f"{level}_all.weir.fst"
    fst = pd.read_csv(path, sep=r"\s+")

    if not REQUIRED_FST.issubset(fst.columns):
        raise ValueError(f"Missing FST columns in {path}")

    n_input = len(fst)

    for col in ("CHROM", "POS", "WEIR_AND_COCKERHAM_FST"):
        fst[col] = pd.to_numeric(fst[col], errors="coerce")

    # Validate genomic positions before excluding nonfinite FST.
    positions_valid = (
        np.isfinite(fst["CHROM"])
        & np.isfinite(fst["POS"])
        & fst["CHROM"].between(1, 29)
        & (fst["CHROM"] % 1 == 0)
        & (fst["POS"] > 0)
        & (fst["POS"] % 1 == 0)
    )

    if not positions_valid.all():
        raise ValueError(f"Invalid genomic positions in {path}")

    if fst[["CHROM", "POS"]].duplicated().any():
        raise ValueError(f"Duplicate genomic positions in {path}")

    fst["CHROM"] = fst["CHROM"].astype(int)
    fst["POS"] = fst["POS"].astype(int)

    fst["FST_RAW"] = fst["WEIR_AND_COCKERHAM_FST"]

    fst = fst.loc[np.isfinite(fst["FST_RAW"])].copy()

    if fst.empty:
        raise ValueError(f"No finite FST values in {path}")

    fst["FST_PLOT"] = fst["FST_RAW"].clip(lower=0)

    threshold = float(fst["FST_PLOT"].quantile(0.999))
    fst["OUTLIER"] = fst["FST_PLOT"] >= threshold

    return fst, n_input, threshold


def chromosome_layout(fst):
    max_positions = fst.groupby("CHROM")["POS"].max()

    missing = set(CHROMOSOMES) - set(max_positions.index)
    if missing:
        raise ValueError(f"Missing chromosomes: {sorted(missing)}")

    offsets = {}
    centres = {}
    running = 0

    for chrom in CHROMOSOMES:
        length = int(max_positions.loc[chrom])
        offsets[chrom] = running
        centres[chrom] = running + length / 2
        running += length

    fst["X"] = fst["POS"] + fst["CHROM"].map(offsets)

    return offsets, centres, running


def annotate_genes(fst, genes, offsets):
    fst["GENE"] = ""
    gene_rows = []

    for g in genes.itertuples(index=False):
        chrom = int(g.chromosome)
        start = int(g.start)
        end = int(g.end)

        mask = (
            (fst["CHROM"] == chrom)
            & fst["POS"].between(start, end)
        )

        if (fst.loc[mask, "GENE"] != "").any():
            raise ValueError(
                f"Overlapping pathway annotations: {g.gene_symbol}"
            )

        fst.loc[mask, "GENE"] = g.gene_symbol

        gene_rows.append({
            "gene": g.gene_symbol,
            "chromosome": chrom,
            "start": start,
            "end": end,
            "fst_snp_count": int(mask.sum()),
            "fst_outlier_count": int(
                (mask & fst["OUTLIER"]).sum()
            ),
            "gene_mid_x": (
                offsets[chrom] + (start + end) / 2
            ),
        })

    return (
        pd.DataFrame(gene_rows)
        .sort_values(["chromosome", "gene_mid_x"])
        .reset_index(drop=True)
    )


def plot_fst(level, fst, threshold, gene_summary, centres, running):
    outliers = fst.loc[fst["OUTLIER"]]
    in_genes = fst.loc[fst["GENE"] != ""]
    gene_outliers = in_genes.loc[in_genes["OUTLIER"]]

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

    for chrom in CHROMOSOMES:
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
    ax.set_ylim(
        0,
        max(1.0, float(fst["FST_PLOT"].max()) * 1.04)
    )
    ax.set_xticks([centres[c] for c in CHROMOSOMES])
    ax.set_xticklabels(
        [str(c) for c in CHROMOSOMES],
        fontsize=8
    )
    ax.set_xlabel("Bovine autosome", fontsize=11)
    ax.set_ylabel(
        r"Per-site $F_{ST}$ (negative estimates shown as zero)",
        fontsize=10
    )
    ax.grid(axis="y", alpha=0.12)
    ax.spines[["top", "right"]].set_visible(False)

    # Preserve the approved H12 version-4 gene-label layout.
    ax_genes.set_xlim(0, running)
    ax_genes.set_ylim(0, 1)
    ax_genes.axis("off")

    for g in gene_summary.itertuples(index=False):
        x = g.gene_mid_x
        label_y = LABEL_ROWS.get(g.gene, 0.88)

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
        f"{level.capitalize()}-level genetic differentiation "
        "across bovine autosomes",
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
            label=(
                "Empirical 99.9th percentile "
                f"({threshold:.3f})"
            )
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

    prefix = f"supplementary_{level}_fst"

    for extension in ("png", "pdf", "svg"):
        kwargs = {"facecolor": "white"}
        if extension == "png":
            kwargs["dpi"] = 300

        fig.savefig(
            OUT / f"{prefix}.{extension}",
            **kwargs
        )

    plt.close(fig)

    return len(outliers), len(in_genes), len(gene_outliers)


def main():
    genes = load_genes()
    OUT.mkdir(parents=True, exist_ok=True)

    summaries = []

    for level in LEVELS:
        fst, n_input, threshold = load_fst(level)

        offsets, centres, running = chromosome_layout(fst)

        gene_summary = annotate_genes(
            fst, genes, offsets
        )

        gene_summary.to_csv(
            OUT / f"{level}_gene_fst_summary.tsv",
            sep="\t",
            index=False
        )

        n_outliers, n_gene_snps, n_gene_outliers = plot_fst(
            level,
            fst,
            threshold,
            gene_summary,
            centres,
            running
        )

        summary = {
            "level": level,
            "input_positions": n_input,
            "finite_positions": len(fst),
            "nonfinite_positions": n_input - len(fst),
            "negative_estimates": int(
                (fst["FST_RAW"] < 0).sum()
            ),
            "empirical_q999": threshold,
            "empirical_outlier_positions": n_outliers,
            "pathway_gene_positions": n_gene_snps,
            "pathway_gene_outliers": n_gene_outliers,
            "genes_with_positions": int(
                (gene_summary["fst_snp_count"] > 0).sum()
            ),
        }

        summaries.append(summary)

        print(f"\n=== {level.upper()} ===")
        for key, value in summary.items():
            print(f"{key}: {value}")

    pd.DataFrame(summaries).to_csv(
        OUT / "H13_fst_validation_summary.tsv",
        sep="\t",
        index=False
    )

    print(f"\nOutputs: {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
