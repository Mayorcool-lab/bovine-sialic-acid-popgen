#!/usr/bin/env python3
"""H14 supplementary figures for pathway SNP MAF and heterozygosity."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("results/H_integration/H14_maf_heterozygosity")
BASE.mkdir(parents=True, exist_ok=True)

maf = pd.read_csv(BASE / "H14_breed_maf.tsv", sep="\t")
het = pd.read_csv(BASE / "H14_heterozygosity_summary.tsv", sep="\t")

het = het.loc[het["breed"] != "ALL"].copy()

assert len(maf) == 4505 * 16
assert len(het) == 16
assert maf.groupby("breed").size().eq(4505).all()
assert set(maf["breed"]) == set(het["breed"])

het["n_samples"] = het["n_samples"].astype(int)

# Preserve consistent alphabetical breed order.
breeds = sorted(het["breed"])
het = het.set_index("breed").loc[breeds].reset_index()

sample_sizes = dict(zip(het["breed"], het["n_samples"]))

def classify_maf(value):
    if pd.isna(value):
        return "Uncallable"
    if value == 0:
        return "Monomorphic"
    if value < 0.05:
        return "Low-frequency"
    return "Common"

maf["category"] = maf["maf"].map(classify_maf)

categories = [
    "Monomorphic",
    "Low-frequency",
    "Common",
    "Uncallable",
]

counts = (
    maf.groupby(["breed", "category"])
    .size()
    .unstack(fill_value=0)
    .reindex(index=breeds, columns=categories, fill_value=0)
)

assert counts.sum(axis=1).eq(4505).all()

percentages = counts / 4505 * 100

labels = [
    f"{breed.replace('_', ' ')} (n={sample_sizes[breed]})"
    + (" *" if sample_sizes[breed] < 5 else "")
    for breed in breeds
]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 300,
})

# Figure S1: MAF distribution
fig, ax = plt.subplots(figsize=(11, 8.5))
y = np.arange(len(breeds))
left = np.zeros(len(breeds))

colors = {
    "Monomorphic": "#6C757D",
    "Low-frequency": "#D9A441",
    "Common": "#3174A1",
    "Uncallable": "#C75C5C",
}

for category in categories:
    values = percentages[category].to_numpy()
    ax.barh(
        y, values, left=left,
        label=category,
        color=colors[category],
        height=0.7,
    )
    left += values

ax.set_yticks(y, labels)
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xlabel("Percentage of pathway SNP records (%)")
ax.set_title(
    "Breed-specific minor allele frequency distribution",
    weight="bold",
    pad=15,
)
ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=4,
    frameon=False,
)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="x", alpha=0.2)
ax.set_axisbelow(True)

fig.text(
    0.5, 0.015,
    "* Small sample (n < 5). "
    "Percentages refer to 4,505 ascertained pathway SNP records.",
    ha="center",
    fontsize=8,
)

fig.tight_layout(rect=(0, 0.06, 1, 1))

for ext in ("png", "pdf", "svg"):
    fig.savefig(
        BASE / f"H14_figure_S1_breed_maf_distribution.{ext}",
        bbox_inches="tight",
    )
plt.close(fig)

# Figure S2: Ho and He
fig, ax = plt.subplots(figsize=(11, 8.5))
y = np.arange(len(breeds))
height = 0.36

ho = het["observed_heterozygosity"].astype(float).to_numpy()
he = het["mean_expected_heterozygosity"].astype(float).to_numpy()

ax.barh(
    y - height / 2, ho,
    height=height,
    label="Observed heterozygosity (Ho)",
    color="#3174A1",
)
ax.barh(
    y + height / 2, he,
    height=height,
    label="Expected heterozygosity (He)",
    color="#D9A441",
)

ax.set_yticks(y, labels)
ax.invert_yaxis()
ax.set_xlabel("Heterozygosity")
ax.set_xlim(0, max(ho.max(), he.max()) * 1.15)
ax.set_title(
    "Observed and expected heterozygosity by breed",
    weight="bold",
    pad=15,
)
ax.legend(
    frameon=False,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=2,
)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="x", alpha=0.2)
ax.set_axisbelow(True)

fig.text(
    0.5, 0.015,
    "* Small sample (n < 5). "
    "Estimates describe ascertained pathway SNP records only.",
    ha="center",
    fontsize=8,
)

fig.tight_layout(rect=(0, 0.05, 1, 1))

for ext in ("png", "pdf", "svg"):
    fig.savefig(
        BASE / f"H14_figure_S2_breed_heterozygosity.{ext}",
        bbox_inches="tight",
    )
plt.close(fig)

# Supplementary summary table
table = het[[
    "breed",
    "n_samples",
    "n_called_genotypes",
    "n_missing_genotypes",
    "call_rate",
    "observed_heterozygosity",
    "mean_expected_heterozygosity",
    "interpretation_caution",
]].copy()

for category in categories:
    table[f"maf_{category.lower().replace('-', '_')}_count"] = (
        counts.loc[breeds, category].to_numpy()
    )
    table[f"maf_{category.lower().replace('-', '_')}_percent"] = (
        percentages.loc[breeds, category].to_numpy()
    )

table.to_csv(
    BASE / "H14_supplementary_breed_summary.tsv",
    sep="\t",
    index=False,
    float_format="%.6f",
)

print("H14 supplementary outputs generated successfully.")
print("Breeds:", len(breeds))
print("SNP records per breed:", 4505)
print("Figures: 2 (PNG, PDF, SVG)")
print("Supplementary tables: 1")
print("Small-sample breeds flagged: Mertolenga, Mirandesa")
print()
for path in sorted(BASE.glob("H14_figure_S*")):
    print(path.name, f"({path.stat().st_size:,} bytes)")
print("H14_supplementary_breed_summary.tsv")
