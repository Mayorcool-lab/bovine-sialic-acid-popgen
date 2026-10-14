#!/usr/bin/env python3
"""Figure 8: validated pathway-level SNP diversity across cattle breeds."""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
H10 = ROOT / "results/H_integration/H10_gene_pi/gene_pi_validation.tsv"
ORIGINAL = ROOT / "data/external/aiyesa_popgen/pi_res_breed.txt"
OUT = ROOT / "results/H_integration/H11_figure8"
OUT.mkdir(parents=True, exist_ok=True)

# H10 supplies validated diversity estimates.
# The original file supplies breed-to-country labels only.
pi = pd.read_csv(H10, sep="\t")
original = pd.read_csv(ORIGINAL, sep="\t")

countries = (
    original.loc[original["Group"] == "Genes", ["Breed", "Country"]]
    .drop_duplicates()
)

assert len(pi) == 16
assert pi["breed"].is_unique
assert countries["Breed"].is_unique

df = pi.merge(
    countries,
    left_on="breed",
    right_on="Breed",
    how="left",
    validate="one_to_one",
)

assert df["Country"].notna().all()
assert df["n_samples"].sum() == 346
assert df["mean_pi_recomputed"].notna().all()

# Holstein-Friesian includes animals sampled in three countries.
# Change its plotting group only; do not modify its H10 estimate.
hf_mask = (
    df["breed"]
    .astype(str)
    .str.replace(r"[\s_-]+", "", regex=True)
    .str.casefold()
    .eq("holsteinfriesian")
)

assert hf_mask.sum() == 1, (
    "Expected exactly one Holstein-Friesian breed row."
)
assert int(df.loc[hf_mask, "n_samples"].iloc[0]) == 15

original_pi = df["mean_pi_recomputed"].copy()
original_n = df["n_samples"].copy()

df.loc[hf_mask, "Country"] = "Multi-country"

assert df["mean_pi_recomputed"].equals(original_pi)
assert df["n_samples"].equals(original_n)
assert len(df) == 16

country_order = [
    "Finland",
    "Netherlands",
    "Portugal",
    "Multi-country",
    "South Africa",
]

assert set(df["Country"]) == set(country_order), (
    f"Unexpected country categories: {sorted(df['Country'].unique())}"
)

df = df.sort_values(["Country", "mean_pi_recomputed", "breed"])
df.to_csv(OUT / "figure8_plot_data.tsv", sep="\t", index=False)

positions = {country: i for i, country in enumerate(country_order)}
palette = plt.get_cmap("Set2")
colors = {
    country: palette(i)
    for i, country in enumerate(country_order)
}

fig, ax = plt.subplots(figsize=(10.5, 6.0))

for country in country_order:
    group = df.loc[df["Country"] == country]
    x = positions[country]

    # A one-breed category is shown as a point, not a boxplot.
    if len(group) > 1:
        ax.boxplot(
            group["mean_pi_recomputed"],
            positions=[x],
            widths=0.42,
            patch_artist=True,
            showfliers=False,
            boxprops={
                "facecolor": colors[country],
                "alpha": 0.35,
                "edgecolor": "#666666",
            },
            medianprops={"color": "#222222", "linewidth": 1.5},
            whiskerprops={"color": "#666666"},
            capprops={"color": "#666666"},
        )

    group = group.reset_index(drop=True)
    count = len(group)

    for j, row in group.iterrows():
        offset = (j - (count - 1) / 2) * 0.12
        px = x + offset
        py = row["mean_pi_recomputed"]
        small = row["n_samples"] < 5

        ax.scatter(
            px,
            py,
            s=72,
            marker="D" if small else "o",
            facecolor="white" if small else colors[country],
            edgecolor="#222222",
            linewidth=1.1,
            zorder=4,
        )

        if small:
            name = row["breed"].replace("_", " ")
            ax.annotate(
                f"{name} (n={row['n_samples']})",
                xy=(px, py),
                xytext=(10, 12 if name == "Mertolenga" else -15),
                textcoords="offset points",
                fontsize=8.5,
                ha="left",
                va="center",
            )

ax.set_xticks(range(len(country_order)))
ax.set_xticklabels(country_order, fontsize=10)

ax.set_ylabel(
    "Mean nucleotide diversity per ascertained SNP (π)",
    fontsize=11,
)
ax.set_xlabel("Country / sampling group", fontsize=11)

ax.set_ylim(0.145, 0.205)
ax.set_xlim(-0.55, len(country_order) - 0.35)
ax.grid(axis="y", alpha=0.20)
ax.set_axisbelow(True)

legend_handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        markerfacecolor="#80b1d3",
        markeredgecolor="#222222",
        markersize=8,
        label="Breed (n ≥ 5)",
    ),
    Line2D(
        [0], [0],
        marker="D",
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="#222222",
        markersize=8,
        label="Small-sample breed (n < 5)",
    ),
]

ax.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.02),
    ncol=2,
    frameon=False,
    fontsize=9,
)

fig.tight_layout(rect=[0, 0.03, 1, 0.95])

for ext in ("png", "pdf", "svg"):
    fig.savefig(
        OUT / f"figure8_pathway_pi.{ext}",
        dpi=300,
        bbox_inches="tight",
    )

plt.close(fig)

print("Figure 8 regenerated successfully")
print(f"Breeds: {len(df)}")
print(f"Animals: {df['n_samples'].sum()}")
print(f"Sampling groups: {', '.join(country_order)}")
print(f"Holstein-Friesian reassigned: {int(hf_mask.sum())} breed")
print(f"Output directory: {OUT}")
