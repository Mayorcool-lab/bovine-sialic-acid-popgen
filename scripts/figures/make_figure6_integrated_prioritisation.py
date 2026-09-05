#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


INPUT = Path(
    "results/figures/source_data/Figure6_integrated_prioritisation.tsv"
)

PNG_OUT = Path(
    "results/figures/Figure6_integrated_prioritisation.png"
)

SVG_OUT = Path(
    "results/figures/Figure6_integrated_prioritisation.svg"
)


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    rows.sort(key=lambda r: int(r["display_order"]))
    return rows


def candidate_label(row):
    return f"{row['SYMBOL']} {row['protein_change']}"


def consequence_short(value):
    mapping = {
        "Missense": "Missense",
        "Missense + splice region": "Missense\n+ splice",
        "Stop gained": "Stop\ngained",
        "Stop lost": "Stop\nlost",
    }
    return mapping.get(value, value)


def sift_short(row):
    prediction = row["SIFT_prediction"].strip()
    score = row["SIFT_score"].strip()

    mapping = {
        "deleterious": "Deleterious",
        "deleterious_low_confidence": "Deleterious\nlow confidence",
        "tolerated": "Tolerated",
        "tolerated_low_confidence": "Tolerated\nlow confidence",
    }

    if not prediction:
        if row["Consequence"] in {"stop_gained", "stop_lost"}:
            return "Not\napplicable"
        return "Unavailable"

    label = mapping.get(prediction, prediction)

    if score:
        return f"{label}\n{score}"

    return label


def h6_short(value):
    mapping = {
        "RESIDUE_LEVEL_DIRECT_OR_MAPPED":
            "Residue\nlevel",
        "RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED":
            "Residue level\nmulti-reference",
        "RESIDUE_LEVEL_SEQUENCE_BASED":
            "Residue level\nsequence-based",
        "PROTEIN_LEVEL_ONLY":
            "Protein\nlevel only",
    }
    return mapping.get(value, value)


def population_short(status, af, breed_n):
    status = status.strip()
    af = af.strip()
    breed_n = breed_n.strip()

    try:
        breed_n_int = int(breed_n)
    except (TypeError, ValueError):
        breed_n_int = 0

    if status == "PRESENT":
        if af == "0.0000":
            return "ALT not observed\n0 breeds"
        return f"AF {af}\n{breed_n_int} breed" + (
            "" if breed_n_int == 1 else "s"
        )

    if status == "NOT_IN_SOURCE":
        return "Not in\nsource"

    if status == "ABSENT_AFTER_QC":
        af_value = af.split()[0] if af else ""
        if af_value:
            return f"QC-limited\nAF {af_value}"
        return "QC-limited"

    return status.replace("_", "\n")


def rate4site_text(row):
    value = row["Rate4Site"].strip()

    if value:
        return value

    if row["Consequence"] == "stop_lost":
        return "Special\ncase"

    return "Not\nassessable"


def priority_short(value):
    if value == "HIGH_PRIORITY":
        return "High\npriority"
    if value == "SECONDARY_PRIORITY":
        return "Secondary\npriority"
    return value


rows = read_rows(INPUT)

assert len(rows) == 15, f"Expected 15 rows, found {len(rows)}"

priority_counts = {}
for row in rows:
    p = row["priority"]
    priority_counts[p] = priority_counts.get(p, 0) + 1

assert priority_counts == {
    "HIGH_PRIORITY": 7,
    "SECONDARY_PRIORITY": 8,
}, f"Unexpected priority split: {priority_counts}"


columns = [
    "Candidate",
    "Consequence",
    "SIFT",
    "Rate4Site",
    "Protein context",
    "South Africa",
    "Europe",
    "Priority",
]

column_widths = [
    1.80,
    1.70,
    2.00,
    1.50,
    2.35,
    1.95,
    1.95,
    1.75,
]

x_positions = [0]
for width in column_widths:
    x_positions.append(x_positions[-1] + width)

total_width = x_positions[-1]

row_height = 0.92
header_height = 1.10
n_rows = len(rows)
total_height = header_height + n_rows * row_height

fig, ax = plt.subplots(figsize=(18.5, 12.5))

ax.set_xlim(0, total_width)
ax.set_ylim(0, total_height)
ax.axis("off")

header_face = "#E6E6E6"
grid_color = "#B8B8B8"

high_priority_face = "#D9EAD3"
secondary_priority_face = "#FFF2CC"

deleterious_face = "#F4CCCC"
tolerated_face = "#D9EAD3"
missing_face = "#E7E6E6"

residue_face = "#D9EAF7"
protein_only_face = "#FCE5CD"

present_face = "#D9EAD3"
not_source_face = "#E7E6E6"
qc_face = "#FCE5CD"
zero_alt_face = "#FFF2CC"

stop_face = "#F4CCCC"
missense_face = "#EAF2F8"

rate_values = [
    float(row["Rate4Site"])
    for row in rows
    if row["Rate4Site"].strip()
]

rate_min = min(rate_values)
rate_max = max(rate_values)

rate_norm = Normalize(
    vmin=rate_min,
    vmax=rate_max,
)

rate_cmap = plt.get_cmap("coolwarm")

for col_idx, title in enumerate(columns):
    x0 = x_positions[col_idx]
    width = column_widths[col_idx]

    ax.add_patch(
        Rectangle(
            (x0, total_height - header_height),
            width,
            header_height,
            facecolor=header_face,
            edgecolor=grid_color,
            linewidth=0.8,
        )
    )

    ax.text(
        x0 + width / 2,
        total_height - header_height / 2,
        title,
        ha="center",
        va="center",
        fontsize=13.5,
        fontweight="bold",
    )

for row_idx, row in enumerate(rows):
    y0 = total_height - header_height - (row_idx + 1) * row_height

    is_high = row["priority"] == "HIGH_PRIORITY"

    candidate_bg = (
        high_priority_face
        if is_high
        else secondary_priority_face
    )

    cell_texts = [
        candidate_label(row),
        consequence_short(row["Consequence_display"]),
        sift_short(row),
        rate4site_text(row),
        h6_short(row["H6_context_level"]),
        population_short(
            row["SA_status"],
            row["SA_AF"],
            row["SA_observed_breed_n"],
        ),
        population_short(
            row["EU_status"],
            row["EU_AF"],
            row["EU_observed_breed_n"],
        ),
        priority_short(row["priority"]),
    ]

    cell_faces = [candidate_bg]

    consequence = row["Consequence"]
    if consequence in {"stop_gained", "stop_lost"}:
        cell_faces.append(stop_face)
    else:
        cell_faces.append(missense_face)

    sift = row["SIFT_prediction"].strip()

    if sift.startswith("deleterious"):
        cell_faces.append(deleterious_face)
    elif sift.startswith("tolerated"):
        cell_faces.append(tolerated_face)
    else:
        cell_faces.append(missing_face)

    rate = row["Rate4Site"].strip()

    if rate:
        cell_faces.append(
            rate_cmap(rate_norm(float(rate)))
        )
    else:
        cell_faces.append(missing_face)

    context = row["H6_context_level"]

    if context == "PROTEIN_LEVEL_ONLY":
        cell_faces.append(protein_only_face)
    else:
        cell_faces.append(residue_face)

    for prefix in ("SA", "EU"):
        status = row[f"{prefix}_status"].strip()
        af = row[f"{prefix}_AF"].strip()

        if status == "NOT_IN_SOURCE":
            cell_faces.append(not_source_face)
        elif status == "ABSENT_AFTER_QC":
            cell_faces.append(qc_face)
        elif status == "PRESENT" and af == "0.0000":
            cell_faces.append(zero_alt_face)
        else:
            cell_faces.append(present_face)

    cell_faces.append(
        high_priority_face
        if is_high
        else secondary_priority_face
    )

    for col_idx, text in enumerate(cell_texts):
        x0 = x_positions[col_idx]
        width = column_widths[col_idx]

        ax.add_patch(
            Rectangle(
                (x0, y0),
                width,
                row_height,
                facecolor=cell_faces[col_idx],
                edgecolor=grid_color,
                linewidth=0.65,
            )
        )

        fontsize = 11.2

        if col_idx == 0:
            fontsize = 12.1

        if col_idx in {2, 4, 5, 6}:
            fontsize = 10.6

        ax.text(
            x0 + width / 2,
            y0 + row_height / 2,
            text,
            ha="center",
            va="center",
            fontsize=fontsize,
            fontweight="bold" if col_idx in {0, 7} else "normal",
            linespacing=1.15,
        )

separator_y = (
    total_height
    - header_height
    - 7 * row_height
)

ax.plot(
    [0, total_width],
    [separator_y, separator_y],
    linewidth=2.2,
    color="#555555",
)

ax.text(
    -0.10,
    total_height
    - header_height
    - (7 * row_height) / 2,
    "HIGH PRIORITY",
    rotation=90,
    ha="center",
    va="center",
    fontsize=11.5,
    fontweight="bold",
)

ax.text(
    -0.10,
    separator_y
    - (8 * row_height) / 2,
    "SECONDARY PRIORITY",
    rotation=90,
    ha="center",
    va="center",
    fontsize=11.5,
    fontweight="bold",
)

sm = ScalarMappable(
    norm=rate_norm,
    cmap=rate_cmap,
)
sm.set_array([])

cbar = fig.colorbar(
    sm,
    ax=ax,
    orientation="horizontal",
    fraction=0.024,
    pad=0.035,
    aspect=48,
)

cbar.set_label(
    "Rate4Site score  |  lower values indicate stronger evolutionary conservation",
    fontsize=11,
)

cbar.ax.tick_params(labelsize=10)

legend_text = (
    "Population cells report cohort AF and number of adequately sampled breeds "
    "with observed ALT. Grey: site not represented in source VCF. "
    "Orange: QC-limited. Yellow: site interrogated but ALT not observed. "
    "SIFT-unavailable and non-assessable evidence remain explicitly missing."
)

fig.text(
    0.5,
    0.018,
    legend_text,
    ha="center",
    va="bottom",
    fontsize=10.2,
    wrap=True,
)

plt.subplots_adjust(
    left=0.055,
    right=0.985,
    top=0.985,
    bottom=0.105,
)

PNG_OUT.parent.mkdir(parents=True, exist_ok=True)

fig.savefig(
    PNG_OUT,
    dpi=600,
    bbox_inches="tight",
)

fig.savefig(
    SVG_OUT,
    bbox_inches="tight",
)

plt.close(fig)

print("Figure 6 generated successfully.")
print()
print(f"Rows: {len(rows)}")
print("High priority: 7")
print("Secondary priority: 8")
print()
print(
    f"Rate4Site plotted range: "
    f"{rate_min:.4f} to {rate_max:.4f}"
)
print()
print(f"Wrote: {PNG_OUT}")
print(f"Wrote: {SVG_OUT}")
