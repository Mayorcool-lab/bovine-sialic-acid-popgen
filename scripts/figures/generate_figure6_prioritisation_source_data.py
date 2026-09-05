#!/usr/bin/env python3

import csv
from pathlib import Path


INPUT = Path(
    "results/H_integration/H_candidate_final_triage_62_FINAL.tsv"
)

OUTDIR = Path("results/figures/source_data")
OUTDIR.mkdir(parents=True, exist_ok=True)

FIGURE6_OUT = OUTDIR / "Figure6_integrated_prioritisation.tsv"
DETAILED_OUT = OUTDIR / "Prioritised_candidates_detailed.tsv"


HIGH_PRIORITY_ORDER = [
    ("NEU3", "F159C"),
    ("ST8SIA4", "I73L"),
    ("ST8SIA2", "R129Q"),
    ("NEU3", "E94K"),
    ("CMAH", "N242S"),
    ("ST3GAL4", "R203S"),
    ("CMAH", "I388M"),
]

SECONDARY_PRIORITY_ORDER = [
    ("CMAH", "P424L"),
    ("CMAH", "F512Y"),
    ("ST8SIA2", "S63L"),
    ("ST8SIA4", "R203Q"),
    ("NEU1", "L275P"),
    ("CMAS", "E128*"),
    ("ST3GAL4", "*396L"),
    ("NANS", "P302A"),
]


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clean(value):
    if value is None:
        return ""
    return value.strip()


def priority_publication_label(status):
    if status == "PRIMARY_SHORTLIST":
        return "HIGH_PRIORITY"
    if status == "RESERVE":
        return "SECONDARY_PRIORITY"
    raise ValueError(f"Unexpected shortlist status: {status}")


def consequence_label(value):
    mapping = {
        "missense_variant": "Missense",
        "missense_variant&splice_region_variant": "Missense + splice region",
        "stop_gained": "Stop gained",
        "stop_lost": "Stop lost",
    }
    return mapping.get(value, value)


def sift_label(row):
    prediction = clean(row["SIFT_prediction"])
    score = clean(row["SIFT_score"])

    if not prediction:
        if row["Consequence"] in {"stop_gained", "stop_lost"}:
            return "Not applicable"
        return "Unavailable"

    label_map = {
        "deleterious": "Deleterious",
        "deleterious_low_confidence": "Deleterious, low confidence",
        "tolerated": "Tolerated",
        "tolerated_low_confidence": "Tolerated, low confidence",
    }

    label = label_map.get(prediction, prediction)

    if score:
        return f"{label} ({score})"

    return label


def rate4site_label(row):
    score = clean(row["H5_rate4site_score"])

    if score:
        return score

    consequence = clean(row["Consequence"])
    coverage = clean(row["H5_coverage_status"])

    if consequence == "stop_lost":
        return "Special case"

    if coverage:
        return "Not assessable"

    return "Not assessable"


def h6_label(value):
    mapping = {
        "RESIDUE_LEVEL_DIRECT_OR_MAPPED": "Residue level",
        "RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED":
            "Residue level, multi-reference",
        "RESIDUE_LEVEL_SEQUENCE_BASED":
            "Residue level, sequence-based",
        "PROTEIN_LEVEL_ONLY": "Protein level only",
    }
    return mapping.get(value, value)


def cohort_label(status, af):
    status = clean(status)
    af = clean(af)

    if status == "PRESENT":
        if af == "0.0000":
            return "Interrogated, ALT not observed"
        if af:
            return f"ALT observed, AF={af}"
        return "ALT observed"

    if status == "NOT_IN_SOURCE":
        return "Not in source"

    if status == "ABSENT_AFTER_QC":
        if af:
            return f"QC-limited, {af}"
        return "QC-limited"

    return status


def breed_label(n, breeds):
    n = clean(n)
    breeds = clean(breeds)

    try:
        n_int = int(n)
    except (TypeError, ValueError):
        n_int = 0

    if n_int == 0:
        return "0"

    if breeds:
        return f"{n_int}: {breeds.replace('_', ' ')}"

    return str(n_int)


def get_key(row):
    return (
        clean(row["SYMBOL"]),
        clean(row["protein_change_standardized"]),
    )


rows = read_rows(INPUT)

assert len(rows) == 62, f"Expected 62 candidates, found {len(rows)}"

status_counts = {}
for row in rows:
    status = row["shortlist_status"]
    status_counts[status] = status_counts.get(status, 0) + 1

assert status_counts == {
    "PRIMARY_SHORTLIST": 7,
    "RESERVE": 8,
    "NOT_SHORTLISTED": 47,
}, f"Unexpected shortlist counts: {status_counts}"

row_lookup = {get_key(row): row for row in rows}

expected_keys = HIGH_PRIORITY_ORDER + SECONDARY_PRIORITY_ORDER

missing = [key for key in expected_keys if key not in row_lookup]
assert not missing, f"Missing prioritised candidates: {missing}"

selected = []

for rank, key in enumerate(HIGH_PRIORITY_ORDER, start=1):
    row = row_lookup[key]
    assert row["shortlist_status"] == "PRIMARY_SHORTLIST", (
        f"{key} expected PRIMARY_SHORTLIST but found "
        f"{row['shortlist_status']}"
    )
    selected.append((rank, row))

offset = len(HIGH_PRIORITY_ORDER)

for rank, key in enumerate(
    SECONDARY_PRIORITY_ORDER,
    start=offset + 1,
):
    row = row_lookup[key]
    assert row["shortlist_status"] == "RESERVE", (
        f"{key} expected RESERVE but found "
        f"{row['shortlist_status']}"
    )
    selected.append((rank, row))

assert len(selected) == 15


figure6_fields = [
    "display_order",
    "candidate",
    "SYMBOL",
    "protein_change",
    "priority",
    "Consequence",
    "Consequence_display",
    "SIFT_prediction",
    "SIFT_score",
    "SIFT_display",
    "Rate4Site",
    "Rate4Site_display",
    "H5_coverage_status",
    "H6_context_level",
    "H6_context_display",
    "SA_status",
    "SA_AF",
    "SA_display",
    "SA_observed_breed_n",
    "SA_observed_breeds",
    "EU_status",
    "EU_AF",
    "EU_display",
    "EU_observed_breed_n",
    "EU_observed_breeds",
]

figure6_rows = []

for display_order, row in selected:
    gene = clean(row["SYMBOL"])
    change = clean(row["protein_change_standardized"])

    figure6_rows.append({
        "display_order": display_order,
        "candidate": f"{gene} {change}",
        "SYMBOL": gene,
        "protein_change": change,
        "priority": priority_publication_label(
            row["shortlist_status"]
        ),
        "Consequence": clean(row["Consequence"]),
        "Consequence_display": consequence_label(
            clean(row["Consequence"])
        ),
        "SIFT_prediction": clean(row["SIFT_prediction"]),
        "SIFT_score": clean(row["SIFT_score"]),
        "SIFT_display": sift_label(row),
        "Rate4Site": clean(row["H5_rate4site_score"]),
        "Rate4Site_display": rate4site_label(row),
        "H5_coverage_status": clean(row["H5_coverage_status"]),
        "H6_context_level": clean(row["H6_context_level"]),
        "H6_context_display": h6_label(
            clean(row["H6_context_level"])
        ),
        "SA_status": clean(row["SA_status"]),
        "SA_AF": clean(row["SA_AF"]),
        "SA_display": cohort_label(
            row["SA_status"],
            row["SA_AF"],
        ),
        "SA_observed_breed_n":
            clean(row["SA_observed_OK_n"]),
        "SA_observed_breeds":
            clean(row["SA_observed_OK_breeds"]).replace("_", " "),
        "EU_status": clean(row["EU_status"]),
        "EU_AF": clean(row["EU_AF"]),
        "EU_display": cohort_label(
            row["EU_status"],
            row["EU_AF"],
        ),
        "EU_observed_breed_n":
            clean(row["EU_observed_OK_n"]),
        "EU_observed_breeds":
            clean(row["EU_observed_OK_breeds"]).replace("_", " "),
    })


with FIGURE6_OUT.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=figure6_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(figure6_rows)


detailed_fields = [
    "candidate",
    "Consequence",
    "SIFT",
    "Rate4Site",
    "Protein_context",
    "South_African_evidence",
    "European_evidence",
    "SA_breeds_with_ALT",
    "EU_breeds_with_ALT",
    "Priority",
    "Biological_rationale",
    "Main_limitation",
]

detailed_rows = []

for _, row in selected:
    gene = clean(row["SYMBOL"])
    change = clean(row["protein_change_standardized"])

    detailed_rows.append({
        "candidate": f"{gene} {change}",
        "Consequence": consequence_label(
            clean(row["Consequence"])
        ),
        "SIFT": sift_label(row),
        "Rate4Site": rate4site_label(row),
        "Protein_context": h6_label(
            clean(row["H6_context_level"])
        ),
        "South_African_evidence": cohort_label(
            row["SA_status"],
            row["SA_AF"],
        ),
        "European_evidence": cohort_label(
            row["EU_status"],
            row["EU_AF"],
        ),
        "SA_breeds_with_ALT": breed_label(
            row["SA_observed_OK_n"],
            row["SA_observed_OK_breeds"],
        ),
        "EU_breeds_with_ALT": breed_label(
            row["EU_observed_OK_n"],
            row["EU_observed_OK_breeds"],
        ),
        "Priority": (
            "High priority"
            if row["shortlist_status"] == "PRIMARY_SHORTLIST"
            else "Secondary priority"
        ),
        "Biological_rationale":
            clean(row["biological_rationale"]),
        "Main_limitation":
            clean(row["main_limitation"]),
    })


with DETAILED_OUT.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=detailed_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(detailed_rows)


print(
    "Figure 6 and detailed prioritised-candidate source data "
    "generated successfully."
)
print()
print(f"Authoritative candidates: {len(rows)}")
print(
    "Frozen split: "
    f"{status_counts['PRIMARY_SHORTLIST']} primary, "
    f"{status_counts['RESERVE']} reserve, "
    f"{status_counts['NOT_SHORTLISTED']} not shortlisted"
)
print()
print("Publication-facing prioritised set:")
print("  High priority: 7")
print("  Secondary priority: 8")
print()
print(f"Figure 6 source rows: {len(figure6_rows)}")
print(f"Detailed prioritised-candidate rows: {len(detailed_rows)}")
print()
print(f"Wrote: {FIGURE6_OUT}")
print(f"Wrote: {DETAILED_OUT}")
