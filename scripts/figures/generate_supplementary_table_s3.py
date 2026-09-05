#!/usr/bin/env python3

import csv
from pathlib import Path


INPUT = Path(
    "results/figures/source_data/Prioritised_candidates_detailed.tsv"
)

OUTPUT = Path(
    "results/figures/source_data/"
    "Supplementary_Table_S3_prioritised_candidates.tsv"
)


EXPECTED_ORDER = [
    "NEU3 F159C",
    "ST8SIA4 I73L",
    "ST8SIA2 R129Q",
    "NEU3 E94K",
    "CMAH N242S",
    "ST3GAL4 R203S",
    "CMAH I388M",
    "CMAH P424L",
    "CMAH F512Y",
    "ST8SIA2 S63L",
    "ST8SIA4 R203Q",
    "NEU1 L275P",
    "CMAS E128*",
    "ST3GAL4 *396L",
    "NANS P302A",
]


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clean(value):
    if value is None:
        return ""
    return value.strip()


def population_summary(row):
    sa = clean(row["South_African_evidence"])
    eu = clean(row["European_evidence"])

    return f"SA: {sa}; EU: {eu}"


def key_interpretation(row):
    rationale = clean(row["Biological_rationale"])
    limitation = clean(row["Main_limitation"])

    if rationale and limitation:
        return f"{rationale} Limitation: {limitation}"

    if rationale:
        return rationale

    if limitation:
        return f"Limitation: {limitation}"

    return ""


rows = read_rows(INPUT)

assert len(rows) == 15, (
    f"Expected 15 prioritised candidates, found {len(rows)}"
)

observed_order = [clean(row["candidate"]) for row in rows]

assert observed_order == EXPECTED_ORDER, (
    "Unexpected candidate order.\n"
    f"Expected: {EXPECTED_ORDER}\n"
    f"Observed: {observed_order}"
)

priority_counts = {}
for row in rows:
    priority = clean(row["Priority"])
    priority_counts[priority] = priority_counts.get(priority, 0) + 1

assert priority_counts == {
    "High priority": 7,
    "Secondary priority": 8,
}, f"Unexpected priority counts: {priority_counts}"


output_fields = [
    "Candidate",
    "Consequence",
    "SIFT",
    "Rate4Site",
    "Population_evidence",
    "Priority",
    "Key_interpretation",
]

output_rows = []

for row in rows:
    output_rows.append({
        "Candidate": clean(row["candidate"]),
        "Consequence": clean(row["Consequence"]),
        "SIFT": clean(row["SIFT"]),
        "Rate4Site": clean(row["Rate4Site"]),
        "Population_evidence": population_summary(row),
        "Priority": clean(row["Priority"]),
        "Key_interpretation": key_interpretation(row),
    })


OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=output_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(output_rows)


print("Supplementary Table S3 manuscript source generated successfully.")
print()
print(f"Rows: {len(output_rows)}")
print(f"High priority: {priority_counts['High priority']}")
print(
    "Secondary priority: "
    f"{priority_counts['Secondary priority']}"
)
print()
print(f"Wrote: {OUTPUT}")
