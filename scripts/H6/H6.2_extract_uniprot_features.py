#!/usr/bin/env python3

import csv
import json
from pathlib import Path

manifest = Path(
    "results/H6/H6.2/processed/H6.2_single_accession_uniprot_targets.tsv"
)
rawdir = Path("results/H6/H6.2/raw/uniprot")
outfile = Path(
    "results/H6/H6.2/processed/H6.2_uniprot_features.tsv"
)

targets = []

with manifest.open() as fh:
    targets = list(csv.DictReader(fh, delimiter="\t"))

fields = [
    "SYMBOL",
    "uniprot_accession",
    "uniprot_entry_type",
    "feature_type",
    "feature_start",
    "feature_end",
    "feature_description",
    "evidence_code",
    "evidence_source",
    "evidence_source_accession",
]

rows = []

for t in targets:
    gene = t["SYMBOL"]
    acc = t["external_accession"]

    path = rawdir / f"{gene}__{acc}.json"

    with path.open() as fh:
        data = json.load(fh)

    entry_type = data.get("entryType", "")

    for feature in data.get("features", []):
        loc = feature.get("location", {})

        start = (
            loc.get("start", {}).get("value", "")
            if isinstance(loc.get("start"), dict)
            else ""
        )

        end = (
            loc.get("end", {}).get("value", "")
            if isinstance(loc.get("end"), dict)
            else ""
        )

        evidences = feature.get("evidences", [])

        if not evidences:
            rows.append({
                "SYMBOL": gene,
                "uniprot_accession": acc,
                "uniprot_entry_type": entry_type,
                "feature_type": feature.get("type", ""),
                "feature_start": start,
                "feature_end": end,
                "feature_description": feature.get("description", ""),
                "evidence_code": "",
                "evidence_source": "",
                "evidence_source_accession": "",
            })
            continue

        for ev in evidences:
            src = ev.get("source") or ""
            src_id = ev.get("id") or ""

            rows.append({
                "SYMBOL": gene,
                "uniprot_accession": acc,
                "uniprot_entry_type": entry_type,
                "feature_type": feature.get("type", ""),
                "feature_start": start,
                "feature_end": end,
                "feature_description": feature.get("description", ""),
                "evidence_code": ev.get("evidenceCode", ""),
                "evidence_source": src,
                "evidence_source_accession": src_id,
            })

with outfile.open("w") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {outfile}")
print(f"Feature-evidence rows: {len(rows)}")
