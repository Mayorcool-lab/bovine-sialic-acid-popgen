#!/usr/bin/env python3

import csv
import json
from pathlib import Path

manifest = Path(
    "results/H6/H6.2/processed/H6.2_single_accession_uniprot_targets.tsv"
)
rawdir = Path("results/H6/H6.2/raw/interpro")
outfile = Path(
    "results/H6/H6.2/processed/H6.2_interpro_intervals.tsv"
)

with manifest.open() as fh:
    targets = list(csv.DictReader(fh, delimiter="\t"))

fields = [
    "SYMBOL",
    "protein_accession",
    "protein_length",
    "source_database",
    "entry_accession",
    "entry_name",
    "entry_type",
    "interval_start",
    "interval_end",
    "dc_status",
    "representative",
    "model",
    "score",
    "subfamily_accession",
    "subfamily_name",
]

rows = []

for t in targets:
    gene = t["SYMBOL"]
    acc = t["external_accession"]

    path = rawdir / f"{gene}__{acc}.json"

    with path.open() as fh:
        data = json.load(fh)

    for result in data.get("results", []):
        meta = result.get("metadata", {}) or {}

        for protein in result.get("proteins", []):
            returned_acc = protein.get("accession", "")
            protein_length = protein.get("protein_length", "")

            if returned_acc.lower() != acc.lower():
                raise RuntimeError(
                    f"{gene}: expected {acc}, found {returned_acc}"
                )

            for loc in protein.get("entry_protein_locations", []):
                subfamily = loc.get("subfamily", {}) or {}

                for frag in loc.get("fragments", []):
                    rows.append({
                        "SYMBOL": gene,
                        "protein_accession": acc,
                        "protein_length": protein_length,
                        "source_database": meta.get("source_database", ""),
                        "entry_accession": meta.get("accession", ""),
                        "entry_name": meta.get("name", ""),
                        "entry_type": meta.get("type", ""),
                        "interval_start": frag.get("start", ""),
                        "interval_end": frag.get("end", ""),
                        "dc_status": frag.get("dc-status", ""),
                        "representative": loc.get("representative", ""),
                        "model": loc.get("model", ""),
                        "score": loc.get("score", ""),
                        "subfamily_accession": subfamily.get("accession", ""),
                        "subfamily_name": subfamily.get("name", ""),
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
print(f"Interval rows: {len(rows)}")
