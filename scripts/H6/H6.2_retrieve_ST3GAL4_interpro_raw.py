#!/usr/bin/env python3

import csv
import json
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

manifest = Path(
    "results/H6/H6.2/processed/H6.2_ST3GAL4_uniprot_targets.tsv"
)
rawdir = Path("results/H6/H6.2/raw/interpro_ST3GAL4")
logfile = Path(
    "results/H6/H6.2/qc/H6.2_ST3GAL4_interpro_retrieval_log.tsv"
)

rawdir.mkdir(parents=True, exist_ok=True)

with manifest.open() as fh:
    targets = list(csv.DictReader(fh, delimiter="\t"))

fields = [
    "SYMBOL",
    "accession",
    "reported_count",
    "results_in_response",
    "retrieved_utc",
    "raw_path",
]

rows = []

for t in targets:
    gene = t["SYMBOL"]
    acc = t["external_accession"]

    url = (
        "https://www.ebi.ac.uk/interpro/api/"
        f"entry/all/protein/uniprot/{acc}/"
    )

    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json"}
    )

    with urllib.request.urlopen(req) as response:
        data = json.load(response)

    rawpath = rawdir / f"{gene}__{acc}.json"

    with rawpath.open("w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")

    reported = data.get("count", "")
    observed = len(data.get("results", []) or [])

    rows.append({
        "SYMBOL": gene,
        "accession": acc,
        "reported_count": reported,
        "results_in_response": observed,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "raw_path": str(rawpath),
    })

with logfile.open("w") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Retrieved records: {len(rows)}")
print(f"Wrote: {logfile}")

for r in rows:
    print(
        r["accession"],
        f'reported={r["reported_count"]}',
        f'observed={r["results_in_response"]}',
        sep="\t"
    )
