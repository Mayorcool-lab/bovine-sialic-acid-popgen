#!/usr/bin/env python3

import csv
import json
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

manifest = Path(
    "results/H6/H6.2/processed/H6.2_ST3GAL4_uniprot_targets.tsv"
)
rawdir = Path("results/H6/H6.2/raw/uniprot_ST3GAL4")
logfile = Path(
    "results/H6/H6.2/qc/H6.2_ST3GAL4_uniprot_retrieval_log.tsv"
)

rawdir.mkdir(parents=True, exist_ok=True)

with manifest.open() as fh:
    targets = list(csv.DictReader(fh, delimiter="\t"))

fields = [
    "SYMBOL",
    "requested_accession",
    "returned_accession",
    "entry_type",
    "sequence_length",
    "entry_version",
    "annotation_update_date",
    "retrieved_utc",
    "raw_path",
]

rows = []

for t in targets:
    gene = t["SYMBOL"]
    acc = t["external_accession"]

    url = f"https://rest.uniprot.org/uniprotkb/{acc}"

    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json"}
    )

    with urllib.request.urlopen(req) as response:
        data = json.load(response)

    returned = data.get("primaryAccession", "")

    if returned != acc:
        raise RuntimeError(
            f"{gene}: requested {acc}, returned {returned}"
        )

    rawpath = rawdir / f"{gene}__{acc}.json"

    with rawpath.open("w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")

    entry_type = data.get("entryType", "")
    seq_length = (data.get("sequence") or {}).get("length", "")
    entry_version = data.get("entryAudit", {}).get("entryVersion", "")
    annotation_date = data.get("entryAudit", {}).get(
        "lastAnnotationUpdateDate", ""
    )

    rows.append({
        "SYMBOL": gene,
        "requested_accession": acc,
        "returned_accession": returned,
        "entry_type": entry_type,
        "sequence_length": seq_length,
        "entry_version": entry_version,
        "annotation_update_date": annotation_date,
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
        r["requested_accession"],
        r["entry_type"],
        f'length={r["sequence_length"]}',
        f'version={r["entry_version"]}',
        f'annotation={r["annotation_update_date"]}',
        sep="\t"
    )
