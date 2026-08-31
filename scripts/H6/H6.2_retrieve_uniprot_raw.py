#!/usr/bin/env python3

import csv
import json
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

manifest = Path(
    "results/H6/H6.2/processed/H6.2_single_accession_uniprot_targets.tsv"
)
outdir = Path("results/H6/H6.2/raw/uniprot")
logfile = Path("results/H6/H6.2/qc/H6.2_uniprot_retrieval_log.tsv")

outdir.mkdir(parents=True, exist_ok=True)
logfile.parent.mkdir(parents=True, exist_ok=True)

rows = []

with manifest.open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gene = r["SYMBOL"]
        accession = r["external_accession"]

        url = f"https://rest.uniprot.org/uniprotkb/{accession}"
        outfile = outdir / f"{gene}__{accession}.json"

        req = urllib.request.Request(
            url,
            headers={"Accept": "application/json"}
        )

        with urllib.request.urlopen(req) as response:
            data = json.load(response)

        with outfile.open("w") as out:
            json.dump(data, out, indent=2)
            out.write("\n")

        returned_accession = data.get("primaryAccession", "")
        sequence_length = data.get("sequence", {}).get("length", "")
        entry_type = data.get("entryType", "")
        entry_version = data.get("entryAudit", {}).get("entryVersion", "")
        annotation_date = data.get("entryAudit", {}).get(
            "lastAnnotationUpdateDate", ""
        )

        if returned_accession != accession:
            raise RuntimeError(
                f"{gene}: requested {accession}, "
                f"received {returned_accession}"
            )

        rows.append({
            "SYMBOL": gene,
            "requested_accession": accession,
            "returned_accession": returned_accession,
            "entry_type": entry_type,
            "sequence_length": sequence_length,
            "entry_version": entry_version,
            "last_annotation_update": annotation_date,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "raw_json": str(outfile),
        })

with logfile.open("w") as fh:
    fields = [
        "SYMBOL",
        "requested_accession",
        "returned_accession",
        "entry_type",
        "sequence_length",
        "entry_version",
        "last_annotation_update",
        "retrieved_utc",
        "raw_json",
    ]

    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Prepared retrieval for {len(rows)} UniProt records")
print(f"Log: {logfile}")
