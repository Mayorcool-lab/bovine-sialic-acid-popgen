#!/usr/bin/env python3

import csv
import json
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

manifest = Path(
    "results/H6/H6.2/processed/H6.2_single_accession_uniprot_targets.tsv"
)
outdir = Path("results/H6/H6.2/raw/interpro")
logfile = Path("results/H6/H6.2/qc/H6.2_interpro_retrieval_log.tsv")

outdir.mkdir(parents=True, exist_ok=True)
logfile.parent.mkdir(parents=True, exist_ok=True)

rows = []

with manifest.open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gene = r["SYMBOL"]
        accession = r["external_accession"]

        url = (
            "https://www.ebi.ac.uk/interpro/api/entry/all/"
            f"protein/uniprot/{accession}/"
        )

        outfile = outdir / f"{gene}__{accession}.json"

        req = urllib.request.Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "DFG-H6.2-functional-annotation"
            }
        )

        with urllib.request.urlopen(req) as response:
            data = json.load(response)

        with outfile.open("w") as out:
            json.dump(data, out, indent=2)
            out.write("\n")

        results = data.get("results", [])
        reported_count = data.get("count", "")

        rows.append({
            "SYMBOL": gene,
            "accession": accession,
            "reported_count": reported_count,
            "results_in_response": len(results),
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "raw_json": str(outfile),
        })

with logfile.open("w") as fh:
    fields = [
        "SYMBOL",
        "accession",
        "reported_count",
        "results_in_response",
        "retrieved_utc",
        "raw_json",
    ]

    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Prepared retrieval for {len(rows)} InterPro records")
print(f"Log: {logfile}")
