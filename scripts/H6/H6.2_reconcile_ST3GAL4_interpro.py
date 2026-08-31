#!/usr/bin/env python3

import csv
import json
from pathlib import Path
from collections import defaultdict

targets = [190, 203, 320, 333, 334]

alndir = Path("results/H6/reference_mapping/pairwise_alignments")
iprdir = Path("results/H6/H6.2/raw/interpro_ST3GAL4")

out_detail = Path(
    "results/H6/H6.2/processed/H6.2_ST3GAL4_interpro_reconciliation_detail.tsv"
)
out_summary = Path(
    "results/H6/H6.2/processed/H6.2_ST3GAL4_interpro_reconciliation_summary.tsv"
)

files = sorted(alndir.glob("ST3GAL4_H5_vs_*.fasta"))

# ------------------------------------------------------------------
# Recover validated H5 -> external coordinates from frozen alignments
# ------------------------------------------------------------------

coord_rows = []

for path in files:
    acc = path.stem.replace("ST3GAL4_H5_vs_", "")

    seqs = []
    header = None
    seq = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

    if header is not None:
        seqs.append((header, "".join(seq)))

    if len(seqs) != 2:
        raise RuntimeError(
            f"{path}: expected 2 aligned sequences, found {len(seqs)}"
        )

    h5_aln = seqs[0][1]
    ext_aln = seqs[1][1]

    if len(h5_aln) != len(ext_aln):
        raise RuntimeError(f"{path}: alignment lengths differ")

    h5_pos = 0
    ext_pos = 0

    for a, b in zip(h5_aln, ext_aln):
        if a != "-":
            h5_pos += 1
        if b != "-":
            ext_pos += 1

        if a != "-" and h5_pos in targets:
            if b == "-":
                raise RuntimeError(
                    f"{acc}: H5 position {h5_pos} maps to gap"
                )

            coord_rows.append({
                "accession": acc,
                "H5_position": h5_pos,
                "H5_residue": a,
                "external_position": ext_pos,
                "external_residue": b,
            })

# ------------------------------------------------------------------
# Intersect coordinates with InterPro/member-database intervals
# ------------------------------------------------------------------

detail = []

for c in coord_rows:
    acc = c["accession"]
    pos = c["external_position"]

    path = iprdir / f"ST3GAL4__{acc}.json"

    with path.open() as fh:
        data = json.load(fh)

    for result in data.get("results", []) or []:
        meta = result.get("metadata", {}) or {}

        db = meta.get("source_database", "")
        entry_acc = meta.get("accession", "")
        name = meta.get("name", "")
        entry_type = meta.get("type", "")

        for protein in result.get("proteins", []) or []:
            for loc in protein.get("entry_protein_locations", []) or []:
                for frag in loc.get("fragments", []) or []:
                    start = frag.get("start")
                    end = frag.get("end")

                    if start is None or end is None:
                        continue

                    if start <= pos <= end:
                        detail.append({
                            **c,
                            "source_database": db,
                            "entry_accession": entry_acc,
                            "entry_name": name,
                            "entry_type": entry_type,
                            "interval_start": start,
                            "interval_end": end,
                        })

fields = [
    "accession",
    "H5_position",
    "H5_residue",
    "external_position",
    "external_residue",
    "source_database",
    "entry_accession",
    "entry_name",
    "entry_type",
    "interval_start",
    "interval_end",
]

with out_detail.open("w") as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(detail)

# ------------------------------------------------------------------
# Summarise cross-reference concordance by candidate + annotation
# ------------------------------------------------------------------

support = defaultdict(set)
names = {}

for r in detail:
    key = (
        r["H5_position"],
        r["H5_residue"],
        r["source_database"],
        r["entry_accession"],
        r["entry_type"],
    )
    support[key].add(r["accession"])
    names[key] = r["entry_name"]

summary = []

for key, accessions in support.items():
    h5pos, h5res, db, entry_acc, entry_type = key

    summary.append({
        "H5_position": h5pos,
        "H5_residue": h5res,
        "source_database": db,
        "entry_accession": entry_acc,
        "entry_name": names[key],
        "entry_type": entry_type,
        "references_overlapping": len(accessions),
        "references_total": 7,
        "all_seven_overlap": "YES" if len(accessions) == 7 else "NO",
        "accessions": ",".join(sorted(accessions)),
    })

summary.sort(
    key=lambda r: (
        int(r["H5_position"]),
        r["source_database"],
        r["entry_accession"],
    )
)

summary_fields = [
    "H5_position",
    "H5_residue",
    "source_database",
    "entry_accession",
    "entry_name",
    "entry_type",
    "references_overlapping",
    "references_total",
    "all_seven_overlap",
    "accessions",
]

with out_summary.open("w") as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(summary)

print(f"Coordinate mappings: {len(coord_rows)}")
print(f"Detailed overlap rows: {len(detail)}")
print(f"Summary rows: {len(summary)}")
print()

print("=== ALL-SEVEN INTERPRO CONTEXT ===")

for r in summary:
    if r["all_seven_overlap"] == "YES":
        print(
            f'{r["H5_residue"]}{r["H5_position"]}\t'
            f'{r["source_database"]}\t'
            f'{r["entry_accession"]}\t'
            f'{r["entry_name"]}'
        )

print()
print(f"Wrote: {out_detail}")
print(f"Wrote: {out_summary}")
