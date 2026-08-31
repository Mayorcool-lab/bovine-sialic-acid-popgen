#!/usr/bin/env python3

import csv
from pathlib import Path

candidates_path = Path(
    "results/H6/H6.2/processed/H6.2_straightforward_candidate_targets.tsv"
)
uniprot_path = Path(
    "results/H6/H6.2/processed/H6.2_uniprot_features.tsv"
)
interpro_path = Path(
    "results/H6/H6.2/processed/H6.2_interpro_intervals.tsv"
)

uniprot_out = Path(
    "results/H6/H6.2/processed/H6.2_candidate_uniprot_overlaps.tsv"
)
interpro_out = Path(
    "results/H6/H6.2/processed/H6.2_candidate_interpro_overlaps.tsv"
)

with candidates_path.open() as fh:
    candidates = list(csv.DictReader(fh, delimiter="\t"))

with uniprot_path.open() as fh:
    uniprot = list(csv.DictReader(fh, delimiter="\t"))

with interpro_path.open() as fh:
    interpro = list(csv.DictReader(fh, delimiter="\t"))

uniprot_fields = [
    "candidate_id",
    "SYMBOL",
    "external_accession",
    "external_position",
    "external_aa",
    "mapping_status",
    "feature_transfer_status",
    "feature_type",
    "feature_start",
    "feature_end",
    "feature_description",
    "evidence_code",
    "evidence_source",
    "evidence_source_accession",
    "uniprot_entry_type",
]

interpro_fields = [
    "candidate_id",
    "SYMBOL",
    "external_accession",
    "external_position",
    "external_aa",
    "mapping_status",
    "feature_transfer_status",
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

uni_rows = []
ipr_rows = []

for c in candidates:
    pos = int(c["external_position"])
    gene = c["SYMBOL"]
    acc = c["external_accession"]

    for f in uniprot:
        if f["SYMBOL"] != gene:
            continue
        if f["uniprot_accession"] != acc:
            continue

        start = int(f["feature_start"])
        end = int(f["feature_end"])

        if start <= pos <= end:
            row = {
                "candidate_id": c["candidate_id"],
                "SYMBOL": gene,
                "external_accession": acc,
                "external_position": c["external_position"],
                "external_aa": c["external_aa"],
                "mapping_status": c["mapping_status"],
                "feature_transfer_status": c["feature_transfer_status"],
            }
            row.update({
                k: f[k]
                for k in [
                    "feature_type",
                    "feature_start",
                    "feature_end",
                    "feature_description",
                    "evidence_code",
                    "evidence_source",
                    "evidence_source_accession",
                    "uniprot_entry_type",
                ]
            })
            uni_rows.append(row)

    for f in interpro:
        if f["SYMBOL"] != gene:
            continue
        if f["protein_accession"] != acc:
            continue

        start = int(f["interval_start"])
        end = int(f["interval_end"])

        if start <= pos <= end:
            row = {
                "candidate_id": c["candidate_id"],
                "SYMBOL": gene,
                "external_accession": acc,
                "external_position": c["external_position"],
                "external_aa": c["external_aa"],
                "mapping_status": c["mapping_status"],
                "feature_transfer_status": c["feature_transfer_status"],
            }
            row.update({
                k: f[k]
                for k in [
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
            })
            ipr_rows.append(row)

with uniprot_out.open("w") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=uniprot_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(uni_rows)

with interpro_out.open("w") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=interpro_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(ipr_rows)

print(f"UniProt overlap rows: {len(uni_rows)}")
print(f"InterPro overlap rows: {len(ipr_rows)}")
print(f"Wrote: {uniprot_out}")
print(f"Wrote: {interpro_out}")
