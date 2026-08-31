#!/usr/bin/env python3

import csv
import json
from pathlib import Path

RAW = Path(
    "results/H6/H6.2/raw/interproscan_CMAH_X2/"
    "CMAH_XP_024839792.1__interproscan.json"
)

MATCH_OUT = Path(
    "results/H6/H6.2/processed/"
    "H6.2_CMAH_RefSeq_X2_interproscan_matches.tsv"
)

AUDIT_OUT = Path(
    "results/H6/H6.2/qc/"
    "H6.2_CMAH_RefSeq_X2_candidate_overlap_audit.tsv"
)

CANDIDATES = [
    ("23:32739285:A>G", 386, "N"),
    ("23:32739292:C>G", 388, "I"),
    ("23:32743918:C>T", 424, "P"),
    ("23:32745866:T>A", 512, "F"),
]

data = json.loads(RAW.read_text())
result = data["results"][0]
sequence = result["sequence"]

match_rows = []

for match in result["matches"]:
    sig = match.get("signature", {})
    entry = sig.get("entry") or {}
    library = (
        sig.get("signatureLibraryRelease", {})
        .get("library", "")
    )

    for loc in match.get("locations", []):
        match_rows.append({
            "reference_protein": "XP_024839792.1",
            "protein_length": len(sequence),
            "library": library,
            "signature_accession": sig.get("accession", ""),
            "signature_name": sig.get("name", "") or "",
            "entry_type": entry.get("type", "") or "",
            "entry_accession": entry.get("accession", "") or "",
            "entry_name": entry.get("name", "") or "",
            "interval_start": loc.get("start", ""),
            "interval_end": loc.get("end", ""),
            "raw_source": str(RAW),
        })

match_fields = [
    "reference_protein",
    "protein_length",
    "library",
    "signature_accession",
    "signature_name",
    "entry_type",
    "entry_accession",
    "entry_name",
    "interval_start",
    "interval_end",
    "raw_source",
]

with MATCH_OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(match_rows)

audit_rows = []

for candidate_id, pos, aa in CANDIDATES:
    observed = sequence[pos - 1]
    if observed != aa:
        raise ValueError(
            f"{candidate_id}: expected {aa} at {pos}, observed {observed}"
        )

    overlaps = [
        row for row in match_rows
        if int(row["interval_start"]) <= pos <= int(row["interval_end"])
    ]

    libraries = sorted({row["library"] for row in overlaps})
    signatures = sorted({
        f'{row["library"]}:{row["signature_accession"]}'
        for row in overlaps
    })

    panther = [row for row in overlaps if row["library"] == "PANTHER"]
    non_panther = [row for row in overlaps if row["library"] != "PANTHER"]

    if non_panther:
        interpretation = "NARROWER_OR_ADDITIONAL_FEATURE_OVERLAP"
    elif panther:
        interpretation = "BROAD_FAMILY_OVERLAP_ONLY"
    else:
        interpretation = "NO_FEATURE_OVERLAP"

    audit_rows.append({
        "candidate_id": candidate_id,
        "protein_position": pos,
        "reference_aa": aa,
        "total_overlap_count": len(overlaps),
        "overlap_libraries": ";".join(libraries),
        "overlap_signatures": ";".join(signatures),
        "panther_overlap": "YES" if panther else "NO",
        "non_panther_overlap": "YES" if non_panther else "NO",
        "interpretation": interpretation,
    })

audit_fields = [
    "candidate_id",
    "protein_position",
    "reference_aa",
    "total_overlap_count",
    "overlap_libraries",
    "overlap_signatures",
    "panther_overlap",
    "non_panther_overlap",
    "interpretation",
]

with AUDIT_OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=audit_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(audit_rows)

print("Raw matches:", len(match_rows))
print("Distinct libraries:", len({r["library"] for r in match_rows}))

for row in audit_rows:
    print(
        row["candidate_id"],
        f'pos={row["protein_position"]}',
        f'overlaps={row["total_overlap_count"]}',
        f'libs={row["overlap_libraries"]}',
        row["interpretation"],
        sep="\t",
    )

if len(match_rows) != 10:
    raise ValueError(f"Expected 10 raw match intervals, found {len(match_rows)}")

if any(
    row["total_overlap_count"] != 1
    or row["overlap_libraries"] != "PANTHER"
    or row["interpretation"] != "BROAD_FAMILY_OVERLAP_ONLY"
    for row in audit_rows
):
    raise ValueError("CMAH X2 candidate-overlap audit did not match validated expectation")

print("CMAH X2 INTERPROSCAN AUDIT: PASS")
