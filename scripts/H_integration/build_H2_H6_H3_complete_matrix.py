#!/usr/bin/env python3
import csv
from pathlib import Path

matrix_path = Path(
    "results/H_integration/H2_H6_62_candidate_evidence_matrix.tsv"
)

h3_summary_path = Path(
    "results/H_integration/H3_62_candidate_breed_summary.tsv"
)

out_path = Path(
    "results/H_integration/H2_H6_H3_62_candidate_complete_evidence_matrix.tsv"
)


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


matrix = read_tsv(matrix_path)
h3 = read_tsv(h3_summary_path)

matrix_map = {r["candidate_id"]: r for r in matrix}
h3_map = {r["candidate_id"]: r for r in h3}

if len(matrix_map) != 62:
    raise RuntimeError(f"Expected 62 matrix candidates, observed {len(matrix_map)}")

if len(h3_map) != 62:
    raise RuntimeError(f"Expected 62 H3 summary candidates, observed {len(h3_map)}")

missing_h3 = sorted(set(matrix_map) - set(h3_map))
extra_h3 = sorted(set(h3_map) - set(matrix_map))

if missing_h3 or extra_h3:
    raise RuntimeError(f"Candidate mismatch: missing_h3={missing_h3}, extra_h3={extra_h3}")

outrows = []

for candidate_id, base in matrix_map.items():
    breed = h3_map[candidate_id]

    if base["SYMBOL"] != breed["SYMBOL"]:
        raise RuntimeError(
            f"{candidate_id}: SYMBOL mismatch "
            f"{base['SYMBOL']} != {breed['SYMBOL']}"
        )

    out = dict(base)

    for field, value in breed.items():
        if field in {"candidate_id", "SYMBOL"}:
            continue

        if field in out:
            raise RuntimeError(
                f"{candidate_id}: duplicate field during H3 join: {field}"
            )

        out[field] = value

    outrows.append(out)

if len(outrows) != 62:
    raise RuntimeError(f"Expected 62 output rows, observed {len(outrows)}")

if len({r["candidate_id"] for r in outrows}) != 62:
    raise RuntimeError("Output candidate IDs are not unique")

fieldnames = list(outrows[0].keys())

ranking_terms = [
    c for c in fieldnames
    if (
        "score" in c.lower()
        or "rank" in c.lower()
        or "priority" in c.lower()
    )
    and "rate4site_score" not in c.lower()
]

with out_path.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(outrows)

print("Output:", out_path)
print("Rows:", len(outrows))
print("Columns:", len(fieldnames))
print("Unique candidates:", len({r["candidate_id"] for r in outrows}))
print("Missing H3 candidates:", len(missing_h3))
print("Extra H3 candidates:", len(extra_h3))
print("Ranking/priority columns added:", ranking_terms)
print("Composite score added:", 0)
print("COMPLETE H2-H6-H3 EVIDENCE MATRIX BUILD: PASS")
