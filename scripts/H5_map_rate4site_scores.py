#!/usr/bin/env python3

import csv
import re
from pathlib import Path

INPUT = Path("results/H5/H5_candidate_site_coverage.tsv")
R4S_DIR = Path("results/H5/rate4site")
OUTPUT = Path("results/H5/H5_62_variant_conservation_evidence.tsv")

EXPECTED = {
    "ASSESSABLE": 57,
    "NOT_ASSESSABLE": 4,
    "STOP_LOST_SPECIAL_CASE": 1,
}


def parse_rate4site(path):
    """Return Rate4Site results indexed by bovine reference position."""

    results = {}

    with path.open() as fh:
        for line in fh:
            fields = line.split()

            # Skip comments, headers and blank lines.
            if not fields or not fields[0].isdigit():
                continue

            # Expected Rate4Site residue row:
            # POS SEQ SCORE QQ-INTERVAL STD MSA_DATA
            if len(fields) < 6:
                raise ValueError(
                    f"Unexpected Rate4Site row in {path}: "
                    f"{line.rstrip()}"
                )

            pos = int(fields[0])
            aa = fields[1]
            score = fields[2]
            qq_interval = fields[3]
            std = fields[4]
            msa_data = fields[5]

            if not re.fullmatch(r"\d+/\d+", msa_data):
                raise ValueError(
                    f"Cannot parse MSA DATA in {path}: "
                    f"{line.rstrip()}"
                )

            if pos in results:
                raise ValueError(
                    f"Duplicate position {pos} in {path}"
                )

            results[pos] = {
                "aa": aa,
                "score": score,
                "qq_interval": qq_interval,
                "std": std,
                "msa_data": msa_data,
            }

    if not results:
        raise ValueError(
            f"No Rate4Site residue rows found in {path}"
        )

    return results


# ------------------------------------------------------------
# Load all production Rate4Site outputs
# ------------------------------------------------------------

score_tables = {}

for path in sorted(R4S_DIR.glob("*_rate4site.txt")):
    name = path.name.replace("_rate4site.txt", "")
    score_tables[name] = parse_rate4site(path)


required_tracks = {
    "CMAH",
    "CMAS",
    "GNE",
    "NANP",
    "NANS",
    "NEU1",
    "NEU3",
    "SLC35A1",
    "ST3GAL1",
    "ST3GAL4",
    "ST6GAL1",
    "ST8SIA2",
    "ST8SIA4",
    "CMAH_RefSeq_X2",
}

missing_tracks = required_tracks - set(score_tables)

if missing_tracks:
    raise SystemExit(
        f"ERROR: missing Rate4Site tracks: "
        f"{sorted(missing_tracks)}"
    )


# ------------------------------------------------------------
# Read validated 62-candidate site coverage table
# ------------------------------------------------------------

with INPUT.open(newline="") as fh:
    reader = csv.DictReader(
        fh,
        delimiter="\t",
    )
    rows = list(reader)


if len(rows) != 62:
    raise SystemExit(
        f"ERROR: expected 62 candidates, found {len(rows)}"
    )


output_fields = list(rows[0].keys()) + [
    "rate4site_status",
    "rate4site_score",
    "rate4site_qq_interval",
    "rate4site_std",
    "rate4site_msa_data",
]


status_counts = {
    "ASSESSABLE": 0,
    "NOT_ASSESSABLE": 0,
    "STOP_LOST_SPECIAL_CASE": 0,
}

errors = 0


# ------------------------------------------------------------
# Map Rate4Site results back to candidate variants
# ------------------------------------------------------------

for row in rows:

    coverage_status = row["coverage_status"]

    if coverage_status not in status_counts:
        raise SystemExit(
            f"ERROR: unexpected coverage status: "
            f"{coverage_status}"
        )

    status_counts[coverage_status] += 1

    # Default: no interpretable Rate4Site evidence.
    row["rate4site_status"] = coverage_status
    row["rate4site_score"] = ""
    row["rate4site_qq_interval"] = ""
    row["rate4site_std"] = ""
    row["rate4site_msa_data"] = ""

    # Only scientifically assessable residue sites receive
    # interpretable conservation evidence.
    if coverage_status != "ASSESSABLE":
        continue


    symbol = row["SYMBOL"]
    coordinate_system = row["coordinate_system"]


    # --------------------------------------------------------
    # Determine correct Rate4Site track
    # --------------------------------------------------------

    if coordinate_system == "RefSeq_X2":

        if symbol != "CMAH":
            raise SystemExit(
                "ERROR: RefSeq_X2 used for non-CMAH row: "
                f"{row}"
            )

        track = "CMAH_RefSeq_X2"


    elif coordinate_system == "Ensembl110":

        track = symbol


    else:
        raise SystemExit(
            f"ERROR: unknown coordinate system: "
            f"{coordinate_system}"
        )


    protein_pos = int(row["protein_position"])


    if track not in score_tables:
        raise SystemExit(
            f"ERROR: missing score table for {track}"
        )


    if protein_pos not in score_tables[track]:
        raise SystemExit(
            f"ERROR: {track} position {protein_pos} "
            "absent from Rate4Site output"
        )


    result = score_tables[track][protein_pos]


    # --------------------------------------------------------
    # Validate reference amino acid
    # --------------------------------------------------------

    if result["aa"] != row["reference_aa"]:

        print(
            "ERROR reference-AA mismatch:",
            symbol,
            row["CHROM"],
            row["POS"],
            f"candidate={row['reference_aa']}",
            f"Rate4Site={result['aa']}",
        )

        errors += 1
        continue


    # --------------------------------------------------------
    # Independently reconcile Rate4Site coverage with
    # our previously validated H5 site-coverage table.
    #
    # H5 table excludes bovine from comparative counts.
    # Rate4Site MSA DATA includes bovine.
    # --------------------------------------------------------

    msa_non_gap, msa_total = map(
        int,
        result["msa_data"].split("/")
    )


    expected_total = (
        int(row["comparative_total"]) + 1
    )

    expected_non_gap = (
        int(row["comparative_non_gap"]) + 1
    )


    if (
        msa_total != expected_total
        or msa_non_gap != expected_non_gap
    ):

        print(
            "ERROR coverage mismatch:",
            symbol,
            row["CHROM"],
            row["POS"],
            f"H5={expected_non_gap}/{expected_total}",
            f"Rate4Site={result['msa_data']}",
        )

        errors += 1
        continue


    # --------------------------------------------------------
    # Validation passed: retain Rate4Site evidence
    # --------------------------------------------------------

    row["rate4site_score"] = result["score"]
    row["rate4site_qq_interval"] = result["qq_interval"]
    row["rate4site_std"] = result["std"]
    row["rate4site_msa_data"] = result["msa_data"]


# ------------------------------------------------------------
# Global validation
# ------------------------------------------------------------

if status_counts != EXPECTED:
    raise SystemExit(
        "ERROR: unexpected candidate classification: "
        f"{status_counts}"
    )


if errors:
    raise SystemExit(
        f"ERROR: validation errors: {errors}"
    )


# ------------------------------------------------------------
# Write authoritative candidate-level H5F table
# ------------------------------------------------------------

with OUTPUT.open(
    "w",
    newline="",
) as fh:

    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=output_fields,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)


scored = sum(
    bool(row["rate4site_score"])
    for row in rows
)


print(f"Candidates: {len(rows)}")
print(
    f"ASSESSABLE: "
    f"{status_counts['ASSESSABLE']}"
)
print(
    f"NOT_ASSESSABLE: "
    f"{status_counts['NOT_ASSESSABLE']}"
)
print(
    "STOP_LOST_SPECIAL_CASE:",
    status_counts["STOP_LOST_SPECIAL_CASE"],
)
print(
    f"Rate4Site scores retained: {scored}"
)
print(
    f"Validation errors: {errors}"
)
print(
    f"Output: {OUTPUT}"
)


if scored != 57:
    raise SystemExit(
        "ERROR: expected 57 retained Rate4Site scores, "
        f"found {scored}"
    )


print(
    "H5F candidate Rate4Site mapping: PASS"
)
