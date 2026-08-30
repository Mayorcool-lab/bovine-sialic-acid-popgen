#!/usr/bin/env python3

import csv
from pathlib import Path

CANDIDATES = Path("results/vep/phase_g_protein_altering_FINAL.tsv")
ALIGN_DIR = Path("results/H5/alignments")
OUT = Path("results/H5/H5_candidate_site_coverage.tsv")

MIN_COMPARATIVE = 6

REFSEQ_CMAH_POSITIONS = {386, 388, 424, 512}
REFSEQ_CMAH_PROTEIN = "XP_024839792.1"

TRANSCRIPT_TO_PROTEIN = {
    "ENSBTAT00000081903": "ENSBTAP00000073225",  # CMAH
    "ENSBTAT00000013373": "ENSBTAP00000013373",  # CMAS
    "ENSBTAT00000001520": "ENSBTAP00000001520",  # GNE
    "ENSBTAT00000000918": "ENSBTAP00000000918",  # NANP
    "ENSBTAT00000025318": "ENSBTAP00000025318",  # NANS
    "ENSBTAT00000067415": "ENSBTAP00000073421",  # NEU1
    "ENSBTAT00000036738": "ENSBTAP00000036593",  # NEU3
    "ENSBTAT00000018680": "ENSBTAP00000018680",  # SLC35A1
    "ENSBTAT00000001538": "ENSBTAP00000001538",  # ST3GAL1
    "ENSBTAT00000075324": "ENSBTAP00000073357",  # ST3GAL4
    "ENSBTAT00000083430": "ENSBTAP00000063809",  # ST6GAL1
    "ENSBTAT00000077891": "ENSBTAP00000057846",  # ST8SIA2
    "ENSBTAT00000028590": "ENSBTAP00000028590",  # ST8SIA4
}


def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq)))

    return records


def get_reference(records, marker):
    refs = [
        (header, seq)
        for header, seq in records
        if marker in header
    ]

    if len(refs) != 1:
        raise RuntimeError(
            f"Expected exactly one reference containing "
            f"{marker}; found {len(refs)}"
        )

    return refs[0]


def map_residue(ref_aln, protein_position):
    residue_number = 0

    for column, aa in enumerate(ref_aln):
        if aa != "-":
            residue_number += 1

            if residue_number == protein_position:
                return column, aa

    return None, None


with open(CANDIDATES) as fh:
    candidates = list(csv.DictReader(fh, delimiter="\t"))

if len(candidates) != 62:
    raise RuntimeError(
        f"Expected 62 Phase G candidates, found {len(candidates)}"
    )


output_rows = []

for row in candidates:

    gene = row["SYMBOL"]
    protein_position = int(row["Protein_position"])
    consequence = row["Consequence"]
    feature = row["Feature"].strip()

    # Four manually rescued CMAH variants use RefSeq X2.
    is_refseq_cmah = (
        gene == "CMAH"
        and not feature
        and protein_position in REFSEQ_CMAH_POSITIONS
    )

    if is_refseq_cmah:
        coordinate_system = "RefSeq_X2"
        reference_transcript = ""
        reference_protein = REFSEQ_CMAH_PROTEIN
        aln_path = ALIGN_DIR / "CMAH_RefSeq_X2_aligned.fasta"
        reference_marker = "RefSeq_X2_REFERENCE"

    else:
        coordinate_system = "Ensembl110"
        reference_transcript = feature

        if reference_transcript not in TRANSCRIPT_TO_PROTEIN:
            raise RuntimeError(
                f"No validated protein mapping for transcript "
                f"{reference_transcript}"
            )

        reference_protein = TRANSCRIPT_TO_PROTEIN[reference_transcript]
        aln_path = ALIGN_DIR / f"{gene}_aligned.fasta"
        reference_marker = "Ensembl110_REFERENCE"

    if not aln_path.exists():
        raise RuntimeError(
            f"Alignment missing for {gene}: {aln_path}"
        )

    records = read_fasta(aln_path)
    ref_header, ref_aln = get_reference(
        records,
        reference_marker
    )

    ref_ungapped = ref_aln.replace("-", "")
    protein_length = len(ref_ungapped)

    # Stop-lost affects the terminal stop codon rather than
    # an encoded reference amino acid.
    if "stop_lost" in consequence:

        if protein_position != protein_length + 1:
            raise RuntimeError(
                f"{gene} {row['CHROM']}:{row['POS']}: "
                f"stop_lost position={protein_position}, "
                f"expected={protein_length + 1}"
            )

        output_rows.append({
            "SYMBOL": gene,
            "CHROM": row["CHROM"],
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "HGVSp": row["HGVSp"],
            "Consequence": consequence,
            "coordinate_system": coordinate_system,
            "reference_transcript": reference_transcript,
            "reference_protein": reference_protein,
            "protein_position": protein_position,
            "reference_aa": "*",
            "alignment_column": "",
            "comparative_total": len(records) - 1,
            "comparative_non_gap": "",
            "coverage_status": "STOP_LOST_SPECIAL_CASE",
        })

        continue

    column, reference_aa = map_residue(
        ref_aln,
        protein_position
    )

    if column is None:
        raise RuntimeError(
            f"{gene} {row['CHROM']}:{row['POS']}: "
            f"cannot map protein position {protein_position}"
        )

    # Validate Phase G reference amino acid where available.
    amino_acids = row["Amino_acids"].strip()

    if amino_acids and "/" in amino_acids:
        expected_reference_aa = amino_acids.split("/")[0]

        if reference_aa != expected_reference_aa:
            raise RuntimeError(
                f"{gene} {row['CHROM']}:{row['POS']}: "
                f"reference AA mismatch: "
                f"PhaseG={expected_reference_aa}, "
                f"alignment={reference_aa}"
            )

    comparative = [
        aln[column]
        for header, aln in records
        if reference_marker not in header
    ]

    non_gap = sum(
        aa != "-"
        for aa in comparative
    )

    status = (
        "ASSESSABLE"
        if non_gap >= MIN_COMPARATIVE
        else "NOT_ASSESSABLE"
    )

    output_rows.append({
        "SYMBOL": gene,
        "CHROM": row["CHROM"],
        "POS": row["POS"],
        "REF": row["REF"],
        "ALT": row["ALT"],
        "HGVSp": row["HGVSp"],
        "Consequence": consequence,
        "coordinate_system": coordinate_system,
        "reference_transcript": reference_transcript,
        "reference_protein": reference_protein,
        "protein_position": protein_position,
        "reference_aa": reference_aa,
        "alignment_column": column + 1,
        "comparative_total": len(comparative),
        "comparative_non_gap": non_gap,
        "coverage_status": status,
    })


if len(output_rows) != 62:
    raise RuntimeError(
        f"Expected 62 output rows, found {len(output_rows)}"
    )


fields = [
    "SYMBOL",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "HGVSp",
    "Consequence",
    "coordinate_system",
    "reference_transcript",
    "reference_protein",
    "protein_position",
    "reference_aa",
    "alignment_column",
    "comparative_total",
    "comparative_non_gap",
    "coverage_status",
]


OUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUT, "w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(output_rows)


from collections import Counter

counts = Counter(
    row["coverage_status"]
    for row in output_rows
)

coordinate_counts = Counter(
    row["coordinate_system"]
    for row in output_rows
)

print("Candidates written:", len(output_rows))

for status in sorted(counts):
    print(f"{status}: {counts[status]}")

print()

for system in sorted(coordinate_counts):
    print(f"{system}: {coordinate_counts[system]}")

print()
print("Output:", OUT)
print("H5 candidate-site mapping: PASS")
