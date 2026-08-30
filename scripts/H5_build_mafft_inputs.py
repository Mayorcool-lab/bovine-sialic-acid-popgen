#!/usr/bin/env python3

import gzip
import re
from pathlib import Path

REFERENCE_FASTA = Path(
    "data/reference/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz"
)

ORTHO_DIR = Path("results/H5/orthologues")
OUT_DIR = Path("results/H5/mafft_inputs")

CMAH_REFSEQ_FASTA = Path(
    "results/H5/refseq/CMAH_XP_024839792.1.fasta"
)
CMAH_REFSEQ_ACCESSION = "XP_024839792.1"
CMAH_REFSEQ_LENGTH = 577

REFERENCE = {
    "CMAH": ("ENSBTAP00000073225", 383),
    "CMAS": ("ENSBTAP00000013373", 434),
    "GNE": ("ENSBTAP00000001520", 722),
    "NANP": ("ENSBTAP00000000918", 248),
    "NANS": ("ENSBTAP00000025318", 359),
    "NEU1": ("ENSBTAP00000073421", 416),
    "NEU3": ("ENSBTAP00000036593", 428),
    "SLC35A1": ("ENSBTAP00000018680", 337),
    "ST3GAL1": ("ENSBTAP00000001538", 339),
    "ST3GAL4": ("ENSBTAP00000073357", 395),
    "ST6GAL1": ("ENSBTAP00000063809", 405),
    "ST8SIA2": ("ENSBTAP00000057846", 375),
    "ST8SIA4": ("ENSBTAP00000028590", 359),
}


def read_fasta(path):
    records = []
    header = None
    seq = []

    opener = gzip.open if str(path).endswith(".gz") else open

    with opener(path, "rt") as fh:
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


def write_fasta(path, records):
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")


# Extract exact release-110 bovine proteins
wanted = {
    protein: (gene, expected_len)
    for gene, (protein, expected_len) in REFERENCE.items()
}

bovine = {}

for header, seq in read_fasta(REFERENCE_FASTA):
    protein = header.split()[0]
    protein = re.sub(r"\.\d+$", "", protein)

    if protein in wanted:
        gene, expected_len = wanted[protein]

        if len(seq) != expected_len:
            raise RuntimeError(
                f"{gene}: expected {expected_len} aa, observed {len(seq)}"
            )

        bovine[gene] = (protein, seq)

if set(bovine) != set(REFERENCE):
    missing = sorted(set(REFERENCE) - set(bovine))
    raise RuntimeError(f"Missing bovine reference proteins: {missing}")

OUT_DIR.mkdir(parents=True, exist_ok=True)

total_comparative = 0

for gene in REFERENCE:
    protein, ref_seq = bovine[gene]

    orthologue_fasta = (
        ORTHO_DIR / f"{gene}_one2one_sequences.fasta"
    )

    comparative = read_fasta(orthologue_fasta)

    output_records = [
        (
            f"{gene}|bos_taurus|{protein}|Ensembl110_REFERENCE",
            ref_seq,
        )
    ] + comparative

    outpath = OUT_DIR / f"{gene}_mafft_input.fasta"

    write_fasta(outpath, output_records)

    total_comparative += len(comparative)

    print(
        f"{gene}\t"
        f"reference=1\t"
        f"comparative={len(comparative)}\t"
        f"total={len(output_records)}"
    )

print()
print("Genes:", len(REFERENCE))
print("Comparative sequences:", total_comparative)
print("Bovine reference sequences:", len(REFERENCE))
print("Total MAFFT input sequences:", total_comparative + len(REFERENCE))


# Build separate CMAH RefSeq-X2 coordinate track.
#
# These four manually rescued CMAH variants use XP_024839792.1
# coordinates and must not be mapped onto the 383-aa Ensembl110
# CMAH reference protein.

refseq_records = read_fasta(CMAH_REFSEQ_FASTA)

if len(refseq_records) != 1:
    raise RuntimeError(
        f"CMAH RefSeq: expected one sequence, "
        f"found {len(refseq_records)}"
    )

refseq_header, refseq_seq = refseq_records[0]
refseq_accession = refseq_header.split()[0]

if refseq_accession != CMAH_REFSEQ_ACCESSION:
    raise RuntimeError(
        f"CMAH RefSeq: expected accession "
        f"{CMAH_REFSEQ_ACCESSION}, observed {refseq_accession}"
    )

if len(refseq_seq) != CMAH_REFSEQ_LENGTH:
    raise RuntimeError(
        f"CMAH RefSeq: expected {CMAH_REFSEQ_LENGTH} aa, "
        f"observed {len(refseq_seq)}"
    )

cmah_comparative = read_fasta(
    ORTHO_DIR / "CMAH_one2one_sequences.fasta"
)

if len(cmah_comparative) != 15:
    raise RuntimeError(
        f"CMAH RefSeq: expected 15 comparative sequences, "
        f"found {len(cmah_comparative)}"
    )

refseq_output_records = [
    (
        f"CMAH|bos_taurus|{CMAH_REFSEQ_ACCESSION}|"
        f"RefSeq_X2_REFERENCE",
        refseq_seq,
    )
] + cmah_comparative

refseq_outpath = (
    OUT_DIR / "CMAH_RefSeq_X2_mafft_input.fasta"
)

write_fasta(
    refseq_outpath,
    refseq_output_records,
)

print()
print(
    "CMAH_RefSeq_X2\t"
    "reference=1\t"
    f"comparative={len(cmah_comparative)}\t"
    f"total={len(refseq_output_records)}"
)

print()
print("Alignment input sets:", len(REFERENCE) + 1)
print(
    "Total sequences including RefSeq track:",
    total_comparative
    + len(REFERENCE)
    + len(refseq_output_records)
)
