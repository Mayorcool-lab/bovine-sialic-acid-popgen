#!/usr/bin/env python3

"""
Phase H6.1: build candidate-to-external-reference mapping evidence.

Authoritative candidate source:
    results/H5/H5_62_variant_conservation_evidence.tsv

Purpose:
    Encode validated protein-reference reconciliation before any
    UniProt/InterPro functional feature transfer.

Important:
    Mapping status describes coordinate/reference validity only.
    It is not evidence that a variant is functional or benign.
"""

import csv
from pathlib import Path

INPUT = Path("results/H5/H5_62_variant_conservation_evidence.tsv")
OUTPUT = Path("results/H6/H6_candidate_reference_mapping.tsv")


# Whole-sequence exact matches established during H6.1.
#
# Because the external sequence is identical to the H5 reference,
# H5 residue position == external residue position.
EXACT_SEQUENCE_MAP = {
    "CMAS": {
        "accession": "Q3SZM5",
        "database": "UniProtKB",
        "evidence": "Exact whole-protein sequence identity to H5 reference",
    },
    "GNE": {
        "accession": "F1MIL9",
        "database": "UniProtKB",
        "evidence": (
            "Official Ensembl Protein-to-UniProtKB mapping plus exact "
            "whole-protein sequence identity to H5 reference"
        ),
    },
    "NANP": {
        "accession": "A7YWM8",
        "database": "UniProtKB",
        "evidence": (
            "Exact whole-protein sequence identity to H5 reference; "
            "A0ABI0NY54 independently has the same exact 248-aa sequence"
        ),
    },
    "NANS": {
        "accession": "F6RGT7",
        "database": "UniProtKB",
        "evidence": "Exact whole-protein sequence identity to H5 reference",
    },
    "SLC35A1": {
        "accession": "F6R148",
        "database": "UniProtKB",
        "evidence": "Exact whole-protein sequence identity to H5 reference",
    },
    "ST6GAL1": {
        "accession": "F1MHF1",
        "database": "UniProtKB",
        "evidence": "Exact whole-protein sequence identity to H5 reference",
    },
}


# Candidate-level coordinate mappings established by MAFFT pairwise
# reconciliation during H6.1.
#
# These proteins are NOT whole-sequence identical to the H5 reference.
# Feature transfer is therefore authorised only at candidate residues
# explicitly reconciled below.
ALIGNMENT_RESIDUE_MAP = {
    "NEU3": {
        "accession": "F1MRJ3",
        "database": "UniProtKB",
        "positions": {
            94: 104,
            136: 146,
            159: 169,
            309: 319,
            356: 366,
            416: 426,
        },
    },
    "ST3GAL1": {
        "accession": "A0AAF6DLT4",
        "database": "UniProtKB",
        "positions": {
            102: 105,
            105: 108,
        },
    },
    "ST8SIA2": {
        "accession": "A0A3Q1LM08",
        "database": "UniProtKB",
        "positions": {
            63: 42,
            129: 108,
        },
    },
    "ST8SIA4": {
        "accession": "F1MVI9",
        "database": "UniProtKB",
        "positions": {
            73: 41,
            198: 166,
            203: 171,
            305: 273,
        },
    },
}


# NEU1 was reconciled against reviewed bovine UniProtKB/Swiss-Prot
# A6BMK7 (415 aa). Eight candidate residues map exactly.
# H5 L275 maps to a gap in A6BMK7 and is therefore explicitly blocked
# from residue-specific UniProt feature transfer.
NEU1_MAP = {
    "accession": "A6BMK7",
    "database": "UniProtKB/Swiss-Prot",
    "positions": {
        8: 8,
        13: 13,
        254: 254,
        294: 293,
        380: 379,
        381: 380,
        382: 381,
        400: 399,
    },
    "unmapped_gap_positions": {275},
}


# CMAH requires two independent coordinate systems.
#
# Ensembl110 ENSBTAP00000073225 (383 aa):
# reconciliation against bovine UniProtKB A0AAA9S174 supports exact
# candidate-residue transfer for positions 107, 210, 232 and 242.
# Positions 381 and 382 occur in a divergent C-terminal isoform context
# and are explicitly blocked from direct residue-specific UniProt transfer.
#
# RefSeq-X2 XP_024839792.1 (577 aa):
# retained as the authoritative reference for N386S, I388M, P424L and
# F512Y. UniParc UPI0000EBE207 independently confirms the exact sequence,
# but no active exact Bos taurus UniProtKB record was identified.
CMAH_MAP = {
    "Ensembl110": {
        "accession": "A0AAA9S174",
        "database": "UniProtKB",
        "positions": {
            107: 107,
            210: 210,
            232: 232,
            242: 242,
        },
        "unresolved_positions": {381, 382},
    },
    "RefSeq_X2": {
        "accession": "XP_024839792.1",
        "database": "RefSeq",
        "uniparc": "UPI0000EBE207",
        "positions": {386, 388, 424, 512},
    },
}


# ST3GAL4 has no single bovine UniProtKB record that is a defensible
# global match to the frozen 395-aa H5 reference. Seven bovine records
# were therefore examined by pairwise MAFFT reconciliation.
#
# Positions 190, 203, 320, 333 and 334 reconcile to the identical
# residue across all seven tested representations.
#
# Positions 4, 16, 43 and 59 are isoform-sensitive: some records map
# exactly while others contain gaps or mismatches. Position 38 has no
# exact match among the seven records. Position 87 is reference-divergent:
# every external record that maps the site carries S rather than H5 R.
#
# Position 396 represents the terminal stop after the 395-aa H5 protein
# and is handled as a stop-lost special case, not an ordinary residue.
ST3GAL4_MAP = {
    "tested_accessions": [
        "A0A3Q1M6Q2",
        "A0AAA9RT99",
        "A0AAA9TIA8",
        "A0AAA9U1T6",
        "A0AAA9TT54",
        "Q70E76",
        "A0AAF6DLN4",
    ],
    "robust_exact_positions": {190, 203, 320, 333, 334},
    "isoform_sensitive_positions": {4, 16, 43, 59},
    "no_exact_match_positions": {38},
    "reference_divergent_positions": {87},
    "stop_lost_positions": {396},
}


def candidate_id(row):
    return (
        f'{row["CHROM"]}:{row["POS"]}:'
        f'{row["REF"]}>{row["ALT"]}'
    )


def main():
    with INPUT.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if len(rows) != 62:
        raise ValueError(f"Expected 62 H5 candidates, found {len(rows)}")

    out = []

    for row in rows:
        gene = row["SYMBOL"]

        record = {
            "candidate_id": candidate_id(row),
            "SYMBOL": gene,
            "CHROM": row["CHROM"],
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "H5_coordinate_system": row["coordinate_system"],
            "H5_reference_protein": row["reference_protein"],
            "H5_position": row["protein_position"],
            "H5_reference_aa": row["reference_aa"],
            "external_database": "",
            "external_accession": "",
            "external_position": "",
            "external_aa": "",
            "mapping_status": "PENDING_RECONCILIATION",
            "feature_transfer_status": "NOT_YET_AUTHORISED",
            "mapping_evidence": "",
        }

        if gene in EXACT_SEQUENCE_MAP:
            m = EXACT_SEQUENCE_MAP[gene]

            record["external_database"] = m["database"]
            record["external_accession"] = m["accession"]
            record["external_position"] = row["protein_position"]
            record["external_aa"] = row["reference_aa"]
            record["mapping_status"] = "EXACT_SEQUENCE_MATCH"
            record["feature_transfer_status"] = "DIRECT_TRANSFER_ALLOWED"
            record["mapping_evidence"] = m["evidence"]

        elif gene == "ST3GAL4":
            h5_pos = int(row["protein_position"])
            tested = ",".join(ST3GAL4_MAP["tested_accessions"])

            record["external_database"] = "UniProtKB"
            record["external_accession"] = tested

            if h5_pos in ST3GAL4_MAP["robust_exact_positions"]:
                record["mapping_status"] = (
                    "ROBUST_MULTI_REFERENCE_RESIDUE_MATCH"
                )
                record["feature_transfer_status"] = (
                    "FEATURE_RECONCILIATION_ALLOWED"
                )
                record["mapping_evidence"] = (
                    "Candidate residue reconciles to the identical "
                    "amino acid across all seven tested bovine "
                    "UniProtKB representations; no single record is "
                    "treated as a globally equivalent protein reference"
                )

            elif h5_pos in ST3GAL4_MAP["isoform_sensitive_positions"]:
                record["mapping_status"] = "ISOFORM_SENSITIVE"
                record["feature_transfer_status"] = (
                    "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "Candidate residue maps exactly in some bovine "
                    "UniProtKB representations but to gap or mismatch "
                    "in others; isoform-specific cherry-picking blocked"
                )

            elif h5_pos in ST3GAL4_MAP["no_exact_match_positions"]:
                record["mapping_status"] = (
                    "NO_EXACT_EXTERNAL_RESIDUE_MATCH"
                )
                record["feature_transfer_status"] = (
                    "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "No exact candidate-residue match was identified "
                    "among the seven tested bovine UniProtKB records"
                )

            elif h5_pos in ST3GAL4_MAP["reference_divergent_positions"]:
                record["mapping_status"] = "REFERENCE_DIVERGENT"
                record["feature_transfer_status"] = (
                    "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "Every tested external record that maps this site "
                    "carries S rather than the H5 reference residue R; "
                    "direct residue-specific transfer blocked"
                )

            elif h5_pos in ST3GAL4_MAP["stop_lost_positions"]:
                record["mapping_status"] = "STOP_LOST_SPECIAL_CASE"
                record["feature_transfer_status"] = (
                    "NO_RESIDUE_LEVEL_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "H5 position 396 is the terminal stop following "
                    "the 395-aa reference protein, not an ordinary "
                    "amino-acid residue"
                )

            else:
                raise ValueError(
                    f"Unexpected ST3GAL4 H5 position {h5_pos}"
                )

        elif gene == "CMAH":
            coord = row["coordinate_system"]
            h5_pos = int(row["protein_position"])

            if coord not in CMAH_MAP:
                raise ValueError(
                    f"Unexpected CMAH coordinate system: {coord}"
                )

            m = CMAH_MAP[coord]

            if coord == "Ensembl110":
                record["external_database"] = m["database"]
                record["external_accession"] = m["accession"]

                if h5_pos in m["positions"]:
                    external_pos = m["positions"][h5_pos]

                    record["external_position"] = str(external_pos)
                    record["external_aa"] = row["reference_aa"]
                    record["mapping_status"] = (
                        "ALIGNMENT_EXACT_RESIDUE_MATCH"
                    )
                    record["feature_transfer_status"] = (
                        "MAPPED_RESIDUE_TRANSFER_ALLOWED"
                    )
                    record["mapping_evidence"] = (
                        "MAFFT reconciliation of Ensembl110 CMAH to "
                        "bovine UniProtKB A0AAA9S174; candidate residue "
                        "maps to identical amino acid"
                    )

                elif h5_pos in m["unresolved_positions"]:
                    record["mapping_status"] = (
                        "UNRESOLVED_ISOFORM_CONTEXT"
                    )
                    record["feature_transfer_status"] = (
                        "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                    )
                    record["mapping_evidence"] = (
                        "Candidate lies in divergent C-terminal region "
                        "between Ensembl110 CMAH and A0AAA9S174; "
                        "residue-specific UniProt transfer blocked"
                    )

                else:
                    raise ValueError(
                        f"Unexpected Ensembl CMAH position {h5_pos}"
                    )

            elif coord == "RefSeq_X2":
                if h5_pos not in m["positions"]:
                    raise ValueError(
                        f"Unexpected RefSeq-X2 CMAH position {h5_pos}"
                    )

                record["external_database"] = m["database"]
                record["external_accession"] = m["accession"]
                record["external_position"] = str(h5_pos)
                record["external_aa"] = row["reference_aa"]
                record["mapping_status"] = (
                    "AUTHORITATIVE_REFSEQ_SEQUENCE"
                )
                record["feature_transfer_status"] = (
                    "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "Candidate defined directly on RefSeq XP_024839792.1; "
                    "exact sequence independently represented by UniParc "
                    f'{m["uniparc"]}; no active exact Bos taurus '
                    "UniProtKB record identified"
                )

        elif gene == "NEU1":
            m = NEU1_MAP
            h5_pos = int(row["protein_position"])

            record["external_database"] = m["database"]
            record["external_accession"] = m["accession"]

            if h5_pos in m["positions"]:
                external_pos = m["positions"][h5_pos]

                record["external_position"] = str(external_pos)
                record["external_aa"] = row["reference_aa"]
                record["mapping_status"] = (
                    "ALIGNMENT_EXACT_RESIDUE_MATCH"
                )
                record["feature_transfer_status"] = (
                    "MAPPED_RESIDUE_TRANSFER_ALLOWED"
                )
                record["mapping_evidence"] = (
                    "MAFFT reconciliation to reviewed bovine "
                    "UniProtKB/Swiss-Prot A6BMK7; candidate residue "
                    "maps to identical amino acid"
                )

            elif h5_pos in m["unmapped_gap_positions"]:
                record["mapping_status"] = "UNMAPPED_GAP"
                record["feature_transfer_status"] = (
                    "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"
                )
                record["mapping_evidence"] = (
                    "H5 residue maps to a gap in reviewed bovine "
                    "UniProtKB/Swiss-Prot A6BMK7; residue-specific "
                    "UniProt feature transfer blocked"
                )

            else:
                raise ValueError(
                    f"Unexpected NEU1 H5 position {h5_pos}"
                )

        elif gene in ALIGNMENT_RESIDUE_MAP:
            m = ALIGNMENT_RESIDUE_MAP[gene]
            h5_pos = int(row["protein_position"])

            if h5_pos not in m["positions"]:
                raise ValueError(
                    f"No validated alignment mapping for "
                    f"{gene} H5 position {h5_pos}"
                )

            external_pos = m["positions"][h5_pos]

            record["external_database"] = m["database"]
            record["external_accession"] = m["accession"]
            record["external_position"] = str(external_pos)
            record["external_aa"] = row["reference_aa"]
            record["mapping_status"] = "ALIGNMENT_EXACT_RESIDUE_MATCH"
            record["feature_transfer_status"] = (
                "MAPPED_RESIDUE_TRANSFER_ALLOWED"
            )
            record["mapping_evidence"] = (
                "MAFFT pairwise reconciliation; candidate residue "
                "mapped to identical amino acid in external reference"
            )

        out.append(record)

    fields = list(out[0].keys())

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(out)

    print("Wrote:", OUTPUT)
    print("Candidates:", len(out))

    from collections import Counter

    status_counts = Counter(r["mapping_status"] for r in out)

    print()
    print("Mapping status counts:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    total_classified = sum(status_counts.values())
    pending = status_counts.get("PENDING_RECONCILIATION", 0)

    print()
    print("Total classified:", total_classified)
    print("PENDING_RECONCILIATION:", pending)

    if total_classified != 62:
        raise ValueError(
            f"Expected 62 classified candidates, found {total_classified}"
        )

    if pending != 0:
        raise ValueError(
            f"{pending} candidates remain pending reconciliation"
        )

    print("Classification completeness: PASS")


if __name__ == "__main__":
    main()
