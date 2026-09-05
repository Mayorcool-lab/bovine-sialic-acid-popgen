#!/usr/bin/env python3

"""
Generate publication-facing source data for Figure 5 and Supplementary Table S2.

Figure 5
--------
Panel A:
    ST3GAL4 multi-reference reconciled candidate context.

Panel B:
    RefSeq-X2 CMAH sequence-based candidate context.

Supplementary Table S2:
    Complete functional-context annotation for all 62 protein-altering
    candidates from the frozen H6.2 evidence table.

This script does not re-derive H6 annotations. It reformats and validates
the frozen H6.2 outputs for publication.
"""

import csv
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[2]

H6_COMPLETE = (
    ROOT
    / "results/H6/H6.2/processed/"
    / "H6.2_complete_62_candidate_functional_context.tsv"
)

ST3GAL4_SOURCE = (
    ROOT
    / "results/H6/H6.2/processed/"
    / "H6.2_ST3GAL4_reconciled_candidate_summary.tsv"
)

CMAH_X2_SOURCE = (
    ROOT
    / "results/H6/H6.2/processed/"
    / "H6.2_CMAH_RefSeq_X2_candidate_overlaps.tsv"
)

OUTDIR = ROOT / "results/figures/source_data"

FIG5A_OUT = OUTDIR / "Figure5_panelA_ST3GAL4_context.tsv"
FIG5B_OUT = OUTDIR / "Figure5_panelB_CMAH_X2_context.tsv"
S2_OUT = OUTDIR / "Supplementary_Table_S2_functional_context.tsv"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------
# Load frozen H6 evidence
# ---------------------------------------------------------------------

complete = read_tsv(H6_COMPLETE)
st3 = read_tsv(ST3GAL4_SOURCE)
cmah = read_tsv(CMAH_X2_SOURCE)

assert len(complete) == 62, (
    f"Expected 62 candidates in H6.2 complete table; found {len(complete)}"
)


# ---------------------------------------------------------------------
# Validate frozen branch/context structure
# ---------------------------------------------------------------------

expected_branches = {
    "STRAIGHTFORWARD": 43,
    "EXCEPTIONAL_OR_SPECIAL_CASE": 10,
    "ST3GAL4_MULTI_REFERENCE_RECONCILED": 5,
    "CMAH_REFSEQ_X2_SEQUENCE_BASED": 4,
}

observed_branches = Counter(r["H6_2_branch"] for r in complete)

assert observed_branches == Counter(expected_branches), (
    "Unexpected H6.2 branch structure:\n"
    f"Observed: {dict(observed_branches)}\n"
    f"Expected: {expected_branches}"
)

expected_context = {
    "RESIDUE_LEVEL_DIRECT_OR_MAPPED": 43,
    "PROTEIN_LEVEL_ONLY": 10,
    "RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED": 5,
    "RESIDUE_LEVEL_SEQUENCE_BASED": 4,
}

observed_context = Counter(r["context_level"] for r in complete)

assert observed_context == Counter(expected_context), (
    "Unexpected context-level structure:\n"
    f"Observed: {dict(observed_context)}\n"
    f"Expected: {expected_context}"
)

straightforward = [
    r for r in complete
    if r["H6_2_branch"] == "STRAIGHTFORWARD"
]

uniprot_counts = Counter(r["uniprot_overlap"] for r in straightforward)

assert uniprot_counts == Counter({"YES": 22, "NO": 21}), (
    "Unexpected UniProt overlap counts in STRAIGHTFORWARD branch:\n"
    f"{dict(uniprot_counts)}"
)


# ---------------------------------------------------------------------
# Protein-change labels
# ---------------------------------------------------------------------
# Publication-facing nomenclature frozen from the candidate evidence.
# Only Figure 5 special-branch candidates require explicit labels here.

protein_changes = {
    # ST3GAL4
    "29:29731343:A>C": "Y190S",
    "29:29732068:C>A": "R203S",
    "29:29733130:A>T": "Y320F",
    "29:29733170:A>G": "I333M",
    "29:29733173:C>A": "H334Q",

    # CMAH RefSeq X2
    "23:32739285:A>G": "N386S",
    "23:32739292:C>G": "I388M",
    "23:32743918:C>T": "P424L",
    "23:32745866:T>A": "F512Y",
}


# ---------------------------------------------------------------------
# Figure 5A: ST3GAL4 multi-reference reconciliation
# ---------------------------------------------------------------------

assert len(st3) == 5, (
    f"Expected 5 reconciled ST3GAL4 candidates; found {len(st3)}"
)

expected_st3_ids = {
    "29:29731343:A>C",
    "29:29732068:C>A",
    "29:29733130:A>T",
    "29:29733170:A>G",
    "29:29733173:C>A",
}

observed_st3_ids = {r["candidate_id"] for r in st3}

assert observed_st3_ids == expected_st3_ids, (
    "Unexpected ST3GAL4 reconciled candidate set:\n"
    f"Observed: {sorted(observed_st3_ids)}"
)

panel_a = []

for r in st3:
    candidate_id = r["candidate_id"]

    assert r["references_reconciled"] == "7"
    assert r["mapping_status"] == "ROBUST_MULTI_REFERENCE_RESIDUE_MATCH"
    assert r["feature_transfer_status"] == "FEATURE_RECONCILIATION_ALLOWED"
    assert r["context_level"] == "RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED"

    panel_a.append({
        "candidate_id": candidate_id,
        "Gene": "ST3GAL4",
        "Protein_change": protein_changes[candidate_id],
        "Protein_position": r["H5_position"],
        "Reference_AA": r["H5_reference_aa"],
        "References_reconciled": r["references_reconciled"],
        "All_seven_overlap_entries": r["all_seven_overlap_entries"],
        "All_seven_overlap_names": r["all_seven_overlap_names"],
        "Partial_overlap_entries": r["partial_overlap_entries"],
        "cd23982_all_seven_overlap": r["cd23982_all_seven_overlap"],
        "Context_level": r["context_level"],
    })

panel_a.sort(key=lambda r: int(r["Protein_position"]))

assert sum(
    r["cd23982_all_seven_overlap"] == "YES"
    for r in panel_a
) == 4

assert [
    r["Protein_change"]
    for r in panel_a
    if r["cd23982_all_seven_overlap"] == "NO"
] == ["Y190S"]

write_tsv(
    FIG5A_OUT,
    panel_a,
    [
        "candidate_id",
        "Gene",
        "Protein_change",
        "Protein_position",
        "Reference_AA",
        "References_reconciled",
        "All_seven_overlap_entries",
        "All_seven_overlap_names",
        "Partial_overlap_entries",
        "cd23982_all_seven_overlap",
        "Context_level",
    ],
)


# ---------------------------------------------------------------------
# Figure 5B: CMAH RefSeq-X2 sequence-based context
# ---------------------------------------------------------------------

assert len(cmah) == 4, (
    f"Expected 4 RefSeq-X2 CMAH candidates; found {len(cmah)}"
)

expected_cmah_ids = {
    "23:32739285:A>G",
    "23:32739292:C>G",
    "23:32743918:C>T",
    "23:32745866:T>A",
}

observed_cmah_ids = {r["candidate_id"] for r in cmah}

assert observed_cmah_ids == expected_cmah_ids, (
    "Unexpected CMAH RefSeq-X2 candidate set:\n"
    f"Observed: {sorted(observed_cmah_ids)}"
)

panel_b = []

for r in cmah:
    candidate_id = r["candidate_id"]

    assert r["reference_protein"] == "XP_024839792.1"
    assert r["annotation_source"] == "InterProScan"
    assert r["interproscan_version"] == "5.78-109.0"
    assert r["library"] == "PANTHER"
    assert r["signature_accession"] == "PTHR46522"
    assert r["entry_accession"] == "IPR027033"
    assert r["interval_start"] == "5"
    assert r["interval_end"] == "576"
    assert r["mapping_status"] == "AUTHORITATIVE_REFSEQ_SEQUENCE"
    assert r["feature_transfer_status"] == "NO_DIRECT_UNIPROT_FEATURE_TRANSFER"

    panel_b.append({
        "candidate_id": candidate_id,
        "Gene": "CMAH",
        "Protein_change": protein_changes[candidate_id],
        "Reference_protein": r["reference_protein"],
        "Protein_position": r["candidate_position"],
        "Reference_AA": r["reference_aa"],
        "Annotation_source": r["annotation_source"],
        "InterProScan_version": r["interproscan_version"],
        "Library": r["library"],
        "Signature_accession": r["signature_accession"],
        "Signature_name": r["signature_name"],
        "InterPro_entry": r["entry_accession"],
        "InterPro_entry_name": r["entry_name"],
        "Family_interval_start": r["interval_start"],
        "Family_interval_end": r["interval_end"],
        "Context_interpretation": (
            "BROAD_FAMILY_OVERLAP_ONLY"
        ),
    })

panel_b.sort(key=lambda r: int(r["Protein_position"]))

write_tsv(
    FIG5B_OUT,
    panel_b,
    [
        "candidate_id",
        "Gene",
        "Protein_change",
        "Reference_protein",
        "Protein_position",
        "Reference_AA",
        "Annotation_source",
        "InterProScan_version",
        "Library",
        "Signature_accession",
        "Signature_name",
        "InterPro_entry",
        "InterPro_entry_name",
        "Family_interval_start",
        "Family_interval_end",
        "Context_interpretation",
    ],
)


# ---------------------------------------------------------------------
# Supplementary Table S2: complete 62-candidate H6 context
# ---------------------------------------------------------------------

s2_fields = [
    "candidate_id",
    "SYMBOL",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "H5_coordinate_system",
    "H5_reference_protein",
    "H5_position",
    "H5_reference_aa",
    "external_database",
    "external_accession",
    "external_position",
    "external_aa",
    "mapping_status",
    "feature_transfer_status",
    "H6_2_branch",
    "context_level",
    "uniprot_overlap",
    "uniprot_feature_types",
    "interpro_overlap",
    "interpro_entry_types",
    "interpro_source_databases",
    "interpro_entry_accessions",
    "st3gal4_references_reconciled",
    "st3gal4_all_seven_overlap_entries",
    "st3gal4_partial_overlap_entries",
    "st3gal4_cd23982_all_seven_overlap",
    "cmah_x2_annotation_source",
    "cmah_x2_interproscan_version",
    "cmah_x2_signature_accession",
    "cmah_x2_signature_name",
    "cmah_x2_entry_accession",
    "cmah_x2_entry_name",
    "cmah_x2_interval",
    "protein_context",
    "context_evidence",
    "residue_specific_interpretation",
]

assert list(complete[0].keys()) == s2_fields, (
    "The frozen H6.2 complete-table schema has changed.\n"
    "Review before regenerating Supplementary Table S2."
)

write_tsv(
    S2_OUT,
    complete,
    s2_fields,
)


# ---------------------------------------------------------------------
# Final validation/report
# ---------------------------------------------------------------------

assert len(panel_a) == 5
assert len(panel_b) == 4
assert len(complete) == 62

print("Figure 5 / Supplementary Table S2 source data generated successfully.")
print()
print("Validated H6.2 branches:")
for key in [
    "STRAIGHTFORWARD",
    "EXCEPTIONAL_OR_SPECIAL_CASE",
    "ST3GAL4_MULTI_REFERENCE_RECONCILED",
    "CMAH_REFSEQ_X2_SEQUENCE_BASED",
]:
    print(f"  {key}: {observed_branches[key]}")

print()
print("Validated context levels:")
for key in [
    "RESIDUE_LEVEL_DIRECT_OR_MAPPED",
    "RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED",
    "RESIDUE_LEVEL_SEQUENCE_BASED",
    "PROTEIN_LEVEL_ONLY",
]:
    print(f"  {key}: {observed_context[key]}")

print()
print("STRAIGHTFORWARD UniProt overlap:")
print(f"  YES: {uniprot_counts['YES']}")
print(f"  NO:  {uniprot_counts['NO']}")

print()
print(f"Figure 5A rows: {len(panel_a)}")
print(f"Figure 5B rows: {len(panel_b)}")
print(f"Supplementary Table S2 rows: {len(complete)}")

print()
print("Outputs:")
print(f"  {FIG5A_OUT.relative_to(ROOT)}")
print(f"  {FIG5B_OUT.relative_to(ROOT)}")
print(f"  {S2_OUT.relative_to(ROOT)}")
