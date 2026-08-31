#!/usr/bin/env python3

import csv
from pathlib import Path

H2 = Path("results/phase_H_62_variant_population_frequencies.tsv")
H4 = Path("results/H4/H4_complete_62_variant_evidence.tsv")
H5 = Path("results/H5/H5_62_variant_conservation_evidence.tsv")
H6 = Path("results/H6/H6.2/processed/H6.2_complete_62_candidate_functional_context.tsv")

OUT = Path("results/H_integration/H2_H6_62_candidate_evidence_matrix.tsv")

KEY = ["CHROM", "POS", "REF", "ALT"]

def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def key(row):
    return tuple(str(row[x]) for x in KEY)

h2 = read_tsv(H2)
h4 = read_tsv(H4)
h5 = read_tsv(H5)
h6 = read_tsv(H6)

maps = {
    "H4": {key(r): r for r in h4},
    "H5": {key(r): r for r in h5},
    "H6": {key(r): r for r in h6},
}

# Explicit fields only: no ranking, no derived score, no evidence collapsing.
fieldnames = [
    # Identity
    "candidate_id",
    "SYMBOL", "CHROM", "POS", "REF", "ALT",
    "HGVSp", "Consequence", "IMPACT",

    # H2 population evidence
    "SA_status", "SA_AF", "SA_N_genotyped", "SA_het", "SA_homalt",
    "EU_status", "EU_AF", "EU_N_genotyped", "EU_het", "EU_homalt",

    # H4 functional prediction
    "SIFT_status", "SIFT_prediction", "SIFT_score",

    # H5 conservation
    "H5_coordinate_system",
    "H5_reference_transcript",
    "H5_reference_protein",
    "H5_protein_position",
    "H5_reference_aa",
    "H5_comparative_total",
    "H5_comparative_non_gap",
    "H5_coverage_status",
    "H5_rate4site_status",
    "H5_rate4site_score",
    "H5_rate4site_qq_interval",
    "H5_rate4site_std",

    # H6 mapping / context confidence
    "H6_external_database",
    "H6_external_accession",
    "H6_external_position",
    "H6_external_aa",
    "H6_mapping_status",
    "H6_feature_transfer_status",
    "H6_branch",
    "H6_context_level",

    # H6 residue/domain evidence
    "H6_uniprot_overlap",
    "H6_uniprot_feature_types",
    "H6_interpro_overlap",
    "H6_interpro_entry_types",
    "H6_interpro_source_databases",
    "H6_interpro_entry_accessions",

    # ST3GAL4 reconciliation
    "H6_st3gal4_references_reconciled",
    "H6_st3gal4_all_seven_overlap_entries",
    "H6_st3gal4_partial_overlap_entries",
    "H6_st3gal4_cd23982_all_seven_overlap",

    # CMAH RefSeq X2 evidence
    "H6_cmah_x2_annotation_source",
    "H6_cmah_x2_interproscan_version",
    "H6_cmah_x2_signature_accession",
    "H6_cmah_x2_signature_name",
    "H6_cmah_x2_entry_accession",
    "H6_cmah_x2_entry_name",
    "H6_cmah_x2_interval",

    # Human-readable context
    "H6_protein_context",
    "H6_context_evidence",
    "H6_residue_specific_interpretation",
]

outrows = []

for r2 in h2:
    k = key(r2)

    if k not in maps["H4"] or k not in maps["H5"] or k not in maps["H6"]:
        raise RuntimeError(f"Missing integration layer for {k}")

    r4 = maps["H4"][k]
    r5 = maps["H5"][k]
    r6 = maps["H6"][k]

    # Frozen identity cross-check.
    symbols = {r2["SYMBOL"], r4["SYMBOL"], r5["SYMBOL"], r6["SYMBOL"]}
    if len(symbols) != 1:
        raise RuntimeError(f"SYMBOL mismatch for {k}: {symbols}")

    if r2["Consequence"] != r4["Consequence"] or r2["Consequence"] != r5["Consequence"]:
        raise RuntimeError(f"Consequence mismatch for {k}")

    out = {
        "candidate_id": r6["candidate_id"],

        "SYMBOL": r2["SYMBOL"],
        "CHROM": r2["CHROM"],
        "POS": r2["POS"],
        "REF": r2["REF"],
        "ALT": r2["ALT"],
        "HGVSp": r4["HGVSp"],
        "Consequence": r2["Consequence"],
        "IMPACT": r4["IMPACT"],

        "SA_status": r2["SA_status"],
        "SA_AF": r2["SA_AF"],
        "SA_N_genotyped": r2["SA_N_genotyped"],
        "SA_het": r2["SA_het"],
        "SA_homalt": r2["SA_homalt"],
        "EU_status": r2["EU_status"],
        "EU_AF": r2["EU_AF"],
        "EU_N_genotyped": r2["EU_N_genotyped"],
        "EU_het": r2["EU_het"],
        "EU_homalt": r2["EU_homalt"],

        "SIFT_status": r4["SIFT_status"],
        "SIFT_prediction": r4["SIFT_prediction"],
        "SIFT_score": r4["SIFT_score"],

        "H5_coordinate_system": r5["coordinate_system"],
        "H5_reference_transcript": r5["reference_transcript"],
        "H5_reference_protein": r5["reference_protein"],
        "H5_protein_position": r5["protein_position"],
        "H5_reference_aa": r5["reference_aa"],
        "H5_comparative_total": r5["comparative_total"],
        "H5_comparative_non_gap": r5["comparative_non_gap"],
        "H5_coverage_status": r5["coverage_status"],
        "H5_rate4site_status": r5["rate4site_status"],
        "H5_rate4site_score": r5["rate4site_score"],
        "H5_rate4site_qq_interval": r5["rate4site_qq_interval"],
        "H5_rate4site_std": r5["rate4site_std"],

        "H6_external_database": r6["external_database"],
        "H6_external_accession": r6["external_accession"],
        "H6_external_position": r6["external_position"],
        "H6_external_aa": r6["external_aa"],
        "H6_mapping_status": r6["mapping_status"],
        "H6_feature_transfer_status": r6["feature_transfer_status"],
        "H6_branch": r6["H6_2_branch"],
        "H6_context_level": r6["context_level"],

        "H6_uniprot_overlap": r6["uniprot_overlap"],
        "H6_uniprot_feature_types": r6["uniprot_feature_types"],
        "H6_interpro_overlap": r6["interpro_overlap"],
        "H6_interpro_entry_types": r6["interpro_entry_types"],
        "H6_interpro_source_databases": r6["interpro_source_databases"],
        "H6_interpro_entry_accessions": r6["interpro_entry_accessions"],

        "H6_st3gal4_references_reconciled": r6["st3gal4_references_reconciled"],
        "H6_st3gal4_all_seven_overlap_entries": r6["st3gal4_all_seven_overlap_entries"],
        "H6_st3gal4_partial_overlap_entries": r6["st3gal4_partial_overlap_entries"],
        "H6_st3gal4_cd23982_all_seven_overlap": r6["st3gal4_cd23982_all_seven_overlap"],

        "H6_cmah_x2_annotation_source": r6["cmah_x2_annotation_source"],
        "H6_cmah_x2_interproscan_version": r6["cmah_x2_interproscan_version"],
        "H6_cmah_x2_signature_accession": r6["cmah_x2_signature_accession"],
        "H6_cmah_x2_signature_name": r6["cmah_x2_signature_name"],
        "H6_cmah_x2_entry_accession": r6["cmah_x2_entry_accession"],
        "H6_cmah_x2_entry_name": r6["cmah_x2_entry_name"],
        "H6_cmah_x2_interval": r6["cmah_x2_interval"],

        "H6_protein_context": r6["protein_context"],
        "H6_context_evidence": r6["context_evidence"],
        "H6_residue_specific_interpretation": r6["residue_specific_interpretation"],
    }

    outrows.append(out)

if len(outrows) != 62:
    raise RuntimeError(f"Expected 62 output rows, got {len(outrows)}")

if len({(r["CHROM"], r["POS"], r["REF"], r["ALT"]) for r in outrows}) != 62:
    raise RuntimeError("Output variant keys are not unique")

OUT.parent.mkdir(parents=True, exist_ok=True)

with OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(outrows)

print(f"Output: {OUT}")
print(f"Rows: {len(outrows)}")
print(f"Columns: {len(fieldnames)}")
print(f"Unique variant keys: {len({(r['CHROM'], r['POS'], r['REF'], r['ALT']) for r in outrows})}")
print("Ranking columns added: 0")
print("Composite scores added: 0")
print("H2-H6 EVIDENCE MATRIX BUILD: PASS")
