#!/usr/bin/env python3

import csv
from pathlib import Path

INFILE = Path(
    "results/H_integration/"
    "H2_H6_H3_62_candidate_complete_evidence_matrix.tsv"
)

OUTFILE = Path(
    "results/H_integration/"
    "H_candidate_prioritisation_review_62.tsv"
)


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


rows = read_tsv(INFILE)

if len(rows) != 62:
    raise RuntimeError(
        f"Expected 62 candidates, observed {len(rows)}"
    )

if len({r["candidate_id"] for r in rows}) != 62:
    raise RuntimeError("Candidate IDs are not unique")


def observation_pattern(r):
    sa = int(r["South_Africa_observed_OK"])
    eu = int(r["Europe_observed_OK"])

    if sa > 0 and eu > 0:
        return "OBSERVED_BOTH"

    if sa > 0 and eu == 0:
        if int(r["Europe_source_unavailable"]) == 12:
            return "SA_OBSERVED_EU_SOURCE_UNAVAILABLE"
        return "OBSERVED_SA_ONLY_WITHIN_AVAILABLE_DATA"

    if sa == 0 and eu > 0:
        if int(r["South_Africa_source_unavailable"]) == 4:
            return "EU_OBSERVED_SA_SOURCE_UNAVAILABLE"
        return "OBSERVED_EU_ONLY_WITHIN_AVAILABLE_DATA"

    if (
        int(r["South_Africa_source_unavailable"]) == 4
        and int(r["Europe_source_unavailable"]) == 12
    ):
        return "NO_COHORT_SOURCE_AVAILABLE"

    return "NO_ALT_OBSERVED_IN_ADEQUATELY_SAMPLED_BREEDS"


def molecular_class(r):
    if r["Consequence"] == "stop_gained":
        return "STOP_GAINED"

    if r["Consequence"] == "stop_lost":
        return "STOP_LOST"

    if r["SIFT_status"] == "SIFT_UNAVAILABLE":
        return "MISSENSE_SIFT_UNAVAILABLE"

    if r["SIFT_prediction"] == "deleterious":
        return "MISSENSE_SIFT_DELETERIOUS"

    if r["SIFT_prediction"] == "deleterious_low_confidence":
        return "MISSENSE_SIFT_DELETERIOUS_LOW_CONFIDENCE"

    if r["SIFT_prediction"] == "tolerated":
        return "MISSENSE_SIFT_TOLERATED"

    if r["SIFT_prediction"] == "tolerated_low_confidence":
        return "MISSENSE_SIFT_TOLERATED_LOW_CONFIDENCE"

    return "MISSENSE_OTHER"


outrows = []

for r in rows:
    out = {
        "candidate_id": r["candidate_id"],
        "SYMBOL": r["SYMBOL"],
        "CHROM": r["CHROM"],
        "POS": r["POS"],
        "REF": r["REF"],
        "ALT": r["ALT"],
        "HGVSp": r["HGVSp"],
        "Consequence": r["Consequence"],
        "IMPACT": r["IMPACT"],

        "molecular_evidence_class": molecular_class(r),

        "SIFT_status": r["SIFT_status"],
        "SIFT_prediction": r["SIFT_prediction"],
        "SIFT_score": r["SIFT_score"],

        "H5_rate4site_status": r["H5_rate4site_status"],
        "H5_rate4site_score": r["H5_rate4site_score"],
        "H5_comparative_non_gap": r["H5_comparative_non_gap"],
        "H5_coordinate_system": r["H5_coordinate_system"],

        "H6_context_level": r["H6_context_level"],
        "H6_branch": r["H6_branch"],
        "H6_mapping_status": r["H6_mapping_status"],
        "H6_feature_transfer_status": r["H6_feature_transfer_status"],
        "H6_protein_context": r["H6_protein_context"],
        "H6_context_evidence": r["H6_context_evidence"],
        "H6_residue_specific_interpretation":
            r["H6_residue_specific_interpretation"],

        "SA_status": r["SA_status"],
        "SA_AF": r["SA_AF"],
        "SA_N_genotyped": r["SA_N_genotyped"],

        "EU_status": r["EU_status"],
        "EU_AF": r["EU_AF"],
        "EU_N_genotyped": r["EU_N_genotyped"],

        "South_Africa_present_OK":
            r["South_Africa_present_OK"],
        "South_Africa_observed_OK":
            r["South_Africa_observed_OK"],
        "South_Africa_observed_OK_breeds":
            r["South_Africa_observed_OK_breeds"],
        "South_Africa_observed_OK_AF":
            r["South_Africa_observed_OK_AF"],
        "South_Africa_source_unavailable":
            r["South_Africa_source_unavailable"],
        "South_Africa_absent_after_qc":
            r["South_Africa_absent_after_qc"],

        "Europe_present_OK":
            r["Europe_present_OK"],
        "Europe_observed_OK":
            r["Europe_observed_OK"],
        "Europe_observed_OK_breeds":
            r["Europe_observed_OK_breeds"],
        "Europe_observed_OK_AF":
            r["Europe_observed_OK_AF"],
        "Europe_source_unavailable":
            r["Europe_source_unavailable"],
        "Europe_absent_after_qc":
            r["Europe_absent_after_qc"],
        "Europe_observed_very_small_n_AF":
            r["Europe_observed_very_small_n_AF"],

        "breed_observation_pattern":
            observation_pattern(r),

        # Deliberately left blank for later scientific review.
        "practical_genotype_phenotype_informativeness": "",
        "evidence_convergence_rationale": "",
        "evidence_limitations_or_contradictions": "",
        "prioritisation_decision": "",
    }

    outrows.append(out)


fieldnames = list(outrows[0].keys())

for forbidden in [
    "composite_score",
    "priority_score",
    "rank",
]:
    if forbidden in fieldnames:
        raise RuntimeError(
            f"Forbidden prioritisation field present: {forbidden}"
        )


with OUTFILE.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(outrows)


print("Output:", OUTFILE)
print("Rows:", len(outrows))
print("Columns:", len(fieldnames))
print(
    "Unique candidates:",
    len({r["candidate_id"] for r in outrows})
)

print(
    "Blank practical-informativeness fields:",
    sum(
        1 for r in outrows
        if r["practical_genotype_phenotype_informativeness"] == ""
    )
)

print(
    "Blank convergence-rationale fields:",
    sum(
        1 for r in outrows
        if r["evidence_convergence_rationale"] == ""
    )
)

print(
    "Blank prioritisation decisions:",
    sum(
        1 for r in outrows
        if r["prioritisation_decision"] == ""
    )
)

print("Composite score added: 0")
print("Numerical rank added: 0")
print("CANDIDATE PRIORITISATION REVIEW SCAFFOLD: PASS")
