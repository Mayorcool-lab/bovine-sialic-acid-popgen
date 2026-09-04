#!/usr/bin/env python3
"""
Figure 2 source data generation.

Derives population- and breed-level frequency data directly from the frozen
Phase H outputs:
  - results/phase_H_62_variant_population_frequencies.tsv
  - results/phase_H_62_variant_breed_frequencies.tsv
  - results/H_integration/H2_H6_H3_62_candidate_complete_evidence_matrix.tsv

Publication-facing variant labels use compact protein notation
(e.g. NEU3 Y136H, CMAS E128*, ST3GAL4 *396L), never slash notation.

Outputs:
  results/figures/source_data/Figure2_panelA_population.tsv
  results/figures/source_data/Figure2_panelB_breed_heatmap.tsv
  results/figures/source_data/Supplementary_Table_S1_candidate_breed_frequencies.tsv
"""

import csv
from collections import Counter

POP_FILE = "results/phase_H_62_variant_population_frequencies.tsv"
BREED_FILE = "results/phase_H_62_variant_breed_frequencies.tsv"
EVIDENCE_FILE = "results/H_integration/H2_H6_H3_62_candidate_complete_evidence_matrix.tsv"

OUT_A = "results/figures/source_data/Figure2_panelA_population.tsv"
OUT_B = "results/figures/source_data/Figure2_panelB_breed_heatmap.tsv"
OUT_S1 = "results/figures/source_data/Supplementary_Table_S1_candidate_breed_frequencies.tsv"


def candidate_id(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}>{alt}"


def compact_protein_change(hgvsp, ref_aa, position):
    """
    Convert repository HGVSp representation to compact protein notation.

    Examples:
      Y/H + Y + 136 -> Y136H
      E/* + E + 128 -> E128*
      */L + * + 396 -> *396L
      N386S          -> N386S
    """
    if "/" in hgvsp:
        left, right = hgvsp.split("/", 1)
        if left != ref_aa:
            raise RuntimeError(
                f"HGVSp/reference-AA mismatch: HGVSp={hgvsp}, ref={ref_aa}, pos={position}"
            )
        return f"{left}{position}{right}"

    # RefSeq-resolved CMAH candidates are already compact.
    return hgvsp


# ---------------------------------------------------------------------
# Build authoritative candidate -> publication label map
# ---------------------------------------------------------------------

label_map = {}
gene_map = {}

with open(EVIDENCE_FILE, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

if len(rows) != 62:
    raise RuntimeError(f"Expected 62 candidates in evidence matrix, got {len(rows)}")

for r in rows:
    cid = r["candidate_id"]
    label = compact_protein_change(
        r["HGVSp"],
        r["H5_reference_aa"],
        r["H5_protein_position"],
    )

    if "/" in label:
        raise RuntimeError(f"Slash nomenclature remained in publication label: {cid} -> {label}")

    label_map[cid] = label
    gene_map[cid] = r["SYMBOL"]

if len(label_map) != 62:
    raise RuntimeError(f"Expected 62 unique candidate labels, got {len(label_map)}")


# ---------------------------------------------------------------------
# Panel A: population-level status and allele frequencies
# ---------------------------------------------------------------------

with open(POP_FILE, newline="") as f:
    pop_rows = list(csv.DictReader(f, delimiter="\t"))

if len(pop_rows) != 62:
    raise RuntimeError(f"Expected 62 population-frequency rows, got {len(pop_rows)}")

status_pairs = Counter((r["SA_status"], r["EU_status"]) for r in pop_rows)

expected_status_pairs = {
    ("PRESENT", "PRESENT"): 31,
    ("PRESENT", "NOT_IN_SOURCE"): 21,
    ("NOT_IN_SOURCE", "PRESENT"): 6,
    ("ABSENT_AFTER_QC", "PRESENT"): 4,
}

if dict(status_pairs) != expected_status_pairs:
    raise RuntimeError(
        f"Unexpected population status distribution: {dict(status_pairs)}"
    )

with open(OUT_A, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "candidate_id",
        "gene",
        "protein_change",
        "variant_label",
        "chromosome",
        "position",
        "REF",
        "ALT",
        "SA_status",
        "SA_AF",
        "SA_N_genotyped",
        "EU_status",
        "EU_AF",
        "EU_N_genotyped",
        "absolute_delta_AF",
    ])

    for r in pop_rows:
        cid = candidate_id(r["CHROM"], r["POS"], r["REF"], r["ALT"])

        if cid not in label_map:
            raise RuntimeError(f"Population candidate missing from evidence matrix: {cid}")

        delta = ""
        if r["SA_status"] == "PRESENT" and r["EU_status"] == "PRESENT":
            delta = f"{abs(float(r['SA_AF']) - float(r['EU_AF'])):.6f}"

        w.writerow([
            cid,
            r["SYMBOL"],
            label_map[cid],
            f"{r['SYMBOL']} {label_map[cid]}",
            r["CHROM"],
            r["POS"],
            r["REF"],
            r["ALT"],
            r["SA_status"],
            r["SA_AF"],
            r["SA_N_genotyped"],
            r["EU_status"],
            r["EU_AF"],
            r["EU_N_genotyped"],
            delta,
        ])

print(f"Saved: {OUT_A} ({len(pop_rows)} candidates)")
print(
    "Panel A status counts:",
    "31 shared;",
    "21 SA-present/EU-source-unavailable;",
    "6 EU-present/SA-source-unavailable;",
    "4 EU-present/SA-sub-threshold"
)


# ---------------------------------------------------------------------
# Panel B + Supplementary Table S1: breed-level frequencies
# ---------------------------------------------------------------------

with open(BREED_FILE, newline="") as f:
    breed_rows = list(csv.DictReader(f, delimiter="\t"))

# Strip CR if present in final field or other text fields.
for r in breed_rows:
    for k, v in r.items():
        if isinstance(v, str):
            r[k] = v.rstrip("\r")

if len(breed_rows) != 62 * 16:
    raise RuntimeError(
        f"Expected 992 candidate-breed rows, got {len(breed_rows)}"
    )

status_counts = Counter(r["STATUS"] for r in breed_rows)
expected_breed_status = {
    "PRESENT": 700,
    "ABSENT_AFTER_QC": 16,
    "NOT_DETECTED_IN_SOURCE": 276,
}

if dict(status_counts) != expected_breed_status:
    raise RuntimeError(
        f"Unexpected breed-level status distribution: {dict(status_counts)}"
    )

nflag_counts = Counter(r["N_FLAG"] for r in breed_rows)
expected_nflags = {
    "OK": 634,
    "VERY_SMALL_N": 82,
    "NA": 276,
}

if dict(nflag_counts) != expected_nflags:
    raise RuntimeError(
        f"Unexpected breed-level N_FLAG distribution: {dict(nflag_counts)}"
    )

breeds = sorted(set(r["BREED"] for r in breed_rows))
if len(breeds) != 16:
    raise RuntimeError(f"Expected 16 breeds, got {len(breeds)}: {breeds}")

# Figure 2B plotting source.
with open(OUT_B, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "candidate_id",
        "gene",
        "protein_change",
        "variant_label",
        "cohort",
        "breed",
        "status",
        "allele_frequency",
        "callable_N",
        "REF_count",
        "ALT_count",
        "sample_size_flag",
    ])

    for r in breed_rows:
        cid = candidate_id(r["CHROM"], r["POS"], r["REF"], r["ALT"])

        if cid not in label_map:
            raise RuntimeError(f"Breed candidate missing from evidence matrix: {cid}")

        w.writerow([
            cid,
            r["SYMBOL"],
            label_map[cid],
            f"{r['SYMBOL']} {label_map[cid]}",
            r["COHORT"],
            r["BREED"],
            r["STATUS"],
            r["AF"],
            r["N_CALLABLE"],
            r["REF_COUNT"],
            r["ALT_COUNT"],
            r["N_FLAG"],
        ])

print(f"Saved: {OUT_B} ({len(breed_rows)} candidate-breed combinations)")


# Publication-facing Supplementary Table S1.
with open(OUT_S1, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow([
        "Gene",
        "Protein_change",
        "Chromosome",
        "Position",
        "REF",
        "ALT",
        "Consequence",
        "Cohort",
        "Breed",
        "Status",
        "Callable_N",
        "REF_count",
        "ALT_count",
        "Allele_frequency",
        "Sample_size_flag",
    ])

    for r in breed_rows:
        cid = candidate_id(r["CHROM"], r["POS"], r["REF"], r["ALT"])

        w.writerow([
            r["SYMBOL"],
            label_map[cid],
            r["CHROM"],
            r["POS"],
            r["REF"],
            r["ALT"],
            r["Consequence"],
            r["COHORT"],
            r["BREED"],
            r["STATUS"],
            r["N_CALLABLE"],
            r["REF_COUNT"],
            r["ALT_COUNT"],
            r["AF"],
            r["N_FLAG"],
        ])

print(f"Saved: {OUT_S1} ({len(breed_rows)} rows)")
print(
    "Breed-level status counts:",
    "700 PRESENT;",
    "16 ABSENT_AFTER_QC;",
    "276 NOT_DETECTED_IN_SOURCE"
)
print(
    "Sample-size flags:",
    "634 OK;",
    "82 VERY_SMALL_N;",
    "276 NA"
)

print("\nFIGURE 2 SOURCE DATA GENERATION: PASS")
