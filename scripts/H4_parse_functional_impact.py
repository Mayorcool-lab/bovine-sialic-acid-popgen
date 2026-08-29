#!/usr/bin/env python3

from collections import Counter
import csv
from pathlib import Path

PHASE_G = Path("results/vep/phase_g_protein_altering_FINAL.tsv")
SIFT_VCF = Path("results/H4/missense_sift_annotated.vcf")

# Load authoritative Phase G candidate set.
with PHASE_G.open() as f:
    phase_g_rows = list(csv.DictReader(f, delimiter="\t"))

assert len(phase_g_rows) == 62, \
    f"Expected 62 Phase G candidates, found {len(phase_g_rows)}"

missense_rows = [
    r for r in phase_g_rows
    if "missense" in r["Consequence"]
]

assert len(missense_rows) == 60, \
    f"Expected 60 missense-containing candidates, found {len(missense_rows)}"

print("Phase G candidates:", len(phase_g_rows))
print("Missense-containing candidates:", len(missense_rows))
print("Non-missense candidates:", len(phase_g_rows) - len(missense_rows))

# Parse H4 VEP/SIFT output by exact genomic allele.
sift_records = {}

with SIFT_VCF.open() as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        key = (fields[0], fields[1], fields[3], fields[4])

        assert key not in sift_records, \
            f"Duplicate H4 SIFT allele encountered: {key}"

        sift_records[key] = fields[7]

# Build exact-allele keys for the 60 Phase G missense-containing candidates.
phase_g_missense_keys = {
    (r["CHROM"], r["POS"], r["REF"], r["ALT"])
    for r in missense_rows
}

sift_keys = set(sift_records)

assert len(sift_records) == 60, \
    f"Expected 60 H4 SIFT VCF records, found {len(sift_records)}"

assert sift_keys == phase_g_missense_keys, (
    "Exact allele reconciliation failed.\n"
    f"Missing from SIFT: {sorted(phase_g_missense_keys - sift_keys)}\n"
    f"Unexpected in SIFT: {sorted(sift_keys - phase_g_missense_keys)}"
)

print("SIFT VCF records:", len(sift_records))
print("Exact Phase G allele matches:", len(sift_keys & phase_g_missense_keys))
print("Missing alleles:", len(phase_g_missense_keys - sift_keys))
print("Unexpected alleles:", len(sift_keys - phase_g_missense_keys))
# Extract SIFT from the exact transcript selected in Phase G.
# This avoids ambiguity where more than one VEP annotation may be marked
# CANONICAL=YES and keeps H4 consistent with Phase G consequence identity.
phase_g_lookup = {
    (r["CHROM"], r["POS"], r["REF"], r["ALT"]): r
    for r in missense_rows
}

transcript_sift = {}

for key, info in sift_records.items():
    phase_g_row = phase_g_lookup[key]
    expected_feature = phase_g_row["Feature"]

    # Manually rescued CMAH variants have no Phase G transcript Feature.
    # In the current Ensembl 110 ARS-UCD1.2 cache, these exact alleles are
    # annotated as intergenic and therefore have no native bovine SIFT score.
    # Do not treat blank Feature-to-blank Feature as transcript matching.
    if not expected_feature:
        assert phase_g_row["SYMBOL"] == "CMAH", \
            f"Unexpected blank Phase G Feature for allele: {key}"
        assert phase_g_row["annotation_source"] != "VEP", \
            f"Blank Feature unexpectedly associated with VEP annotation: {key}"

        transcript_sift[key] = ("SIFT_UNAVAILABLE", "", "")
        continue

    csq_value = None

    for item in info.split(";"):
        if item.startswith("CSQ="):
            csq_value = item[4:]
            break

    assert csq_value is not None, \
        f"No CSQ annotation found for allele: {key}"

    matching_annotations = []

    for annotation in csq_value.split(","):
        x = annotation.split("|")
        x += [""] * (25 - len(x))

        # CSQ schema:
        # field 7 = Feature, field 25 = SIFT (1-based numbering)
        if x[6] == expected_feature:
            matching_annotations.append(x)

    assert len(matching_annotations) == 1, (
        f"Expected exactly one H4 annotation matching Phase G transcript "
        f"{expected_feature} for allele {key}; "
        f"found {len(matching_annotations)}"
    )

    x = matching_annotations[0]
    sift = x[24]

    if not sift:
        transcript_sift[key] = ("SIFT_UNAVAILABLE", "", "")
        continue

    if "(" in sift and sift.endswith(")"):
        prediction, score = sift[:-1].split("(", 1)
    else:
        prediction, score = sift, ""

    transcript_sift[key] = ("SIFT_AVAILABLE", prediction, score)

# Summarise transcript-matched SIFT coverage and predictions.
status_counts = {}
prediction_counts = {}

for status, prediction, score in transcript_sift.values():
    status_counts[status] = status_counts.get(status, 0) + 1

    if prediction:
        prediction_counts[prediction] = \
            prediction_counts.get(prediction, 0) + 1

# These assertions freeze the validated H4.2 result.
assert status_counts == {
    "SIFT_AVAILABLE": 56,
    "SIFT_UNAVAILABLE": 4,
}, f"Unexpected SIFT coverage: {status_counts}"

assert prediction_counts == {
    "deleterious": 6,
    "deleterious_low_confidence": 6,
    "tolerated": 38,
    "tolerated_low_confidence": 6,
}, f"Unexpected SIFT prediction distribution: {prediction_counts}"

expected_unavailable = {
    ("23", "32739285", "A", "G"),
    ("23", "32739292", "C", "G"),
    ("23", "32743918", "C", "T"),
    ("23", "32745866", "T", "A"),
}

observed_unavailable = {
    key
    for key, (status, prediction, score) in transcript_sift.items()
    if status == "SIFT_UNAVAILABLE"
}

assert observed_unavailable == expected_unavailable, (
    "Unexpected SIFT-unavailable allele set.\n"
    f"Observed: {sorted(observed_unavailable)}\n"
    f"Expected: {sorted(expected_unavailable)}"
)

print("Transcript-matched SIFT status:")
for status in sorted(status_counts):
    print(f"  {status}: {status_counts[status]}")

print("Transcript-matched SIFT predictions:")
for prediction in sorted(prediction_counts):
    print(f"  {prediction}: {prediction_counts[prediction]}")

print("SIFT-unavailable alleles:")
for key in sorted(observed_unavailable):
    r = phase_g_lookup[key]
    print(
        f"  {r['SYMBOL']} "
        f"{r['CHROM']}:{r['POS']}:{r['REF']}>{r['ALT']} "
        f"Feature={r['Feature']}"
    )


# Build the complete 62-variant H4 evidence table: 60 missense-containing
# variants (with SIFT results) + 2 HIGH-impact stop_gained/stop_lost
# variants (SIFT not applicable by design, not a coverage gap).
final_rows = []
for r in phase_g_rows:
    key = (r["CHROM"], r["POS"], r["REF"], r["ALT"])
    if key in phase_g_lookup:
        status, pred, score = transcript_sift[key]
        sift_status = "SIFT_AVAILABLE" if pred else status
    else:
        sift_status, pred, score = "NOT_APPLICABLE_HIGH_IMPACT", "", ""

    final_rows.append({
        "SYMBOL": r["SYMBOL"], "CHROM": r["CHROM"], "POS": r["POS"],
        "REF": r["REF"], "ALT": r["ALT"],
        "HGVSp": r["HGVSp"] or r["Amino_acids"],
        "Consequence": r["Consequence"], "IMPACT": r["IMPACT"],
        "SIFT_status": sift_status, "SIFT_prediction": pred, "SIFT_score": score
    })

assert len(final_rows) == 62, f"Expected 62 final rows, got {len(final_rows)}"

FINAL_TABLE = Path("results/H4/H4_complete_62_variant_evidence.tsv")
final_columns = ["SYMBOL","CHROM","POS","REF","ALT","HGVSp","Consequence","IMPACT","SIFT_status","SIFT_prediction","SIFT_score"]
with FINAL_TABLE.open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=final_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(final_rows)

print()
print("=== Complete 62-variant H4 evidence table ===")
print(Counter(r["SIFT_status"] for r in final_rows))
print(f"Saved: {FINAL_TABLE}")
