#!/usr/bin/env python3
"""
Phase G8: combine the clean 13-gene VEP-assigned protein-altering
variants with the manually-annotated CMAH rescued-region coding
variants (G9) into the final, definitive protein-altering candidate
table for the corrected 10,261-variant dataset.

IMPORTANT: the manual rescue loop filters on Consequence exactly like
the VEP loop does. I515I (synonymous) must NOT appear in this table,
only the 4 manually-rescued MISSENSE variants (N386S, I388M, P424L,
F512Y) belong here.
"""

import csv
from collections import Counter

interesting = []

with open("results/vep/phase_g_final_target_genes_corrected.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        c = row["Consequence"]
        if "missense_variant" in c or "stop_gained" in c or "stop_lost" in c:
            row["annotation_source"] = "VEP"
            interesting.append(row)

with open("results/vep/cmah_rescued_region_override.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["annotation_source"] == "manual_RefSeq_rescue_G9":
            c = row["Consequence"]
            if "missense_variant" in c or "stop_gained" in c or "stop_lost" in c:
                interesting.append(row)

print("=== FINAL PROTEIN-ALTERING / STOP VARIANTS (corrected, complete dataset) ===")
print("Total:", len(interesting))
print()
print("=== BY GENE ===")
counts = Counter(r["SYMBOL"] for r in interesting)
for gene, n in counts.most_common():
    print(f"{gene}\t{n}")

print()
print("=== BY SOURCE ===")
source_counts = Counter(r["annotation_source"] for r in interesting)
for src, n in source_counts.most_common():
    print(f"{src}\t{n}")

columns = ["CHROM","POS","REF","ALT","SYMBOL","Gene","Feature","Consequence","IMPACT",
           "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","Protein_position","Amino_acids",
           "Codons","Existing_variation","CANONICAL","annotation_source"]

with open("results/vep/phase_g_protein_altering_FINAL.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    for row in interesting:
        clean_row = {col: row.get(col, "") for col in columns}
        writer.writerow(clean_row)

print()
print(f"Saved: results/vep/phase_g_protein_altering_FINAL.tsv ({len(interesting)} variants)")

# Self-check against known facts
manual_missense_expected = {32739285, 32739292, 32743918, 32745866}
manual_in_output = {int(r["POS"]) for r in interesting if r["annotation_source"] == "manual_RefSeq_rescue_G9"}
print()
print("=== SELF-CHECK ===")
print(f"Expected manual rescue positions (4 missense): {sorted(manual_missense_expected)}")
print(f"Actual manual rescue positions in output: {sorted(manual_in_output)}")
if manual_in_output == manual_missense_expected:
    print("PASS: manual rescue set matches expected 4 missense variants exactly.")
else:
    print("FAIL: mismatch detected, review output before trusting this table.")
