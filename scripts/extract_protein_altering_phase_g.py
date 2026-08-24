#!/usr/bin/env python3
"""
Phase G: extract final protein-altering/stop variants from the corrected,
target-gene-filtered VEP annotation. Produces the definitive 58-variant
candidate set used in the preliminary manuscript.
"""

import csv
from collections import Counter

interesting = []

with open("results/vep/phase_g_final_target_genes.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        c = row["Consequence"]
        if "missense_variant" in c or "stop_gained" in c or "stop_lost" in c:
            interesting.append(row)

print("=== PROTEIN-ALTERING / STOP VARIANTS (corrected, filtered dataset) ===")
print("Total:", len(interesting))

print("\n=== BY GENE ===")
counts = Counter(r["SYMBOL"] for r in interesting)
for gene, n in counts.most_common():
    print(f"{gene}\t{n}")

print("\n=== INDIVIDUAL VARIANTS ===")
for r in interesting:
    print(r["CHROM"], r["POS"], f'{r["REF"]}>{r["ALT"]}', r["SYMBOL"],
          r["Consequence"], r["IMPACT"], r["HGVSp"], r["Amino_acids"], sep="\t")

with open("results/vep/phase_g_protein_altering_final.tsv", "w", newline="") as out:
    if interesting:
        writer = csv.DictWriter(out, fieldnames=interesting[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(interesting)

print(f"\nSaved: results/vep/phase_g_protein_altering_final.tsv ({len(interesting)} variants)")
