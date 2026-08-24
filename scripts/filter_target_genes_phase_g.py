#!/usr/bin/env python3
"""
Phase G: restrict annotated variants to the 13 real target genes.
A small number of variants (144) were annotated to off-target genes
(KIRREL3, SLC44A4, RARS2, 5S_rRNA) due to genomic proximity to target
genes within the 2kb upstream extraction window. This script separates
those out explicitly rather than silently discarding them.
"""

import csv
from collections import Counter

target_genes = {"GNE","NANS","NANP","CMAS","SLC35A1","CMAH","ST8SIA2","ST8SIA4",
                 "ST6GAL1","ST3GAL1","ST3GAL4","NEU1","NEU3"}

kept = []
excluded = []

with open("results/vep/phase_g_canonical_annotations_fixed.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["SYMBOL"] in target_genes:
            kept.append(row)
        else:
            excluded.append(row)

print(f"Kept (target genes): {len(kept)}")
print(f"Excluded (off-target/no-symbol): {len(excluded)}")
print()
print("Excluded breakdown:")
for gene, n in Counter(r["SYMBOL"] or "[NO SYMBOL]" for r in excluded).most_common():
    print(f"  {gene}: {n}")

columns = list(kept[0].keys())
with open("results/vep/phase_g_final_target_genes.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(kept)

with open("results/vep/phase_g_excluded_offtarget.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(excluded)

print()
print("Saved: results/vep/phase_g_final_target_genes.tsv (the real dataset)")
print("Saved: results/vep/phase_g_excluded_offtarget.tsv (documented exclusions)")
