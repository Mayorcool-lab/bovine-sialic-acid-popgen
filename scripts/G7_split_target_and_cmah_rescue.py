#!/usr/bin/env python3
"""
Phase G7: split the corrected canonical annotation table into two
deliberately separate outputs:

1. A clean 13-target-gene table, based purely on VEP's own SYMBOL
   assignment.
2. A CMAH rescued-region override table, capturing all variants in the
   RefSeq-validated CMAH extension (23:32,731,973-32,753,209) that VEP
   could not correctly attribute (Ensembl's gene model ends at
   32,731,972). This preserves 5 manually-annotated coding variants
   (N386S, I388M, P424L, F512Y, I515I; see scripts/
   annotate_cmah_refseq_rescued.py) plus 388 correctly-intronic CMAH
   variants that would otherwise be silently discarded by a naive
   SYMBOL-based filter.
"""

import csv
from collections import Counter

target_genes = {"GNE","NANS","NANP","CMAS","SLC35A1","CMAH","ST8SIA2","ST8SIA4",
                 "ST6GAL1","ST3GAL1","ST3GAL4","NEU1","NEU3"}

manually_annotated = {
    32739285: {"consequence": "missense_variant", "aa_change": "N386S", "codon": 386, "cdna": 1157},
    32739292: {"consequence": "missense_variant", "aa_change": "I388M", "codon": 388, "cdna": 1164},
    32743918: {"consequence": "missense_variant", "aa_change": "P424L", "codon": 424, "cdna": 1271},
    32745866: {"consequence": "missense_variant", "aa_change": "F512Y", "codon": 512, "cdna": 1535},
    32750628: {"consequence": "synonymous_variant", "aa_change": "I515I", "codon": 515, "cdna": 1545},
}

# --- Part 1: clean 13-gene VEP-assigned table ---
kept = []
excluded = []

with open("results/vep/phase_g_canonical_annotations_corrected.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["SYMBOL"] in target_genes:
            kept.append(row)
        else:
            excluded.append(row)

print(f"Kept (target genes, VEP-assigned): {len(kept)}")
print(f"Excluded (off-target or blank symbol): {len(excluded)}")
print()
print("Excluded breakdown:")
for gene, n in Counter(r["SYMBOL"] or "[NO SYMBOL]" for r in excluded).most_common():
    print(f"  {gene}: {n}")

columns = list(kept[0].keys())
with open("results/vep/phase_g_final_target_genes_corrected.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(kept)

with open("results/vep/phase_g_excluded_offtarget_corrected.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(excluded)

print()
print("Saved: results/vep/phase_g_final_target_genes_corrected.tsv")
print("Saved: results/vep/phase_g_excluded_offtarget_corrected.tsv")

# --- Part 2: CMAH rescued-region override table ---
rows = []
with open("results/vep/phase_g_excluded_offtarget_corrected.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        chrom, pos = row["CHROM"], int(row["POS"])
        if chrom == "23" and 32731972 < pos <= 32753209:
            pos_int = int(row["POS"])
            if pos_int in manually_annotated:
                m = manually_annotated[pos_int]
                row["SYMBOL"] = "CMAH"
                row["Gene"] = "ENSBTAG00000003892"
                row["Consequence"] = m["consequence"]
                row["HGVSp"] = m["aa_change"]
                row["HGVSc"] = f"c.{m['cdna']}"
                row["Protein_position"] = str(m["codon"])
                row["annotation_source"] = "manual_RefSeq_rescue_G9"
            else:
                row["SYMBOL"] = "CMAH"
                row["Gene"] = "ENSBTAG00000003892"
                row["Consequence"] = "intron_variant"
                row["annotation_source"] = "RefSeq_rescue_region_noncoding"
            rows.append(row)

columns2 = list(rows[0].keys())
with open("results/vep/cmah_rescued_region_override.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns2, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

coding = sum(1 for r in rows if r["annotation_source"] == "manual_RefSeq_rescue_G9")
noncoding = sum(1 for r in rows if r["annotation_source"] == "RefSeq_rescue_region_noncoding")
print()
print(f"CMAH rescued-region variants recovered: {len(rows)}")
print(f"  Manually annotated coding (G9): {coding}")
print(f"  Correctly attributed non-coding: {noncoding}")
print(f"Saved: results/vep/cmah_rescued_region_override.tsv")
