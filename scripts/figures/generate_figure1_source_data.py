#!/usr/bin/env python3
"""
Figure 1 source data generation.

Re-derives Panel B, C, and D data directly from the frozen Phase G outputs:
  - data/reference/candidate_gene_coordinates_provisional.tsv
  - results/vep/phase_g_final_target_genes_corrected.tsv
  - results/vep/phase_g_protein_altering_FINAL.tsv

Panel A (pathway schematic) uses no repository data; tier/gene assignments
are fixed by manuscript Table 1 and are hardcoded in the plotting script.

Outputs:
  results/figures/source_data/Figure1_panelB_gene_landscape.tsv
  results/figures/source_data/Figure1_panelC_consequences.tsv
  results/figures/source_data/Figure1_panelD_candidates.tsv
"""

import csv
from collections import Counter, defaultdict

COORD_FILE = "data/reference/candidate_gene_coordinates_provisional.tsv"
TARGET_FILE = "results/vep/phase_g_final_target_genes_corrected.tsv"
CANDIDATE_FILE = "results/vep/phase_g_protein_altering_FINAL.tsv"

OUT_B = "results/figures/source_data/Figure1_panelB_gene_landscape.tsv"
OUT_C = "results/figures/source_data/Figure1_panelC_consequences.tsv"
OUT_D = "results/figures/source_data/Figure1_panelD_candidates.tsv"

UPSTREAM_BP = 2000

lengths = {}
with open(COORD_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        lengths[r["gene_symbol"]] = (
            int(r["end"]) - int(r["start"]) + 1 + UPSTREAM_BP
        )

variants = Counter()
with open(TARGET_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        variants[r["SYMBOL"]] += 1

candidates = Counter()
with open(CANDIDATE_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        candidates[r["SYMBOL"]] += 1

if len(lengths) != 13:
    raise RuntimeError(f"Expected 13 genes in coordinate file, got {len(lengths)}")
if sum(variants.values()) != 9724:
    raise RuntimeError(f"Expected 9724 target-gene variants, got {sum(variants.values())}")
if sum(candidates.values()) != 62:
    raise RuntimeError(f"Expected 62 protein-altering candidates, got {sum(candidates.values())}")

with open(OUT_B, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["gene", "analysed_kb", "target_gene_variants", "variants_per_kb", "protein_altering_candidates"])
    for g in sorted(lengths):
        kb = lengths[g] / 1000
        w.writerow([g, f"{kb:.3f}", variants[g], f"{variants[g]/kb:.3f}", candidates[g]])

print(f"Saved: {OUT_B} ({len(lengths)} genes)")

groups = Counter()
with open(TARGET_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        c = r["Consequence"]
        if "missense_variant" in c or c in {"stop_gained", "stop_lost"}:
            groups["Protein-altering"] += 1
        elif "synonymous_variant" in c:
            groups["Synonymous"] += 1
        elif "splice" in c:
            groups["Splice-associated"] += 1
        elif "5_prime_UTR_variant" in c or "3_prime_UTR_variant" in c:
            groups["UTR"] += 1
        elif "upstream_gene_variant" in c or "downstream_gene_variant" in c:
            groups["Upstream/downstream"] += 1
        elif "intron_variant" in c:
            groups["Intronic"] += 1
        else:
            groups["Other"] += 1

if sum(groups.values()) != 9724:
    raise RuntimeError(f"Panel C total mismatch: {sum(groups.values())} != 9724")
if groups.get("Protein-altering", 0) != 58:
    raise RuntimeError(f"Expected 58 VEP-native protein-altering, got {groups.get('Protein-altering', 0)}")

with open(OUT_C, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["consequence_group", "count"])
    for k, v in groups.items():
        w.writerow([k, v])

print(f"Saved: {OUT_C} (total {sum(groups.values())}, scope: 9,724 VEP-annotated target-gene variants)")

CONSEQUENCE_MAP = {
    "missense_variant": "missense",
    "missense_variant&splice_region_variant": "missense_splice",
    "stop_gained": "stop_gained",
    "stop_lost": "stop_lost",
}

breakdown = defaultdict(Counter)
with open(CANDIDATE_FILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        ctype = CONSEQUENCE_MAP.get(r["Consequence"])
        if ctype is None:
            raise RuntimeError(f"Unrecognised consequence in candidate file: {r['Consequence']}")
        breakdown[r["SYMBOL"]][ctype] += 1

total_check = sum(sum(v.values()) for v in breakdown.values())
if total_check != 62:
    raise RuntimeError(f"Panel D total mismatch: {total_check} != 62")

with open(OUT_D, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["gene", "missense", "missense_splice", "stop_gained", "stop_lost", "total"])
    for g in sorted(breakdown):
        row = breakdown[g]
        total = sum(row.values())
        w.writerow([g, row.get("missense", 0), row.get("missense_splice", 0),
                    row.get("stop_gained", 0), row.get("stop_lost", 0), total])

print(f"Saved: {OUT_D} (total {total_check} candidates across {len(breakdown)} genes)")
print("\nFIGURE 1 SOURCE DATA GENERATION: PASS")
