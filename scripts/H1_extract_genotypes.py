#!/usr/bin/env python3
"""
Phase H, checkpoint H1: extract full per-animal genotypes for all 62
protein-altering candidates, both cohorts, using the corrected,
QC-filtered VCFs already validated in H0.
"""

import csv
import subprocess

candidates = []
with open("results/vep/phase_g_protein_altering_FINAL.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates.append({
            "CHROM": row["CHROM"], "POS": row["POS"], "REF": row["REF"],
            "ALT": row["ALT"], "SYMBOL": row["SYMBOL"],
            "Consequence": row["Consequence"], "HGVSp": row["HGVSp"]
        })

sa_dir = "data/processed/qc_filtered_corrected/south_africa"
eu_dir = "data/processed/qc_filtered_corrected/europe"

def get_genotypes(vcf_path, chrom, pos):
    cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
           "-f", "[%SAMPLE=%GT\t]\n", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if not result.stdout.strip():
        return {}
    genos = {}
    for entry in result.stdout.strip().split("\t"):
        entry = entry.strip()
        if not entry or "=" not in entry:
            continue
        sample, gt = entry.split("=")
        genos[sample] = gt
    return genos

rows = []
for c in candidates:
    sa_vcf = f"{sa_dir}/{c['SYMBOL']}_qc.vcf.gz"
    eu_vcf = f"{eu_dir}/{c['SYMBOL']}_qc.vcf.gz"
    sa_genos = get_genotypes(sa_vcf, c["CHROM"], c["POS"])
    eu_genos = get_genotypes(eu_vcf, c["CHROM"], c["POS"])

    for sample, gt in sa_genos.items():
        rows.append({"CHROM": c["CHROM"], "POS": c["POS"], "REF": c["REF"],
                      "ALT": c["ALT"], "SYMBOL": c["SYMBOL"],
                      "Consequence": c["Consequence"], "HGVSp": c["HGVSp"],
                      "COHORT": "South_Africa", "SAMPLE": sample, "GT": gt})
    for sample, gt in eu_genos.items():
        rows.append({"CHROM": c["CHROM"], "POS": c["POS"], "REF": c["REF"],
                      "ALT": c["ALT"], "SYMBOL": c["SYMBOL"],
                      "Consequence": c["Consequence"], "HGVSp": c["HGVSp"],
                      "COHORT": "Europe", "SAMPLE": sample, "GT": gt})

columns = ["CHROM","POS","REF","ALT","SYMBOL","Consequence","HGVSp","COHORT","SAMPLE","GT"]
with open("results/phase_H_genotypes_raw.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Total genotype records extracted: {len(rows)}")
print(f"Expected: 62 variants x 346 animals = {62*346}")
sa_count = sum(1 for r in rows if r["COHORT"] == "South_Africa")
eu_count = sum(1 for r in rows if r["COHORT"] == "Europe")
print(f"South Africa records: {sa_count} (expected up to 62 x 120 = {62*120})")
print(f"Europe records: {eu_count} (expected up to 62 x 226 = {62*226})")
print(f"Saved: results/phase_H_genotypes_raw.tsv")
