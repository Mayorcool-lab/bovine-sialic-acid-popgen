#!/usr/bin/env python3
"""
Phase H, checkpoint H0: validate that all 62 protein-altering candidates
from phase_g_protein_altering_FINAL.tsv genuinely exist in the QC-filtered
cohort VCFs (data/processed/qc_filtered_corrected/), matching exactly on
CHROM, POS, REF, ALT. This must pass before any genotype extraction happens.
"""

import csv
import subprocess

candidates = []
with open("results/vep/phase_g_protein_altering_FINAL.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates.append((row["CHROM"], row["POS"], row["REF"], row["ALT"], row["SYMBOL"]))

print(f"Total candidates to validate: {len(candidates)}")

def check_in_vcf(vcf_path, chrom, pos, ref, alt):
    cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
           "-f", "%CHROM\t%POS\t%REF\t%ALT\n", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        c, p, r, a = line.split("\t")
        if c == chrom and p == pos and r == ref and alt in a.split(","):
            return True
    return False

sa_dir = "data/processed/qc_filtered_corrected/south_africa"
eu_dir = "data/processed/qc_filtered_corrected/europe"

results = []
for chrom, pos, ref, alt, symbol in candidates:
    sa_vcf = f"{sa_dir}/{symbol}_qc.vcf.gz"
    eu_vcf = f"{eu_dir}/{symbol}_qc.vcf.gz"
    in_sa = check_in_vcf(sa_vcf, chrom, pos, ref, alt)
    in_eu = check_in_vcf(eu_vcf, chrom, pos, ref, alt)
    results.append((symbol, chrom, pos, ref, alt, in_sa, in_eu))

print()
print("=== VALIDATION RESULTS ===")
fail_count = 0
for symbol, chrom, pos, ref, alt, in_sa, in_eu in results:
    status = "OK" if (in_sa or in_eu) else "FAIL - NOT FOUND IN EITHER COHORT"
    if not (in_sa or in_eu):
        fail_count += 1
    print(f"{symbol}\t{chrom}:{pos}\t{ref}>{alt}\tSA={in_sa}\tEU={in_eu}\t{status}")

print()
if fail_count == 0:
    print(f"PASS: all {len(candidates)} candidates confirmed present in at least one cohort VCF.")
else:
    print(f"FAIL: {fail_count} candidates NOT found in either cohort VCF. Do not proceed to H1 until resolved.")
