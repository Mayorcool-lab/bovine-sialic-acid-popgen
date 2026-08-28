#!/usr/bin/env python3
"""
Phase H, checkpoint H1d: for the 27 NOT_IN_SOURCE candidate x cohort
cases, determine whether the absence reflects genuine invariant/
reference status (real information) or simply an absence of variant
calling at that position (missing information, not zero).

Method: query the cohort's FULL gene-region normalized VCF (not just
the exact candidate position) to see whether the position falls WITHIN
a region that had genotype calls nearby (suggesting real coverage/
calling occurred, and this specific position was genuinely invariant),
versus a region with no nearby calls at all (suggesting a coverage or
calling gap, not a confirmed-invariant site).
"""

import csv
import subprocess

# Load the H1b classification to get the 27 NOT_IN_SOURCE cases
not_in_source = []
with open("results/phase_H_cohort_absence_classification.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["STATUS"] == "NOT_IN_SOURCE":
            not_in_source.append(row)

print(f"Investigating {len(not_in_source)} NOT_IN_SOURCE cases")
print()

for r in not_in_source:
    gene, chrom, pos, cohort = r["SYMBOL"], r["CHROM"], r["POS"], r["COHORT"]
    pos_int = int(pos)

    if cohort == "South_Africa":
        if gene == "CMAH":
            vcf = "data/processed/corrected_cmah/south_africa/CMAH_norm.vcf.gz"
        else:
            vcf = f"data/processed/normalized/{gene}_norm.vcf.gz"
    else:
        if gene == "CMAH":
            vcf = "data/processed/corrected_cmah/europe/CMAH_norm.vcf.gz"
        else:
            vcf = f"data/processed/european/normalized/{gene}_norm.vcf.gz"

    # Check for any called variants within 5kb of this position
    window_start = max(1, pos_int - 5000)
    window_end = pos_int + 5000
    cmd = ["bcftools", "query", "-r", f"{chrom}:{window_start}-{window_end}",
           "-f", "%POS\n", vcf]
    result = subprocess.run(cmd, capture_output=True, text=True)
    nearby_positions = [p for p in result.stdout.strip().split("\n") if p]

    nearest_dist = None
    if nearby_positions:
        nearest_dist = min(abs(int(p) - pos_int) for p in nearby_positions)

    status = "NEARBY_CALLS_EXIST" if nearby_positions else "NO_NEARBY_CALLS_AT_ALL"
    print(f"{gene}\t{chrom}:{pos}\t{cohort}\tnearby_calls={len(nearby_positions)}\tnearest_dist={nearest_dist}\t{status}")
