#!/usr/bin/env python3
"""
Phase H, checkpoint H1b: for each of the 62 candidates, classify why a
variant record might be absent in one cohort's QC-filtered VCF. Checks
three possible upstream stages to distinguish:
  PRESENT           - exact record exists in the final QC VCF
  ABSENT_AFTER_QC   - record existed in an earlier pipeline stage
                      (normalized/biallelic) but was removed by QC
  NOT_IN_SOURCE     - record never existed even at the normalized stage,
                      i.e. genuinely invariant/not called in that cohort
"""

import csv
import subprocess

candidates = []
with open("results/vep/phase_g_protein_altering_FINAL.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates.append({
            "CHROM": row["CHROM"], "POS": row["POS"], "REF": row["REF"],
            "ALT": row["ALT"], "SYMBOL": row["SYMBOL"]
        })

def record_exists(vcf_path, chrom, pos, ref, alt):
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

results = []
for c in candidates:
    for cohort, qc_dir, norm_dir in [
        ("South_Africa", "data/processed/qc_filtered_corrected/south_africa",
         "data/processed/normalized"),
        ("Europe", "data/processed/qc_filtered_corrected/europe",
         "data/processed/european/normalized")
    ]:
        # CMAH uses a different normalized path (corrected extraction)
        if c["SYMBOL"] == "CMAH":
            norm_path = f"data/processed/corrected_cmah/{'south_africa' if cohort=='South_Africa' else 'europe'}/CMAH_norm.vcf.gz"
        else:
            norm_path = f"{norm_dir}/{c['SYMBOL']}_norm.vcf.gz"

        qc_path = f"{qc_dir}/{c['SYMBOL']}_qc.vcf.gz"

        in_qc = record_exists(qc_path, c["CHROM"], c["POS"], c["REF"], c["ALT"])
        if in_qc:
            status = "PRESENT"
        else:
            in_norm = record_exists(norm_path, c["CHROM"], c["POS"], c["REF"], c["ALT"])
            status = "ABSENT_AFTER_QC" if in_norm else "NOT_IN_SOURCE"

        results.append({
            "SYMBOL": c["SYMBOL"], "CHROM": c["CHROM"], "POS": c["POS"],
            "REF": c["REF"], "ALT": c["ALT"], "COHORT": cohort, "STATUS": status
        })

columns = ["SYMBOL","CHROM","POS","REF","ALT","COHORT","STATUS"]
with open("results/phase_H_cohort_absence_classification.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

from collections import Counter
print("=== Classification summary ===")
for status, n in Counter(r["STATUS"] for r in results).most_common():
    print(f"{status}: {n}")

print()
print("=== ABSENT_AFTER_QC entries (need closest attention) ===")
for r in results:
    if r["STATUS"] == "ABSENT_AFTER_QC":
        print(f"{r['SYMBOL']}\t{r['CHROM']}:{r['POS']}\t{r['COHORT']}")

print()
print("Saved: results/phase_H_cohort_absence_classification.tsv")
