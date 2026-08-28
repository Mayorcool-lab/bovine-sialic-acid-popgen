#!/usr/bin/env python3
"""
Phase H, checkpoint H2: calculate South Africa vs. Europe allele
frequencies for all 62 protein-altering candidates, using the verified
three-way classification from H1b-H1d:
  PRESENT                - use genotype data from the final QC VCF
  ABSENT_AFTER_QC         - use the pre-site-QC normalized/biallelic VCF,
                            apply DP>=4/GQ>=20 genotype masking directly,
                            calculate AF from only the callable animals
  NOT_DETECTED_IN_SOURCE  - report as not detected, NOT as AF=0.000
"""

import csv
import subprocess
from collections import Counter

candidates = []
with open("results/vep/phase_g_protein_altering_FINAL.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates.append({
            "CHROM": row["CHROM"], "POS": row["POS"], "REF": row["REF"],
            "ALT": row["ALT"], "SYMBOL": row["SYMBOL"],
            "Consequence": row["Consequence"], "HGVSp": row["HGVSp"]
        })

classification = {}
with open("results/phase_H_cohort_absence_classification.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        key = (row["SYMBOL"], row["CHROM"], row["POS"], row["COHORT"])
        classification[key] = row["STATUS"]

sa_qc_dir = "data/processed/qc_filtered_corrected/south_africa"
eu_qc_dir = "data/processed/qc_filtered_corrected/europe"

def norm_path(gene, cohort):
    if gene == "CMAH":
        sub = "south_africa" if cohort == "South_Africa" else "europe"
        return f"data/processed/corrected_cmah/{sub}/CMAH_norm.vcf.gz"
    if cohort == "South_Africa":
        return f"data/processed/normalized/{gene}_norm.vcf.gz"
    return f"data/processed/european/normalized/{gene}_norm.vcf.gz"

def parse_gt_alleles(gt):
    gt = gt.replace("|", "/")
    if gt in ("./.", "."):
        return None
    return gt.split("/")

def get_genotypes_simple(vcf_path, chrom, pos, ref, alt):
    cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
           "-f", "%REF\t%ALT\n", vcf_path]
    check = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in check.stdout.strip().splitlines():
        r, a = line.split("\t")
        if r == ref and alt in a.split(","):
            break
    else:
        raise RuntimeError(f"REF/ALT mismatch or not found at {chrom}:{pos} in {vcf_path}: expected {ref}>{alt}")

    cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
           "-f", "[%SAMPLE=%GT\t]\n", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    genos = {}
    for entry in result.stdout.strip().split("\t"):
        entry = entry.strip()
        if not entry or "=" not in entry:
            continue
        sample, gt = entry.split("=")
        genos[sample] = gt
    return genos

def get_genotypes_with_qc_mask(vcf_path, chrom, pos, ref, alt):
    check_cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
                 "-f", "%REF\t%ALT\n", vcf_path]
    check = subprocess.run(check_cmd, capture_output=True, text=True, check=True)
    for line in check.stdout.strip().splitlines():
        r, a = line.split("\t")
        if r == ref and alt in a.split(","):
            break
    else:
        raise RuntimeError(f"REF/ALT mismatch or not found at {chrom}:{pos} in {vcf_path}: expected {ref}>{alt}")

    cmd = ["bcftools", "query", "-r", f"{chrom}:{pos}-{pos}",
           "-f", "[%SAMPLE\t%GT\t%DP\t%GQ\n]", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    genos = {}
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        sample, gt, dp, gq = line.split("\t")
        try:
            dp_n = int(dp)
        except ValueError:
            dp_n = None
        try:
            gq_n = int(gq)
        except ValueError:
            gq_n = None
        if dp_n is not None and dp_n >= 4 and gq_n is not None and gq_n >= 20 and gt not in ("./.", ".|."):
            genos[sample] = gt
    return genos

def calc_freq(genos):
    valid = {s: g for s, g in genos.items() if g not in ("./.", ".|.")}
    n = len(valid)
    if n == 0:
        return None, 0, 0, 0, 0
    ref_count = 0
    alt_count = 0
    het = hom_alt = hom_ref = 0
    for gt in valid.values():
        alleles = parse_gt_alleles(gt)
        if alleles is None:
            continue
        ref_count += alleles.count("0")
        alt_count += alleles.count("1")
        if alleles == ["0", "0"]:
            hom_ref += 1
        elif alleles == ["1", "1"]:
            hom_alt += 1
        elif set(alleles) == {"0", "1"}:
            het += 1
    total = ref_count + alt_count
    af = alt_count / total if total > 0 else None
    return af, n, het, hom_alt, hom_ref

rows = []
for c in candidates:
    key_sa = (c["SYMBOL"], c["CHROM"], c["POS"], "South_Africa")
    key_eu = (c["SYMBOL"], c["CHROM"], c["POS"], "Europe")
    status_sa = classification.get(key_sa)
    status_eu = classification.get(key_eu)

    allowed = {"PRESENT", "ABSENT_AFTER_QC", "NOT_IN_SOURCE"}
    if status_sa not in allowed:
        raise RuntimeError(f"Invalid/missing SA classification for {key_sa}: {status_sa}")
    if status_eu not in allowed:
        raise RuntimeError(f"Invalid/missing EU classification for {key_eu}: {status_eu}")

    if status_sa == "NOT_IN_SOURCE":
        sa_af, sa_n, sa_het, sa_homalt, sa_homref = (None, 0, 0, 0, 0)
        sa_report = "NOT_DETECTED_IN_SOURCE"
    elif status_sa == "ABSENT_AFTER_QC":
        vcf = norm_path(c["SYMBOL"], "South_Africa")
        genos = get_genotypes_with_qc_mask(vcf, c["CHROM"], c["POS"], c["REF"], c["ALT"])
        sa_af, sa_n, sa_het, sa_homalt, sa_homref = calc_freq(genos)
        sa_report = (f"{sa_af:.4f}" if sa_af is not None else "NA") + " (SUB-THRESHOLD CALL RATE)"
    else:
        vcf = f"{sa_qc_dir}/{c['SYMBOL']}_qc.vcf.gz"
        genos = get_genotypes_simple(vcf, c["CHROM"], c["POS"], c["REF"], c["ALT"])
        sa_af, sa_n, sa_het, sa_homalt, sa_homref = calc_freq(genos)
        sa_report = f"{sa_af:.4f}" if sa_af is not None else "NA"

    if status_eu == "NOT_IN_SOURCE":
        eu_af, eu_n, eu_het, eu_homalt, eu_homref = (None, 0, 0, 0, 0)
        eu_report = "NOT_DETECTED_IN_SOURCE"
    elif status_eu == "ABSENT_AFTER_QC":
        vcf = norm_path(c["SYMBOL"], "Europe")
        genos = get_genotypes_with_qc_mask(vcf, c["CHROM"], c["POS"], c["REF"], c["ALT"])
        eu_af, eu_n, eu_het, eu_homalt, eu_homref = calc_freq(genos)
        eu_report = (f"{eu_af:.4f}" if eu_af is not None else "NA") + " (SUB-THRESHOLD CALL RATE)"
    else:
        vcf = f"{eu_qc_dir}/{c['SYMBOL']}_qc.vcf.gz"
        genos = get_genotypes_simple(vcf, c["CHROM"], c["POS"], c["REF"], c["ALT"])
        eu_af, eu_n, eu_het, eu_homalt, eu_homref = calc_freq(genos)
        eu_report = f"{eu_af:.4f}" if eu_af is not None else "NA"

    rows.append({
        "SYMBOL": c["SYMBOL"], "CHROM": c["CHROM"], "POS": c["POS"],
        "REF": c["REF"], "ALT": c["ALT"], "Consequence": c["Consequence"],
        "HGVSp": c["HGVSp"],
        "SA_status": status_sa, "SA_AF": sa_report, "SA_N_genotyped": sa_n,
        "SA_het": sa_het, "SA_homalt": sa_homalt,
        "EU_status": status_eu, "EU_AF": eu_report, "EU_N_genotyped": eu_n,
        "EU_het": eu_het, "EU_homalt": eu_homalt,
    })

columns = ["SYMBOL","CHROM","POS","REF","ALT","Consequence","HGVSp",
           "SA_status","SA_AF","SA_N_genotyped","SA_het","SA_homalt",
           "EU_status","EU_AF","EU_N_genotyped","EU_het","EU_homalt"]

with open("results/phase_H_62_variant_population_frequencies.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Total variants processed: {len(rows)}")
print()
print("=== Status summary ===")
print("SA:", Counter(r["SA_status"] for r in rows))
print("EU:", Counter(r["EU_status"] for r in rows))
print()
print("=== Spot-check: ABSENT_AFTER_QC sites (should match H1c: 107/107/105/98) ===")
for r in rows:
    if r["SA_status"] == "ABSENT_AFTER_QC" or r["EU_status"] == "ABSENT_AFTER_QC":
        print(f"{r['SYMBOL']}\t{r['CHROM']}:{r['POS']}\tSA_N={r['SA_N_genotyped']}\tSA_AF={r['SA_AF']}\tEU_N={r['EU_N_genotyped']}\tEU_AF={r['EU_AF']}")

print()
print("Saved: results/phase_H_62_variant_population_frequencies.tsv")
