#!/usr/bin/env python3

import csv
import subprocess
from collections import Counter

# Definitive H3 breed mapping:
# South Africa = 120 South African animals only.
# Europe = verified 226-animal working cohort only.
sample_breed = {}

with open("data/metadata/combined_metadata_409.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["country"] == "South Africa":
            sample = row["vcf_id"]
            if sample in sample_breed:
                raise RuntimeError(f"Duplicate South African sample ID: {sample}")
            sample_breed[sample] = (row["breed"], "South_Africa")

with open("data/processed/qc_characterization/europe_226_breed_mapping.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sample = row["vcf_id"]
        if sample in sample_breed:
            raise RuntimeError(f"Sample occurs in both cohorts: {sample}")
        sample_breed[sample] = (row["breed"], "Europe")

sa = Counter(
    breed for breed, cohort in sample_breed.values()
    if cohort == "South_Africa"
)

eu = Counter(
    breed for breed, cohort in sample_breed.values()
    if cohort == "Europe"
)

assert sum(sa.values()) == 120
assert sum(eu.values()) == 226
assert len(sample_breed) == 346

print("H3 metadata PASS")
print("SA:", dict(sorted(sa.items())))
print("EU:", dict(sorted(eu.items())))
print("TOTAL:", len(sample_breed))

# Load the definitive 62 protein-altering candidates
candidates = []
with open("results/vep/phase_g_protein_altering_FINAL.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates.append({
            "CHROM": row["CHROM"],
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "SYMBOL": row["SYMBOL"],
            "Consequence": row["Consequence"],
            "HGVSp": row["HGVSp"]
        })

# Load verified H1b cohort classification
classification = {}
with open("results/phase_H_cohort_absence_classification.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        key = (row["SYMBOL"], row["CHROM"], row["POS"], row["COHORT"])
        classification[key] = row["STATUS"]

allowed_statuses = {"PRESENT", "ABSENT_AFTER_QC", "NOT_IN_SOURCE"}

# Validate complete 62 x 2 cohort classification
for c in candidates:
    for cohort in ("South_Africa", "Europe"):
        key = (c["SYMBOL"], c["CHROM"], c["POS"], cohort)
        status = classification.get(key)

        if status not in allowed_statuses:
            raise RuntimeError(
                f"Missing/invalid classification for "
                f"{c['SYMBOL']} {c['CHROM']}:{c['POS']} {cohort}: {status}"
            )

print("\nH3 candidate/classification validation PASS")
print("Candidates:", len(candidates))
print("Expected candidate/cohort combinations:", len(candidates) * 2)

sa_qc_dir = "data/processed/qc_filtered_corrected/south_africa"
eu_qc_dir = "data/processed/qc_filtered_corrected/europe"

def norm_path(gene, cohort):
    if gene == "CMAH":
        sub = "south_africa" if cohort == "South_Africa" else "europe"
        return f"data/processed/corrected_cmah/{sub}/CMAH_norm.vcf.gz"
    if cohort == "South_Africa":
        return f"data/processed/normalized/{gene}_norm.vcf.gz"
    return f"data/processed/european/normalized/{gene}_norm.vcf.gz"

def verify_ref_alt(vcf_path, chrom, pos, ref, alt):
    cmd = [
        "bcftools", "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f", "%REF\t%ALT\n",
        vcf_path
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        r, a = line.split("\t")
        if r == ref and alt in a.split(","):
            return True

    return False

def get_genotypes(vcf_path, chrom, pos, ref, alt, with_qc_mask=False):
    if not verify_ref_alt(vcf_path, chrom, pos, ref, alt):
        raise RuntimeError(
            f"REF/ALT mismatch or site not found at "
            f"{chrom}:{pos} in {vcf_path}: expected {ref}>{alt}"
        )

    if with_qc_mask:
        cmd = [
            "bcftools", "query",
            "-r", f"{chrom}:{pos}-{pos}",
            "-f", "[%SAMPLE\t%GT\t%DP\t%GQ\n]",
            vcf_path
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        genos = {}
        for line in result.stdout.strip().splitlines():
            if not line:
                continue

            sample, gt, dp, gq = line.split("\t")

            try:
                dp_n = int(dp)
                gq_n = int(gq)
            except ValueError:
                continue

            if dp_n >= 4 and gq_n >= 20 and gt not in ("./.", ".|."):
                genos[sample] = gt

        return genos

    cmd = [
        "bcftools", "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f", "[%SAMPLE=%GT\t]\n",
        vcf_path
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    genos = {}
    for entry in result.stdout.strip().split("\t"):
        entry = entry.strip()
        if not entry or "=" not in entry:
            continue

        sample, gt = entry.split("=", 1)
        genos[sample] = gt

    return genos

print("\nH3 genotype functions loaded PASS")

SA_BREEDS = [
    "Afrikaner",
    "Bonsmara",
    "Nguni",
    "Tuli",
]

EU_BREEDS = [
    "Western_Finncattle",
    "Northern_Finncattle",
    "Eastern_Finncattle",
    "Holstein_Friesian",
    "Deep_Red",
    "Dutch_Belted",
    "Dutch_Friesian",
    "Meuse_Rhine_Yssel",
    "Groningen_White_Headed",
    "Barrosa",
    "Mirandesa",
    "Mertolenga",
]

def small_n_flag(n_callable):
    """
    Flag unstable breed-level estimates caused by very small callable N.
    """
    if n_callable < 5:
        return "VERY_SMALL_N"
    if n_callable < 10:
        return "SMALL_N"
    return "OK"

print("\nH3 breed definitions loaded PASS")

def calculate_breed_af(genos, breed, cohort):
    """
    Calculate allele frequency among callable animals belonging to one breed.
    Returns callable N, REF count, ALT count, and ALT allele frequency.
    """
    n_callable = 0
    ref_count = 0
    alt_count = 0

    for sample, gt in genos.items():
        meta = sample_breed.get(sample)

        if meta is None:
            continue

        sample_breed_name, sample_cohort = meta

        if sample_breed_name != breed or sample_cohort != cohort:
            continue

        if gt in ("./.", ".|.", "."):
            continue

        alleles = gt.replace("|", "/").split("/")

        if any(a not in ("0", "1") for a in alleles):
            raise RuntimeError(
                f"Unexpected genotype {gt} for sample {sample}"
            )

        n_callable += 1
        ref_count += alleles.count("0")
        alt_count += alleles.count("1")

    total_alleles = ref_count + alt_count
    af = alt_count / total_alleles if total_alleles else None

    return n_callable, ref_count, alt_count, af


print("H3 breed-frequency calculator loaded PASS")

production_rows = []

for candidate in candidates:
    gene = candidate["SYMBOL"]
    chrom = candidate["CHROM"]
    pos = candidate["POS"]
    ref = candidate["REF"]
    alt = candidate["ALT"]

    for cohort, breeds in [
        ("South_Africa", SA_BREEDS),
        ("Europe", EU_BREEDS),
    ]:
        status = classification[(gene, chrom, pos, cohort)]

        if status == "NOT_IN_SOURCE":
            for breed in breeds:
                production_rows.append({
                    "SYMBOL": gene,
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "Consequence": candidate["Consequence"],
                    "HGVSp": candidate["HGVSp"],
                    "COHORT": cohort,
                    "BREED": breed,
                    "STATUS": "NOT_DETECTED_IN_SOURCE",
                    "N_CALLABLE": "",
                    "REF_COUNT": "",
                    "ALT_COUNT": "",
                    "AF": "",
                    "N_FLAG": "NA",
                })
            continue

        if status == "PRESENT":
            vcf_path = (
                f"{sa_qc_dir}/{gene}_qc.vcf.gz"
                if cohort == "South_Africa"
                else f"{eu_qc_dir}/{gene}_qc.vcf.gz"
            )
            genos = get_genotypes(
                vcf_path,
                chrom,
                pos,
                ref,
                alt,
                with_qc_mask=False,
            )

        elif status == "ABSENT_AFTER_QC":
            genos = get_genotypes(
                norm_path(gene, cohort),
                chrom,
                pos,
                ref,
                alt,
                with_qc_mask=True,
            )

        else:
            raise RuntimeError(
                f"Unexpected classification {status} for "
                f"{gene} {chrom}:{pos} {cohort}"
            )

        for breed in breeds:
            n_callable, ref_count, alt_count, af = calculate_breed_af(
                genos,
                breed,
                cohort,
            )

            production_rows.append({
                "SYMBOL": gene,
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "Consequence": candidate["Consequence"],
                "HGVSp": candidate["HGVSp"],
                "COHORT": cohort,
                "BREED": breed,
                "STATUS": status,
                "N_CALLABLE": n_callable,
                "REF_COUNT": ref_count,
                "ALT_COUNT": alt_count,
                "AF": "" if af is None else f"{af:.6f}",
                "N_FLAG": small_n_flag(n_callable),
            })

print("\nH3 production row generation PASS")
print("Rows generated:", len(production_rows))

# ------------------------------------------------------------
# H3 -> H2 reconciliation audit
# ------------------------------------------------------------

h2 = {}

with open("results/phase_H_62_variant_population_frequencies.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        key = (
            row["SYMBOL"],
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"],
        )
        h2[key] = row

if len(h2) != 62:
    raise RuntimeError(
        f"Expected 62 H2 variants, found {len(h2)}"
    )

reconciliation_pass = 0

for candidate in candidates:
    key = (
        candidate["SYMBOL"],
        candidate["CHROM"],
        candidate["POS"],
        candidate["REF"],
        candidate["ALT"],
    )

    h2_row = h2[key]

    for cohort, prefix in [
        ("South_Africa", "SA"),
        ("Europe", "EU"),
    ]:
        rows = [
            r for r in production_rows
            if (
                r["SYMBOL"],
                r["CHROM"],
                r["POS"],
                r["REF"],
                r["ALT"],
            ) == key
            and r["COHORT"] == cohort
        ]

        expected_breeds = 4 if cohort == "South_Africa" else 12

        if len(rows) != expected_breeds:
            raise RuntimeError(
                f"{key} {cohort}: expected {expected_breeds} breed rows, "
                f"found {len(rows)}"
            )

        h2_status = h2_row[f"{prefix}_status"]
        h2_af = h2_row[f"{prefix}_AF"]
        h2_n = int(h2_row[f"{prefix}_N_genotyped"])

        if h2_status == "NOT_IN_SOURCE":
            if h2_af != "NOT_DETECTED_IN_SOURCE":
                raise RuntimeError(
                    f"{key} {cohort}: H2 NOT_IN_SOURCE has unexpected AF {h2_af}"
                )

            if h2_n != 0:
                raise RuntimeError(
                    f"{key} {cohort}: H2 NOT_IN_SOURCE has N={h2_n}, expected 0"
                )

            if any(r["STATUS"] != "NOT_DETECTED_IN_SOURCE" for r in rows):
                raise RuntimeError(
                    f"{key} {cohort}: H3 failed to preserve NOT_IN_SOURCE"
                )

            if any(r["AF"] != "" or r["N_CALLABLE"] != "" for r in rows):
                raise RuntimeError(
                    f"{key} {cohort}: H3 assigned frequency/counts to NOT_IN_SOURCE"
                )

            reconciliation_pass += 1
            continue

        if any(r["STATUS"] != h2_status for r in rows):
            raise RuntimeError(
                f"{key} {cohort}: H3 status differs from H2 {h2_status}"
            )

        breed_n = sum(int(r["N_CALLABLE"]) for r in rows)
        breed_ref = sum(int(r["REF_COUNT"]) for r in rows)
        breed_alt = sum(int(r["ALT_COUNT"]) for r in rows)

        if breed_n != h2_n:
            raise RuntimeError(
                f"{key} {cohort}: breed N={breed_n}, H2 N={h2_n}"
            )

        total_alleles = breed_ref + breed_alt

        if total_alleles != 2 * breed_n:
            raise RuntimeError(
                f"{key} {cohort}: allele count {total_alleles} "
                f"does not equal 2*N ({2 * breed_n})"
            )

        breed_af = breed_alt / total_alleles if total_alleles else None

        if breed_af is None:
            raise RuntimeError(
                f"{key} {cohort}: no callable alleles despite H2 status {h2_status}"
            )

        h2_af_numeric = h2_af.split()[0]

        if f"{breed_af:.4f}" != h2_af_numeric:
            raise RuntimeError(
                f"{key} {cohort}: reconstructed AF={breed_af:.4f}, "
                f"H2 AF={h2_af}"
            )

        reconciliation_pass += 1

print("\nH3 -> H2 reconciliation PASS")
print("Variant/cohort combinations reconciled:", reconciliation_pass)

# ---------------------------------------------------------------------------
# Write authoritative H3 breed-level frequency table only after reconciliation
# ---------------------------------------------------------------------------

output_path = "results/phase_H_62_variant_breed_frequencies.tsv"

fieldnames = [
    "SYMBOL",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Consequence",
    "HGVSp",
    "COHORT",
    "BREED",
    "STATUS",
    "N_CALLABLE",
    "REF_COUNT",
    "ALT_COUNT",
    "AF",
    "N_FLAG",
]

with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(production_rows)

print("\nH3 output written PASS")
print("Rows written:", len(production_rows))
print("Saved:", output_path)
