#!/usr/bin/env python3
"""H14: Pooled and breed-specific MAF for bovine pathway SNPs."""

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

VCF = Path("data/external/aiyesa_popgen/mega_shared.vcf.gz")
MANIFEST = Path(
    "results/H_integration/H9_candidate_breed_fst/"
    "breed_sample_manifest_346.tsv"
)
OUT = Path("results/H_integration/H14_maf_heterozygosity")

OUT.mkdir(parents=True, exist_ok=True)

with MANIFEST.open(newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))

sample_breed = {row["vcf_id"]: row["breed"] for row in manifest}

if len(manifest) != 346 or len(sample_breed) != 346:
    raise SystemExit("ERROR: Expected 346 unique manifest samples")

breed_sizes = Counter(sample_breed.values())

if len(breed_sizes) != 16:
    raise SystemExit("ERROR: Expected 16 breeds")

pooled_rows = []
breed_rows = []
missing_format_counts = Counter()

with gzip.open(VCF, "rt") as handle:
    samples = None

    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]

            if len(samples) != 346 or set(samples) != set(sample_breed):
                raise SystemExit("ERROR: VCF and manifest samples differ")

            if len(samples) != len(set(samples)):
                raise SystemExit("ERROR: Duplicate VCF sample IDs")

            continue

        fields = line.rstrip("\n").split("\t")

        chrom, pos, snp_id, ref, alt = fields[:5]

        if (
            len(ref) != 1
            or len(alt) != 1
            or "," in alt
            or ref not in "ACGT"
            or alt not in "ACGT"
        ):
            raise SystemExit(
                f"ERROR: Unexpected SNP alleles at {chrom}:{pos}"
            )

        fmt = fields[8].split(":")

        if "GT" not in fmt:
            raise SystemExit(
                f"ERROR: Missing GT field at {chrom}:{pos}"
            )

        gt_index = fmt.index("GT")

        if len(fields[9:]) != len(samples):
            raise SystemExit(
                f"ERROR: Sample count mismatch at {chrom}:{pos}"
            )

        counts = defaultdict(lambda: [0, 0, 0])

        for sample, entry in zip(samples, fields[9:]):
            values = entry.split(":")

            if gt_index >= len(values):
                raise SystemExit(
                    f"ERROR: Missing GT value at {chrom}:{pos}"
                )

            raw_gt = values[gt_index]
            genotype = raw_gt.replace("|", "/").split("/")

            # Missing calls may appear as ".", "./.", ".|.",
            # or partially missing diploid genotypes.
            if "." in genotype:
                missing_format_counts[raw_gt] += 1
                continue

            if len(genotype) != 2:
                raise SystemExit(
                    f"ERROR: Non-diploid called genotype "
                    f"{raw_gt} at {chrom}:{pos}"
                )

            if any(a not in ("0", "1") for a in genotype):
                raise SystemExit(
                    f"ERROR: Unexpected genotype "
                    f"{raw_gt} at {chrom}:{pos}"
                )

            breed = sample_breed[sample]

            counts[breed][0] += genotype.count("0")
            counts[breed][1] += genotype.count("1")
            counts[breed][2] += 1

        def frequency_row(breed, ref_count, alt_count, called, n):
            allele_number = ref_count + alt_count

            if allele_number != 2 * called:
                raise SystemExit(
                    "ERROR: Allele-count inconsistency"
                )

            if called > n:
                raise SystemExit(
                    "ERROR: Called samples exceed breed size"
                )

            if allele_number:
                af = alt_count / allele_number
                maf = min(af, 1 - af)
                polymorphic = int(
                    ref_count > 0 and alt_count > 0
                )
            else:
                af = None
                maf = None
                polymorphic = "NA"

            return {
                "chrom": chrom,
                "pos": pos,
                "id": snp_id,
                "ref": ref,
                "alt": alt,
                "breed": breed,
                "n_samples": n,
                "n_called": called,
                "n_missing": n - called,
                "ref_allele_count": ref_count,
                "alt_allele_count": alt_count,
                "allele_number": allele_number,
                "alt_frequency": (
                    f"{af:.10f}" if af is not None else "NA"
                ),
                "maf": (
                    f"{maf:.10f}" if maf is not None else "NA"
                ),
                "polymorphic": polymorphic,
            }

        pooled_ref = sum(v[0] for v in counts.values())
        pooled_alt = sum(v[1] for v in counts.values())
        pooled_called = sum(v[2] for v in counts.values())

        pooled_rows.append(
            frequency_row(
                "ALL",
                pooled_ref,
                pooled_alt,
                pooled_called,
                len(samples),
            )
        )

        for breed in sorted(breed_sizes):
            ref_count, alt_count, called = counts[breed]

            breed_rows.append(
                frequency_row(
                    breed,
                    ref_count,
                    alt_count,
                    called,
                    breed_sizes[breed],
                )
            )

if len(pooled_rows) != 4505:
    raise SystemExit(
        f"ERROR: Expected 4505 SNPs, found {len(pooled_rows)}"
    )

if len(breed_rows) != 4505 * 16:
    raise SystemExit("ERROR: Incomplete breed-SNP table")

total_missing = sum(
    row["n_missing"] for row in pooled_rows
)

if total_missing != 74353:
    raise SystemExit(
        f"ERROR: Missing genotype total {total_missing} "
        "does not match independently verified 74353"
    )

polymorphic = sum(
    int(row["polymorphic"]) for row in pooled_rows
)

zero_alt = sum(
    row["alt_allele_count"] == 0
    for row in pooled_rows
)

fixed_alt = sum(
    row["alt_allele_count"] == row["allele_number"]
    for row in pooled_rows
)

if (polymorphic, zero_alt, fixed_alt) != (4298, 205, 2):
    raise SystemExit(
        "ERROR: Pooled allele classifications do not "
        "match independently verified bcftools results"
    )

fields = list(pooled_rows[0])

for filename, rows in (
    ("H14_pooled_maf.tsv", pooled_rows),
    ("H14_breed_maf.tsv", breed_rows),
):
    with (OUT / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

print("H14 MAF calculation complete")
print("Samples:", len(samples))
print("Breeds:", len(breed_sizes))
print("Pooled SNP rows:", len(pooled_rows))
print("Breed-SNP rows:", len(breed_rows))
print("Missing genotype calls:", total_missing)
print("Polymorphic pooled SNPs:", polymorphic)
print("ALT absent:", zero_alt)
print("ALT fixed:", fixed_alt)
print("Missing genotype formats:")
for gt, count in sorted(missing_format_counts.items()):
    print(f"  {gt}: {count}")
print("Output directory:", OUT)
