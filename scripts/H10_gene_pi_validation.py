#!/usr/bin/env python3
"""Independent breed-level variant-site pi validation for H10."""
import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path("data/external/aiyesa_popgen")
MANIFEST = Path(
    "results/H_integration/H9_candidate_breed_fst/"
    "breed_sample_manifest_346.tsv"
)
OUT = Path("results/H_integration/H10_gene_pi")
OUT.mkdir(parents=True, exist_ok=True)

with MANIFEST.open() as f:
    manifest = list(csv.DictReader(f, delimiter="\t"))

sample_breed = {r["vcf_id"]: r["breed"] for r in manifest}
if len(sample_breed) != 346:
    raise SystemExit("ERROR: Manifest does not contain 346 unique IDs")

with (ROOT / "pi_res_breed.txt").open() as f:
    reported = {
        r["Breed"]: r
        for r in csv.DictReader(f, delimiter="\t")
        if r["Group"] == "Genes"
    }

groups = defaultdict(list)
sums = defaultdict(float)
valid = defaultdict(int)
polymorphic = defaultdict(int)
missing_gt = defaultdict(int)

with gzip.open(ROOT / "mega_shared.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            ids = line.rstrip("\n").split("\t")[9:]
            if len(ids) != 346 or set(ids) != set(sample_breed):
                raise SystemExit("ERROR: VCF and manifest IDs differ")
            for i, sid in enumerate(ids):
                groups[sample_breed[sid]].append(i)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields[4].split(",")) != 1:
            raise SystemExit("ERROR: Multiallelic record encountered")

        gt_index = fields[8].split(":").index("GT")
        genotypes = [
            x.split(":")[gt_index]
            for x in fields[9:]
        ]

        for breed, indices in groups.items():
            counts = [0, 0]
            for i in indices:
                gt = genotypes[i].replace("|", "/").split("/")
                for allele in gt:
                    if allele == ".":
                        missing_gt[breed] += 1
                    elif allele in ("0", "1"):
                        counts[int(allele)] += 1
                    else:
                        raise SystemExit(
                            f"ERROR: Unexpected allele {allele}"
                        )

            n = sum(counts)
            if n < 2:
                continue

            pi = (
                n / (n - 1)
                * (1 - sum((c / n) ** 2 for c in counts))
            )
            sums[breed] += pi
            valid[breed] += 1
            if min(counts) > 0:
                polymorphic[breed] += 1

rows = []
for breed in sorted(groups):
    ref = reported.get(breed)
    if ref is None:
        raise SystemExit(f"ERROR: Missing Aiyesa result for {breed}")

    estimate = sums[breed] / valid[breed]
    reference = float(ref["mean_pi"])
    rows.append({
        "breed": breed,
        "n_samples": len(groups[breed]),
        "aiyesa_n_samples": int(ref["n_samples"]),
        "valid_sites": valid[breed],
        "polymorphic_sites": polymorphic[breed],
        "aiyesa_n_sites": int(ref["n_sites"]),
        "mean_pi_recomputed": f"{estimate:.8f}",
        "mean_pi_aiyesa": f"{reference:.8f}",
        "difference": f"{estimate-reference:+.8f}",
        "missing_allele_calls": missing_gt[breed],
    })

output = OUT / "gene_pi_validation.tsv"
with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("Breed | n | Valid | Polymorphic | Aiyesa sites | Our pi | Aiyesa pi | Difference")
for r in rows:
    print(
        f"{r['breed']} | {r['n_samples']} | "
        f"{r['valid_sites']} | {r['polymorphic_sites']} | "
        f"{r['aiyesa_n_sites']} | "
        f"{r['mean_pi_recomputed']} | "
        f"{r['mean_pi_aiyesa']} | {r['difference']}"
    )

print("\nSample count mismatches:")
print([
    r["breed"] for r in rows
    if r["n_samples"] != r["aiyesa_n_samples"]
])

print("\nSaved:", output)
