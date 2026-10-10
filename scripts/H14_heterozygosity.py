#!/usr/bin/env python3
"""H14: Observed and expected heterozygosity at pathway SNP records."""

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


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def variant_key(row):
    return (
        row["chrom"],
        row["pos"],
        row["ref"],
        row["alt"],
        row["breed"],
    )


manifest = read_tsv(MANIFEST)
sample_breed = {r["vcf_id"]: r["breed"] for r in manifest}
breed_sizes = Counter(sample_breed.values())

assert len(manifest) == len(sample_breed) == 346
assert len(breed_sizes) == 16

pooled = read_tsv(OUT / "H14_pooled_maf.tsv")
breed_maf = read_tsv(OUT / "H14_breed_maf.tsv")

assert len(pooled) == 4505
assert len(breed_maf) == 4505 * 16

maf_lookup = {
    variant_key(r): r
    for r in pooled + breed_maf
}

assert len(maf_lookup) == 4505 * 17, (
    "Non-unique variant identifiers"
)

groups = ["ALL"] + sorted(breed_sizes)
sizes = {"ALL": 346, **breed_sizes}

site_rows = []
totals = defaultdict(Counter)
seen_keys = set()

with gzip.open(VCF, "rt") as handle:
    samples = None

    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]

            assert len(samples) == len(set(samples)) == 346
            assert set(samples) == set(sample_breed)
            continue

        f = line.rstrip("\n").split("\t")

        chrom, pos, ref, alt = f[0], f[1], f[3], f[4]

        assert len(ref) == len(alt) == 1
        assert ref in "ACGT" and alt in "ACGT"
        assert "," not in alt

        fmt = f[8].split(":")
        assert "GT" in fmt

        gt_index = fmt.index("GT")

        assert len(f[9:]) == len(samples)

        counts = defaultdict(Counter)

        for sample, entry in zip(samples, f[9:]):
            values = entry.split(":")

            assert gt_index < len(values)

            raw_gt = values[gt_index]
            gt = raw_gt.replace("|", "/").split("/")

            if "." in gt:
                continue

            assert len(gt) == 2, (
                chrom, pos, sample, raw_gt
            )

            assert all(a in ("0", "1") for a in gt), (
                chrom, pos, sample, raw_gt
            )

            breed = sample_breed[sample]
            category = gt.count("1")

            for group in ("ALL", breed):
                counts[group][category] += 1

        for group in groups:
            c = counts[group]

            hom_ref = c[0]
            het = c[1]
            hom_alt = c[2]

            called = hom_ref + het + hom_alt
            missing = sizes[group] - called

            key = (chrom, pos, ref, alt, group)

            assert key in maf_lookup, (
                "Missing MAF record", key
            )

            assert key not in seen_keys, (
                "Repeated variant-group record", key
            )

            seen_keys.add(key)
            m = maf_lookup[key]

            assert called == int(m["n_called"]), key
            assert missing == int(m["n_missing"]), key

            assert 2 * hom_ref + het == int(
                m["ref_allele_count"]
            ), key

            assert 2 * hom_alt + het == int(
                m["alt_allele_count"]
            ), key

            if called:
                p = (2 * hom_alt + het) / (2 * called)
                ho = het / called
                he = 2 * p * (1 - p)
            else:
                ho = None
                he = None

            site_rows.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "breed": group,
                "n_samples": sizes[group],
                "n_called": called,
                "n_missing": missing,
                "n_hom_ref": hom_ref,
                "n_het": het,
                "n_hom_alt": hom_alt,
                "ho": (
                    f"{ho:.10f}" if ho is not None else "NA"
                ),
                "he": (
                    f"{he:.10f}" if he is not None else "NA"
                ),
            })

            totals[group]["called"] += called
            totals[group]["missing"] += missing
            totals[group]["het"] += het
            totals[group]["sites_callable"] += int(called > 0)

            if he is not None:
                totals[group]["sum_he"] += he

assert len(site_rows) == 4505 * 17
assert seen_keys == set(maf_lookup)

assert totals["ALL"]["missing"] == 74353
assert totals["ALL"]["called"] == 1484377

summary_rows = []

for group in groups:
    t = totals[group]

    called = t["called"]
    callable_sites = t["sites_callable"]

    assert called + t["missing"] == sizes[group] * 4505

    summary_rows.append({
        "breed": group,
        "n_samples": sizes[group],
        "n_pathway_snp_records": 4505,
        "n_callable_snp_records": callable_sites,
        "n_called_genotypes": called,
        "n_missing_genotypes": t["missing"],
        "n_heterozygous_genotypes": t["het"],
        "call_rate": (
            f"{called / (sizes[group] * 4505):.10f}"
        ),
        "observed_heterozygosity": (
            f"{t['het'] / called:.10f}"
            if called else "NA"
        ),
        "mean_expected_heterozygosity": (
            f"{t['sum_he'] / callable_sites:.10f}"
            if callable_sites else "NA"
        ),
        "interpretation_caution": (
            "small_sample_n_lt_5"
            if sizes[group] < 5 else ""
        ),
    })


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


write_tsv(
    OUT / "H14_site_heterozygosity.tsv",
    site_rows,
)

write_tsv(
    OUT / "H14_heterozygosity_summary.tsv",
    summary_rows,
)

print("H14 heterozygosity analysis complete")
print("Site-group records:", len(site_rows))
print("Summary groups:", len(summary_rows))
print("All genotype counts match validated MAF tables")
print("All 4,505 SNP records retained")
print()
print("Group\tN\tHo\tMean_He\tCall_rate")

for row in summary_rows:
    print(
        f"{row['breed']}\t"
        f"{row['n_samples']}\t"
        f"{row['observed_heterozygosity']}\t"
        f"{row['mean_expected_heterozygosity']}\t"
        f"{row['call_rate']}"
    )
