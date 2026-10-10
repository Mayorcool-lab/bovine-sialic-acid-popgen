#!/usr/bin/env python3
"""H8: Six pairwise country-level Weir-Cockerham FST comparisons."""

from pathlib import Path
from itertools import combinations
from collections import Counter
import math
import subprocess
import tempfile

import pandas as pd

from H7_candidate_regional_fst import (
    ROOT, INPUT, VCFS, EXPECTED,
    samples, extract, stats, independent_wc,
)

OUT = ROOT / "results/H_integration/H8_candidate_country_fst"

COUNTRIES = (
    "Finland",
    "Netherlands",
    "Portugal",
    "South Africa",
)

EXPECTED_COUNTRIES = {
    "Finland": 80,
    "Netherlands": 120,
    "Portugal": 26,
    "South Africa": 120,
}

COHORT = {
    "Finland": "europe",
    "Netherlands": "europe",
    "Portugal": "europe",
    "South Africa": "south_africa",
}


def load_manifest():
    corrected = pd.read_csv(
        ROOT / "Aiyesa_handoff_2026/europe_226_breed_mapping.tsv",
        sep="\t", dtype=str,
    )

    original = pd.read_csv(
        ROOT / "data/metadata/sample_metadata_european.tsv",
        sep="\t", dtype=str,
    )

    african = pd.read_csv(
        ROOT / "data/metadata/sample_metadata_south_africa.tsv",
        sep="\t", dtype=str,
    )

    europe = corrected.merge(
        original[["ena_sample_accession", "country"]],
        on="ena_sample_accession",
        how="left",
        validate="one_to_one",
    )

    if europe["country"].isna().any():
        raise ValueError("Unmatched European country assignments")

    europe["cohort"] = "europe"
    african["cohort"] = "south_africa"

    columns = [
        "vcf_id", "ena_sample_accession",
        "breed", "country", "cohort",
    ]

    manifest = pd.concat(
        [europe[columns], african[columns]],
        ignore_index=True,
    )

    if len(manifest) != 346:
        raise ValueError("Expected exactly 346 metadata samples")

    if manifest["vcf_id"].duplicated().any():
        raise ValueError("Duplicate VCF sample IDs in manifest")

    observed = manifest["country"].value_counts().to_dict()

    if observed != EXPECTED_COUNTRIES:
        raise ValueError(f"Country counts differ: {observed}")

    if (
        manifest.loc[
            manifest["cohort"] == "south_africa", "country"
        ] != "South Africa"
    ).any():
        raise ValueError("South African cohort/country mismatch")

    return manifest


def main():
    candidates = pd.read_csv(
        INPUT, sep="\t", dtype={"CHROM": str}
    )

    if len(candidates) != 62:
        raise ValueError("Expected 62 candidate variants")

    required = {
        "candidate_id", "SYMBOL",
        "CHROM", "POS", "REF", "ALT",
    }

    if not required.issubset(candidates.columns):
        raise ValueError("Candidate table missing required columns")

    if candidates.duplicated(
        ["CHROM", "POS", "REF", "ALT"]
    ).any():
        raise ValueError("Duplicate candidate variants")

    manifest = load_manifest()
    genes = sorted(candidates["SYMBOL"].unique())

    paths = {}
    cohort_ids = {}

    for cohort, expected_n in EXPECTED.items():
        for gene in genes:
            path = VCFS / cohort / f"{gene}_qc.vcf.gz"

            if not path.is_file():
                raise FileNotFoundError(path)

            ids = samples(path)

            if len(ids) != expected_n or len(set(ids)) != expected_n:
                raise ValueError(
                    f"Invalid VCF sample list: {cohort}/{gene}"
                )

            if cohort not in cohort_ids:
                cohort_ids[cohort] = ids
            elif ids != cohort_ids[cohort]:
                raise ValueError(
                    f"Different sample order: {cohort}/{gene}"
                )

            paths[(cohort, gene)] = path

    if set(cohort_ids["europe"]) & set(
        cohort_ids["south_africa"]
    ):
        raise ValueError("Overlapping cohort sample IDs")

    country_ids = {}

    for country in COUNTRIES:
        cohort = COHORT[country]

        assigned = set(
            manifest.loc[
                manifest["country"] == country, "vcf_id"
            ]
        )

        ids = [
            sid for sid in cohort_ids[cohort]
            if sid in assigned
        ]

        if len(ids) != EXPECTED_COUNTRIES[country]:
            raise ValueError(
                f"VCF country membership mismatch: {country}"
            )

        country_ids[country] = ids

    all_vcf_ids = set(
        cohort_ids["europe"] + cohort_ids["south_africa"]
    )

    if all_vcf_ids != set(manifest["vcf_id"]):
        raise ValueError("VCF IDs and manifest do not match")

    # Extract every candidate once from each cohort.
    extracted = {}

    for _, candidate in candidates.iterrows():
        cid = str(candidate["candidate_id"])
        gene = str(candidate["SYMBOL"])
        chrom = str(candidate["CHROM"])
        pos = int(candidate["POS"])
        ref = str(candidate["REF"])
        alt = str(candidate["ALT"])

        extracted[cid] = {}

        for cohort in EXPECTED:
            calls = extract(
                paths[(cohort, gene)],
                chrom, pos, ref, alt,
                EXPECTED[cohort],
            )

            extracted[cid][cohort] = (
                None if calls is None
                else dict(zip(cohort_ids[cohort], calls))
            )

    results = []

    with tempfile.TemporaryDirectory(
        prefix="H8_country_fst_"
    ) as directory:

        tmp = Path(directory)

        for country_a, country_b in combinations(COUNTRIES, 2):
            pair_name = (
                country_a.replace(" ", "_")
                + "__"
                + country_b.replace(" ", "_")
            )

            ids_a = country_ids[country_a]
            ids_b = country_ids[country_b]

            (tmp / "a.samples").write_text(
                "\n".join(ids_a) + "\n"
            )

            (tmp / "b.samples").write_text(
                "\n".join(ids_b) + "\n"
            )

            eligible = []

            for _, candidate in candidates.iterrows():
                cid = str(candidate["candidate_id"])
                chrom = str(candidate["CHROM"])
                pos = int(candidate["POS"])

                cohort_a = COHORT[country_a]
                cohort_b = COHORT[country_b]

                source_a = extracted[cid][cohort_a]
                source_b = extracted[cid][cohort_b]

                calls_a = (
                    None if source_a is None
                    else [source_a[sid] for sid in ids_a]
                )

                calls_b = (
                    None if source_b is None
                    else [source_b[sid] for sid in ids_b]
                )

                st_a = stats(calls_a)
                st_b = stats(calls_b)

                row = {
                    "candidate_id": cid,
                    "SYMBOL": str(candidate["SYMBOL"]),
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": str(candidate["REF"]),
                    "ALT": str(candidate["ALT"]),
                    "country_a": country_a,
                    "country_b": country_b,
                    "comparison": pair_name,
                    "country_a_n": len(ids_a),
                    "country_b_n": len(ids_b),
                    "country_a_called": st_a["called"],
                    "country_b_called": st_b["called"],
                    "country_a_alt_af": st_a["alt_af"],
                    "country_b_alt_af": st_b["alt_af"],
                    "country_a_hom_ref": st_a["hom_ref"],
                    "country_a_het": st_a["het"],
                    "country_a_hom_alt": st_a["hom_alt"],
                    "country_b_hom_ref": st_b["hom_ref"],
                    "country_b_het": st_b["het"],
                    "country_b_hom_alt": st_b["hom_alt"],
                    "fst_wc": None,
                    "fst_wc_independent": None,
                }

                if source_a is None:
                    row["fst_status"] = "absent_country_a_cohort_vcf"
                elif source_b is None:
                    row["fst_status"] = "absent_country_b_cohort_vcf"
                elif st_a["called"] < 2 or st_b["called"] < 2:
                    row["fst_status"] = "insufficient_called_genotypes"
                else:
                    row["fst_status"] = "pending"
                    eligible.append(
                        (row, calls_a + calls_b, st_a, st_b)
                    )

                results.append(row)

            eligible.sort(
                key=lambda x: (
                    int(x[0]["CHROM"])
                    if x[0]["CHROM"].isdigit()
                    else 10**9,
                    x[0]["POS"],
                )
            )

            if not eligible:
                continue

            vcf_path = tmp / "pair.vcf"

            with vcf_path.open("w") as f:
                f.write("##fileformat=VCFv4.2\n")

                for chrom in sorted({
                    item[0]["CHROM"] for item in eligible
                }):
                    f.write(f"##contig=<ID={chrom}>\n")

                f.write(
                    '##FORMAT=<ID=GT,Number=1,Type=String,'
                    'Description="Diploid genotype">\n'
                )

                f.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
                    "\tINFO\tFORMAT\t"
                    + "\t".join(ids_a + ids_b)
                    + "\n"
                )

                for row, calls, _, _ in eligible:
                    f.write("\t".join([
                        row["CHROM"],
                        str(row["POS"]),
                        row["candidate_id"],
                        row["REF"],
                        row["ALT"],
                        ".", "PASS", ".", "GT",
                        *calls,
                    ]) + "\n")

            prefix = tmp / "pair_result"

            process = subprocess.run(
                [
                    "vcftools",
                    "--vcf", str(vcf_path),
                    "--weir-fst-pop", str(tmp / "a.samples"),
                    "--weir-fst-pop", str(tmp / "b.samples"),
                    "--out", str(prefix),
                ],
                capture_output=True,
                text=True,
            )

            if process.returncode:
                raise RuntimeError(
                    f"{pair_name} failed:\n"
                    + process.stdout + process.stderr
                )

            fst_file = Path(str(prefix) + ".weir.fst")

            fst_df = pd.read_csv(
                fst_file,
                sep=r"\s+",
                dtype={"CHROM": str},
            )

            lookup = {}

            for record in fst_df.itertuples(index=False):
                key = (str(record.CHROM), int(record.POS))

                if key in lookup:
                    raise ValueError(
                        f"Duplicate VCFtools position: {pair_name} {key}"
                    )

                lookup[key] = float(
                    record.WEIR_AND_COCKERHAM_FST
                )

            expected_keys = {
                (row["CHROM"], row["POS"])
                for row, _, _, _ in eligible
            }

            unexpected = set(lookup) - expected_keys

            if unexpected:
                raise ValueError(
                    f"Unexpected positions: {pair_name} {unexpected}"
                )

            for row, _, st_a, st_b in eligible:
                key = (row["CHROM"], row["POS"])
                value = lookup.get(key, float("nan"))
                independent = independent_wc(st_a, st_b)

                if not math.isfinite(value):
                    row["fst_status"] = "not_estimable_vcftools"

                    if math.isfinite(independent):
                        raise ValueError(
                            f"VCFtools missing finite estimate: "
                            f"{pair_name} {key}"
                        )

                    continue

                if not math.isfinite(independent):
                    raise ValueError(
                        f"Independent estimate undefined: "
                        f"{pair_name} {key}"
                    )

                if abs(value - independent) > 1e-5:
                    raise ValueError(
                        f"Estimator disagreement: {pair_name} "
                        f"{key}: {value} vs {independent}"
                    )

                row["fst_wc"] = value
                row["fst_wc_independent"] = independent
                row["fst_status"] = "estimated"

            print(
                f"{pair_name}: "
                f"{sum(x[0]['fst_status'] == 'estimated' for x in eligible)} "
                f"finite estimates"
            )

    result = pd.DataFrame(results)

    if len(result) != 62 * 6:
        raise ValueError("Expected 372 candidate-pair rows")

    if result.duplicated(
        ["candidate_id", "comparison"]
    ).any():
        raise ValueError("Duplicate candidate-pair records")

    if (result["fst_status"] == "pending").any():
        raise ValueError("Unresolved pending records")

    OUT.mkdir(parents=True, exist_ok=True)

    result.to_csv(
        OUT / "candidate_country_pairwise_fst_62.tsv",
        sep="\t", index=False,
        na_rep="NA", float_format="%.9f",
    )

    manifest.to_csv(
        OUT / "country_sample_manifest_346.tsv",
        sep="\t", index=False,
    )

    summary = [
        "H8: Candidate-level country pairwise FST",
        "Samples: Finland=80, Netherlands=120, "
        "Portugal=26, South Africa=120",
        "Candidates: 62",
        "Country pairs: 6",
        f"Candidate-pair records: {len(result)}",
        "",
        "Status counts:",
        *[
            f"{status}: {count}"
            for status, count in
            result["fst_status"].value_counts().items()
        ],
        "",
        "Estimated counts by comparison:",
        *[
            f"{pair}: {count}"
            for pair, count in (
                result.loc[
                    result["fst_status"] == "estimated"
                ].groupby("comparison").size().items()
            )
        ],
        "",
        "Finite estimates independently checked within 1e-5.",
        "Absent cohort variant records never imputed as reference.",
        "Raw negative FST estimates retained.",
        "Country groups may be confounded by breed composition.",
    ]

    (OUT / "candidate_country_fst_run_summary.txt").write_text(
        "\n".join(summary) + "\n"
    )

    print("\n" + "\n".join(summary))
    print("\nPASS: H8 completed")


if __name__ == "__main__":
    main()
