#!/usr/bin/env python3
"""H7: Reproducible regional Weir-Cockerham FST for 62 candidates."""

from pathlib import Path
from collections import Counter
import math
import subprocess
import tempfile

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / (
    "results/H6/H6.2/processed/"
    "H6.2_complete_62_candidate_functional_context.tsv"
)
VCFS = ROOT / "Aiyesa_handoff_2026/qc_filtered_corrected"
OUT = ROOT / "results/H_integration/H7_candidate_regional_fst"

EXPECTED = {"europe": 226, "south_africa": 120}
NEU3 = "15:54265634:T>C"


def command(args):
    return subprocess.run(
        args, text=True, capture_output=True, check=True
    ).stdout


def samples(vcf):
    return command(["bcftools", "query", "-l", str(vcf)]).splitlines()


def normalize(gt):
    if gt in (".", "./.", ".|."):
        return "./."
    gt = gt.replace("|", "/")
    if gt == "1/0":
        gt = "0/1"
    if gt not in ("0/0", "0/1", "1/1"):
        raise ValueError(f"Unsupported genotype: {gt}")
    return gt


def extract(vcf, chrom, pos, ref, alt, n):
    text = command([
        "bcftools", "query",
        "-r", f"{chrom}:{pos}-{pos}",
        "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
        str(vcf),
    ])

    matching = []

    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue

        key = (
            str(fields[0]),
            int(fields[1]),
            fields[2],
            fields[3],
        )

        if key == (chrom, pos, ref, alt):
            matching.append(fields[4:])

    if not matching:
        return None

    if len(matching) != 1:
        raise ValueError(f"Duplicate exact variant: {chrom}:{pos}")

    if len(matching[0]) != n:
        raise ValueError(f"Genotype/sample mismatch: {chrom}:{pos}")

    return [normalize(gt) for gt in matching[0]]


def stats(calls):
    if calls is None:
        return dict.fromkeys(
            ("called", "missing", "hom_ref", "het", "hom_alt", "alt_af")
        )

    c = Counter(calls)
    n = c["0/0"] + c["0/1"] + c["1/1"]

    return {
        "called": n,
        "missing": c["./."],
        "hom_ref": c["0/0"],
        "het": c["0/1"],
        "hom_alt": c["1/1"],
        "alt_af": (c["0/1"] + 2*c["1/1"]) / (2*n) if n else None,
    }


def independent_wc(left, right):
    """Diploid biallelic Weir-Cockerham 1984 single-site estimator."""
    import numpy as np

    counts = np.array([
        [left["hom_ref"], left["het"], left["hom_alt"]],
        [right["hom_ref"], right["het"], right["hom_alt"]],
    ], dtype=float)

    n = counts.sum(axis=1)
    if (n < 2).any():
        return float("nan")

    r = len(n)
    p = (counts[:, 1] + 2*counts[:, 2]) / (2*n)
    h = counts[:, 1] / n
    nb = n.mean()
    pb = (n*p).sum() / n.sum()
    hb = (n*h).sum() / n.sum()
    nc = (n.sum() - (n*n).sum()/n.sum()) / (r-1)
    s2 = (n*(p-pb)**2).sum() / ((r-1)*nb)

    a = (nb/nc) * (
        s2 - (pb*(1-pb) - ((r-1)/r)*s2 - hb/4)/(nb-1)
    )
    b = (nb/(nb-1)) * (
        pb*(1-pb) - ((r-1)/r)*s2
        - ((2*nb-1)/(4*nb))*hb
    )
    c = hb/2
    denominator = a+b+c

    return a/denominator if denominator != 0 else float("nan")


def main():
    df = pd.read_csv(INPUT, sep="\t", dtype={"CHROM": str})
    required = ["candidate_id", "SYMBOL", "CHROM", "POS", "REF", "ALT"]

    if len(df) != 62 or not set(required).issubset(df.columns):
        raise ValueError("Candidate table failed schema/row-count validation")

    if df.duplicated(["CHROM", "POS", "REF", "ALT"]).any():
        raise ValueError("Duplicate candidate variants")

    genes = sorted(df["SYMBOL"].unique())
    paths = {}
    cohort_ids = {}

    for cohort, expected in EXPECTED.items():
        for gene in genes:
            path = VCFS / cohort / f"{gene}_qc.vcf.gz"
            if not path.is_file():
                raise FileNotFoundError(path)

            ids = samples(path)
            if len(ids) != expected or len(set(ids)) != expected:
                raise ValueError(f"Invalid sample list: {cohort}/{gene}")

            if cohort not in cohort_ids:
                cohort_ids[cohort] = ids
            elif ids != cohort_ids[cohort]:
                raise ValueError(f"Sample order differs: {cohort}/{gene}")

            paths[(cohort, gene)] = path

    eu_ids = cohort_ids["europe"]
    sa_ids = cohort_ids["south_africa"]

    if set(eu_ids) & set(sa_ids):
        raise ValueError("Overlapping sample IDs")

    rows = []
    eligible = []

    for _, candidate in df.iterrows():
        gene = str(candidate["SYMBOL"])
        chrom = str(candidate["CHROM"])
        pos = int(candidate["POS"])
        ref = str(candidate["REF"])
        alt = str(candidate["ALT"])
        cid = str(candidate["candidate_id"])

        calls = {
            cohort: extract(
                paths[(cohort, gene)], chrom, pos, ref, alt,
                EXPECTED[cohort]
            )
            for cohort in EXPECTED
        }

        eu = stats(calls["europe"])
        sa = stats(calls["south_africa"])

        row = {
            "candidate_id": cid,
            "SYMBOL": gene,
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "europe_record_present": calls["europe"] is not None,
            "south_africa_record_present": calls["south_africa"] is not None,
        }

        for cohort, st in (("europe", eu), ("south_africa", sa)):
            row.update({f"{cohort}_{k}": v for k, v in st.items()})

        row["regional_fst_wc"] = None
        row["regional_fst_wc_independent"] = None

        if calls["europe"] is None:
            row["fst_status"] = "absent_europe_record"
        elif calls["south_africa"] is None:
            row["fst_status"] = "absent_south_africa_record"
        elif eu["called"] < 2 or sa["called"] < 2:
            row["fst_status"] = "insufficient_called_genotypes"
        else:
            row["fst_status"] = "pending"
            eligible.append((row, calls["europe"] + calls["south_africa"], eu, sa))

        rows.append(row)

    eligible.sort(
        key=lambda item: (
            int(item[0]["CHROM"]) if item[0]["CHROM"].isdigit() else 10**9,
            item[0]["POS"],
        )
    )

    with tempfile.TemporaryDirectory(prefix="H7_fst_") as directory:
        tmp = Path(directory)

        for cohort, ids in cohort_ids.items():
            (tmp / f"{cohort}.samples").write_text(
                "\n".join(ids) + "\n"
            )

        vcf = tmp / "candidate_genotypes.vcf"

        with vcf.open("w") as f:
            f.write("##fileformat=VCFv4.2\n")

            for chrom in sorted({x[0]["CHROM"] for x in eligible}):
                f.write(f"##contig=<ID={chrom}>\n")

            f.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,'
                'Description="Diploid genotype">\n'
            )
            f.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(eu_ids + sa_ids) + "\n"
            )

            for row, calls, _, _ in eligible:
                f.write("\t".join([
                    row["CHROM"], str(row["POS"]),
                    row["candidate_id"], row["REF"], row["ALT"],
                    ".", "PASS", ".", "GT", *calls
                ]) + "\n")

        prefix = tmp / "regional"

        process = subprocess.run([
            "vcftools",
            "--vcf", str(vcf),
            "--weir-fst-pop", str(tmp / "europe.samples"),
            "--weir-fst-pop", str(tmp / "south_africa.samples"),
            "--out", str(prefix),
        ], text=True, capture_output=True)

        if process.returncode:
            raise RuntimeError(process.stdout + process.stderr)

        fst_file = Path(str(prefix) + ".weir.fst")

        # CRITICAL FIX: preserve chromosome identifiers as strings.
        fst_df = pd.read_csv(
            fst_file,
            sep=r"\s+",
            dtype={"CHROM": str},
        )

        lookup = {}

        for record in fst_df.itertuples(index=False):
            key = (str(record.CHROM), int(record.POS))
            if key in lookup:
                raise ValueError(f"Duplicate VCFtools output: {key}")
            lookup[key] = float(record.WEIR_AND_COCKERHAM_FST)

        expected_keys = {
            (row["CHROM"], row["POS"])
            for row, _, _, _ in eligible
        }

        unexpected = set(lookup) - expected_keys

        if unexpected:
            raise ValueError(
                f"Unexpected VCFtools positions: {sorted(unexpected)[:10]}"
            )

        for row, _, eu, sa in eligible:
            key = (row["CHROM"], row["POS"])
            value = lookup.get(key, float("nan"))
            independent = independent_wc(eu, sa)

            if not math.isfinite(value):
                row["fst_status"] = "not_estimable_vcftools"
                continue

            if not math.isfinite(independent):
                raise ValueError(f"Independent FST undefined: {key}")

            if abs(value - independent) > 1e-5:
                raise ValueError(
                    f"FST disagreement at {key}: "
                    f"VCFtools={value}, independent={independent}"
                )

            row["regional_fst_wc"] = value
            row["regional_fst_wc_independent"] = independent
            row["fst_status"] = "estimated"

    result = pd.DataFrame(rows)

    pilot = result[result["candidate_id"] == NEU3]

    if len(pilot) != 1:
        raise ValueError("NEU3 pilot missing")

    pilot = pilot.iloc[0]

    if (
        pilot["fst_status"] != "estimated"
        or abs(float(pilot["regional_fst_wc"]) - 0.369494) > 1e-5
        or int(pilot["europe_called"]) != 219
        or int(pilot["south_africa_called"]) != 116
    ):
        raise ValueError("NEU3 pilot validation failed")

    OUT.mkdir(parents=True, exist_ok=True)

    output = OUT / "candidate_regional_fst_62.tsv"
    result.to_csv(
        output, sep="\t", index=False,
        na_rep="NA", float_format="%.9f"
    )

    counts = result["fst_status"].value_counts()

    summary = [
        "H7: Candidate regional Weir-Cockerham FST",
        f"Total candidates: {len(result)}",
        f"VCFtools eligible: {len(eligible)}",
        f"NEU3 validated FST: {float(pilot['regional_fst_wc']):.6f}",
        "",
        "Status counts:",
        *[f"{status}: {count}" for status, count in counts.items()],
        "",
        "Absent variant-only VCF records were not imputed as reference.",
        "FST values are raw, including negative estimates.",
        "Every finite FST estimate was checked independently.",
    ]

    (OUT / "candidate_regional_fst_run_summary.txt").write_text(
        "\n".join(summary) + "\n"
    )

    print("\n".join(summary))
    print(f"\nOutput: {output}")
    print("\nPASS: H7 completed and NEU3 independently validated")


if __name__ == "__main__":
    main()
