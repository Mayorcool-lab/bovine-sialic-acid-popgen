#!/usr/bin/env python3
"""Generate Supplementary Table S4 from validated H9 breed-level FST results."""

from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

SOURCE = (
    ROOT / "results/H_integration/H9_candidate_breed_fst/"
    "candidate_breed_pairwise_fst_62.tsv"
)

MANIFEST = (
    ROOT / "results/H_integration/H9_candidate_breed_fst/"
    "breed_sample_manifest_346.tsv"
)

OUTPUT = (
    ROOT / "results/figures/source_data/"
    "Supplementary_Table_S4_candidate_breed_fst.tsv"
)

EXPECTED_STATUSES = {
    "estimated": 2028,
    "not_estimable_vcftools": 2463,
    "absent_breed_a_cohort_vcf": 2168,
    "absent_breed_b_cohort_vcf": 766,
    "insufficient_called_genotypes": 15,
}

def main():
    df = pd.read_csv(SOURCE, sep="\t", dtype={"CHROM": str})
    manifest = pd.read_csv(MANIFEST, sep="\t")

    assert len(manifest) == 346
    assert manifest["vcf_id"].is_unique

    breeds = sorted(manifest["breed"].unique())
    assert len(breeds) == 16

    assert len(df) == 7440
    assert df["candidate_id"].nunique() == 62

    expected_pairs = set(combinations(breeds, 2))
    actual_pairs = set(zip(df["breed_a"], df["breed_b"]))
    assert actual_pairs == expected_pairs

    assert not df.duplicated(
        ["candidate_id", "breed_a", "breed_b"]
    ).any()

    assert df.groupby("candidate_id").size().eq(120).all()
    assert df.groupby(["breed_a", "breed_b"]).size().eq(62).all()

    counts = df["fst_status"].value_counts().to_dict()
    assert counts == EXPECTED_STATUSES, counts

    estimated = df["fst_status"].eq("estimated")
    fst = pd.to_numeric(df["fst_wc"], errors="coerce")
    independent = pd.to_numeric(
        df["fst_wc_independent"], errors="coerce"
    )

    assert np.isfinite(fst[estimated]).all()
    assert np.isfinite(independent[estimated]).all()
    assert (
        np.abs(fst[estimated] - independent[estimated]) <= 1e-5
    ).all()
    assert not np.isfinite(fst[~estimated]).any()

    breed_sizes = manifest["breed"].value_counts()

    assert df["breed_a_n"].eq(
        df["breed_a"].map(breed_sizes)
    ).all()
    assert df["breed_b_n"].eq(
        df["breed_b"].map(breed_sizes)
    ).all()

    expected_small = (
        (df["breed_a_n"] < 5) |
        (df["breed_b_n"] < 5)
    )
    actual_small = (
        df["small_breed_pair"].astype(str).str.lower()
        .isin(["true", "1", "yes"])
    )
    assert expected_small.eq(actual_small).all()

    for side in ("a", "b"):
        called = df[f"breed_{side}_called"]
        hom_ref = df[f"breed_{side}_hom_ref"]
        het = df[f"breed_{side}_het"]
        hom_alt = df[f"breed_{side}_hom_alt"]
        af = df[f"breed_{side}_alt_af"]

        present = (
            called.notna() & hom_ref.notna() &
            het.notna() & hom_alt.notna()
        )

        assert (
            called[present] ==
            (hom_ref + het + hom_alt)[present]
        ).all()

        valid = present & (called > 0) & af.notna()
        expected_af = (
            (het[valid] + 2 * hom_alt[valid]) /
            (2 * called[valid])
        )

        assert (
            np.abs(af[valid] - expected_af) <= 1e-5
        ).all()

    undefined = df["fst_status"].eq("not_estimable_vcftools")
    assert (df.loc[undefined, "breed_a_called"] >= 2).all()
    assert (df.loc[undefined, "breed_b_called"] >= 2).all()
    assert (df.loc[undefined, "breed_a_alt_af"] == 0).all()
    assert (df.loc[undefined, "breed_b_alt_af"] == 0).all()

    # Retain all original H9 fields and numerical values.
    table = df.copy()

    table = table.sort_values(
        ["SYMBOL", "CHROM", "POS", "breed_a", "breed_b"],
        kind="stable"
    ).reset_index(drop=True)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(
        OUTPUT,
        sep="\t",
        index=False,
        na_rep="NA",
        float_format="%.9f"
    )

    # Verify the generated file.
    saved = pd.read_csv(OUTPUT, sep="\t", dtype={"CHROM": str})

    assert len(saved) == 7440
    assert saved["fst_status"].value_counts().to_dict() == EXPECTED_STATUSES

    print("PASS: Supplementary Table S4 generated")
    print("Output:", OUTPUT.relative_to(ROOT))
    print("Rows:", len(saved))
    print("Candidates:", saved["candidate_id"].nunique())
    print("Breed pairs:", len(expected_pairs))
    print("Finite FST estimates:", int(estimated.sum()))
    print("Undefined/reference-fixed:", int(undefined.sum()))
    print("Small-breed comparisons:", int(expected_small.sum()))
    print("Original H9 values preserved.")

if __name__ == "__main__":
    main()
