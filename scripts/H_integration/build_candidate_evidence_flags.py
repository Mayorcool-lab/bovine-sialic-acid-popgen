import csv
from pathlib import Path
from collections import Counter

src = Path(
    "results/H_integration/"
    "H_candidate_prioritisation_review_62.tsv"
)

out = Path(
    "results/H_integration/"
    "H_candidate_evidence_flags_62.tsv"
)

with src.open() as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if len(rows) != 62:
    raise RuntimeError(f"Expected 62 rows, observed {len(rows)}")


def yesno(x):
    return "YES" if x else "NO"


outrows = []

for r in rows:

    sift_del = r["SIFT_prediction"] == "deleterious"

    sift_del_lc = (
        r["SIFT_prediction"] ==
        "deleterious_low_confidence"
    )

    sift_unavailable = (
        r["SIFT_status"] == "SIFT_UNAVAILABLE"
    )

    high_impact = r["Consequence"] in {
        "stop_gained",
        "stop_lost",
    }

    r4s_assessable = (
        r["H5_rate4site_status"] == "ASSESSABLE"
    )

    if r4s_assessable:
        r4s = float(r["H5_rate4site_score"])
        negative_r4s = r4s < 0
    else:
        negative_r4s = False

    residue_context = (
        r["H6_context_level"].startswith(
            "RESIDUE_LEVEL"
        )
    )

    sa_obs = int(r["South_Africa_observed_OK"])
    eu_obs = int(r["Europe_observed_OK"])

    observed_sa = sa_obs > 0
    observed_eu = eu_obs > 0
    observed_both = observed_sa and observed_eu
    observed_any = observed_sa or observed_eu

    outrow = {
        "candidate_id": r["candidate_id"],
        "SYMBOL": r["SYMBOL"],
        "HGVSp": r["HGVSp"],
        "Consequence": r["Consequence"],

        "SIFT_prediction": r["SIFT_prediction"],
        "SIFT_score": r["SIFT_score"],
        "flag_SIFT_deleterious":
            yesno(sift_del),
        "flag_SIFT_deleterious_low_confidence":
            yesno(sift_del_lc),
        "flag_SIFT_unavailable":
            yesno(sift_unavailable),

        "flag_high_impact_consequence":
            yesno(high_impact),

        "H5_rate4site_status":
            r["H5_rate4site_status"],
        "H5_rate4site_score":
            r["H5_rate4site_score"],
        "flag_negative_Rate4Site":
            (
                yesno(negative_r4s)
                if r4s_assessable
                else "NOT_ASSESSABLE"
            ),

        "H6_context_level":
            r["H6_context_level"],
        "flag_residue_level_context":
            yesno(residue_context),

        "SA_AF": r["SA_AF"],
        "EU_AF": r["EU_AF"],

        "SA_observed_OK_breeds":
            r["South_Africa_observed_OK"],
        "EU_observed_OK_breeds":
            r["Europe_observed_OK"],

        "flag_ALT_observed_SA_OK":
            yesno(observed_sa),
        "flag_ALT_observed_EU_OK":
            yesno(observed_eu),
        "flag_ALT_observed_both":
            yesno(observed_both),
        "flag_ALT_observed_any_OK":
            yesno(observed_any),

        "breed_observation_pattern":
            r["breed_observation_pattern"],
    }

    outrows.append(outrow)


fieldnames = list(outrows[0].keys())

with out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(outrows)

print("Output:", out)
print("Rows:", len(outrows))
print(
    "Unique candidates:",
    len({r['candidate_id'] for r in outrows})
)

print("\n=== FACTUAL FLAG COUNTS ===")

for field in [
    "flag_SIFT_deleterious",
    "flag_SIFT_deleterious_low_confidence",
    "flag_SIFT_unavailable",
    "flag_high_impact_consequence",
    "flag_negative_Rate4Site",
    "flag_residue_level_context",
    "flag_ALT_observed_SA_OK",
    "flag_ALT_observed_EU_OK",
    "flag_ALT_observed_both",
    "flag_ALT_observed_any_OK",
]:
    print(
        field,
        dict(sorted(
            Counter(r[field] for r in outrows).items()
        ))
    )

print("\n=== KEY CONVERGENCE COMBINATIONS ===")

def count_where(fn):
    return sum(1 for r in outrows if fn(r))

print(
    "SIFT deleterious + negative R4S:",
    count_where(
        lambda r:
        r["flag_SIFT_deleterious"] == "YES"
        and r["flag_negative_Rate4Site"] == "YES"
    )
)

print(
    "SIFT deleterious + negative R4S + ALT observed:",
    count_where(
        lambda r:
        r["flag_SIFT_deleterious"] == "YES"
        and r["flag_negative_Rate4Site"] == "YES"
        and r["flag_ALT_observed_any_OK"] == "YES"
    )
)

print(
    "SIFT deleterious + negative R4S + residue context + ALT observed:",
    count_where(
        lambda r:
        r["flag_SIFT_deleterious"] == "YES"
        and r["flag_negative_Rate4Site"] == "YES"
        and r["flag_residue_level_context"] == "YES"
        and r["flag_ALT_observed_any_OK"] == "YES"
    )
)

print(
    "Negative R4S + ALT observed both populations:",
    count_where(
        lambda r:
        r["flag_negative_Rate4Site"] == "YES"
        and r["flag_ALT_observed_both"] == "YES"
    )
)

print(
    "High-impact + ALT observed adequately sampled:",
    count_where(
        lambda r:
        r["flag_high_impact_consequence"] == "YES"
        and r["flag_ALT_observed_any_OK"] == "YES"
    )
)

print("\nComposite score added: 0")
print("Ranking added: 0")
print("Prioritisation decision added: 0")
print("FACTUAL EVIDENCE-FLAG TABLE: PASS")
