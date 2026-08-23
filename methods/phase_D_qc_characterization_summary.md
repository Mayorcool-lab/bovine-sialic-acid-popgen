# Phase D: QC Characterization Summary

Real distributions from the actual dataset, both cohorts, 13 genes,
9,684 SNP sites (South Africa) / 6,842 SNP sites (Europe), before any
filtering applied.

| Metric | South Africa (n=120) | Europe (n=226) |
|---|---|---|
| Ts/Tv ratio | 2.48 | 2.52 |
| Per-animal missingness range | 0.01% - 3.8% | 0.01% - 1.05% |
| Median genotype depth (DP) | ~8-9x | ~9x |
| Median genotype quality (GQ) | ~48-53 | ~49-54 |
| DP < 4 rate (per gene) | 3-8% | 3-8% |
| GQ < 20 rate (per gene) | 1-5% | 1-2% |

## Interpretation

- Ts/Tv ratios (2.48, 2.52) both consistent with real coding-region variant
  data, not artifacts (genome-wide cattle expectation ~2.0-2.1; higher in
  coding regions is expected and reassuring).
- Per-animal missingness is low in both cohorts; no animal is an outlier
  requiring exclusion on this basis alone.
- DP and GQ distributions are consistent between cohorts, supporting a
  shared QC threshold rather than cohort-specific rules.

## Site missingness findings

Quantified distribution (not just top-N extremes):

| Metric | South Africa (n=9,684 sites) | Europe (n=6,842 sites) |
|---|---|---|
| Median F_MISS | 0.0 | 0.0 |
| 95th percentile | 0.84% | 0.22% |
| 99th percentile | 1.27% | 1.12% |
| Maximum | 83.2% | 79.4% |
| Sites >1% missing | 189 (1.95%) | 83 (1.21%) |
| Sites >5% missing | 19 (0.20%) | 19 (0.28%) |
| Sites >10% missing | 11 (0.11%) | 13 (0.19%) |
| Sites >20% missing | 3 (0.03%) | 8 (0.12%) |

Worst individual sites, both cohorts independently:

| Position | Gene | SA missing | EU missing |
|---|---|---|---|
| 23:32,674,499 | CMAH (~2.2kb into gene) | 83.2% | 79.4% |
| 7:100,067,228 | ST8SIA4 | 29.0% | 39.5% |
| 7:100,057,2xx cluster | ST8SIA4 | 10-13% | 22-29% |
| 14:7,879,4xx | ST3GAL1 | (not in SA top 10) | 21-23% |
| 29:29,736-748k | ST3GAL4 | 11-12% | (not in EU top 10) |

## Interpretation (site missingness)

- Site missingness is genuinely low across the large majority of sites in
  both cohorts: median F_MISS is 0.0, and 95% of sites have well under 1%
  missing genotypes. High-missingness sites are a small, quantified minority,
  not a general dataset-wide problem: only 0.11% (SA) and 0.19% (EU) of
  sites exceed 10% missingness.
- Position 23:32,674,499 (CMAH) is the most severe outlier in both cohorts
  independently (83.2% / 79.4% missing). The same site performing poorly
  in both independently-processed cohorts is a real, checkable observation.
  A specific technical cause (e.g. mappability or sequence complexity at
  this position) is a plausible hypothesis, not yet confirmed; this should
  be treated as a candidate region-specific technical issue pending further
  investigation of sequence context, not stated as an established mechanism.
- Missingness clusters in specific sub-regions of ST8SIA4, ST3GAL1, and
  ST3GAL4 in both cohorts. This consistent pattern across the two cohorts
  supports a real, region-specific cause rather than cohort-specific error,
  though the underlying mechanism has not been established.

## Implication for Phase E

- Missingness will not be a binding constraint at any reasonable threshold.
  A site call-rate threshold of >=90% (<=10% missing) would remove 11 sites
  in South Africa and 13 in Europe, together well under 0.2% of all sites
  in either cohort. This is a well-justified, minimally invasive threshold
  based on the actual quantified distribution, not an assumed convention.
- DP and GQ are the more meaningful levers for threshold selection. Given
  median DP ~8-9x and 3-8% of genotypes below DP=4, a DP>=4 threshold
  removes a modest, defensible minority of genotypes without being overly
  aggressive. Full per-gene breakdown in data/processed/qc_characterization/
  and results/qc/dp_gq_summary.tsv.

## Heterozygosity findings (breed-stratified)

Pooled-cohort heterozygosity (all breeds combined per continent) showed
elevated mean F (~+0.11-0.14 in both cohorts), consistent with the Wahlund
effect: pooling genetically distinct breeds inflates apparent homozygosity
relative to a single-population Hardy-Weinberg expectation.

Breed-stratified recalculation:

**Europe (10 breeds, n>=15 each; Mertolenga n=2 and Mirandesa n=4 excluded
as too small for stable estimates):**
All 10 breeds show mean F near zero (range -0.024 to +0.079), confirming
the pooled elevated F was a population-structure artifact, not real excess
inbreeding, once animals are compared only within their true breed.

**South Africa (4 breeds, n=20-40 each):**
Bonsmara (0.012) and Nguni (0.061) normalized similarly to the European
breeds. Afrikaner (0.156) and Tuli (0.129) remained elevated even after
breed-stratification. This may reflect genuine breed-specific effective
population size differences (plausible for Afrikaner, a composite breed
with documented smaller Ne) or panel-size noise (13-gene panel is a small,
non-random marker set, not ideal for precise inbreeding estimation).
Not treated as a data quality problem; noted for correct interpretation
if these breeds show notable patterns in later population-genetic analyses.

Individual-animal outlier check (pooled, 3 SD threshold): only 1 of 346
animals flagged. No evidence of individual sample-quality problems.

## Phase D: COMPLETE
All planned characterization done: DP, GQ, Ts/Tv, per-animal missingness,
per-site missingness (quantified distribution), site QUAL, heterozygosity
(pooled and breed-stratified). Ready for Phase E threshold justification.
