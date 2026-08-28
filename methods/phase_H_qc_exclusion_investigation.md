# Phase H: QC-Excluded Candidate Variant Investigation

## Why this investigation happened

Phase H (H1) initially extracted genotypes for all 62 protein-altering
candidates and found fewer records than the theoretical maximum
(15,506 instead of 21,452). H1b classified every candidate x cohort
combination into three categories:

| Status | Count |
|---|---|
| PRESENT (in final QC VCF) | 93 |
| NOT_IN_SOURCE (absent even before QC) | 27 |
| ABSENT_AFTER_QC (existed, removed by QC) | 4 |

Total: 93+27+4 = 124 = 62 variants x 2 cohorts. Complete, no loss.

This document covers the 4 ABSENT_AFTER_QC cases, all in South Africa.

## The 4 ABSENT_AFTER_QC candidates, fully explained

All four fail the Phase E site call-rate threshold (>=90%), confirmed
by direct per-genotype DP/GQ audit against the actual QC rule:

| Gene | Position | Callable/120 | Call rate | Genotype pattern (callable) |
|---|---|---|---|---|
| ST3GAL4 | 29:29,729,953 | 107 | 89.17% | 41 0/0, 46 0/1, 20 1/1 (common variant) |
| CMAS | 5:88,162,466 | 107 | 89.17% | 105 0/0, 2 0/1 (rare) |
| NANS | 8:63,212,470 | 105 | 87.50% | 96 0/0, 9 0/1 (uncommon) |
| CMAH (F512Y) | 23:32,745,866 | 98 | 81.67% | see below |

All fail by a narrow margin (81.7-89.2%), not a large gap from the 90%
threshold. This is a real, correctly-applied QC exclusion, not an
error, but it removes real, biologically meaningful variation. ST3GAL4
in particular is a common variant (66/107 callable animals are
carriers) lost purely on a borderline call-rate technicality.

## CMAH F512Y: full accounting (verified, no file discrepancy)

Initial audit of two different files appeared to show different carrier
counts (10 vs. 16), raising a concern of a possible file inconsistency.
Direct verification (comparing all 120 raw genotype:DP:GQ values between
`CMAH_norm.vcf.gz` and `CMAH_biallelic_snps.vcf.gz` at this exact
position) confirmed ZERO mismatches: both files contain identical raw
data. There was no file discrepancy.

The two numbers instead reflected two different, both correct,
quantities:

- **16 raw non-reference genotype calls** (15 heterozygous, 1
  homozygous-alt) before any DP/GQ filtering
- **10 of those 16 pass DP>=4 and GQ>=20** (the remaining 6, including
  the sole 1/1 call, fail genotype-level QC: e.g. AFR21 DP=1/GQ=0,
  AFR17 GQ=12, AFR32 GQ=11, Bonsmara4 GQ=11, Nguni21 GQ=0, NGI31 GQ=11)
- Across all 120 animals: 98 callable genotypes total (site call rate
  81.67%), below the 90% threshold, hence exclusion from the final QC
  VCF

**Correct summary statement:** F512Y has 16 raw non-reference genotype
calls in South Africa. After the predefined DP>=4/GQ>=20 genotype
filters, 6 of these fail QC (including the only homozygous-alt call),
leaving 10 reliable heterozygous carriers among 98/120 total callable
genotypes (81.67% site call rate), below the 90% threshold. F512Y
therefore did not pass Phase E QC in South Africa; it passed in Europe.

No correction to G9's biological conclusion is needed (F512Y fails QC
in South Africa, passes in Europe, is unchanged and correct). If G9's
original text stated a specific carrier count, it should be checked
against this fuller accounting before being cited further, but the core
finding stands.

## Governing rule established for Phase H genotype matrix construction

Absence of a variant record in a cohort's final QC-filtered VCF is NOT
automatically equivalent to 0/0 (homozygous reference) for every animal
in that cohort. Three distinct situations exist, and must be recorded
as distinct in the master genotype/frequency table, not silently
collapsed:

1. **PRESENT** — genuine data exists, use it directly.
2. **ABSENT_AFTER_QC** — real variation exists but did not meet the
   pre-defined QC threshold as a complete site. Report using whatever
   partial, callable data exists (e.g. ST3GAL4: 107/120 callable),
   explicitly labeled as sub-threshold and reported alongside its real
   call rate, not silently omitted and not treated as a full 120-animal
   0/0 result.
3. **NOT_IN_SOURCE** — genuinely absent even before QC (either truly
   invariant in that cohort, or excluded at an earlier pipeline stage
   than QC). Still not automatically 0/0 without further verification;
   treatment of these 27 cases is the next open item.

## Still open

- Define the specific rule for how the 27 NOT_IN_SOURCE cases should be
  encoded in the final Phase H genotype/frequency matrix.
- Decide manuscript-level treatment of the 4 ABSENT_AFTER_QC variants:
  report with explicit call-rate caveat (recommended, especially for
  ST3GAL4 given its real, common carrier frequency) vs. exclude from
  cross-cohort frequency comparison entirely.

## H1d CORRECTION: NOT_IN_SOURCE cases cannot be confirmed as 0/0

An earlier version of this section concluded that nearby variant-
calling activity (confirmed for all 27/27 cases, most within 15-50bp)
justified encoding these as confirmed 0/0 genotypes. This conclusion
was too strong and has been corrected.

**The actual, defensible interpretation:** nearby calls demonstrate the
surrounding region was actively covered and variant-called, arguing
against a gross regional data-absence problem. They do NOT establish
reference-confidence at the exact candidate position itself. Variant-
only VCFs emit no record whatsoever for positions where no alternate
allele was called, this is structurally different from a reference-
confidence VCF or gVCF, which would explicitly confirm 0/0 with its own
quality metrics. Absence of a record here means "no ALT detected," not
"confirmed REF in every animal." Treating it as the latter would assert
more certainty than the data supports.

**Establishing true site-level reference confidence would require
BAM/gVCF-level information or re-calling, which is disproportionate to
what this preliminary manuscript needs.**

## FINAL Phase H genotype matrix construction rule (H0-H1d, corrected)

| Classification | Treatment in final matrix |
|---|---|
| PRESENT (93 cases) | Use genotype data directly from QC VCF |
| ABSENT_AFTER_QC (4 cases) | Use partial callable data with explicit call-rate caveat; not treated as full-cohort 0/0 |
| NOT_IN_SOURCE (27 cases) | Encoded as `NOT_DETECTED_IN_SOURCE`, NOT as 0/0. Valid to report as "not detected as a variant in that cohort's source data." NOT valid to report as "allele frequency = 0.000" in that cohort, since exact-site reference confidence was not established. |

H1 (genotype extraction) is now finalized using this corrected rule.
No H1e planned; this level of resolution is sufficient and appropriately
scoped for the preliminary manuscript.
