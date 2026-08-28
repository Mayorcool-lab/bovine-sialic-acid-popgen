# Phase H3: Breed-level allele frequency validation

## Objective

H3 calculates breed-level allele frequencies for all 62 protein-altering
candidate variants across the verified South African and European cattle
cohorts, while preserving the H1/H2 three-way source-status logic and
independently reconciling breed-level counts back to the authoritative H2
cohort-level results.

## Definitive metadata scope

South Africa was restricted to animals in the combined metadata table with
`country == "South Africa"`:

- Afrikaner: 40
- Bonsmara: 20
- Nguni: 40
- Tuli: 20

Total South Africa: 120 animals.

Europe used only the previously verified 226-animal mapping:

- Western_Finncattle: 25
- Northern_Finncattle: 25
- Eastern_Finncattle: 25
- Holstein_Friesian: 15
- Deep_Red: 24
- Dutch_Belted: 23
- Dutch_Friesian: 24
- Meuse_Rhine_Yssel: 23
- Groningen_White_Headed: 21
- Barrosa: 15
- Mirandesa: 4
- Mertolenga: 2

Total Europe: 226 animals.

Grand total: 346 animals.

This avoids leakage of the 63 unresolved additional Portuguese samples.

## Candidate and source-status validation

The authoritative 62-candidate set was loaded from:

`results/vep/phase_g_protein_altering_FINAL.tsv`

For every candidate in both cohorts, H3 required an existing H1/H2
classification:

1. `PRESENT`
2. `ABSENT_AFTER_QC`
3. `NOT_IN_SOURCE`

All 62 x 2 = 124 candidate-by-cohort combinations were present and valid.

## Genotype extraction rules

### PRESENT

Genotypes were extracted from the final QC-filtered VCF for the relevant
gene and cohort.

### ABSENT_AFTER_QC

Genotypes were extracted from the normalized pre-site-QC VCF and filtered
at genotype level using the established thresholds:

- DP >= 4
- GQ >= 20
- missing genotypes excluded

This preserves partial callable evidence without pretending the site passed
the full >=90% site call-rate threshold.

### NOT_IN_SOURCE

No genotype frequency was inferred.

These rows are encoded as:

`NOT_DETECTED_IN_SOURCE`

with blank callable count, allele count, and AF fields.

They are not encoded as AF=0.0000 because exact-site reference confidence
was not established from variant-only source VCFs.

## REF/ALT identity safeguard

Before genotype extraction, H3 explicitly verifies the candidate REF/ALT
identity at the queried genomic position.

A hard RuntimeError is raised if the candidate allele does not match the VCF
record.

This prevents silent use of a different variant occurring at the same
position.

## Breed-level calculation

For each callable genotype:

- `0/0` contributes two REF alleles
- `0/1` or `1/0` contributes one REF and one ALT allele
- `1/1` contributes two ALT alleles
- missing genotypes are excluded

Unexpected allele values outside the expected 0/1 candidate representation
trigger a hard RuntimeError.

The output contains one row per candidate x breed:

62 variants x 16 breeds = 992 rows.

Output:

`results/phase_H_62_variant_breed_frequencies.tsv`

## H3 -> H2 reconciliation

Before the H3 table is written, breed-level counts are reconstructed to the
cohort level and compared against the authoritative H2 table.

For every candidate and cohort:

- breed callable N values must sum to H2 `N_genotyped`
- REF + ALT allele counts must equal 2 x callable N
- reconstructed breed-level AF must equal H2 AF to four decimal places
- source-status labels must agree with H2
- `NOT_IN_SOURCE` must remain frequency-undefined

The four H2 `ABSENT_AFTER_QC` AF values retain their
`(SUB-THRESHOLD CALL RATE)` annotation; H3 compares only the numeric AF
component while preserving the annotation in H2.

Final reconciliation:

124 / 124 candidate-by-cohort combinations PASS.

## H2 correction discovered through H3 audit

The first H3 reconciliation attempt detected a discrepancy for
ST3GAL1 14:7943747 in South Africa:

- H3 breed-summed callable N = 114
- H2 callable N = 115

Direct VCF inspection showed one genotype encoded as single-field missing
`.` in addition to five `./.` genotypes.

This revealed that H2 had excluded `./.` and `.|.` but not `.` when defining
callable N.

A global audit identified 8 affected candidate/cohort combinations.

H2 was corrected to treat `.`, `./.`, and `.|.` as missing. After rerunning
H2, only `N_genotyped` changed in those 8 rows; AF, heterozygote counts,
homozygous-ALT counts, source classifications, and the four H1c
ABSENT_AFTER_QC spot checks were unchanged.

After this controlled correction, H3 reconciled 124/124 combinations.

## Structural audit of final H3 table

Final output:

- total rows: 992
- unique variants: 62
- South Africa rows: 248
- Europe rows: 744
- integrity errors: 0

Status counts:

- PRESENT: 700
- ABSENT_AFTER_QC: 16
- NOT_DETECTED_IN_SOURCE: 276

All frequency-bearing rows satisfied:

- REF_COUNT + ALT_COUNT = 2 x N_CALLABLE
- AF within [0,1]

No `NOT_DETECTED_IN_SOURCE` row contained frequency or genotype-count data.

## Small-sample safeguard

Breed-level rows were flagged by callable sample size:

- `VERY_SMALL_N`: N < 5
- `SMALL_N`: N < 10
- `OK`: N >= 10
- `NA`: no frequency because variant was NOT_IN_SOURCE

Observed:

- OK: 634 rows
- VERY_SMALL_N: 82 rows
- NA: 276 rows
- SMALL_N: 0 rows

All 82 VERY_SMALL_N rows were confined to the two known tiny European breed
groups:

- Mirandesa: 41 rows, callable N=3-4
- Mertolenga: 41 rows, callable N=1-2

No larger breed unexpectedly fell below N=10.

These tiny-breed frequency estimates are retained for completeness but
should not be treated as strong evidence of breed-specific enrichment.

## Frequency plausibility audit

Among adequately sampled (`N_FLAG=OK`) breed/variant observations:

- usable rows: 634
- polymorphic rows: 227
- reference-fixed rows: 407
- ALT-fixed rows: 0
- impossible AF values: 0

Several high-frequency ALT observations were observed in adequately sampled
breeds, including variants in NEU3, ST3GAL4, SLC35A1, NEU1, and CMAH.
These are biologically interesting patterns for later prioritisation but are
not, by themselves, sufficient for candidate ranking.

Functional consequence, conservation, broader population contrast, and
genome-wide differentiation/selection evidence must be integrated before
final prioritisation.

## H3: CLOSED

H3 passes:

- definitive 346-animal metadata validation
- 62-candidate / 124 candidate-cohort classification validation
- explicit REF/ALT identity checking
- correct PRESENT / ABSENT_AFTER_QC / NOT_IN_SOURCE handling
- 992-row production generation
- 124/124 H3 -> H2 reconciliation
- zero structural integrity errors
- explicit small-N flagging
- frequency plausibility checks

Authoritative output:

`results/phase_H_62_variant_breed_frequencies.tsv`
