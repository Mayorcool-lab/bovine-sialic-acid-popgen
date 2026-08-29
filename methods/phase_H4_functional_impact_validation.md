# Phase H4 — Functional-impact annotation and validation

## Objective

Phase H4 adds functional-impact evidence to the 62 protein-altering/stop
variants defined in the authoritative Phase G candidate table:

`results/vep/phase_g_protein_altering_FINAL.tsv`

Phase G remains authoritative for variant identity, gene assignment,
consequence, protein change, and manual CMAH rescue. H4 adds native
Bos taurus SIFT evidence where technically available.

No structural modelling, molecular dynamics, docking, or other new
structure-based analysis is performed in H4.

## H4.1 Consequence classes

The 62 candidates comprise:

- 58 `missense_variant`
- 2 `missense_variant&splice_region_variant`
- 1 `stop_gained`
- 1 `stop_lost`

Therefore:

- 60 variants contain a missense consequence and are eligible for SIFT.
- 2 HIGH-impact stop variants are not evaluated by SIFT.

The splice-region labels on the two missense+splice-region variants are
retained as VEP consequence annotations and are not interpreted as proof
of altered splicing.

## H4.2 Native bovine SIFT analysis

SIFT was generated with Ensembl VEP using the locally installed bovine
cache and ARS-UCD1.2 reference.

Environment/cache information established during H4:

- Ensembl VEP: 110.1
- Ensembl API: 110
- Species: Bos taurus
- Assembly: ARS-UCD1.2
- Cache: Ensembl 110 ARS-UCD1.2
- Genebuild: 2018-12
- dbSNP: 150
- Native bovine SIFT data available in the cache

Input:

`results/H4/missense_candidates.vcf`

VEP/SIFT output:

`results/H4/missense_sift_annotated.vcf`

The VEP run used native bovine SIFT prediction and score (`--sift b`).

## Exact allele reconciliation

All 60 missense-containing Phase G candidates were represented exactly
once in the H4 VEP/SIFT input/output by:

`CHROM + POS + REF + ALT`

Validation result:

- Phase G missense-containing candidates: 60
- H4 VEP/SIFT records: 60
- Exact allele matches: 60
- Missing alleles: 0
- Unexpected alleles: 0

## Transcript-selection rule

Initial exploratory parsing attempted to use VEP `CANONICAL=YES`.

This assumption was rejected because at least one variant
(NEU1 23:27494657 A>G) had more than one transcript marked
`CANONICAL=YES`.

The authoritative Phase G transcript (`Feature`) is therefore used for
H4 SIFT extraction.

For every Phase G missense candidate with a non-empty `Feature`, the H4
parser requires exactly one VEP annotation matching that transcript ID.

This preserves transcript consistency between Phase G consequence
annotation and H4 functional-impact prediction.

## CMAH manual-rescue exception

Four manually rescued CMAH missense variants have no Phase G `Feature`
because the Ensembl 110 ARS-UCD1.2 gene model does not represent the
relevant RefSeq-supported CMAH coding extension correctly.

The four variants are:

- CMAH 23:32739285 A>G — N386S
- CMAH 23:32739292 C>G — I388M
- CMAH 23:32743918 C>T — P424L
- CMAH 23:32745866 T>A — F512Y

In the fresh H4 VEP run, all four exact alleles are annotated as
`intergenic_variant` with blank transcript Feature and blank SIFT.

Therefore these variants are explicitly classified as:

`SIFT_UNAVAILABLE`

They are not assigned substitute transcripts, human-orthologue scores,
or fabricated SIFT values.

Phase G manual RefSeq rescue remains authoritative for their coding
identity and protein changes.

## SIFT results

Among the 60 missense-containing variants:

- SIFT available: 56
- SIFT unavailable: 4

Prediction distribution among the 56 scoreable variants:

- deleterious: 6
- deleterious_low_confidence: 6
- tolerated: 38
- tolerated_low_confidence: 6

Low-confidence categories are preserved explicitly and are not collapsed
into the corresponding high-confidence categories.

The two non-missense HIGH-impact candidates are:

- ST3GAL4 29:29738209 G>T — stop_lost — HIGH
- CMAS 5:88162466 C>A — stop_gained — HIGH

For these variants SIFT is recorded as:

`NOT_APPLICABLE_HIGH_IMPACT`

This represents non-applicability of the missense predictor, not missing
coverage.

## CMAH IMPACT metadata correction discovered during H4

During validation of the final H4 table, the four manually rescued CMAH
missense variants were found to have:

`Consequence=missense_variant`
but
`IMPACT=MODIFIER`

Inspection of `scripts/G7_split_target_and_cmah_rescue.py` showed that
the manual RefSeq rescue replaced the coding consequence but did not
replace the original VEP `IMPACT` value inherited from the intergenic
annotation.

The same issue affected the manually rescued synonymous CMAH variant
I515I.

The correction was made at the upstream G7 source:

- manually rescued `missense_variant` -> `IMPACT=MODERATE`
- manually rescued `synonymous_variant` -> `IMPACT=LOW`

G7 was regenerated and validated:

- N386S — missense_variant — MODERATE
- I388M — missense_variant — MODERATE
- P424L — missense_variant — MODERATE
- F512Y — missense_variant — MODERATE
- I515I — synonymous_variant — LOW

G8 was then regenerated.

The authoritative Phase G protein-altering table remained at 62
candidates with the same gene distribution and the same four manually
rescued CMAH missense variants.

H4 was subsequently regenerated from the corrected Phase G table.

## Final H4 evidence table

Authoritative H4 output:

`results/H4/H4_complete_62_variant_evidence.tsv`

Final SIFT-status accounting:

- `SIFT_AVAILABLE`: 56
- `SIFT_UNAVAILABLE`: 4
- `NOT_APPLICABLE_HIGH_IMPACT`: 2

Total: 62 variants.

Final integrity audit:

- Phase G rows: 62
- H4 rows: 62
- Exact Phase G/H4 allele-set match: PASS
- SYMBOL agreement: PASS
- Consequence agreement: PASS
- IMPACT agreement: PASS
- SIFT-status logic: PASS
- Integrity errors: 0

## Interpretation boundary

SIFT predictions are computational functional-impact predictions and are
not experimental evidence of altered protein function.

`deleterious` must not be equated with disease-causing or pathogenic in
the context of natural bovine population variation.

HIGH-impact stop consequences indicate stronger predicted molecular
consequence but do not by themselves establish phenotypic effect.

Phase H4 provides one evidence layer for later integration with
conservation, known functional/domain context, population frequencies,
breed-level architecture, and genome-wide evolutionary analyses.
