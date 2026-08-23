# sialic_acid_popgen — Full Analysis Workflow

## Scientific question

Determine how naturally occurring SNP variation across the bovine sialic-acid
metabolic network differs among European and South African cattle populations,
which variants are predicted to affect biological function, and whether pathway
loci show unusual population differentiation or signatures of selection relative
to the genome-wide background.

## Discovery cohorts

**South African cohort**: 120 cattle (Afrikaner 40, Bonsmara 20, Nguni 40, Tuli 20)
**European cohort**: 289 cattle, 12 breed groups
**Total**: 409 cattle

## Candidate network (13 provisional genes)

GNE, NANS, NANP, CMAS, SLC35A1, CMAH, ST8SIA2, ST8SIA4, ST6GAL1, ST3GAL1,
ST3GAL4, NEU1, NEU3

## Core rule: SNPs vs indels

Primary scientific analysis = SNPs, especially high-quality biallelic SNPs.
Indels are retained in the normalized archive for completeness but are not part
of the main population-genomic analysis unless an obviously high-impact
functional indel emerges. Normalize everything first, preserve a master
normalized file (SNPs + indels), then derive the SNP dataset from it.

```
RAW VCF -> NORMALIZED MASTER -> SNP analysis dataset (MAIN STUDY)
                              -> indels retained (secondary/archive only)
```

---

## PHASE A — Source-data verification

**A1. South African source: DONE.** Confirmed ARS-UCD1.2, unphased, 240 animals
in source VCF, 120 South African target animals identified via `optibov_id`/
`vcf_id`, regional random access confirmed working, FreeBayes source reports
v1.0.2 (header), papers' methods text states v1.3.1.

**A2. European source: DONE, not open.** Confirmed earlier in this project:
identical header (`##source=freeBayes v1.0.2`, `##fileDate=20240226`) to the
African file. Same discrepancy present identically in both datasets, consistent
with same underlying pipeline, not evidence of divergence.

**A3. Correct project documentation.** Update `dataset_manifest.tsv`,
`methods/17e_harmonisation_strategy.md`, and the checklist to record the actual
header-confirmed FreeBayes v1.0.2 alongside the papers' stated v1.3.1, noting
the discrepancy explicitly rather than picking one as authoritative.

---

## PHASE B — Targeted 13-gene extraction

## PHASE B — Targeted 13-gene extraction

B1. South African extraction: DONE. 120 animals, 13 genes, 10,765 raw records.
B2. European extraction: DONE. 226 verified animals, 13 genes, 7,628 raw
records. (63 Portuguese animals pending, see notes/STATUS.md)

PHASE B COMPLETE: 346 total animals, 18,393 total raw regional records
across both cohorts, before normalization/QC.
---

## PHASE C — Normalization

Both cohorts undergo identical processing, per gene.

**C1.** `bcftools norm` against ARS-UCD1.2: left-aligns indels, standardizes
representation, decomposes multiallelic records, checks alleles against
reference.

**C2.** Preserve normalized master file (SNPs + indels) as the reproducible
archive.

**C3.** Derive the primary SNP dataset (biallelic SNPs) from the normalized
master. This becomes the main analytical dataset.

---

## PHASE D — QC characterization before filtering

Before deciding thresholds, inspect actual distributions for both populations:
per-sample missingness, per-site missingness, DP distribution, GQ distribution,
allele frequency distribution, heterozygosity, transition/transversion ratio,
SNP count per gene, sample-level genotype completeness. Prevents blindly
copying arbitrary QC thresholds.

---

## PHASE E — Apply common QC

One harmonized QC framework for both datasets, informed by Phase D's actual
distributions, not assumed in advance. Genotype QC (min depth, min GQ) and
site QC (call rate/missingness, SNP quality, biallelic status, possibly MAF
depending on analysis). MAF filtering is analysis-specific, not universal,
see Phase F.

---

## PHASE F — Analysis-specific SNP datasets

**F1. Functional SNP dataset.** Question: which SNPs may alter biological
function? Retain sufficiently reliable SNPs including rare variants: coding,
missense, synonymous, splice, UTR, proximal regulatory, biologically relevant
intronic. No strong MAF filter, a rare missense SNP could be the key result.

---

## PHASE G — Functional annotation

Run VEP against ARS-UCD1.2. For each SNP: gene, transcript, consequence, codon
change, amino-acid change, position, impact category, existing variant ID.
Apply the validated cattle-applicable tools from the V1 specification for
missense SNPs. Goal: thousands of SNPs down to a manageable candidate set.

---

## PHASE H — Population-level allele analysis

For every pathway SNP: allele frequency per breed, per South Africa overall,
per Europe overall, private/near-private alleles, rare vs common, breed-
enriched alleles, continental frequency differences. This is where the
biological story starts appearing, but does not yet prove selection.

---

## PHASE I — Genome-wide population-genomics dataset

**Note on scale: this is substantially larger than any step completed so far.**
Separate from the 13-gene functional dataset, PCA/ADMIXTURE/genome-background
FST/iHS/XP-EHH cannot use candidate genes alone. Requires processing the full
African (~75GB) and European (~60GB) genomes, chromosome-by-chromosome rather
than downloading everything at once. This is Aiyesa's territory per his
established role (population genetics execution), not a solo extension of the
13-gene work.

```
African/European whole genome -> normalize/QC chromosome by chromosome
                                -> harmonize common representation
                                -> 409-animal population-genomics dataset
```

---

## PHASE J — PCA

Genome-wide, quality-controlled, LD-pruned SNPs. Question: what is the actual
genetic structure of these 409 cattle? Do breeds cluster as expected, do
African/European populations separate strongly, is there admixture, are there
outliers, do breed labels agree with genomic ancestry? Comes before final
population comparisons.

## PHASE K — ADMIXTURE

LD-pruned genome-wide SNP set. Evaluate ancestry structure across plausible K
values. PCA + ADMIXTURE inform 17f.

## PHASE L — Freeze population comparisons (17f)

Only after PCA/ADMIXTURE. Likely hierarchy: within South Africa (Nguni vs
Afrikaner, Nguni vs Bonsmara, Tuli vs Afrikaner, etc.), within Europe
(biologically justified breed comparisons), Africa vs Europe (broader,
interpreted cautiously). Aiyesa and Idowu advise; final study-design decision
remains yours.

## PHASE M — Genome-wide FST

Sliding windows, V1 default 50kb window/20kb step (cited, cattle literature
precedent). Determine where the 13 pathway genes fall in the genome-wide
empirical distribution (e.g. "CMAH FST percentile = 98.7%" is meaningful;
"CMAH FST = 0.31" alone is not).

## PHASE N — Nucleotide diversity (π)

Genome-wide diversity. Is diversity unusually reduced/elevated around
candidate genes? Interpret alongside population history.

## PHASE O — Tajima's D

Genome-wide sliding-window. Strongly negative = excess rare variants,
potential recent sweep/expansion. Strongly positive = excess intermediate-
frequency variants, balancing selection/structure. Interpret against
genome-wide distribution.

## PHASE P — Phasing

Both datasets currently unphased (confirmed via VCF headers). Phase genome-
wide/chromosome-scale data using Beagle or agreed method before any haplotype
analysis.

## PHASE Q — iHS

Phased chromosome-scale data. Is an allele embedded in an unusually long
haplotype within a population, suggesting recent positive selection? Compute
genome/chromosome-wide first, then examine candidate loci.

## PHASE R — XP-EHH

Phased chromosome-scale data, between populations. Has a haplotype approached
fixation in one population but not another? Compute genome/chromosome-wide
first, then interrogate candidate loci.

## PHASE S — Integrate selection statistics

Do not declare selection from one extreme statistic alone. Strongest evidence:
convergence of high FST + unusual Tajima's D + reduced π + extreme iHS +
extreme XP-EHH at the same locus. More independent signals converging =
stronger inference.
