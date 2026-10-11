# H15: PCA and genome-wide FST verification

## Scope

This audit assesses the reproducibility and manuscript consistency of
collaborator-supplied genome-wide population-genomic results. It does
not recompute the original genome-wide PCA or FST estimates.

## Source data

- data/external/aiyesa_popgen/breed_all.weir.fst
- data/external/aiyesa_popgen/country_all.weir.fst
- data/external/aiyesa_popgen/region_all.weir.fst
- data/external/aiyesa_popgen/pop_gen_script.R
- data/reference/candidate_gene_coordinates_provisional.tsv

The original collaborator files remain unchanged.

## Genome-wide FST verification

All three FST tables contain 950,219 unique, identically ordered
autosomal genomic positions across chromosomes 1–29.

| Metric | Breed | Country | Region |
| --- | ---: | ---: | ---: |
| Input positions | 950219 | 950219 | 950219 |
| Finite FST estimates | 949858 | 949880 | 949880 |
| Nonfinite estimates | 361 | 339 | 339 |
| Negative estimates | 49143 | 86001 | 182777 |
| Raw arithmetic mean FST | 0.08591516 | 0.05914626 | 0.06208302 |
| Empirical 99.9th percentile | 0.45172629 | 0.47501850 | 0.57975208 |
| Positions meeting threshold | 950 | 950 | 950 |

The empirical thresholds were calculated after clipping negative
FST estimates to zero for descriptive plotting. Raw negative
estimates were preserved in the source data.

These thresholds are descriptive, not formal selection tests.
The arithmetic means are not genome-wide weighted FST estimates.

## Pathway-gene positions

Using the 13 gene intervals documented in
data/reference/candidate_gene_coordinates_provisional.tsv,
213 genomic positions overlap 12 of the 13 pathway genes.

NANP contains no overlapping position among these genome-wide
FST markers. This does not establish that NANP lacks variation.

All 213 pathway positions have finite FST estimates at breed,
country and regional levels. None exceeds its respective
genome-wide empirical 99.9th-percentile threshold.

The maximum pathway FST values are:

| Level | Maximum FST | Gene | Position |
| --- | ---: | --- | --- |
| Breed | 0.389977 | ST3GAL1 | chr14:7932660 |
| Country | 0.423442 | ST3GAL1 | chr14:7932660 |
| Region | 0.537383 | ST3GAL1 | chr14:7932660 |

The three grouping schemes are related and do not constitute
independent replication of this genomic signal.

The manuscript previously reported 220 pathway positions.
The independently reproducible count using the documented
gene intervals is 213. The provenance of the earlier 220
count could not be established.

## Population metadata

The independently validated pathway cohort contains 346 animals,
representing 16 breeds and four countries.

All 346 records were reconciled against the original metadata
using ENA sample accessions, with no breed or country conflicts.

Aiyesa's original R script documents sample matching and
population grouping. However, the final genome-wide PLINK/VCFtools
sample manifest was not supplied, so the exact genome-wide
analysis cohort cannot be independently reconstructed from the
handover files.

## PCA verification

The original R script reads pca_matched.eigenvec and
pca_matched.eigenval, but these files were not supplied.

Consequently, the manuscript's approximate PC1 and PC2 variance
percentages (25% and 11%) remain collaborator-reported and
cannot be independently verified from the available files.

The original genome-wide genotype dataset required to reproduce
PCA was also unavailable.

## Pairwise FST

The supplied breed, country and regional per-site FST files
do not contain pairwise population FST estimates.

The manuscript's pairwise FST value of approximately 0.2
requires confirmation from the original pairwise output.

## Plotting-script review

The original R script loads region_all.weir.fst into a variable
named breed. Its first Manhattan plot is titled breed-level
despite using regional data; its subsequent pathway-annotated
plot is correctly titled region-level.

This is a plotting-label inconsistency, not evidence that the
original FST estimates are incorrect.

The curated H12/H13 figures use the appropriately named FST files.

## Remaining limitations

- Original genome-wide genotype QC and LD pruning cannot be
  independently reproduced from the supplied pathway-only VCF.
- The original genome-wide sample manifest is unavailable.
- Original PCA eigenvectors and eigenvalues are unavailable.
- Original pairwise FST estimates are unavailable.
- The FST files do not independently identify their reference
  assembly. The project context supports ARS-UCD1.2, but
  coordinate compatibility is not conclusively established.
- The original gene_coord.txt used in Aiyesa's R script
  was not supplied.

## Manuscript decision

Use 220 pathway-associated positions, defined by the documented
13 gene bodies plus strand-aware 2-kb upstream extensions.
The 213-position count applies to gene bodies alone.
Retain the independently verified FST thresholds and pathway maxima.

Treat PCA variance percentages and pairwise FST estimates as
collaborator-reported pending supporting outputs.

Do not interpret high FST values alone as evidence of selection.

## Addendum: reconciliation of 213 versus 220 pathway SNP positions

A subsequent independent audit resolved the difference between the
213 gene-body SNP positions originally documented above and the
220 pathway-associated positions reported in the manuscript.

Using the 13 gene coordinates and strand annotations in
data/reference/candidate_gene_coordinates_provisional.tsv,
the following overlap counts were independently reproduced
against the 950,219 genome-wide FST positions:

- Gene bodies only: 213 positions.
- Gene bodies plus 2 kb upstream, strand-aware: 220 positions.
- Gene bodies plus 2 kb on both sides: 230 positions.

The seven additional upstream positions are:

| Gene | Chromosome | Position | Breed FST |
| --- | ---: | ---: | ---: |
| CMAS | 5 | 88170608 | 0.135474 |
| NEU1 | 23 | 27492663 | 0.0418101 |
| NEU1 | 23 | 27493079 | 0.0617763 |
| NEU1 | 23 | 27494455 | 0.101666 |
| NEU3 | 15 | 54251750 | 0.112011 |
| ST8SIA2 | 21 | 14597504 | 0.0284394 |
| ST8SIA4 | 7 | 100139496 | 0.114097 |

For positive-strand genes, the upstream interval extends
2,000 bp before the annotated gene start. For negative-strand
genes, it extends 2,000 bp beyond the annotated gene end.

Using the 220-position definition, all three FST scans retain
their previously verified pathway maxima:

- Breed: 0.389977.
- Country: 0.423442.
- Region: 0.537383.

All three maxima occur at ST3GAL1, chr14:7932660.

No pathway position meets the corresponding genome-wide
empirical 99.9th-percentile threshold in any scan.

Manuscript decision: retain 220 pathway-associated positions
because the manuscript defines the pathway intervals as
gene bodies plus 2 kb upstream.

The earlier 213-position result remains correct specifically
for gene bodies without the upstream extension.

This reconciliation does not resolve the separately documented
uncertainty concerning the original FST reference assembly.
