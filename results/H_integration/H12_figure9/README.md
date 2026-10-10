# H12: Figure 9, breed-level per-site FST

## Purpose

Reproduce the breed-level FST Manhattan plot using Aiyesa's original
breed-level FST estimates, with annotations for the 13 bovine
sialic-acid pathway genes.

## Input files

- `data/external/aiyesa_popgen/breed_all.weir.fst`
- `data/reference/candidate_gene_coordinates_provisional.tsv`

The plotting script is `scripts/H12_figure9_breed_fst.py`.

Aiyesa's original R script and source datasets remain unchanged.

## Methods

The input contains 950,219 autosomal SNP positions across chromosomes
1–29. Nonfinite FST estimates were excluded.

Negative Weir and Cockerham FST estimates were retained in the raw
data but displayed as zero, following the original plotting convention.

The empirical 99.9th percentile of the displayed FST distribution
was used to identify high-FST positions. This is a descriptive
outlier threshold, not a formal statistical test for selection.

Cumulative chromosome positions were calculated using the maximum
observed SNP position on each chromosome, following Aiyesa's
original plotting approach.

Gene annotations use Ensembl release 110 coordinates on the
ARS-UCD1.2 bovine reference assembly.

## Results

- Input SNP positions: 950,219
- Finite FST estimates: 949,858
- Nonfinite estimates excluded: 361
- Negative FST estimates: 49,143
- Empirical 99.9th percentile: 0.451726
- Positions meeting the empirical threshold: 950
- FST positions within the 13 pathway genes: 213
- Pathway genes containing at least one FST position: 12 of 13
- Pathway-gene SNPs meeting the empirical threshold: 0
- NANP: no FST positions within its annotated gene boundaries

## Interpretation and limitations

No observed SNP within the 13 annotated pathway genes met the
empirical 99.9th-percentile threshold.

This does not establish an absence of population differentiation
or selection at these genes.

The FST files do not independently document their reference assembly.
Aiyesa's shared VCF specifies ARS-UCD1.2, but only 135 of its
4,504 unique variant positions overlap the breed-level FST table.
Therefore, coordinate compatibility is supported by the available
project context but is not independently established for the FST files.

The FST estimates and their upstream sample grouping and filtering
were supplied by the collaborator and were not recomputed here.

## Approved figure

- `figure9_breed_fst_v4.png`: 300 dpi
- `figure9_breed_fst_v4.pdf`: publication-quality export
- `figure9_breed_fst_v4.svg`: vector export
- `figure9_gene_fst_summary.tsv`: gene-level SNP counts

Figure version 4 uses separate horizontal gene labels and short
position ticks to avoid annotation overlap.
