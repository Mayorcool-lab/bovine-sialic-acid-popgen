# H13: Figure 9, genetic differentiation at three population-grouping levels

## Purpose

Generate the main manuscript Figure 9, comparing per-site genetic
differentiation across bovine autosomes using three alternative
population-grouping levels:

- Panel A: breed
- Panel B: country
- Panel C: region

All three panels use the same 950,219 genomic SNP positions.

## Data provenance

Original FST estimates were supplied by collaborator Leke V. Aiyesa:

- data/external/aiyesa_popgen/breed_all.weir.fst
- data/external/aiyesa_popgen/country_all.weir.fst
- data/external/aiyesa_popgen/region_all.weir.fst

The original pop_gen_script.R constructs sample groupings using
breed, country and region metadata.

The precise upstream FST calculation commands and reference
assembly of the supplied FST tables have not been independently
verified.

Gene annotations use Ensembl release 110 coordinates on the
Bos taurus ARS-UCD1.2 reference assembly.

## Methods

Each dataset contains 950,219 unique genomic positions across
autosomes 1-29. Genomic positions match exactly across datasets.

Nonfinite FST estimates were excluded. Negative Weir and
Cockerham FST estimates were preserved in the input data but
displayed as zero.

For each grouping level, the empirical 99.9th percentile of
the displayed FST values was calculated separately.

The threshold is descriptive and is not a formal statistical
test for selection.

Cumulative chromosome coordinates follow the maximum observed
SNP position per chromosome, consistent with the H12 workflow.

The 13 pathway genes are annotated using provisional gene
coordinates. The same chromosome positions and gene annotations
are used across all three panels.

## Validation results

| Metric | Breed | Country | Region |
| --- | ---: | ---: | ---: |
| Input SNP positions | 950219 | 950219 | 950219 |
| Finite FST estimates | 949858 | 949880 | 949880 |
| Nonfinite estimates | 361 | 339 | 339 |
| Negative estimates | 49143 | 86001 | 182777 |
| Empirical 99.9th percentile | 0.451726 | 0.475019 | 0.579752 |
| Empirical outlier positions | 950 | 950 | 950 |
| Pathway-gene SNP positions | 213 | 213 | 213 |
| Pathway-gene empirical outliers | 0 | 0 | 0 |
| Genes with at least one FST position | 12 | 12 | 12 |

NANP has no observed SNP positions within its annotated gene
boundaries in these FST tables.

## Interpretation and limitations

None of the observed pathway-gene SNP positions reaches the
group-specific empirical 99.9th-percentile threshold.

This does not establish an absence of population differentiation
or selection at the pathway genes.

Breed, country and region analyses represent alternative groupings
of the underlying data, not independent datasets.

The source FST files do not independently document their reference
assembly. Coordinate compatibility therefore remains provisional.

The FST estimates were supplied by the collaborator and were not
recomputed in H13.

## Reproducibility

From the repository root, in the sialic_plot environment:

    python scripts/H13_supplementary_fst.py
    python scripts/H13_figure9_three_level_fst.py

## Main manuscript outputs

- figure9_breed_country_region_fst.png (300 dpi)
- figure9_breed_country_region_fst.pdf
- figure9_breed_country_region_fst.svg

## Supporting outputs

- H13_fst_validation_summary.tsv
- country_gene_fst_summary.tsv
- region_gene_fst_summary.tsv

The original H12 breed-level figure remains preserved under
results/H_integration/H12_figure9/ and Git commit 22d7e31.

The final manuscript Figure 9 consists of three vertically
stacked panels with aligned chromosomes, shared gene annotations,
and separate empirical thresholds.
