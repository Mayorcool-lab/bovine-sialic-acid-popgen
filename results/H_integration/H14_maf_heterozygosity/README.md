# H14: Pathway SNP Minor Allele Frequency and Heterozygosity

## 1. Purpose

This analysis characterizes minor allele frequency (MAF), observed
heterozygosity (Ho), and expected heterozygosity (He) across 16 cattle
breeds represented in the shared bovine sialic-acid pathway SNP dataset.

The analysis is descriptive and does not constitute a genome-wide
diversity analysis or a formal test of selection, inbreeding, or
Hardy-Weinberg equilibrium.

## 2. Input data

Genotype source:
`data/external/aiyesa_popgen/mega_shared.vcf.gz`

Authoritative sample-to-breed mapping:
`results/H_integration/H9_candidate_breed_fst/breed_sample_manifest_346.tsv`

Dataset characteristics:

- 346 diploid animals
- 16 cattle breeds
- 4,505 biallelic SNP records
- 4,504 unique genomic positions
- Two separate biallelic records at chr23:32681432 (T>A and T>G)
- 1,558,730 possible genotype calls

The canonical SNP identity is chromosome, position, REF, and ALT.
Records at the same position with different ALT alleles are retained
separately.

The input VCF is a pathway-specific, variant-ascertained dataset.
Its SNP records must not be interpreted as a complete inventory of
all genomic positions or as genome-wide markers.

## 3. Scripts

Run from the repository root, using the project Python environment:

1. `python scripts/H14_maf_analysis.py`
2. `python scripts/H14_heterozygosity.py`
3. `python scripts/H14_supplementary_figures.py`

The scripts calculate allele frequencies, genotype heterozygosity,
breed-level summaries, and supplementary figures.

The original external genotype dataset is not modified.

## 4. Genotype quality and completeness

The genotype dataset was examined for genotype completeness and
available DP/GQ values.

Observed characteristics:

- Complete genotype calls: 1,484,377
- Missing or partially missing calls: 74,353
- Overall genotype call rate: 95.2299%
- Complete calls satisfied DP >= 4 and GQ >= 20
- All SNP records satisfied a pooled site call rate >= 90%

These are verification results for the supplied VCF, not evidence
that the H14 scripts independently performed the original upstream
genotype-masking and site-filtering procedures.

No MAF threshold was applied to remove SNP records from H14.

## 5. MAF calculation

For each callable biallelic SNP record:

- Allele number (AN) = 2 * number of called diploid genotypes
- Alternative allele frequency = ALT allele count / AN
- MAF = min(ALT allele frequency, 1 - ALT allele frequency)

Missing or partially missing genotypes do not contribute to
allele counts.

Pooled MAF was calculated across all 346 animals.
Breed-specific MAF was calculated independently for each breed.

If AN = 0, allele frequency and MAF are undefined.

The pooled dataset contains:

- 4,298 polymorphic SNP records
- 207 monomorphic SNP records in the analyzed cohort
- 918 records with MAF < 0.01
- 1,125 records with 0.01 <= MAF < 0.05
- 2,462 records with MAF >= 0.05

The pooled mean MAF is approximately 0.134577.

## 6. Heterozygosity calculation

Observed heterozygosity (Ho) is the proportion of complete
genotype calls that are heterozygous.

Breed-level and pooled Ho are calculated as:

Ho = total heterozygous genotype calls / total called genotypes

Expected heterozygosity at a callable biallelic SNP record is:

He = 2p(1-p)

where p is the alternative allele frequency within the relevant
breed or pooled population.

Summary He is the unweighted mean of per-record He across
callable SNP records.

Ho and He therefore use different aggregation schemes:
Ho is weighted by the number of called genotypes, whereas
summary He is averaged across callable SNP records.

Pooled results:

- Heterozygous genotype calls: 247,776
- Ho = 0.1669225540
- Mean He = 0.1876212249

## 7. Output files

Allele frequency tables:

- `H14_pooled_maf.tsv`
- `H14_breed_maf.tsv`

Heterozygosity tables:

- `H14_site_heterozygosity.tsv`
- `H14_heterozygosity_summary.tsv`

Supplementary summary:

- `H14_supplementary_breed_summary.tsv`

Supplementary figures, each supplied in PNG, PDF, and SVG:

- `H14_figure_S1_breed_maf_distribution`
- `H14_figure_S2_breed_heterozygosity`

Figure S1 summarizes breed-level monomorphic, low-frequency,
common, and uncallable SNP-record proportions.

Figure S2 compares observed and expected heterozygosity by breed.

The S1 and S2 labels are provisional until supplementary figure
numbering is reconciled with the complete manuscript.

## 8. Independent consistency checks

The following checks passed:

- All 346 sample assignments reconstruct the pooled cohort.
- Breed-level called genotypes sum to 1,484,377.
- Breed-level missing genotypes sum to 74,353.
- Breed-level heterozygous genotypes sum to 247,776.
- Pooled Ho is reproduced directly from genotype counts.
- Every one of the 4,505 pooled SNP records has 16 breed entries.
- Breed-level called and missing genotype counts reconstruct
  pooled counts at every SNP record.
- Breed-level REF and ALT allele counts reconstruct pooled
  allele counts at every SNP record.
- AN = 2 * called diploid genotype count.
- REF allele count + ALT allele count = AN.
- Alternative allele frequencies and MAF values agree with
  independently recalculated values.

These checks validate internal numerical consistency. They do
not independently establish the accuracy of the original
sequencing, variant calling, or upstream genotype quality control.

## 9. Interpretation limitations

1. Results describe ascertained SNP records within the selected
   sialic-acid pathway regions, not genome-wide genetic diversity.
2. SNP-based Ho and He are not interchangeable with per-nucleotide
   diversity estimates.
3. Mertolenga (n=2) and Mirandesa (n=4) have very small samples;
   their estimates require particular caution.
4. Sample-size differences affect the probability of detecting
   low-frequency and polymorphic variants within breeds.
5. Pooled He may reflect differentiation among breeds and should
   not be interpreted as within-breed expected heterozygosity.
6. Differences between Ho and He alone do not establish
   inbreeding, selection, or Hardy-Weinberg disequilibrium.
7. The source VCF and the original collaborator-generated
   population-genomic datasets are preserved unchanged.

## 10. Manuscript integration

These H14 results provide pathway-specific supplementary
descriptions of allele-frequency distributions and heterozygosity.

They should be reconciled with the manuscript's genome-wide
diversity analysis and supplementary figure numbering before
submission.

The analysis does not replace the genome-wide PCA, diversity,
differentiation, or allele-frequency-spectrum analyses contributed
by Leke V. Aiyesa.
