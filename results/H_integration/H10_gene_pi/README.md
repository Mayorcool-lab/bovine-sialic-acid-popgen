# H10: Independent validation of breed-level nucleotide diversity

## Dataset and provenance

The analysis uses:
- data/external/aiyesa_popgen/mega_shared.vcf.gz
- data/external/aiyesa_popgen/pi_res_breed.txt
- results/H_integration/H9_candidate_breed_fst/breed_sample_manifest_346.tsv

The shared VCF contains 346 animals from 16 breeds and 4,505
biallelic SNP records.

## Estimator

The independent analysis is implemented in:
scripts/H10_gene_pi_validation.py

For each breed and SNP, nucleotide diversity is calculated as:

    pi = n / (n - 1) * (1 - p_ref^2 - p_alt^2)

Here, n is the number of called alleles, and p_ref and p_alt
are their observed frequencies.

Sites with fewer than two called alleles are excluded.
Sites monomorphic within a breed are retained and contribute zero.

The breed estimate is the arithmetic mean of per-site diversity
across eligible positions in the shared SNP panel.

## Independent validation

All 16 breeds were recalculated.

The absolute differences from Aiyesa's reported gene-level
mean pi values were below 0.0006.

This establishes close numerical agreement, not exact
reproduction of Aiyesa's undocumented estimator.

Aiyesa's reported site counts differ from the eligible-site
counts obtained independently. A 90% within-breed call-rate
filter does not explain those differences.

The original denominator selection remains unresolved.
Subsequent figures and text should use the independently
calculated values.

## Cohort discrepancy

The verified 346-animal manifest contains 23 Dutch Belted
animals. Aiyesa's gene-level result reports 22.

The independent calculation uses all 23 manifest animals.
The reason for the difference in Aiyesa's run is unknown.

## Small sample sizes

Mertolenga has 2 animals and Mirandesa has 4 animals.

Both remain in the complete H10 results. They are excluded
from the planned primary Figure 8 comparison because their
small sample sizes limit reliable breed-level interpretation.

## Interpretation and limitations

These estimates measure mean nucleotide diversity per
ascertained SNP position in the shared 4,505-SNP panel.

They are not nucleotide diversity per callable base across
complete gene regions or the genome.

The variant-only VCF does not provide the callable reference
positions needed for those calculations.

SNP ascertainment, missing genotypes and unequal breed
sample sizes may affect comparisons.

Country and regional summaries are descriptive comparisons
of breed-level estimates, not statistical tests of
population differences.

## Output

results/H_integration/H10_gene_pi/gene_pi_validation.tsv

## Downstream use

Task 3 will rebuild Figure 8 using the independently
calculated gene-level values. Genome-wide pi results
from the unverified 408-animal analysis will not be used.
