# H7: Candidate-level regional Weir–Cockerham FST

## Objective
Estimate per-SNP differentiation between European and South African cattle cohorts for the 62 protein-altering candidate variants.

## Inputs
- H6.2 authoritative 62-candidate table
- Corrected gene-specific European VCFs (226 animals)
- Corrected gene-specific South African VCFs (120 animals)

## Methods
- Exact CHROM/POS/REF/ALT matching
- Genotype extraction using BCFtools 1.9
- Weir–Cockerham FST estimation using VCFtools 0.1.16
- Independent mathematical verification of every finite estimate
- Temporary normalization of missing genotype notation
- No reference-genotype imputation for absent VCF records

## Results
- 62 candidates assessed
- 31 candidates represented in both cohort VCFs
- 30 finite regional FST estimates
- 1 undefined estimate
- 21 candidates absent from the European variant-only VCF
- 10 candidates absent from the South African variant-only VCF

## Undefined CMAH estimate
Variant 23:32727639:G>T had 219 callable European and
117 callable South African genotypes.

All callable genotypes were homozygous reference.
Both ALT allele frequencies were zero.

The single-site FST denominator is therefore zero.
This result is undefined, not FST = 0.

## Independent validation
NEU3 15:54265634:T>C:
- Europe ALT AF: 0.344749
- South Africa ALT AF: 0.827586
- Regional Weir–Cockerham FST: 0.369494

All 30 finite estimates agreed with the independent
Weir–Cockerham implementation within 1e-5.

Eight small negative estimates were retained without truncation.

## Interpretation limitations
- Absent records in variant-only VCFs do not establish
  homozygous-reference genotypes.
- Estimates describe differentiation at individual candidate SNPs.
- The arithmetic mean of candidate-level FST values is not
  a genome-wide FST estimate.
- Differentiation alone does not demonstrate selection or causality.
- Country-level and breed-level analyses remain separate tasks.

## Outputs
- candidate_regional_fst_62.tsv
- candidate_regional_fst_run_summary.txt

## Reproduction
Activate the sialic_acid_popgen Conda environment, then run:

python3 scripts/H7_candidate_regional_fst.py
