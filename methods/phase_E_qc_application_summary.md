# Phase E: QC Application Summary

Applied thresholds (justified by Phase D characterization):
- Genotype-level: DP >= 4, GQ >= 20 (masked to missing if failing, not
  removed as sites)
- Site-level: call rate >= 90% (<=10% missing) AFTER genotype masking

## Important finding: raw vs. true missingness

Phase D's raw site-missingness characterization (before genotype QC)
showed only 0.11-0.19% of sites exceeding 10% missing. After applying
DP/GQ genotype-level masking, TRUE site missingness is substantially
higher for several genes: e.g. CMAS lost 49% of sites (281->143),
ST8SIA4 lost 30% (1590->1109), CMAH lost 22% (1349->1051).

This is not a pipeline error. It reflects that low-quality individual
genotype calls (failing DP/GQ) were being counted as "present" in the
raw missingness check, masking the true unreliability of certain sites.
Once genotype quality is properly enforced, sites with genotype failures
concentrated among many animals (rather than evenly spread) correctly
fail the site-level call-rate threshold. Verified directly on CMAS:
individual genotype masking rate was 11.6% (3,907/33,720), but this was
unevenly distributed - some sites had 0% genotype failure, others had
22-49% - causing a nonlinear relationship between genotype-level masking
rate and site-level failure rate.

## Final QC-passed SNP counts

| Gene | South Africa | Europe |
|---|---|---|
| CMAH | 1051 | 1202 |
| GNE | 331 | 285 |
| NANS | 219 | 247 |
| NANP | 74 | 48 |
| CMAS | 143 | 150 |
| SLC35A1 | 206 | 188 |
| ST8SIA2 | 1360 | 883 |
| ST8SIA4 | 1109 | 811 |
| ST6GAL1 | 824 | 530 |
| ST3GAL1 | 1254 | 1135 |
| ST3GAL4 | 699 | 673 |
| NEU1 | 90 | 103 |
| NEU3 | 192 | 184 |
| **Total** | **7552** | **6439** |

## Independent verification and cohort-level pattern

Full site retention machine-verified directly against source files
(results/qc/phase_e_site_retention.tsv), exact match to all numbers above.

Retention rate differs systematically by cohort: South Africa retained
51-86% per gene, whereas Europe retained 88-97%, with no overlap in
per-gene retention ranges. Overall, South Africa retained 7,552/9,684
sites (77.98%) and Europe retained 6,439/6,842 sites (94.11%).

The reason for this cohort-level difference is not attributed to sample
size alone. Although the smaller South African cohort makes each failed
genotype contribute a larger discrete increment to site missingness,
the magnitude and consistency of the observed retention difference
indicate differences in the distribution of genotype-level DP/GQ
failures across sites. Because the two cohorts originate from distinct
source datasets/sequencing efforts, this pattern is treated as a
technical characteristic of the datasets rather than as a biological
difference between African and European cattle. Identical QC thresholds
were applied to both cohorts.

## Phase E: COMPLETE
QC-filtered, analysis-ready SNP datasets for both cohorts, all 13 genes.
Ready for Phase F/G: functional annotation (VEP).
