17e: Harmonisation / QC Strategy

Status: V1, provisional.

Steps
Normalize variants:
bcftools norm -f Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz -m -both -O z -o african_normalized.vcf.gz african_raw.vcf.gz
Restrict to biallelic SNPs:
bcftools view -v snps -m2 -M2 -O z -o african_biallelic.vcf.gz african_normalized.vcf.gz
Check REF/ALT against reference:
bcftools norm --check-ref e -f Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz african_biallelic.vcf.gz
Apply QC thresholds:
vcftools --gzvcf african_biallelic.vcf.gz --minDP 4 --max-missing 0.9 --recode --recode-INFO-all --out african_qc
Merge African and European after both go through steps 1-4:
bcftools merge -O z -o merged_african_european.vcf.gz african_qc.recode.vcf.gz european_qc.recode.vcf.gz
Verify sample count after merge, expect 409:
bcftools query -l merged_african_european.vcf.gz | wc -l
