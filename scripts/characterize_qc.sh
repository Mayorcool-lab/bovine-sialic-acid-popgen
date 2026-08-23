#!/bin/bash
set -euo pipefail

OUTDIR="data/processed/qc_characterization"
mkdir -p "$OUTDIR"

for cohort in south_africa europe; do
    VCF="$OUTDIR/${cohort}_all_genes_merged.vcf.gz"
    echo ">>> $cohort"

    vcftools --gzvcf "$VCF" --missing-indv --out "$OUTDIR/${cohort}_missing_indv" 2>/dev/null
    vcftools --gzvcf "$VCF" --missing-site --out "$OUTDIR/${cohort}_missing_site" 2>/dev/null
    vcftools --gzvcf "$VCF" --het --out "$OUTDIR/${cohort}_heterozygosity" 2>/dev/null
    vcftools --gzvcf "$VCF" --freq --out "$OUTDIR/${cohort}_allele_freq" 2>/dev/null
    bcftools stats "$VCF" | grep "^TSTV" > "$OUTDIR/${cohort}_tstv.txt"

    echo "  missing-indv: $(wc -l < ${OUTDIR}/${cohort}_missing_indv.imiss) lines"
    echo "  missing-site: $(wc -l < ${OUTDIR}/${cohort}_missing_site.lmiss) lines"
    echo "  heterozygosity: $(wc -l < ${OUTDIR}/${cohort}_heterozygosity.het) lines"
    echo "  ts/tv: $(cat ${OUTDIR}/${cohort}_tstv.txt)"
done

echo "Done. Outputs in $OUTDIR/"
