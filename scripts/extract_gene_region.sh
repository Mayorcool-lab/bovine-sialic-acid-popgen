#!/bin/bash
set -euo pipefail
# Usage: ./extract_gene_region.sh <gene_symbol> <vcf_url> <local_csi> <sample_list> <output_dir>
# Reads coordinates from data/reference/candidate_gene_coordinates_provisional.tsv
# Adds strand-aware 2kb upstream buffer, extracts region, subsets to given samples

GENE=$1
VCF_URL=$2
LOCAL_CSI=$3
SAMPLE_LIST=$4
OUTDIR=$5
UPSTREAM_BP=2000

COORDS=$(grep -P "^${GENE}\t" data/reference/candidate_gene_coordinates_provisional.tsv)

if [ -z "$COORDS" ]; then
    echo "ERROR: Gene $GENE not found in coordinates file"
    exit 1
fi

CHR=$(echo "$COORDS" | cut -f3)
START=$(echo "$COORDS" | cut -f4)
END=$(echo "$COORDS" | cut -f5)
STRAND=$(echo "$COORDS" | cut -f6)

if [ "$STRAND" == "+" ]; then
    EXT_START=$((START - UPSTREAM_BP))
    EXT_END=$END
else
    EXT_START=$START
    EXT_END=$((END + UPSTREAM_BP))
fi
if [ "$EXT_START" -lt 1 ]; then
    EXT_START=1
fi

echo "=== $GENE ==="
echo "Strand: $STRAND | Window: ${CHR}:${EXT_START}-${EXT_END}"

mkdir -p "$OUTDIR"

# Step 1: extract region from remote VCF using local index
bcftools view -r "${CHR}:${EXT_START}-${EXT_END}" --regions-overlap 0 \
  "${VCF_URL}##idx##${LOCAL_CSI}" -O z -o "${OUTDIR}/${GENE}_full.vcf.gz"
bcftools index "${OUTDIR}/${GENE}_full.vcf.gz"

# Step 2: subset to target samples
bcftools view -S "$SAMPLE_LIST" "${OUTDIR}/${GENE}_full.vcf.gz" -O z -o "${OUTDIR}/${GENE}_subset.vcf.gz"
bcftools index "${OUTDIR}/${GENE}_subset.vcf.gz"

# Step 3: report
SAMPLES=$(bcftools query -l "${OUTDIR}/${GENE}_subset.vcf.gz" | wc -l)
VARIANTS=$(bcftools view -H "${OUTDIR}/${GENE}_subset.vcf.gz" | wc -l)
echo "Result: $SAMPLES samples, $VARIANTS variants -> ${OUTDIR}/${GENE}_subset.vcf.gz"

# Cleanup intermediate full-cohort file
rm "${OUTDIR}/${GENE}_full.vcf.gz" "${OUTDIR}/${GENE}_full.vcf.gz.csi"
