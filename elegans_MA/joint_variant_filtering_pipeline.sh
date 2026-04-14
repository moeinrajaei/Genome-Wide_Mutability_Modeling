#!/bin/bash
#SBATCH --job-name=j_s10_s12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=m.rifat@ufl.edu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --time=480:00:00
#SBATCH --output=j_s10_s12_%j.out
#SBATCH --account=juannanzhou --qos=juannanzhou

set -euo pipefail

# =============================================================================
# SETTINGS — edit here only
# =============================================================================
WORKDIR="/orange/juannanzhou/Rifat_CB_raw/J"
REF="c_elegans.PRJNA13758.WS292.genomic.fa"
BASE="J.bp.GVCFs"          # base name for all output files

cd "$WORKDIR"
module load gatk/
module load bcftools/

JAVA_MEM="-Xmx8g"

# Step 9 (GenotypeGVCFs) was run separately.
# This script assumes J.bp.GVCFs.vcf already exists in $WORKDIR.
JOINT="${BASE}.vcf"
[[ -f "$JOINT" ]] || { echo "ERROR: Joint VCF not found: $JOINT" >&2; exit 1; }
echo "[INFO] Joint VCF found: $JOINT — skipping Step 9."

# =============================================================================
# STEP 10: Split SNPs and indels
# Produces: .snp.vcf, .indel.vcf  <- intermediate, deleted later
# =============================================================================
echo "[STEP 10] Splitting SNPs and indels..."

SNP="${BASE}.snp.vcf"
INDEL="${BASE}.indel.vcf"

gatk --java-options "$JAVA_MEM" SelectVariants \
    -R "$REF" -V "$JOINT" \
    --select-type-to-include SNP \
    -O "$SNP"

gatk --java-options "$JAVA_MEM" SelectVariants \
    -R "$REF" -V "$JOINT" \
    --select-type-to-include INDEL \
    -O "$INDEL"

echo "[STEP 10] Done: $SNP, $INDEL"

# =============================================================================
# STEP 11: DP filtering (3x and 10x) for each variant type
# For each type x depth: apply filter -> nocall -> remove filtered sites
# All intermediate files deleted at end
# =============================================================================
# NOTE: run_dp_filter writes its output path to a temp file rather than
# using $() command substitution. $() captures ALL stdout including GATK
# verbose progress logs, which silently swallows diagnostic output and
# can mix log lines into the captured variable value.
# =============================================================================

run_dp_filter() {
    local INPUT="$1"
    local DP="$2"
    local TAG="$3"
    local RETFILE="$4"    # temp file to write the output path into

    local FILTERED="${BASE}.${TAG}.vcf"
    local NOCALL="${BASE}.${TAG}.nocall.vcf"
    local CLEAN="${BASE}.${TAG}.filtered.vcf"

    echo "[STEP 11] DP${DP} filter on $INPUT..."

    gatk --java-options "$JAVA_MEM" VariantFiltration \
        -R "$REF" -V "$INPUT" \
        --genotype-filter-name "DP" \
        --genotype-filter-expression "DP < ${DP}.0" \
        -O "$FILTERED"

    gatk --java-options "$JAVA_MEM" SelectVariants \
        -V "$FILTERED" \
        --set-filtered-gt-to-nocall \
        -O "$NOCALL"

    gatk --java-options "$JAVA_MEM" SelectVariants \
        -R "$REF" -V "$NOCALL" \
        --exclude-filtered true \
        --exclude-non-variants true \
        -O "$CLEAN"

    # Write output path to temp file instead of stdout capture
    echo "$CLEAN" > "$RETFILE"
}

run_dp_filter "$SNP"   3  "snp.3x"    .snp_3x_out
run_dp_filter "$INDEL" 3  "indel.3x"  .indel_3x_out
run_dp_filter "$SNP"   10 "snp.10x"   .snp_10x_out
run_dp_filter "$INDEL" 10 "indel.10x" .indel_10x_out

SNP_3X_FILTERED=$(cat .snp_3x_out)
INDEL_3X_FILTERED=$(cat .indel_3x_out)
SNP_10X_FILTERED=$(cat .snp_10x_out)
INDEL_10X_FILTERED=$(cat .indel_10x_out)
rm -f .snp_3x_out .indel_3x_out .snp_10x_out .indel_10x_out

echo "[STEP 11] Done."

# =============================================================================
# STEP 12: Biallelic filtering
# Produces 4 final files <- KEPT
# =============================================================================
echo "[STEP 12] Applying biallelic filter..."

SNP_3X_FINAL="${BASE}.snp.3x.filtered.biallelic.vcf"
INDEL_3X_FINAL="${BASE}.indel.3x.filtered.biallelic.vcf"
SNP_10X_FINAL="${BASE}.snp.10x.filtered.biallelic.vcf"
INDEL_10X_FINAL="${BASE}.indel.10x.filtered.biallelic.vcf"

bcftools view --max-alleles 2 "$SNP_3X_FILTERED"   -o "$SNP_3X_FINAL"
bcftools view --max-alleles 2 "$INDEL_3X_FILTERED"  -o "$INDEL_3X_FINAL"
bcftools view --max-alleles 2 "$SNP_10X_FILTERED"   -o "$SNP_10X_FINAL"
bcftools view --max-alleles 2 "$INDEL_10X_FILTERED" -o "$INDEL_10X_FINAL"

echo "[STEP 12] Done."

# =============================================================================
# CLEANUP: Remove all intermediate files, keep joint VCF and 4 biallelic VCFs
# =============================================================================
echo "[CLEANUP] Removing intermediate files..."

rm -f \
    "${BASE}.snp.vcf"              "${BASE}.snp.vcf.idx" \
    "${BASE}.indel.vcf"            "${BASE}.indel.vcf.idx" \
    "${BASE}.snp.3x.vcf"           "${BASE}.snp.3x.vcf.idx" \
    "${BASE}.snp.3x.nocall.vcf"    "${BASE}.snp.3x.nocall.vcf.idx" \
    "${BASE}.snp.3x.filtered.vcf"  "${BASE}.snp.3x.filtered.vcf.idx" \
    "${BASE}.indel.3x.vcf"         "${BASE}.indel.3x.vcf.idx" \
    "${BASE}.indel.3x.nocall.vcf"  "${BASE}.indel.3x.nocall.vcf.idx" \
    "${BASE}.indel.3x.filtered.vcf"  "${BASE}.indel.3x.filtered.vcf.idx" \
    "${BASE}.snp.10x.vcf"           "${BASE}.snp.10x.vcf.idx" \
    "${BASE}.snp.10x.nocall.vcf"    "${BASE}.snp.10x.nocall.vcf.idx" \
    "${BASE}.snp.10x.filtered.vcf"  "${BASE}.snp.10x.filtered.vcf.idx" \
    "${BASE}.indel.10x.vcf"         "${BASE}.indel.10x.vcf.idx" \
    "${BASE}.indel.10x.nocall.vcf"  "${BASE}.indel.10x.nocall.vcf.idx" \
    "${BASE}.indel.10x.filtered.vcf"  "${BASE}.indel.10x.filtered.vcf.idx"

echo "[CLEANUP] Done."

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo " Pipeline complete. Files retained:"
echo "  [joint]     $JOINT"
echo "  [biallelic] $SNP_3X_FINAL"
echo "  [biallelic] $INDEL_3X_FINAL"
echo "  [biallelic] $SNP_10X_FINAL"
echo "  [biallelic] $INDEL_10X_FINAL"
echo "============================================================"
