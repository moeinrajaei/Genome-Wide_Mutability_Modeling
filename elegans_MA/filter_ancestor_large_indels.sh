#!/usr/bin/env bash

set -euo pipefail



module load bcftools



PATTERN="*mutation.vcf"

OUT_SUFFIX="filtered.vcf.gz"



# Accept phased or unphased hom-ref as "ancestor clean"

ACCEPT_PHASED=1   # set to 0 if you only want 0/0



shopt -s nullglob

files=( $PATTERN )



if (( ${#files[@]} == 0 )); then

  echo "[ERROR] No files matched: $PATTERN"

  exit 1

fi



for vcf in "${files[@]}"; do

  echo "============================================================"

  echo "[INFO] Input: $vcf"



  # Compress with bcftools (avoids broken bgzip) + index

  gz="${vcf}.gz"

  if [[ ! -f "$gz" ]]; then

    echo "[INFO] Compressing -> $gz"

    bcftools view -Oz -o "$gz" "$vcf"

  fi



  if [[ ! -f "${gz}.csi" && ! -f "${gz}.tbi" ]]; then

    echo "[INFO] Indexing -> $gz"

    bcftools index -f "$gz"

  fi



  # Read all sample names into an array (0-based indices)

  mapfile -t all_samples < <(bcftools query -l "$gz")



  # Find ancestor sample indices

  anc_idx=()

  anc_names=()

  for i in "${!all_samples[@]}"; do

    s="${all_samples[$i]}"

    if [[ "$s" == *BaerAncestor* || "$s" == XZG0* ]]; then

      anc_idx+=("$i")

      anc_names+=("$s")

    fi

  done



  if (( ${#anc_idx[@]} == 0 )); then

    echo "[WARN] No BaerAncestor* or XZG0* found in $vcf"

    out="${vcf%.vcf}.${OUT_SUFFIX}"

    cp -f "$gz" "$out"

    bcftools index -f "$out"

    continue

  fi



  echo "[INFO] Ancestors found (index -> name):"

  for j in "${!anc_idx[@]}"; do

    echo "  - ${anc_idx[$j]} -> ${anc_names[$j]}"

  done



  # Build expression using numeric indices (most compatible)

  expr=""

  for idx in "${anc_idx[@]}"; do

    if (( ACCEPT_PHASED == 1 )); then

      expr+="(GT[$idx]=\"0/0\" || GT[$idx]=\"0|0\") && "

    else

      expr+="GT[$idx]=\"0/0\" && "

    fi

  done

  expr="${expr% && }"



  out="${vcf%.vcf}.${OUT_SUFFIX}"

  echo "[INFO] Output: $out"

  echo "[INFO] Expr: $expr"


bcftools view -i "${expr} && (TYPE!=\"indel\" || abs(strlen(ALT)-strlen(REF))<=20)" -Oz -o "$out" "$gz"

bcftools index -f "$out"



  before=$(bcftools view -H "$gz"  | wc -l)

  after=$(bcftools view -H "$out" | wc -l)

  echo "[INFO] Sites before: $before"

  echo "[INFO] Sites after : $after"

done



echo "============================================================"

echo "[DONE] Ancestor-only filtering complete."


