#!/bin/sh
salloc -t 2:00:00 --mem=8G

cd /blue/baer/mutability
mkdir GC_content
cd GC_content


module load bcftools
module load bedtools

### checking the number of variants:
grep -c '^[^#]' *.vcf > checking_number_variants.txt


###GR_Data
mkdir GR_Data
cd GR_Data

GENOME_FA="/blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa"
FAI="${GENOME_FA}.fai"   


###10bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

for i in *.step1.bed; do
    awk '$2-=6' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=5' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10bp.txt}"
done

rm *.final.bed

### 100bp:
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 51 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-51; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=50' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100bp.txt}"
done

rm *.final.bed

### 1k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done


# step2: subtract 501 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-501; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 500 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 500, then cap at max (BED end is 0-based, half-open)
          end = $3 + 500;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_1k_bp.txt}"
done

rm *.final.bed

### 10k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 5001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-5001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 5000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 5000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 5000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10k_bp.txt}"
done

rm *.final.bed

### 100k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 50001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-50001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 50000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 50000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 50000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

# keep tabs clean
for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
    bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100k_bp.txt}"
done

rm *.final.bed

# Merge all matching files:
paste -d $'\t' $(printf "%s\n" *mev.snp.final.GC* | sort -V) > mev.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *N2.snp.final.GC* | sort -V) > N2.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *PB306.snp.final.GC* | sort -V) > PB306.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *mev.indel.final.GC* | sort -V) > mev.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *N2.indel.final.GC* | sort -V) > N2.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *PB306.indel.final.GC* | sort -V) > PB306.indel.final.GC.tsv

rm *.txt
rm *.bed


###C_Strains data
cd ..
mkdir C_Strains
cd C_Strains

GENOME_FA="/orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa"
FAI="${GENOME_FA}.fai"   

echo "$GENOME_FA"
echo "$FAI"

###10bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

for i in *.step1.bed; do
    awk '$2-=6' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=5' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10bp.txt}"
done

rm *.final.bed

### 100bp:
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 51 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-51; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=50' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100bp.txt}"
done

rm *.final.bed

### 1k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done


# step2: subtract 501 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-501; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 500 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 500, then cap at max (BED end is 0-based, half-open)
          end = $3 + 500;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_1k_bp.txt}"
done

rm *.final.bed

### 10k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 5001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-5001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 5000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 5000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 5000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10k_bp.txt}"
done

rm *.final.bed

### 100k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 50001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-50001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 50000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 50000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 50000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

# keep tabs clean
for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
    bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100k_bp.txt}"
done

rm *.final.bed

# Merge all matching files:
paste -d $'\t' $(printf "%s\n" *C.snp.final.GC* | sort -V) > C.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *C.indel.final.GC* | sort -V) > C.indel.final.GC.tsv


rm *.txt
rm *.bed



###E_G_J_M_XZ7_XZ8_C_Strains data mapping to WS292.
cd ..
mkdir E_G_J_M_XZ7_XZ8_C_Strains
cd E_G_J_M_XZ7_XZ8_C_Strains

GENOME_FA="/orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa"
FAI="${GENOME_FA}.fai"   

cp -r -p /orange/baer/Raw_vcf_data/Mutation_Files/*.3x.mutation.filtered.vcf.gz .

gunzip *.gz

mv G.bp.GVCFs.snp.3x.mutation.filtered.vcf G.snp.vcf
mv G.bp.GVCFs.indel.3x.mutation.filtered.vcf G.indel.vcf
mv E.bp.GVCFs.snp.3x.mutation.filtered.vcf E.snp.vcf
mv E.bp.GVCFs.indel.3x.mutation.filtered.vcf E.indel.vcf
mv C.bp.GVCFs.snp.3x.mutation.filtered.vcf C.snp.vcf
mv C.bp.GVCFs.indel.3x.mutation.filtered.vcf C.indel.vcf
mv XZ8.bp.GVCFs.snp.3x.mutation.filtered.vcf XZ8.snp.vcf
mv XZ8.bp.GVCFs.indel.3x.mutation.filtered.vcf XZ8.indel.vcf
mv XZ7.bp.GVCFs.snp.3x.mutation.filtered.vcf XZ7.snp.vcf
mv XZ7.bp.GVCFs.indel.3x.mutation.filtered.vcf XZ7.indel.vcf
mv M.bp.GVCFs.snp.3x.mutation.filtered.vcf M.snp.vcf
mv M.bp.GVCFs.indel.3x.mutation.filtered.vcf M.indel.vcf
mv J.bp.GVCFs.snp.3x.mutation.filtered.vcf J.snp.vcf
mv J.bp.GVCFs.indel.3x.mutation.filtered.vcf J.indel.vcf

###10bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

for i in *.step1.bed; do
    awk '$2-=6' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=5' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10bp.txt}"
done

rm *.final.bed

### 100bp:
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 51 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-51; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk '$3+=50' "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100bp.txt}"
done

rm *.final.bed

### 1k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done


# step2: subtract 501 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-501; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 500 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 500, then cap at max (BED end is 0-based, half-open)
          end = $3 + 500;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_1k_bp.txt}"
done

rm *.final.bed

### 10k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 5001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-5001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 5000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 5000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 5000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10k_bp.txt}"
done

rm *.final.bed

### 100k_bp
## Prepare the genomic_region_file:
for i in *.vcf;do
bcftools query -f '%CHROM %POS\n' "$i" > "${i/%.vcf/.bed}"
done

for i in *.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.bed/.tab.bed}"
done

for i in *.tab.bed; do
cut -f2 "$i" > "${i/%.tab.bed/.tab.txt}" 
done

for i in *.tab.bed; do
    file=`basename $i .bed`
    paste "$i" "$file.txt" > "${i/%.tab.bed/.step1.bed}"
done

# step2: subtract 50001 from col2, but set to 0 if result < 0
for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-50001; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2.bed}"
done

# step3: add 50000 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  # read chrom sizes from .fai (field2)
        {
          max = sz[$1];
          if (max == "" ) {                          # chrom not found in .fai
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 50000, then cap at max (BED end is 0-based, half-open)
          end = $3 + 50000;
          if (end > max) end = max;
          if (end < 0) end = 0;                      # paranoia guard
          $3 = end;

          # ensure start < end; if not, make a 1bp interval at the end cap
          if ($2 >= $3) $2 = ($3 > 0 ? $3 - 1 : 0);

          print
        }' "$FAI" "$i" > "${i/%.step2.bed/.step3.bed}"
done

# keep tabs clean
for i in *.step3.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3.bed/.final.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed
rm *.indel.bed

## Calculate GC-content per region:
for i in *.final.bed; do
    bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100k_bp.txt}"
done

rm *.final.bed

# Merge all matching files:
paste -d $'\t' $(printf "%s\n" *E.snp.final.GC* | sort -V) > E.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *E.indel.final.GC* | sort -V) > E.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *G.snp.final.GC* | sort -V) > G.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *G.indel.final.GC* | sort -V) > G.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *C.snp.final.GC* | sort -V) > C.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *C.indel.final.GC* | sort -V) > C.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *J.snp.final.GC* | sort -V) > J.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *J.indel.final.GC* | sort -V) > J.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *M.snp.final.GC* | sort -V) > M.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *M.indel.final.GC* | sort -V) > M.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *XZ7.snp.final.GC* | sort -V) > XZ7.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *XZ7.indel.final.GC* | sort -V) > XZ7.indel.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *XZ8.snp.final.GC* | sort -V) > XZ8.snp.final.GC.tsv
paste -d $'\t' $(printf "%s\n" *XZ8.indel.final.GC* | sort -V) > XZ8.indel.final.GC.tsv


rm *.txt


#End
