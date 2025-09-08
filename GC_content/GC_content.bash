#!/bin/sh
salloc -t 2:00:00 --mem=8G

cd /vast/palmer/scratch/townsend/mr2683/project/test_7_20_2025

mkdir GC_content
cd GC_content
cp -r -p ../*.snp.vcf .

module load BCFtools
module load BEDTools

### checking the number of variants:
grep -c '^[^#]' *.vcf > checking_number_variants.txt

### Prepare the genomic_region_file:
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

rm *.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.step3.bed
rm *.snp.bed

### Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi ../c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.bed/.GC.txt}"
done

rm *.final.bed

#End
