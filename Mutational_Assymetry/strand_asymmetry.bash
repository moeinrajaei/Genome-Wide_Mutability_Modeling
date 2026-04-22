#!/bin/sh
salloc -t 2:00:00 --mem=8G

# + is forward
# - is reverse

module load bedops
module load bcftools
module load bedtools
module load vcftools

cd /blue/baer/mutability
mkdir strand_asymmetry
cd strand_asymmetry

### converting GFF3 to BED format
gff2bed < c_elegans.PRJNA13758.WS263.annotations.gff3 > c_elegans.PRJNA13758.WS263.annotations.gff3.bed
gff2bed < c_elegans.PRJNA13758.WS292.annotations.gff3 > c_elegans.PRJNA13758.WS292.annotations.gff3.bed

###GR_Data
mkdir GR_Data
cd GR_Data

cp -r -p /blue/baer/mutability/motif/GR_Data/*.vcf .

### extract list of positions:
 
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
    awk '$2-=1' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step2.bed/.final.bed}"
done

for i in *.final.bed; do
bedtools intersect -a ../c_elegans.PRJNA13758.WS263.annotations.gff3.bed -b "$i" > "${i/%.final.bed/.gff3.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.snp.bed
rm *.indel.bed

# filters lines where the seventh column is WormBase:
for i in *.gff3.bed; do
cat "$i"  | awk '$7=="WormBase" {print $0}' > "${i/%.gff3.bed/.wormbase.bed}"
done

for i in *.wormbase.bed; do
  awk '$8=="gene" {print $1, $2, $3, $6, $7, $8}' "$i" > "${i/%.wormbase.bed/.wormbase.final.bed}"
done

#remove the duplicates:
for i in *.wormbase.final.bed; do
awk '!seen[$3]++' "$i"  > "${i/%.wormbase.final.bed/.wormbase.unique.bed}"
done

rm *.final.bed 
rm *.wormbase.bed
rm *.gff3.bed

###E_G_J_M_XZ7_XZ8_C_Strains data mapping to WS292.
cd ..
mkdir E_G_J_M_XZ7_XZ8_C_Strains
cd E_G_J_M_XZ7_XZ8_C_Strains

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

### extract list of positions:
 
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
    awk '$2-=1' "$i" > "${i/%.step1.bed/.step2.bed}"
done

for i in *.step2.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step2.bed/.final.bed}"
done

for i in *.final.bed; do
bedtools intersect -a ../c_elegans.PRJNA13758.WS292.annotations.gff3.bed -b "$i" > "${i/%.final.bed/.gff3.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed
rm *.snp.bed
rm *.indel.bed

# filters lines where the seventh column is WormBase:
for i in *.gff3.bed; do
cat "$i"  | awk '$7=="WormBase" {print $0}' > "${i/%.gff3.bed/.wormbase.bed}"
done

for i in *.wormbase.bed; do
  awk '$8=="gene" {print $1, $2, $3, $6, $7, $8}' "$i" > "${i/%.wormbase.bed/.wormbase.final.bed}"
done

#remove the duplicates:
for i in *.wormbase.final.bed; do
awk '!seen[$3]++' "$i"  > "${i/%.wormbase.final.bed/.wormbase.unique.bed}"
done

rm *.final.bed 
rm *.wormbase.bed
rm *.gff3.bed

#End

