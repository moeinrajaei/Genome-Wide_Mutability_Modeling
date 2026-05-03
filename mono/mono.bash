#!/bin/sh
salloc -t 2:00:00 --mem=8G

module load vcftools


### extracting the variants(SNVs) from mononucleotides run:

# Preparing the data:

cd /blue/baer/mutability
mkdir mono
cd mono

cp -r -p /blue/baer/mutability/motif/GR_Data/*.vcf .

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


cp -r -p /blue/baer/share/moeinraja/reference_WS263/mono.5bp.final.bed .

for i in *.vcf; do
vcftools --vcf "$i" --extract-FORMAT-info GT --bed mono.5bp.final.bed --out "${i/%.vcf/.vcf.mono.5bp}";
done

for i in *.vcf.mono.5bp.GT.FORMAT; do 
cut -f1,2 "$i" > "${i/%.GT.FORMAT/.GT.FORMAT.bed}"
done

for i in *.mono.5bp.GT.FORMAT.bed; do
    x=`basename $i .mono.5bp.GT.FORMAT.bed`
    vcftools --vcf "$x" --recode --out "${i/%.mono.5bp.GT.FORMAT.bed/.mono.5bp}" --positions "$i" --recode-INFO-all
done

rm *.vcf.mono.5bp.GT.FORMAT
rm *.vcf.mono.5bp.GT.FORMAT.bed
rm *.vcf.mono.5bp.log

#End







