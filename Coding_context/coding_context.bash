#!/bin/sh
salloc -t 2:00:00 --mem=8G

cd /blue/baer/mutability
mkdir coding_context 
cd coding_context

module load java/23

wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip -o snpEff_latest_core.zip
cd snpEff

###annotations
#Preparing data:
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

#using snpEff

for i in ../*.vcf; do
java -Xmx8g -jar snpEff.jar -stats "${i/%.vcf/.ex1.html}" -v WBcel235.99 "$i" > "${i/%.vcf/.snpEff.vcf}"
done

for i in ../*.snpEff.vcf; do
java -Xmx2g -jar SnpSift.jar extractFields "$i" \
"CHROM" "POS" "REF" "ALT" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" \
"ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_C" "ANN[0].HGVS_P" "ANN[0].GENEID" \
> "${i/%.snpEff.vcf/.highest_impact_snpEff.tsv}"
done

cd ..
rm *.snpEff.vcf
rm *.genes.txt

#using vep
mkdir vep
cd vep

module load vep/113
module load bcftools

for i in ../*.vcf; do
vep -i "$i" -o "${i/%.vcf/.vep.vcf}" --species caenorhabditis_elegans --vcf --cache   --flag_pick
done

# best representative annotation based on: canonical transcript, transcript support level, biotype, consequence severity
for i in ../*.vep.vcf; do
bcftools +split-vep "$i" -s pick -f '%CHROM\t%POS\t%SYMBOL\t%Consequence\t%IMPACT\n' > "${i/%.vep.vcf/.best_annotation_vep.tsv}"
done

# The closest to “highest impact” in the sense of: stop_gained > frameshift > missense > synonymous > intronic, etc.
for i in ../*.vep.vcf; do
bcftools +split-vep "$i" -s worst -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\n' > "${i/%.vep.vcf/.worst_annotation_vep.tsv}"
done

rm ../*.vep.vcf

# Example: A variant may have:
#transcript A: synonymous_variant (LOW)
#transcript B: stop_gained (HIGH)
#Then: 
#-s pick may return transcript A (if canonical)
#-s worst will return transcript B

#check the reference genome, using the above vep command the references genome is WS282. 
#113_WBcel235 = Ensembl release 113 + WormBase WS282 annotation
#Genome assembly: WBcel235 / ce11
#Gene annotation version: WormBase WS282

#End
