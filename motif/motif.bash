#!/bin/sh
salloc -t 2:00:00 --mem=8G

cd /blue/baer/mutability/motif
mkdir GR_Data
cd GR_Data
module load bcftools
module load bedtools

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

##3bp_middle

for i in *.step1.bed; do
    awk '$2-=2' "$i" > "${i/%.step1.bed/.step2_3bp_middle.bed}"
done

for i in *.step2_3bp_middle.bed; do
    awk '$3+=1' "$i" > "${i/%.step2_3bp_middle.bed/.step3_3bp_middle.bed}"
done

for i in *.step3_3bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_middle.bed/.final_3bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_middle.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_3bp_middle.bed/.3bp_middle.txt}"
done


##5bp_middle

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_5bp_middle.bed}"
done

for i in *.step2_5bp_middle.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_5bp_middle.bed/.step3_5bp_middle.bed}"
done

for i in *.step3_5bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_5bp_middle.bed/.final_5bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_5bp_middle.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_5bp_middle.bed/.5bp_middle.txt}"
done


##7bp_middle

for i in *.step1.bed; do
    awk '$2-=4' "$i" > "${i/%.step1.bed/.step2_7bp_middle.bed}"
done

for i in *.step2_7bp_middle.bed; do
    awk '$3+=3' "$i" > "${i/%.step2_7bp_middle.bed/.step3_7bp_middle.bed}"
done

for i in *.step3_7bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_7bp_middle.bed/.final_7bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_7bp_middle.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_7bp_middle.bed/.7bp_middle.txt}"
done


##9bp_middle

for i in *.step1.bed; do
    awk '$2-=5' "$i" > "${i/%.step1.bed/.step2_9bp_middle.bed}"
done

for i in *.step2_9bp_middle.bed; do
    awk '$3+=4' "$i" > "${i/%.step2_9bp_middle.bed/.step3_9bp_middle.bed}"
done

for i in *.step3_9bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_9bp_middle.bed/.final_9bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_9bp_middle.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_9bp_middle.bed/.9bp_middle.txt}"
done


##99bp_middle

for i in *.step1.bed; do
    awk '$2-=50' "$i" > "${i/%.step1.bed/.step2_99bp_middle.bed}"
done

for i in *.step2_99bp_middle.bed; do
    awk '$3+=49' "$i" > "${i/%.step2_99bp_middle.bed/.step3_99bp_middle.bed}"
done

for i in *.step3_99bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_99bp_middle.bed/.final_99bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_99bp_middle.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_99bp_middle.bed/.99bp_middle.txt}"
done


##3bp_first

for i in *.step1.bed; do
    awk '$2-=1' "$i" > "${i/%.step1.bed/.step2_3bp_first.bed}"
done

for i in *.step2_3bp_first.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_3bp_first.bed/.step3_3bp_first.bed}"
done

for i in *.step3_3bp_first.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_first.bed/.final_3bp_first.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_first.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_3bp_first.bed/.3bp_first.txt}"
done


##3bp_last

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_3bp_last.bed}"
done

for i in *.step2_3bp_last.bed; do
    awk '$3+=0' "$i" > "${i/%.step2_3bp_last.bed/.step3_3bp_last.bed}"
done

for i in *.step3_3bp_last.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_last.bed/.final_3bp_last.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_last.bed; do
    bedtools getfasta -bedOut -fi /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa -bed "$i" > "${i/%.final_3bp_last.bed/.3bp_last.txt}"
done


##merging the files:
for prefix in mev N2 PB306; do
  for type in snp indel; do

    files=(
      ${prefix}.${type}.3bp_first.txt
      ${prefix}.${type}.3bp_last.txt
      ${prefix}.${type}.3bp_middle.txt
      ${prefix}.${type}.5bp_middle.txt
      ${prefix}.${type}.7bp_middle.txt
      ${prefix}.${type}.9bp_middle.txt
      ${prefix}.${type}.99bp_middle.txt
    )

    # skip if first file doesn't exist (prevents errors)
    [[ -f "${files[0]}" ]] || continue

    # build header
    header=""
    for f in "${files[@]}"; do
      ncol=$(awk '{print NF; exit}' "$f")
      base=$(basename "$f" .txt)
      for ((i=1; i<=ncol; i++)); do
        header+="${base}_col${i}\t"
      done
    done
    header=${header%$'\t'}

    # write merged file
    {
      echo -e "$header"
      paste "${files[@]}"
    } > ${prefix}.${type}.all_merged.with_header.txt

  done
done


rm *.bed
rm *.3bp_last.txt
rm *.3bp_middle.txt
rm *.3bp_first.txt
rm *.5bp_middle.txt
rm *.7bp_middle.txt
rm *.9bp_middle.txt
rm *.99bp_middle.txt
rm *.tab.txt


###C_Strains data

cd ..
cd C_Strains

mv C.bp.GVCFs.indel.3x.mutation.filtered.vcf C.indel.vcf
mv C.bp.GVCFs.snp.3x.mutation.filtered.vcf C.snp.vcf
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

##3bp_middle

for i in *.step1.bed; do
    awk '$2-=2' "$i" > "${i/%.step1.bed/.step2_3bp_middle.bed}"
done

for i in *.step2_3bp_middle.bed; do
    awk '$3+=1' "$i" > "${i/%.step2_3bp_middle.bed/.step3_3bp_middle.bed}"
done

for i in *.step3_3bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_middle.bed/.final_3bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_3bp_middle.bed/.3bp_middle.txt}"
done


##5bp_middle

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_5bp_middle.bed}"
done

for i in *.step2_5bp_middle.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_5bp_middle.bed/.step3_5bp_middle.bed}"
done

for i in *.step3_5bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_5bp_middle.bed/.final_5bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_5bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_5bp_middle.bed/.5bp_middle.txt}"
done


##7bp_middle

for i in *.step1.bed; do
    awk '$2-=4' "$i" > "${i/%.step1.bed/.step2_7bp_middle.bed}"
done

for i in *.step2_7bp_middle.bed; do
    awk '$3+=3' "$i" > "${i/%.step2_7bp_middle.bed/.step3_7bp_middle.bed}"
done

for i in *.step3_7bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_7bp_middle.bed/.final_7bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_7bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_7bp_middle.bed/.7bp_middle.txt}"
done


##9bp_middle

for i in *.step1.bed; do
    awk '$2-=5' "$i" > "${i/%.step1.bed/.step2_9bp_middle.bed}"
done

for i in *.step2_9bp_middle.bed; do
    awk '$3+=4' "$i" > "${i/%.step2_9bp_middle.bed/.step3_9bp_middle.bed}"
done

for i in *.step3_9bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_9bp_middle.bed/.final_9bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_9bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_9bp_middle.bed/.9bp_middle.txt}"
done


##99bp_middle

for i in *.step1.bed; do
    awk '$2-=50' "$i" > "${i/%.step1.bed/.step2_99bp_middle.bed}"
done

for i in *.step2_99bp_middle.bed; do
    awk '$3+=49' "$i" > "${i/%.step2_99bp_middle.bed/.step3_99bp_middle.bed}"
done

for i in *.step3_99bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_99bp_middle.bed/.final_99bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_99bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_99bp_middle.bed/.99bp_middle.txt}"
done


##3bp_first

for i in *.step1.bed; do
    awk '$2-=1' "$i" > "${i/%.step1.bed/.step2_3bp_first.bed}"
done

for i in *.step2_3bp_first.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_3bp_first.bed/.step3_3bp_first.bed}"
done

for i in *.step3_3bp_first.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_first.bed/.final_3bp_first.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_first.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_3bp_first.bed/.3bp_first.txt}"
done


##3bp_last

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_3bp_last.bed}"
done

for i in *.step2_3bp_last.bed; do
    awk '$3+=0' "$i" > "${i/%.step2_3bp_last.bed/.step3_3bp_last.bed}"
done

for i in *.step3_3bp_last.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_last.bed/.final_3bp_last.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_last.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa -bed "$i" > "${i/%.final_3bp_last.bed/.3bp_last.txt}"
done


##merging the files:
for prefix in C; do
  for type in snp indel; do

    files=(
      ${prefix}.${type}.3bp_first.txt
      ${prefix}.${type}.3bp_last.txt
      ${prefix}.${type}.3bp_middle.txt
      ${prefix}.${type}.5bp_middle.txt
      ${prefix}.${type}.7bp_middle.txt
      ${prefix}.${type}.9bp_middle.txt
      ${prefix}.${type}.99bp_middle.txt
    )

    # skip if first file doesn't exist (prevents errors)
    [[ -f "${files[0]}" ]] || continue

    # build header
    header=""
    for f in "${files[@]}"; do
      ncol=$(awk '{print NF; exit}' "$f")
      base=$(basename "$f" .txt)
      for ((i=1; i<=ncol; i++)); do
        header+="${base}_col${i}\t"
      done
    done
    header=${header%$'\t'}

    # write merged file
    {
      echo -e "$header"
      paste "${files[@]}"
    } > ${prefix}.${type}.all_merged.with_header.txt

  done
done


rm *.bed
rm *.3bp_last.txt
rm *.3bp_middle.txt
rm *.3bp_first.txt
rm *.5bp_middle.txt
rm *.7bp_middle.txt
rm *.9bp_middle.txt
rm *.99bp_middle.txt
rm *.tab.txt	


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

##3bp_middle

for i in *.step1.bed; do
    awk '$2-=2' "$i" > "${i/%.step1.bed/.step2_3bp_middle.bed}"
done

for i in *.step2_3bp_middle.bed; do
    awk '$3+=1' "$i" > "${i/%.step2_3bp_middle.bed/.step3_3bp_middle.bed}"
done

for i in *.step3_3bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_middle.bed/.final_3bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_3bp_middle.bed/.3bp_middle.txt}"
done


##5bp_middle

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_5bp_middle.bed}"
done

for i in *.step2_5bp_middle.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_5bp_middle.bed/.step3_5bp_middle.bed}"
done

for i in *.step3_5bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_5bp_middle.bed/.final_5bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_5bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_5bp_middle.bed/.5bp_middle.txt}"
done


##7bp_middle

for i in *.step1.bed; do
    awk '$2-=4' "$i" > "${i/%.step1.bed/.step2_7bp_middle.bed}"
done

for i in *.step2_7bp_middle.bed; do
    awk '$3+=3' "$i" > "${i/%.step2_7bp_middle.bed/.step3_7bp_middle.bed}"
done

for i in *.step3_7bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_7bp_middle.bed/.final_7bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_7bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_7bp_middle.bed/.7bp_middle.txt}"
done


##9bp_middle

for i in *.step1.bed; do
    awk '$2-=5' "$i" > "${i/%.step1.bed/.step2_9bp_middle.bed}"
done

for i in *.step2_9bp_middle.bed; do
    awk '$3+=4' "$i" > "${i/%.step2_9bp_middle.bed/.step3_9bp_middle.bed}"
done

for i in *.step3_9bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_9bp_middle.bed/.final_9bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_9bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_9bp_middle.bed/.9bp_middle.txt}"
done


##99bp_middle

for i in *.step1.bed; do
    awk 'BEGIN{OFS="\t"} { tmp=$2-50; if (tmp<0) tmp=0; $2=tmp; print }' "$i" > "${i/%.step1.bed/.step2_99bp_middle.bed}"
done

for i in *.step2_99bp_middle.bed; do
    awk '$3+=49' "$i" > "${i/%.step2_99bp_middle.bed/.step3_99bp_middle.bed}"
done

for i in *.step3_99bp_middle.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_99bp_middle.bed/.final_99bp_middle.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_99bp_middle.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_99bp_middle.bed/.99bp_middle.txt}"
done


##3bp_first

for i in *.step1.bed; do
    awk '$2-=1' "$i" > "${i/%.step1.bed/.step2_3bp_first.bed}"
done

for i in *.step2_3bp_first.bed; do
    awk '$3+=2' "$i" > "${i/%.step2_3bp_first.bed/.step3_3bp_first.bed}"
done

for i in *.step3_3bp_first.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_first.bed/.final_3bp_first.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_first.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_3bp_first.bed/.3bp_first.txt}"
done


##3bp_last

for i in *.step1.bed; do
    awk '$2-=3' "$i" > "${i/%.step1.bed/.step2_3bp_last.bed}"
done

for i in *.step2_3bp_last.bed; do
    awk '$3+=0' "$i" > "${i/%.step2_3bp_last.bed/.step3_3bp_last.bed}"
done

for i in *.step3_3bp_last.bed; do
    awk -v OFS='\t' '{ $1=$1; print }' "$i" > "${i/%.step3_3bp_last.bed/.final_3bp_last.bed}"
done

#extract the sequence contex based on the bed file:
for i in *.final_3bp_last.bed; do
    bedtools getfasta -bedOut -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa -bed "$i" > "${i/%.final_3bp_last.bed/.3bp_last.txt}"
done
 
	  
##merging the files:
for prefix in E G J C M XZ7 XZ8; do
  for type in snp indel; do

    files=(
      ${prefix}.${type}.3bp_first.txt
      ${prefix}.${type}.3bp_last.txt
      ${prefix}.${type}.3bp_middle.txt
      ${prefix}.${type}.5bp_middle.txt
      ${prefix}.${type}.7bp_middle.txt
      ${prefix}.${type}.9bp_middle.txt
      ${prefix}.${type}.99bp_middle.txt
    )

    # skip if first file doesn't exist (prevents errors)
    [[ -f "${files[0]}" ]] || continue

    # build header
    header=""
    for f in "${files[@]}"; do
      ncol=$(awk '{print NF; exit}' "$f")
      base=$(basename "$f" .txt)
      for ((i=1; i<=ncol; i++)); do
        header+="${base}_col${i}\t"
      done
    done
    header=${header%$'\t'}

    # write merged file
    {
      echo -e "$header"
      paste "${files[@]}"
    } > ${prefix}.${type}.all_merged.with_header.txt

  done
done


rm *.bed
rm *.3bp_last.txt
rm *.3bp_middle.txt
rm *.3bp_first.txt
rm *.5bp_middle.txt
rm *.7bp_middle.txt
rm *.9bp_middle.txt
rm *.99bp_middle.txt
rm *.tab.txt	



#End


