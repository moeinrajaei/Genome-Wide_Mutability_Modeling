#!/bin/sh
salloc -t 2:00:00 --mem=8G

###GC_Content_nonmutated
cd /blue/baer/mutability
mkdir nonmutated
cd nonmutated
mkdir GC_content
cd GC_content

cp -r -p ../random_1000000_non_mutant_sites.txt .

awk 'NR>1 {print $2"\t"$3}' random_1000000_non_mutant_sites.txt | sort -k1,1 -k2,2n > random_1000000_non_mutant_sites.bed

module load bcftools
module load bedtools

GENOME_FA="/orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa"
FAI="${GENOME_FA}.fai"   

###10bp
## Prepare the genomic_region_file:
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

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10bp.txt}"
done

rm *.final.bed

### 100bp:
## Prepare the genomic_region_file:
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

# step3: add 50 to col3; cap (and fix) to chromosome size from .fai
for i in *.step2.bed; do
    awk -v OFS='\t' '
        NR==FNR { sz[$1]=$2; next }                  
        {
          max = sz[$1];
          if (max == "" ) {                          
            printf("WARNING: no size for %s\n", $1) > "/dev/stderr";
            next;
          }

          # clamp start into [0, max-1]
          if ($2 < 0) $2 = 0;
          if ($2 >= max) $2 = (max>0 ? max-1 : 0);

          # extend end by 50, then cap at max (BED end is 0-based, half-open)
          end = $3 + 50;
          if (end > max) end = max;
          if (end < 0) end = 0;                     
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

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100bp.txt}"
done

rm *.final.bed

### 1k_bp
## Prepare the genomic_region_file:
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

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_1k_bp.txt}"
done

rm *.final.bed

### 10k_bp
## Prepare the genomic_region_file:
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

## Calculate GC-content per region:
for i in *.final.bed; do
bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_10k_bp.txt}"
done

rm *.final.bed

### 100k_bp
## Prepare the genomic_region_file:
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

## Calculate GC-content per region:
for i in *.final.bed; do
    bedtools nuc -fi "$GENOME_FA" -bed "$i" > "${i/%.bed/.GC_100k_bp.txt}"
done

rm *.final.bed

# Merge all matching files:
paste -d $'\t' $(printf "%s\n" *random_1000000_non_mutant_sites.final.GC* | sort -V) > random_1000000_non_mutant_sites.final.GC.tsv

#End





### motif
cd /blue/baer/mutability/nonmutated/
mkdir motif
cd motif

cp -r -p ../random_1000000_non_mutant_sites.txt .

awk 'NR>1 {print $2"\t"$3}' random_1000000_non_mutant_sites.txt | sort -k1,1 -k2,2n > random_1000000_non_mutant_sites.bed


### extract list of positions:
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


# create genome size file from fasta index
module load samtools
cut -f1,2 /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa.fai > c_elegans.genome

# clip BED intervals to chromosome boundaries
for i in *.final_99bp_middle.bed; do
    bedtools slop \
        -i "$i" \
        -g c_elegans.genome \
        -b 0 > "${i/%.final_99bp_middle.bed/.clipped.bed}"

    bedtools getfasta \
        -bedOut \
        -fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa \
        -bed "${i/%.final_99bp_middle.bed/.clipped.bed}" \
        > "${i/%.final_99bp_middle.bed/.99bp_middle.txt}"
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
for prefix in random_1000000_non_mutant_sites; do
    files=(
      ${prefix}.3bp_first.txt
      ${prefix}.3bp_last.txt
      ${prefix}.3bp_middle.txt
      ${prefix}.5bp_middle.txt
      ${prefix}.7bp_middle.txt
      ${prefix}.9bp_middle.txt
      ${prefix}.99bp_middle.txt
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


###strand_asymmetry


# + is forward
# - is reverse

module load bedops
module load bcftools
module load bedtools
module load vcftools

cd /blue/baer/mutability/nonmutated
mkdir strand_asymmetry
cd strand_asymmetry

#Preparing the data:
cp -r -p /blue/baer/mutability/strand_asymmetry/c_elegans.PRJNA13758.WS292.annotations.gff3 .
cp -r -p ../random_1000000_non_mutant_sites.txt .

awk 'NR>1 {print $2"\t"$3}' random_1000000_non_mutant_sites.txt | sort -k1,1 -k2,2n > random_1000000_non_mutant_sites.bed

# extract list of positions:
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

### converting GFF3 to BED format
gff2bed < c_elegans.PRJNA13758.WS292.annotations.gff3 > c_elegans.PRJNA13758.WS292.annotations.gff3.bed


for i in *.final.bed; do
bedtools intersect -a c_elegans.PRJNA13758.WS292.annotations.gff3.bed -b "$i" > "${i/%.final.bed/.gff3.bed}"
done

rm *.tab.txt
rm *.tab.bed
rm *.step1.bed
rm *.step2.bed

# filters lines where the seventh column is WormBase:
for i in *sites.gff3.bed; do
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




### Extract reference alleles

# 1. Make BED file: CHROM start end
awk 'NR>1 {print $2"\t"$3-1"\t"$3}' \
random_1000000_non_mutant_sites.txt | \
sort -k1,1 -k2,2n > positions.bed

# 2. Extract reference base
bedtools getfasta \
-fi /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa \
-bed positions.bed \
-tab > ref_alleles.txt

# 3. Convert to: CHROM POS REF
awk -F'[:\t-]' 'BEGIN{OFS="\t"} {print $1, $3, $4}' \
ref_alleles.txt > nonmutant_with_ref.txt

# 4. Sort by C. elegans chromosome order and position
awk '
BEGIN{
  chr["I"]=1; chr["II"]=2; chr["III"]=3;
  chr["IV"]=4; chr["V"]=5; chr["X"]=6; chr["MtDNA"]=7;
}
{
  print chr[$1]"\t"$0
}' nonmutant_with_ref.txt | \
sort -k1,1n -k3,3n | \
cut -f2- > nonmutant_with_ref.sorted.txt

# 5. Add header
(echo -e "CHROM\tPOS\tREF"; cat nonmutant_with_ref.sorted.txt) \
> nonmutant_with_ref.sorted.header.txt

rm nonmutant_with_ref.sorted.txt nonmutant_with_ref.txt ref_alleles.txt positions.bed


# generate fake ALT allele:
awk '
BEGIN{
    srand(123)
    OFS="\t"
}
NR==1{
    print $0, "ALT"
    next
}
{
    ref=$3

    if(ref=="A"){
        split("C G T", alt, " ")
    } else if(ref=="C"){
        split("A G T", alt, " ")
    } else if(ref=="G"){
        split("A C T", alt, " ")
    } else if(ref=="T"){
        split("A C G", alt, " ")
    } else {
        next
    }

    print $0, alt[int(rand()*3)+1]
}' nonmutant_with_ref.sorted.header.txt > nonmutant_with_fake_mutations.txt

#End 


 
### extracting the variants from mononucleotides run:

# Preparing the data:

cd /blue/baer/mutability/nonmutated/
mkdir mono
cd mono

cp -r -p ../random_1000000_non_mutant_sites.txt .

awk 'NR>1 {print $2"\t"$3}' random_1000000_non_mutant_sites.txt | sort -k1,1 -k2,2n > random_1000000_non_mutant_sites.bed

cp -r -p /blue/baer/share/moeinraja/reference_WS263/mono.5bp.final.bed .

module load bedtools

cp -r -p ../nonmutant_with_fake_mutations.txt .

awk 'NR>1{print $1"\t"$2-1"\t"$2}' nonmutant_with_fake_mutations.txt > nonmutant.bed

bedtools intersect -u -a nonmutant.bed -b mono.5bp.final.bed > mono_hits.bed

awk '
BEGIN{OFS="\t"}
NR==FNR{
    hit[$1":"$3]=1
    next
}
FNR==1{
    print $0, "mono"
    next
}
{
    key=$1":"$2
    print $0, (key in hit ? 1 : 0)
}
' mono_hits.bed nonmutant_with_fake_mutations.txt \
> nonmutant_with_fake_mutations_mono.txt


#End

###coding_context

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
cp -r -p ../mono/nonmutant_with_fake_mutations_mono.txt .


### 1. Convert fake mutation table to VCF

input="nonmutant_with_fake_mutations_mono.txt"
vcf="nonmutant_fake_mutations.vcf"

(
echo "##fileformat=VCFv4.2"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,".",$3,$4,".","PASS","."}' "$input"
) > "$vcf"


### 2. Run snpEff

java -Xmx8g -jar snpEff/snpEff.jar -stats nonmutant_fake_mutations.snpEff.html \
-v WBcel235.99 "$vcf" > nonmutant_fake_mutations.snpEff.vcf


### 3. Extract annotation fields

java -Xmx2g -jar snpEff/SnpSift.jar extractFields nonmutant_fake_mutations.snpEff.vcf \
"CHROM" "POS" "REF" "ALT" \
"ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" \
"ANN[0].GENE" "ANN[0].FEATUREID" \
"ANN[0].HGVS_C" "ANN[0].HGVS_P" "ANN[0].GENEID" \
> nonmutant_fake_mutations.highest_impact_snpEff.tsv

#keep  mono column in the final snpEff table:

awk 'BEGIN{OFS="\t"} NR==FNR{
  if(NR>1) mono[$1":"$2]=$5
  next
}
FNR==1{
  print $0,"mono"
  next
}
{
  print $0,mono[$1":"$2]
}' nonmutant_with_fake_mutations_mono.txt \
nonmutant_fake_mutations.highest_impact_snpEff.tsv \
> nonmutant_fake_mutations.highest_impact_snpEff.with_mono.tsv

#End

