
salloc -t 2:00:00 --mem=8G

module load BEDOPS
module load BCFtools
module load BEDTools
module load VCFtools

cd /vast/palmer/scratch/townsend/mr2683/project/test_7_20_2025

### converting GFF3 to BED format
gff2bed < c_elegans.PRJNA13758.WS263.annotations.gff3 > c_elegans.PRJNA13758.WS263.annotations.gff3.bed

mkdir strand_asymmetry
cd strand_asymmetry

### extract list of positions:
cp -r -p ../*.vcf .
 
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

# + is forward
# - is reverse

# extract _forward_bed:
for i in *.wormbase.unique.bed; do
  awk '$4=="+" {print $1, $3}' "$i" > "${i/%.wormbase.unique.bed/.forward.bed}"
done

#extract _reverse.bed
for i in *.wormbase.unique.bed; do
  awk '$4=="-" {print $1, $3}' "$i" > "${i/%.wormbase.unique.bed/.reverse.bed}"
done

vcftools --vcf rare.snp.vcf --recode --out rare_reverse --positions rare.snp.reverse.bed --recode-INFO-all
vcftools --vcf mev.snp.vcf --recode --out mev_reverse --positions mev.snp.reverse.bed --recode-INFO-all
vcftools --vcf N2.snp.vcf --recode --out N2_reverse --positions N2.snp.reverse.bed --recode-INFO-all
vcftools --vcf PB306.snp.vcf --recode --out PB306_reverse --positions PB306.snp.reverse.bed --recode-INFO-all

vcftools --vcf rare.snp.vcf --recode --out rare_reverse --positions rare.snp.reverse.bed --recode-INFO-all
vcftools --vcf mev.snp.vcf --recode --out mev_forward --positions mev.snp.forward.bed --recode-INFO-all
vcftools --vcf N2.snp.vcf --recode --out N2_forward --positions N2.snp.forward.bed --recode-INFO-all
vcftools --vcf PB306.snp.vcf --recode --out PB306_forward --positions PB306.snp.forward.bed --recode-INFO-all

rm *.unique.bed 
rm *.final.bed 
rm *.wormbase.bed
rm *.gff3.bed

### base substitution type analysis for each of *.recode.vcf files(line 70-160):
mkdir whole
cd whole
cp -r -p ../*.recode.vcf .

## rare
mkdir rare
cd rare 
cp -r -p ../rare_reverse.recode.vcf .
cp -r -p ../rare_forward.recode.vcf .

#rare_reverse 
mkdir reverse
cd reverse 

for sample in `bcftools query -l ../rare_reverse.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../rare_reverse.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > rare_reverse_ts_tv.txt

paste *.type.txt > rare_reverse_BS_type.txt

ls *.vcf > rare.reverse.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#rare_forward
cd ..
mkdir forward
cd forward


for sample in `bcftools query -l ../rare_forward.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../rare_forward.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > rare_forward_ts_tv.txt

paste *.type.txt > rare_forward_BS_type.txt

ls *.vcf > rare.forward.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

##mev
cd ../..

mkdir mev
cd mev 
cp -r -p ../mev_reverse.recode.vcf .
cp -r -p ../mev_forward.recode.vcf .

#mev_reverse
mkdir reverse
cd reverse

for sample in `bcftools query -l ../mev_reverse.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../mev_reverse.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > mev_reverse_ts_tv.txt

paste *.type.txt > mev_reverse_BS_type.txt

ls *.vcf > mev.reverse.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#mev_forward:
cd ..
mkdir forward
cd forward


for sample in `bcftools query -l ../mev_forward.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../mev_forward.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > mev_forward_ts_tv.txt

paste *.type.txt > mev_forward_BS_type.txt

ls *.vcf > mev.forward.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#N2
cd ../..

mkdir N2
cd N2 
cp -r -p ../N2_reverse.recode.vcf .
cp -r -p ../N2_forward.recode.vcf .

#N2_reverse
mkdir reverse
cd reverse

for sample in `bcftools query -l ../N2_reverse.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../N2_reverse.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > N2_reverse_ts_tv.txt

paste *.type.txt > N2_reverse_BS_type.txt

ls *.vcf > N2.reverse.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#N2_forward:
cd ..
mkdir forward
cd forward

for sample in `bcftools query -l ../N2_forward.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../N2_forward.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > N2_forward_ts_tv.txt

paste *.type.txt > N2_forward_BS_type.txt

ls *.vcf > N2.forward.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#PB306
cd ../..

mkdir PB306
cd PB306 
cp -r -p ../PB306_reverse.recode.vcf .
cp -r -p ../PB306_forward.recode.vcf .

#PB306_reverse
mkdir reverse
cd reverse

for sample in `bcftools query -l ../PB306_reverse.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../PB306_reverse.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > PB306_reverse_ts_tv.txt

paste *.type.txt > PB306_reverse_BS_type.txt

ls *.vcf > PB306.reverse.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#PB306_forward:
cd ..
mkdir forward
cd forward


for sample in `bcftools query -l ../PB306_forward.recode.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../PB306_forward.recode.vcf
done

for i in *.vcf; do
bcftools stats "$i" > "${i/%.vcf/.stats.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^TSTV | cut -f 3- > "${i/%.stats.txt/.ts_tv.txt}"
done

for i in *.stats.txt; do
cat "$i" | grep ^ST | cut -f 3- > "${i/%.stats.txt/.type.txt}"
done

cut -f -3 *.ts_tv.txt > PB306_forward_ts_tv.txt

paste *.type.txt > PB306_forward_BS_type.txt

ls *.vcf > PB306.forward.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

###MR_Here
# then the files were copied to an excel file named "Base-sub spectra.xlsx"

