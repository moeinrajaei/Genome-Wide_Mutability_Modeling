#!/bin/sh
salloc -t 2:00:00 --mem=8G

cd /blue/baer/mutability
mkdir coding_context 
cd coding_context

module load java/23

wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip
unzip -o snpEff_latest_core.zip
cd snpEff



# Building database WS295:
vi snpEff.config

# C_elegans genome, version 295
C.elegans.295.genome : C_elegans
# Use mitochondrial codon table for C. elegans mitochondrial DNA
C.elegans.295.MtDNA.codonTable : Invertebrate_Mitochondrial

ESC
:wq
ENTER

mkdir -p data/C.elegans.295
cd data/C.elegans.295/

cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA275000.WS295.genomic.fa .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA275000.WS295.genomic.fa.fai .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA275000.WS295.genomic.dict .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA275000.WS295.annotations.gff3.gz .

mv c_elegans.PRJNA275000.WS295.annotations.gff3.gz genes.gff.gz
mv c_elegans.PRJNA275000.WS295.genomic.fa sequences.fa
mv c_elegans.PRJNA275000.WS295.genomic.fa.fai sequences.fa.fai
mv c_elegans.PRJNA275000.WS295.genomic.dict sequences.dict

cd ../..


# Building database WS292:
vi snpEff.config

# C_elegans genome, version 292
C.elegans.292.genome : C_elegans
# Use mitochondrial codon table for C. elegans mitochondrial DNA
C.elegans.292.MtDNA.codonTable : Invertebrate_Mitochondrial

ESC
:wq
ENTER

mkdir -p data/C.elegans.292
cd data/C.elegans.292/

cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.fa.fai .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/c_elegans.PRJNA13758.WS292.genomic.dict .

mv c_elegans.PRJNA13758.WS292.annotations.gff3.gz genes.gff.gz
mv c_elegans.PRJNA13758.WS292.genomic.fa sequences.fa
mv c_elegans.PRJNA13758.WS292.genomic.fa.fai sequences.fa.fai
mv c_elegans.PRJNA13758.WS292.genomic.dict sequences.dict

cd ../..


# Building database WS263:
vi snpEff.config

# C_elegans genome, version 263
C.elegans.263.genome : C_elegans
# Use mitochondrial codon table for C. elegans mitochondrial DNA
C.elegans.263.MtDNA.codonTable : Invertebrate_Mitochondrial

ESC
:wq
ENTER

mkdir -p data/C.elegans.263
cd data/C.elegans.263/

cp -r -p /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa .
cp -r -p /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.fa.fai .
cp -r -p /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.genomic.dict .
cp -r -p /blue/baer/share/moeinraja/reference_WS263/c_elegans.PRJNA13758.WS263.annotations.gff3.gz .

mv c_elegans.PRJNA13758.WS263.annotations.gff3.gz genes.gff.gz
mv c_elegans.PRJNA13758.WS263.genomic.fa sequences.fa
mv c_elegans.PRJNA13758.WS263.genomic.fa.fai sequences.fa.fai
mv c_elegans.PRJNA13758.WS263.genomic.dict sequences.dict

cd ../..


# Building a database CB4856:
vi snpEff.config

# C_elegans genome, version CB4856
C.elegans.CB4856.genome : C_elegans
# Use mitochondrial codon table for C. elegans mitochondrial DNA
C.elegans.CB4856.MtDNA.codonTable : Invertebrate_Mitochondrial

ESC
:wq
ENTER

mkdir -p data/C.elegans.CB4856
cd data/C.elegans.CB4856/

cp -r -p /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.fa.fai .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/CB4856_genomic.dict .
cp -r -p /orange/baer/Raw_vcf_data/Index_files/CB4856.genomic.gff3 .

gzip CB4856.genomic.gff3

mv CB4856.genomic.gff3.gz genes.gff.gz
mv CB4856_genomic.fa sequences.fa
mv CB4856_genomic.fa.fai sequences.fa.fai
mv CB4856_genomic.dict sequences.dict

cd ../..


java -jar snpEff.jar build -gff3 -v C.elegans.CB4856 -noCheckCds -noCheckProtein
java -jar snpEff.jar build -gff3 -v C.elegans.263 -noCheckCds -noCheckProtein
java -jar snpEff.jar build -gff3 -v C.elegans.292 -noCheckCds -noCheckProtein
java -jar snpEff.jar build -gff3 -v C.elegans.295 -noCheckCds -noCheckProtein

#annotations

for i in ../*.vcf; do
java -Xmx8g -jar snpEff.jar -stats "${i/%.vcf/.ex1.html}" -v WBcel235.99 "$i" > "${i/%.vcf/.ann.vcf}"
done


java -Xmx8g -jar snpEff.jar -stats C.snp.html -v C.elegans.CB4856 ../C.snp.vcf > C.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats C.indel.html -v C.elegans.CB4856 ../C.indel.vcf > C.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats PB306.snp.html -v C.elegans.263 ../PB306.snp.vcf > PB306.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats PB306.indel.html -v C.elegans.263 ../PB306.indel.vcf > PB306.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats N2.snp.html -v C.elegans.263 ../N2.snp.vcf > N2.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats N2.indel.html -v C.elegans.263 ../N2.indel.vcf > N2.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats mev.snp.html -v C.elegans.263 ../mev.snp.vcf > mev.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats mev.indel.html -v C.elegans.263 ../mev.indel.vcf > mev.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats E.snp.html -v C.elegans.292 ../E.snp.vcf > E.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats E.indel.html -v C.elegans.292 ../E.indel.vcf > E.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats G.snp.html -v C.elegans.292 ../G.snp.vcf > G.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats G.indel.html -v C.elegans.292 ../G.indel.vcf > G.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats J.snp.html -v C.elegans.295 ../J.snp.vcf > J.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats J.indel.html -v C.elegans.295 ../J.indel.vcf > J.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats M.snp.html -v C.elegans.295 ../M.snp.vcf > M.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats M.indel.html -v C.elegans.295 ../M.indel.vcf > M.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats XZ7.snp.html -v C.elegans.295 ../XZ7.snp.vcf > XZ7.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats XZ7.indel.html -v C.elegans.295 ../XZ7.indel.vcf > XZ7.indel.ann.vcf

java -Xmx8g -jar snpEff.jar -stats XZ8.snp.html -v C.elegans.295 ../XZ8.snp.vcf > XZ8.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats XZ8.indel.html -v C.elegans.295 ../XZ8.indel.vcf > XZ8.indel.ann.vcf



for i in *.ann.vcf; do
java -Xmx2g -jar SnpSift.jar extractFields -s "," -e "." "$i" CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" > "${i/%.ann.vcf/.ann.txt}"
done

#End


WBcel235.99

java -Xmx8g -jar snpEff.jar -stats N2.snp.html -v C.elegans.263 ../N2.snp.vcf > N2.snp.ann.vcf
java -Xmx8g -jar snpEff.jar -stats N2.snp.html -v  -o gatk C.elegans.263 ../N2.snp.vcf > N2.snp.ann.best.vcf


java -Xmx2g -jar SnpSift.jar extractFields -s "," -e "." N2.snp.ann.best.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].HGVS_P" > N2.snp.ann.best.ann.txt

# how to extract the one with highest IMPACT:



java -Xmx8g -jar snpEff.jar -stats J.snp.html -o gatk -v WBcel235.99 ../J.snp.vcf > J.snp.ann.vcf



java -Xmx2g -jar SnpSift.jar extractFields -s "," -e "." J.snp.ann.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].HGVS_P" > J.snp.ann.txt

java -Xmx8g -jar snpEff.jar -stats PB306.snp.html -v WBcel235.75 ../PB306.snp.vcf > PB306.snp.WBcel235.75.ann.vcf
java -Xmx2g -jar SnpSift.jar extractFields -s "," -e "." PB306.snp.WBcel235.99.ann.vcf CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].HGVS_P" > PB306.snp.WBcel235.75.ann.txt




mkdir vep
cd vep

module load vep/113

vep -i J.snp.vcf -o J.vep.txt --vcf --species caenorhabditis_elegans --cache  --offline

vep -i mev.snp.vcf.gz -o output.vep.vcf --vcf --fasta sequences.fa --gff genes.gff --species caenorhabditis_elegans --offline


vep -i mev.snp.vcf.gz -o output.vep.vcf --vcf --fasta sequences.fa --gff genes.gff --species caenorhabditis_elegans --offline --force_overwrite --fork 4  --everything
