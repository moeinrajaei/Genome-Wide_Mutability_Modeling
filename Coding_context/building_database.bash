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


#End