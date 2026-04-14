#!/bin/bash

module load bcftools

# Loop through all .vcf files in the current directory

for vcf_file in *filtered.biallelic.vcf

do

    # Define the output file name by replacing the original filename's extension

    output_file="${vcf_file%filtered.biallelic.vcf}mutation.vcf"

    

    # Run bcftools command on each file

    bcftools view -i 'COUNT(GT="alt")=1 && COUNT(GT="het")=0 && COUNT(GT="./.")<=10' "$vcf_file" -Oz -o "$output_file"


    echo "Processed $vcf_file -> $output_file"

done


