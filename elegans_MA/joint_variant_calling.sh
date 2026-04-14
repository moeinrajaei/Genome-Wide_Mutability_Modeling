#!/bin/sh
#SBATCH --job-name=s9_c_dummy_snp # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=120gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=s9_c_dummy_snp_%j.out # Standard output and error log
#SBATCH --account=juannanzhou --qos=juannanzhou-b
cd /orange/juannanzhou/Rifat_CB_raw/C/dummy_snp
### step 9: Call variants jointly using GenotypeGVCFs in GATK:
module load gatk/
gatk --java-options "-Xmx60g" GenotypeGVCFs -R elegans_snp_dummy.fa -V gendb://C_dummy_snp -O C_dummy_snp.bp.GVCFs.vcf
