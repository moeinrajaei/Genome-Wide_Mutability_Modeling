#!/bin/sh
#SBATCH --job-name=s7j # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=s7j_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer
cd /orange/juannanzhou/Rifat_CB_raw/J

### step 7: Call variants per-sample using HaplotypeCaller (in BP_RESOLUTION mode):
module load gatk/
module load samtools/
module load picard/
for i in *.markduplicates.bam; do

# define output file name
out="${i/%.markduplicates.bam/.bp.g.vcf.gz}"

# skip if output already exists
 if [ -f "$out" ]; then
        echo "Skipping $i — $out already exists."
        continue
    fi

    echo "Processing $i ..."

java -Xmx6g -jar /apps/gatk/4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller \
-R c_elegans.PRJNA13758.WS292.genomic.fa \
-I "$i" \
-O "$out" \
-ERC BP_RESOLUTION
done
