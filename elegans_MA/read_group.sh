#! /bin/sh
#SBATCH --job-name=s5j # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=6 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=480:00:00 # Time limit hrs:min:sec
#SBATCH --output=s5j_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer 
cd /orange/juannanzhou/Rifat_CB_raw/J
### step 5: Add read group information using picard:
module load picard/2.25.5/
for i in *.PP.bam; do
x=`basename $i .PP.bam`
picard AddOrReplaceReadGroups \
-I "$i" \
-O "${i/%.PP.bam/.PP.RG.bam}" \
-RGID 300_MA \
-RGLB "$x" \
-RGPL illumina \
-RGPU 300_MA \
-RGSM "$x"
done
