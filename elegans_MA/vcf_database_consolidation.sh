#!/bin/sh
#SBATCH --job-name=s8_j_dummy_indel # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=m.rifat@ufl.edu # Where to send mail
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --ntasks=1 # Run a single task
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=s8_j_dummy_snp_%j.out # Standard output and error log
#SBATCH --account=baer --qos=baer-b
cd /orange/juannanzhou/Rifat_CB_raw/J/dummy_indel
### step 8: consolidating the vcf.gz files using GenomicsDBImport:
module load gatk/
gatk --java-options "-Xmx4g" GenomicsDBImport \
-V H7GWFDSX7_J301_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J302_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J312_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J314_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J317_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J318_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J322_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J325_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J330_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J332_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J336_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J337_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J339_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J347_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J353_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J356_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J357_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J361_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J366_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J368_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J372_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J373_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J376_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J382_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J390_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J391_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J393_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J394_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J396_POOLRET74.bp.g.vcf.gz \
-V H7GWFDSX7_J397_POOLRET74.bp.g.vcf.gz \
-V J300_POOLRET73_S475_L003.bp.g.vcf.gz \
-V J303_POOLRET72_S208_L002.bp.g.vcf.gz \
-V J304_POOLRET73_S436_L003.bp.g.vcf.gz \
-V J305_POOLRET73_S529_L003.bp.g.vcf.gz \
-V J306_POOLRET73_S418_L003.bp.g.vcf.gz \
-V J307_POOLRET73_S496_L003.bp.g.vcf.gz \
-V J308_POOLRET73_S554_L003.bp.g.vcf.gz \
-V J309_POOLRET73_S440_L003.bp.g.vcf.gz \
-V J310_POOLRET73_S526_L003.bp.g.vcf.gz \
-V J311_POOLRET72_S217_L002.bp.g.vcf.gz \
-V J313_POOLRET73_S557_L003.bp.g.vcf.gz \
-V J315_POOLRET73_S455_L003.bp.g.vcf.gz \
-V J316_POOLRET72_S311_L002.bp.g.vcf.gz \
-V J319_POOLRET72_S272_L002.bp.g.vcf.gz \
-V J320_POOLRET73_S515_L003.bp.g.vcf.gz \
-V J321_POOLRET72_S328_L002.bp.g.vcf.gz \
-V J323_POOLRET73_S550_L003.bp.g.vcf.gz \
-V J324_POOLRET73_S439_L003.bp.g.vcf.gz \
-V J326_POOLRET73_S402_L003.bp.g.vcf.gz \
-V J327_POOLRET73_S400_L003.bp.g.vcf.gz \
-V J328_POOLRET72_S204_L002.bp.g.vcf.gz \
-V J329_POOLRET72_S267_L002.bp.g.vcf.gz \
-V J331_POOLRET72_S375_L002.bp.g.vcf.gz \
-V J333_POOLRET73_S456_L003.bp.g.vcf.gz \
-V J334_POOLRET72_S327_L002.bp.g.vcf.gz \
-V J335_POOLRET72_S262_L002.bp.g.vcf.gz \
-V J338_POOLRET72_S309_L002.bp.g.vcf.gz \
-V J340_POOLRET72_S314_L002.bp.g.vcf.gz \
-V J341_POOLRET72_S268_L002.bp.g.vcf.gz \
-V J342_POOLRET72_S302_L002.bp.g.vcf.gz \
-V J343_POOLRET73_S499_L003.bp.g.vcf.gz \
-V J344_POOLRET72_S249_L002.bp.g.vcf.gz \
-V J345_POOLRET73_S491_L003.bp.g.vcf.gz \
-V J346_POOLRET72_S255_L002.bp.g.vcf.gz \
-V J348_POOLRET73_S416_L003.bp.g.vcf.gz \
-V J349_POOLRET72_S293_L002.bp.g.vcf.gz \
-V J350_POOLRET73_S433_L003.bp.g.vcf.gz \
-V J351_POOLRET73_S486_L003.bp.g.vcf.gz \
-V J352_POOLRET72_S287_L002.bp.g.vcf.gz \
-V J354_POOLRET73_S565_L003.bp.g.vcf.gz \
-V J355_POOLRET73_S443_L003.bp.g.vcf.gz \
-V J358_POOLRET72_S256_L002.bp.g.vcf.gz \
-V J359_POOLRET73_S477_L003.bp.g.vcf.gz \
-V J360_POOLRET73_S466_L003.bp.g.vcf.gz \
-V J362_POOLRET72_S373_L002.bp.g.vcf.gz \
-V J363_POOLRET73_S419_L003.bp.g.vcf.gz \
-V J364_POOLRET73_S520_L003.bp.g.vcf.gz \
-V J365_POOLRET72_S226_L002.bp.g.vcf.gz \
-V J367_POOLRET72_S334_L002.bp.g.vcf.gz \
-V J369_POOLRET72_S347_L002.bp.g.vcf.gz \
-V J370_POOLRET73_S501_L003.bp.g.vcf.gz \
-V J371_POOLRET73_S415_L003.bp.g.vcf.gz \
-V J374_POOLRET73_S556_L003.bp.g.vcf.gz \
-V J375_POOLRET72_S253_L002.bp.g.vcf.gz \
-V J377_POOLRET72_S285_L002.bp.g.vcf.gz \
-V J378_POOLRET72_S279_L002.bp.g.vcf.gz \
-V J379_POOLRET72_S260_L002.bp.g.vcf.gz \
-V J380_POOLRET73_S461_L003.bp.g.vcf.gz \
-V J381_POOLRET72_S238_L002.bp.g.vcf.gz \
-V J383_POOLRET73_S391_L003.bp.g.vcf.gz \
-V J384_POOLRET72_S335_L002.bp.g.vcf.gz \
-V J385_POOLRET72_S235_L002.bp.g.vcf.gz \
-V J386_POOLRET72_S321_L002.bp.g.vcf.gz \
-V J387_POOLRET73_S459_L003.bp.g.vcf.gz \
-V J388_POOLRET73_S417_L003.bp.g.vcf.gz \
-V J389_POOLRET72_S209_L002.bp.g.vcf.gz \
-V J392_POOLRET72_S318_L002.bp.g.vcf.gz \
-V J395_POOLRET73_S429_L003.bp.g.vcf.gz \
-V J398_POOLRET73_S505_L003.bp.g.vcf.gz \
-V J399_POOLRET73_S553_L003.bp.g.vcf.gz \
-V H7GWFDSX7_JU360_BaerAncestor_POOLRET74_a.bp.g.vcf.gz \
-V H7GWFDSX7_JU360_BaerAncestor_POOLRET74_b.bp.g.vcf.gz \
-V H7GWFDSX7_JU360_BaerAncestor_POOLRET74_c.bp.g.vcf.gz \
--genomicsdb-workspace-path J_dummy_indel \
--intervals I --intervals II --intervals III --intervals IV --intervals V --intervals X --intervals MtDNA
