cd /vast/palmer/scratch/townsend/mr2683/project/test_7_20_2025


mkdir GR_data
cd GR_Data

cp -r -p ../*.vcf .

## rare
mkdir rare
cd rare 
cp -r -p ../rare.snp.vcf .

for sample in `bcftools query -l ../rare.snp.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../rare.snp.vcf
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

cut -f -3 *.ts_tv.txt > rare_ts_tv.txt

paste *.type.txt > rare_BS_type.txt

ls *.vcf > rare.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

##mev
cd ..

mkdir mev
cd mev 
cp -r -p ../mev.snp.vcf .

for sample in `bcftools query -l ../mev.snp.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../mev.snp.vcf
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

cut -f -3 *.ts_tv.txt > mev_ts_tv.txt

paste *.type.txt > mev_BS_type.txt

ls *.vcf > mev.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#N2
cd ..

mkdir N2
cd N2 
cp -r -p ../N2.snp.vcf .

for sample in `bcftools query -l ../N2.snp.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../N2.snp.vcf
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

cut -f -3 *.ts_tv.txt > N2_ts_tv.txt

paste *.type.txt > N2_BS_type.txt

ls *.vcf > N2.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#PB306
cd ..

mkdir PB306
cd PB306 
cp -r -p ../PB306.snp.vcf .

for sample in `bcftools query -l ../PB306.snp.vcf`; do
    bcftools view -c1 -Ov -s $sample -o $sample.vcf ../PB306.snp.vcf
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

cut -f -3 *.ts_tv.txt > PB306_ts_tv.txt

paste *.type.txt > PB306_BS_type.txt

ls *.vcf > PB306.txt

rm *.type.txt
rm *.ts_tv.txt
rm *.stats.txt

#End

#make the graph and check with GR_paper.