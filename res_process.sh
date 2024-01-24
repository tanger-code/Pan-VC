#!/bin/bash

# ********************************************************************
bcftools=/usr/local/bin/bcftools
bedtools=/data/software/bedtools2/bin/bedtools
python=/data/home/tangen/.conda/envs/tee01/bin/python
fa_file=/data/home/tangen/data/HPRC_with_HG002/HG002_linear_fasta/GCA_chr15.fa
# ********************************************************************

cp ../vcf_rmDup* ./
mv sample_all.txt sample_all.vcf
sed -i 's#1\/2#0\/1#g' sample_all.vcf
$python vcf_rmDup_1.py
bgzip -f sample_all_rmDup.vcf
$bcftools norm -f $fa_file -m -both sample_all_rmDup.vcf.gz -O z -o sample_all_rmDup_norm.vcf.gz -cw
$bcftools sort sample_all_rmDup_norm.vcf.gz -O z -o sample_all_rmDup_norm_sort.vcf.gz
$bcftools index -t sample_all_rmDup_norm_sort.vcf.gz

$bcftools view -v indels sample_all_rmDup_norm_sort.vcf.gz -o sample_indel.vcf
$python vcf_rmDup_2.py
bgzip -f sample_indel_rmDup.vcf
mv sample_indel_rmDup.vcf.gz sample_indels.vcf.gz
$bcftools index -t sample_indels.vcf.gz

$bcftools view -v snps sample_all_rmDup_norm_sort.vcf.gz -o z_snps.vcf
$bcftools view -v mnps sample_all_rmDup_norm_sort.vcf.gz -o z_mnps.vcf
$bcftools norm -a z_mnps.vcf -o z_mnp2snp.vcf
$bcftools concat z_snps.vcf z_mnp2snp.vcf -o z_total_snps.vcf
$python vcf_rmDup_snp.py
bgzip -f sample_snps.vcf
$bcftools sort sample_snps.vcf.gz -O z -o sample_snps.vcf.gz
$bcftools index -t sample_snps.vcf.gz

rm -f sample_all_rmDup.vcf.gz
rm -f sample_all_rmDup_norm.vcf.gz
rm -f sample_all_rmDup_norm_sort.vcf.gz
rm -f sample_all_rmDup_norm_sort.vcf.gz.tbi
rm -f sample_indel.vcf
rm -f z_mnp2snp.vcf
rm -f z_mnps.vcf
rm -f z_snps.vcf
rm -f z_total_snps.vcf
rm -f vcf_rmDup*


# bcftools view -v snps sample_all_rmDup_norm_sort.vcf.gz -o z_snps.vcf
# bcftools view -v mnps sample_all_rmDup_norm_sort.vcf.gz -o z_mnps.vcf
# sed -i 's#1\/2#0\/1#g' z_mnps.vcf
# bcftools norm -a z_mnps.vcf -o z_mnp2snp.vcf
# bcftools concat z_snps.vcf z_mnp2snp.vcf -o z_total_snps.vcf
# python vcf_rmDup_snp.py
# bgzip -f sample_snp.vcf
