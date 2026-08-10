#get snps used by diem
awk -F " " 'NR > 2 {print $1 "\t" $2 "\t" $3}' ../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output.bed > Snps_used_by_diem.bed
#get snps used by diem within cactus regions
bedtools intersect -a Snps_used_by_diem.bed -b cactus4_12_merged.bed -u > Snps_used_by_diem_in_cactus_regions.bed

#create a vcf with these snps
bcftools index ../../include_lgra_all_indvs/Include_Lgra_all_chr_all_sampled_diem.vcf.gz

#all indvs
bcftools view -R Snps_used_by_diem_in_cactus_regions.bed -S lcaslcor_samples.list -Oz -m2 -M2 -v snps -o All_indvs_LcasLcorLtimLgra_snps_used_by_diem_in_cactus_biallelic.vcf.gz \
../../include_lgra_all_indvs/Include_Lgra_all_chr_all_sampled_diem.vcf.gz --threads 60

#only high coverage
bcftools view -R Snps_used_by_diem_in_cactus_regions.bed -S high_cov_samples.list -Oz -m2 -M2 -v snps \
-o High_cov_indvs_LcasLcorLtimLgra_snps_used_by_diem_in_cactus_biallelic.vcf.gz \
../../include_lgra_all_indvs/Include_Lgra_all_chr_all_sampled_diem.vcf.gz --threads 60

#run r code
Rscript snp_heatmap.r \
All_indvs_LcasLcor_snps_used_by_diem_in_cactus_biallelic.vcf.gz \
All_indvs_LcasLcor_heatmap.pdf \
sample_order.txt


Rscript snp_heatmap.r \
High_cov_indvs_LcasLcor_snps_used_by_diem_in_cactus_biallelic.vcf.gz \
High_cov_indvs_LcasLcor_heatmap.pdf \
sample_order_high_cov.txt
