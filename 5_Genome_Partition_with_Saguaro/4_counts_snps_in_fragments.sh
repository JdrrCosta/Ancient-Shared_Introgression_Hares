#get vcf with cactus4 regions
bcftools view -R cactus_iter12_cactus4_windows.bed -Oz -o Saguaro_vcf_iter12_cactus4_regions.vcf.gz \
Hares_saguaro.vcf.gz --threads 60

#get snps
bcftools view -v snps Saguaro_vcf_iter12_cactus4_regions.vcf.gz -Oz \
-o Saguaro_vcf_iter12_cactus4_regions_snps.vcf.gz --threads 60

#count snps
bedtools intersect -a cactus_iter12_cactus4_windows.bed -b Saguaro_vcf_iter12_cactus4_regions_snps.vcf.gz \
-c > Snp_counts_per_fragment_cactus4.bed
