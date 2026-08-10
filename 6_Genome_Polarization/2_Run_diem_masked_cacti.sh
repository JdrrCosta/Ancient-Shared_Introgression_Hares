###MASKED CACTUS REGIONS##############

#get cactus4 and 12 sites
cat cactus_iter12_cactus4_windows.bed \
saguaro_iter12_cactus12_windows.bed \
| sort -k1,1 -k2,2n | bedtools merge > sites_masked.bed

#run Diem with masked cacti
python Run_diem_extra_plots_with_Masking_and_cactus_masking.py \
No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz \
--cores 40

#run Diem for diffrent DI tresholds 
python Run_diem_DI_threholds.py \
../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz \
--cores 40

#calculate HI for each ltim per snp
#no threshold
python calculate_average_HI_per_indv.py \
-i 5species_masked_cactus_no_auto_polarity_DI.txt \
-o HI_masked_cacti_no_auto_no_threshold.txt
#-40
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_cacti_mask_polarity_DI_40.txt \
-o HI_masked_cacti_no_auto_DI_40.txt
#-30
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_cacti_mask_polarity_DI_30.txt \
-o HI_masked_cacti_no_auto_DI_30.txt
#-20
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_cacti_mask_polarity_DI_20.txt \
-o HI_masked_cacti_no_auto_DI_20.txt
#-15
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_cacti_mask_polarity_DI_15.txt \
-o HI_masked_cacti_no_auto_DI_15.txt
#-10
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_cacti_mask_polarity_DI_10.txt \
-o HI_masked_cacti_no_auto_DI_10.txt


