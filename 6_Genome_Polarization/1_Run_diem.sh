#Filter bcf for polarization
bcftools view -S samples.list -Oz \
-o No_lame_5_indv_all_chr_diem.vcf.gz Hares_filt.vcf.gz \
--threads 80

#run Diem
python Run_diem_extra_plots_with_Masking.py \
No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz \
--cores 40

#extract sites that group LcasLcor e Leur
awk -F " " '{print $1,$2,$3,$10,$12}' \
../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output.bed \
> No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_snp_polarity.txt

#get the sites that group lcas ltim lcor 
python get_sites_grouping_lcaslcorltimleur_vs_lgra.py

#get bed file with these sites
awk -F " " 'NR>1 {print $1 "\t" $2 "\t" $3}' filtered_sites_group_LcasLcorLeur.txt > Autopomotphic_Lgra_sites.bed

#get file wih sites and polarities
awk -F " " '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $12 "\t" $13}' \
../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output.bed \
> No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI.txt

#filter
awk -F " " '$6 >= -20 {print $0}' \
No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI.txt \
> No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI_20.txt

#remove autopomorphies of Lgra
bedtools subtract \
-a No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI_20.txt \
-b Autopomotphic_Lgra_sites.bed \
> 5species_masked_Leur_no_auto_DI_20.bed

#calulate averge HI for Ltim population per snp 
python calculate_average_HI_snp.py

#based on that input calculate the average HI for the boostrap regions and cactus regions
python calc_HI_regions.py 

#run for diffrent Di thresholds
python Run_diem_DI_threholds.py \
../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz \
--cores 40

#calculate HI for each ltim per snp
#no threshold
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI.txt \
-o HI_masked_leur_no_auto_no_threshold.txt
#-40
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI_40.txt \
-o HI_masked_leur_no_auto_DI_40.txt
#-30
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI_30.txt \
-o HI_masked_leur_no_auto_DI_30.txt
#-20
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI_20.txt \
-o HI_masked_leur_no_auto_DI_20.txt
#-15
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI_15.txt \
-o HI_masked_leur_no_auto_DI_15.txt
#-10
python calculate_average_HI_per_indv.py \
-i 5species_masked_Leur_no_auto_polarity_DI_10.txt \
-o HI_masked_leur_no_auto_DI_10.txt
