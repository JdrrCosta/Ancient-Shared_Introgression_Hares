###apply DI and autossome filter
#get file wih sites and polarities
awk -F " " '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $12 "\t" $13}' \
../No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output.bed \
> No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI.txt

#filter
awk -F " " '$6 >= -18 && $6 <= -12 {print $0}' \
No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI.txt \
> No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI_18_to_12.txt

#remove autopomorphies of Lgra
bedtools subtract \
-a No_lame_5_indv_all_chr_masked_leur_diem.vcf.gz.diem_output_snp_polarity_DI_18_to_12.txt \
-b Autopomotphic_Lgra_sites.bed \
> 5species_masked_Leur_no_auto_DI_18_to_12.bed


####calcukate the amount of sites with high DI
#in cactus4
bedtools intersect -a ../5species_masked_Leur_no_auto_DI_18_to_12.bed \
-b cactus_iter12_cactus4_windows.bed -u | wc -l \
> cactus4_empirical_sites_with_DI_18_to_12.txt


#for replicates
for i in {1..2000}; do
    bedtools intersect \
        -a ../5species_masked_Leur_no_auto_DI_18_to_12.bed \
        -b circular_boostrap/bootstrap_cactus4_rep_${i}.bed \
        -u | wc -l
done > cactus4_boostrap_sites_with_DI_18_to_12.txt


#in cactus12
bedtools intersect -a ../5species_masked_Leur_no_auto_DI_18_to_12.bed \
-b saguaro_iter12_cactus12_windows.bed -u | wc -l \
> cactus12_empirical_sites_with_DI_18_to_12.txt



#for replicates
for i in {1..2000}; do
    bedtools intersect \
        -a ../5species_masked_Leur_no_auto_DI_18_to_12.bed \
        -b circular_boostrap/bootstrap_rep_${i}.bed \
        -u | wc -l
done > cactus12_boostrap_sites_with_DI_18_to_12.txt
