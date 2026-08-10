#divide the vcf according to desired cactus
bcftools filter -R cactus_iter12_cactus4_windows.bed -Oz \
-o Saguaro_vcf_iter12_cactus4_regions.vcf.gz \
High_cov_mono_filt.vcf.gz --threads 60

#parse vcf
python parseVCF.py \
-i Saguaro_vcf_iter12_cactus4_regions_mono.vcf.gz \
| gzip > Saguaro_vcf_iter12_cactus4_regions.geno.gz


#divide vcf in boostrap windows
seq 1 100 | xargs -n 1 -P 10 sh -c 'bcftools filter \
-R ../circular_boostrap/bootstrap_cactus4_rep_$0.bed -Oz \
-o /mnt/jcosta/A_Lcas_Lcor/cactus4_boostrap/Circular_bootstrap_$0_cactus4_regions.vcf.gz \
High_cov_mono_filt.vcf.gz --threads 10'

#parse vcf
seq 1 100 | xargs -n 1 -P 10 sh -c 'python parseVCF.py \
-i /mnt/jcosta/A_Lcas_Lcor/cactus4_boostrap/Circular_bootstrap_$0_cactus4_regions.vcf.gz \
| gzip > Circular_bootstrap_$0_cactus4_regions.geno.gz'


#Lcas and Ltim
#cactus4
python dxy_calc_homemade.py --geno Saguaro_vcf_iter12_cactus4_regions.geno.gz --popfile pop.info --pop1 Lcas \
--pop2 Ltim > Saguaro_iter12_cactus4_Lcas_Ltim.dxy

#circular boostarp
seq 1 100 | xargs -n 1 -P 10 sh -c 'python dxy_calc_homemade.py --geno Circular_bootstrap_$0_cactus4_regions.geno.gz \
--popfile pop.info --pop1 Lcas --pop2 Ltim >> Circular_bootstrap_cactus4_Lcas_Ltim.dxy'

#Lcor and Ltim
#cactus4
python dxy_calc_homemade.py \
--geno Saguaro_vcf_iter12_cactus4_regions.geno.gz \
--popfile pop.info --pop1 Lcor --pop2 Ltim \
> Saguaro_iter12_cactus4_Lcor_Ltim.dxy

#circular boostarp
seq 1 100 | xargs -n 1 -P 30 sh -c 'python dxy_calc_homemade.py \
--geno Circular_bootstrap_$0_cactus4_regions.geno.gz \
--popfile pop.info --pop1 Lcor --pop2 Ltim \
>> Circular_bootstrap_cactus4_Lcor_Ltim.dxy'
