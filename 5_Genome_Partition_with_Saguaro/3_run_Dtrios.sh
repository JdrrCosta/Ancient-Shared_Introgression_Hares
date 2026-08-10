#Filter vcf to get only numerical chr and high-cov samples
bcftools view -S samples.list -Oz -o Hares_high_cov_filt_Dsuit.vcf.gz \
Hares_filt.vcf.gz --threads 30

#run Dsuite Dtrios
~/A_Lcas_Lcor_paper/software/Dsuite/Build/Dsuite Dtrios \
--ABBAclustering Hares_high_cov_filt_Dsuit.vcf.gz SETS_trios.txt

