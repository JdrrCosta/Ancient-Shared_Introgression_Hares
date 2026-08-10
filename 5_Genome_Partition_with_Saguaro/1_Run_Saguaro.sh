#filter vcf to remove monomorphic sites.
#code to create the vcf with monos is on the dxy page
bcftools view -O z -o Hares_filt.vcf.gz -i 'MAF>0.009' -m2 Hares_mono_filt.vcf.gz --threads 60

#restrict analysis to just some individuals
bcftools view -Oz -o Hares_saguaro.vcf.gz \
-s Lcas1,Lcas2,Lcas3,Lcas4,Lcas5,Lcor6,Lcor7,Lcor8,Lcor9,Lcor10,Lgra1,Lgra10,Ltim1,Ltim2,Ltim4,Leur3,Leur7,Leur10,Lame1 \
Hares_filt.vcf.gz --threads 60

#run script to turn into binary
./VCF2HMMFeature -i Hares_saguaro.vcf -o Saguaro_input

#run saguaro
~/A_Lcas_Lcor_paper/software/saguarogw-code/Saguaro -f Saguaro_input -o Saguaro_output_iter_12 -iter 12'
