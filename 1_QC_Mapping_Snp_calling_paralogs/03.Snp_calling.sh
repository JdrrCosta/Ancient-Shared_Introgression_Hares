#Snp calling
#list bam files
ls ./5_realign/*.realigned.bam > 6_Ltim_bam_files

#SNP call
bcftools mpileup -a AD,DP,SP -Ou -f ~/Z_refs/GCA_033115175.1_mLepEur2.pri_genomic.fna --bam-list 6_Ltim_bam.files \
--threads 15 | bcftools call -f GQ,GP -m -Ob --threads 15 -o ./6_call/Ltim_call_Leur_ref.bcf.gz

## To prepare our VCF for querying we next index it using tabix
tabix -p bcf Ltim_call_Leur_ref.bcf.gz

## How many unfiltered variants?
bcftools view -H Ltim_call_Leur_ref.bcf.gz | wc -l > Ltim_call_Leur_ref_unfilt.var

#change chr names
bcftools annotate --rename-chrs European_ref_chr.renamemap -Ob -o Ltim_call_Leur_ref_rename.bcf.gz Ltim_call_Leur_ref.bcf.gz --threads 30
