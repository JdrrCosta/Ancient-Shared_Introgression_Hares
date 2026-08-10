#Depth per sample in the filtered vcf
vcftools --gzvcf Hares_filt.vcf.gz --depth --out Samples

#genotype quality per sample
bcftools query -f '[%SAMPLE\t%GQ\n]' Hares_filt.vcf.gz \
| awk '$2 != "." {sum[$1] += $2; n[$1]++} END {for (s in sum) print s, sum[s]/n[s]}' \
> GQ_per_sample.txt

#NgsParalog
#transform the bam files into mpileup format. The -q and -Q flags serve to not delete any reads based on coverage the only filters were to 
#remove dulicates and umppaed reads 
#this will output the likelihood ratio of mismapping reads covering each site.
#Collumn 5 of the output will be the mismapping Lrs which are asymptotically distributed as a 50/50 mixture of chi-suqare(1 d.f) and 
#chi-square (0 d.f) under the null of no mapping problems 
# I am using only the mapping quality filter because i did not use any minind or mincov filters for the S
seq 1 23 | xargs -n 1 -P 3 sh -c 'samtools mpileup -b Lcas_bam.filelist -q 0 -Q 0 --ff UNMAP,DUP -r chr$0 | \
~/A_Lcas_Lcor_paper/software/ngsParalog/ngsParalog calcLR -infile - -outfile Lcas_ngsparalog_chr$0.lr -minQ 30 -minind 11 -mincov 1'

#get coverage per site file
samtools depth -f Lcas_bam.filelist -q 30 | awk '{
    count = 0;
    sum = 0;
    for (i = 3; i <= NF; i++) {
        if ($i >= 1) { count++ }  # count samples with coverage >=1
        sum += $i;
    }
    if (count >= 11) {           # require at least 11 individuals
        avg = sum / (NF - 2);
        print $1, $2, avg;
    }
}' > Lcas_ngsparalog_depth.depth

#separate depth file by chr
for i in {1..23}
do 
  awk -v chr="chr$i" '$1 == chr' Lcas_ngsparalog_depth.depth > Lcas_ngsparalog_depth_chr$i.depth
done

#new filter for depth file so it has the same excat sites as the lr file
for i in {1..23}
do
awk 'NR==FNR {a[$2]; next} ($2 in a)' Lcas_ngsparalog_chr$i.lr Lcas_ngsparalog_depth_chr$i.depth > Lcas_ngsparalog_depth_mutual_chr$i.depth
awk 'NR==FNR {a[$2]; next} ($2 in a)' Lcas_ngsparalog_depth_chr$i.depth Lcas_ngsparalog_chr$i.lr > Lcas_ngsparalog_mutual_chr$i.lr
done

#run dupHMM per chr
seq 1 23 | xargs -n 1 -P 5 sh -c 'Rscript dupHMM.R --lrfile Lcas_ngsparalog_mutual_chr$0.lr --outfile Lcas_ngsparlog_chr$0.regions \
--covfile Lcas_ngsparalog_depth_mutual_chr$0.depth --lrquantile 0.98'

#Compare circular boostrap values to empirical values
#Empirical value of percentage of overlap between cactus regions and poorly mapped regions
bedtools intersect -a cactus_iter12_cactus4_windows.bed \
-b Lcor_ngsparlog_all.regions.rf -wo \
| awk '{sum += $NF} END {print sum}'

#calculte percentage
echo "scale=2; (10200997 / 418006609) * 100" 

#Boostrap values of percentage of overlap between cactus regions and poorly mapped regions
for i in {1..2000}; do   bedtools intersect \
-a ../circular_boostrap/bootstrap_cactus4_rep_$i.bed \
-b Lcor_ngsparlog_all.regions.rf -wo \
| awk '{sum += $NF} END {print sum}' \
>> Number_overlaps_cactus4_ngsParalog_Lcor.txt; done

#calculte percentage
awk '{print ($1 / 418006609 * 100)}' Number_overlaps_cactus4_ngsParalog_Lcor.txt \
> Percentage_overlaps_cactus4_bootstrap_ngsParalog_Lcor.txt

