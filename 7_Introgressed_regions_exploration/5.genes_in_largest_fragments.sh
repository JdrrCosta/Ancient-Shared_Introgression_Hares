#overlap largest fragments and genes
#cactus4
bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
-b Cactus4_top_0.01_percent_fragments.tsv | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_cactus4_largest_fragments.txt

#cactus12
bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
-b Cactus12_top_0.01_percent_fragments.tsv | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_cactus12_largest_fragments.txt

#genes in void regions
#sort file
sort -k1,1V -k2,2n cactus4_12_merged.bed \
> cactus4_12_merged.sorted.bed
#get uncovered regions
bedtools complement \
    -i cactus4_12_merged.sorted.bed \
    -g genome.txt > Cactus4_12_uncovered_regions.bed
    
    
#intersect these with regions that were analysed by saguaro
#first get the saguaro analysed regions
grep "cactus" LocalTrees.out \
| awk 'NR > 52 {gsub(":", "", $2); print $2 "\t" $3 "\t" $5}' \
> Saguaro_analysed_regions.bed

#second merge adajacent fragments
bedtools merge -i Saguaro_analysed_regions.bed \
> Saguaro_analysed_regions_merged.bed

#thrid intersect with void regions
bedtools intersect -a Cactus4_12_uncovered_regions.bed \
-b Saguaro_analysed_regions_merged.bed \
> Cactus4_12_void_regions_saguaro.bed 

#size of uncovered regions
awk -F " " '{print $1 "\t" $2 "\t" $3 "\t" $3-$2}' \
Cactus4_12_void_regions_saguaro.bed | sort -nr -k4 


bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
-b Cactus4_12_top_0.01_percent_uncovered_fragments_saguaro.tsv | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_uncovered_largest_fragments_saguaro.txt



#Count number of basepairs in each chromossome of
#all regions analised by saguaro
awk '
{
    bp[$1] += ($3 - $2)
}
END {
    for (chr in bp)
        print chr "\t" bp[chr]
}
' Saguaro_analysed_regions_merged.bed | sort -k2,2nr

#regions void of cactus4 and 12
awk '
{   
    bp[$1] += ($3 - $2)
}
END {
    for (chr in bp)
        print chr "\t" bp[chr]
}
' Cactus4_12_void_regions_saguaro.bed | sort -k2,2nr
