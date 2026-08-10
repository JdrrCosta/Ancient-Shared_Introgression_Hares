#extract bed file of repeats from the output of repreat masker
awk -F '\t' '$5 ~ /^chr/ {print $5 "\t" $6 "\t" $7 "\t" $10}' \
Repeat_masker_Leur_out.tab > Repeat_masker_Leur_out.bed

#keep only autossomes
grep -v -E "chrY|chrX" Repeat_masker_Leur_out.bed \
> Repeat_masker_Leur_out_num.bed

#sort and merge overlapping regions
sort -k1,1V -k2,2n Repeat_masker_Leur_out_num.bed \
| bedtools merge > Repeat_masker_Leur_out_num_merged_sorted.bed

#overlap with file of uncovered regions
bedtools intersect -a Cactus4_12_uncovered_regions.sorted.bed \
-b Repeat_masker_Leur_out_num_merged_sorted.bed -wo \
| awk '{sum+=$NF} END{print sum}'

#bp in void regions
awk '{sum += $3 - $2} END {print sum}' \
Cactus4_12_uncovered_regions.sorted.bed

#Do the same for cactus regions
#overlap with file of cactus regions
bedtools intersect -a cactus4_12_merged.sorted.bed \
-b Repeat_masker_Leur_out_num_merged_sorted.bed -wo \
| awk '{sum+=$NF} END{print sum}'

#bp in cactus regions
awk '{sum += $3 - $2} END {print sum}' \
cactus4_12_merged.sorted.bed

#Same for the top 0.1% of void regions
bedtools intersect -a Cactus4_12_top_0.01_percent_uncovered_fragments_saguaro.tsv \
-b Repeat_masker_Leur_out_num_merged_sorted.bed -wo \
| awk '{sum+=$NF} END{print sum}'


#bp in top 0.1% of void regions
awk '{sum += $3 - $2} END {print sum}' \
Cactus4_12_top_0.01_percent_uncovered_fragments_saguaro.tsv
