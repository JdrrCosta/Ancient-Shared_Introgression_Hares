#start position +-10bp file
awk '{start = ($2 - 10 < 0 ? 0 : $2 - 10); print $1"\t"start"\t"$2}' \
cactus_iter12_cactus4_windows.bed \
> cactus4_strPos_minus_10bp.bed

awk -F " " '{print $1"\t"($2)"\t"($2 + 11)}' \
cactus_iter12_cactus4_windows.bed \
> cactus4_strPos_plus_10bp.bed

#end position +-10bb file
awk '{start = ($3 - 10 < 0 ? 0 : $3 - 10); print $1"\t"start"\t"$3}' \
cactus_iter12_cactus4_windows.bed \
> cactus4_endPos_minus_10bp.bed

awk -F " " '{print $1"\t"($3)"\t"($3 + 10)}' \
cactus_iter12_cactus4_windows.bed \
> cactus4_endPos_plus_10bp.bed

#merge
cat cactus4_strPos_minus_10bp.bed cactus4_endPos_plus_10bp.bed \
| sort -k1,1 -k2,2n > Negative_cactus4_10bp.bed

cat cactus4_strPos_plus_10bp.bed cactus4_endPos_minus_10bp.bed \
| sort -k1,1 -k2,2n > Positive_cactus4_10bp.bed


#get average polarization per fragment
sort -k1,1V -k2,2n Positive_cactus4_10bp.bed \
> Positive_cactus4_10bp.sort.bed

sort -k1,1V -k2,2n Negative_cactus4_10bp.bed \
> Negative_cactus4_10bp.sort.bed


#overlap with HI file
#negative file
bedtools intersect \
-a Ltim_averages_masked_leur_no_auto_DI_20.txt \
-b Negative_cactus4_10bp.bed -wb | \
awk -F " " '{print $1,$2,$3,$4, ($2-$7)}' \
> Cactus4_Negative_positions_cline_HI_10bp_no_auto_DI_20.txt

#positive file
bedtools intersect \
-a Ltim_averages_masked_leur_no_auto_DI_20.txt \
-b Positive_cactus4_10bp.bed -wb | \
awk -F " " '{print $1,$2,$3,$4, ($2-$6)}' \
> Cactus4_Positive_positions_cline_HI_10bp_no_auto_DI_20.txt 
