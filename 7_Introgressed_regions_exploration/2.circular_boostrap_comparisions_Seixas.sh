#sort bed files
sort -k1,1 -k2,2n -k3,3n All_genes_no_loc.bed \
> All_genes_no_loc.sorted.bed

sort -k1,1 -k2,2n -k3,3n Seixas_genes_filt.bed \
> Seixas_genes_filt.sorted.bed

#overlaps all genes
bedtools intersect -a saguaro_iter12_cactus4_windows.bed \
-b All_genes_no_loc.bed -wo \
| awk '{sum += $NF} END {print sum}' \
> Cactus4_bp_overlaps_all_genes.txt


#overlaps seixas genes
bedtools intersect -a saguaro_iter12_cactus4_windows.bed \
-b Seixas_genes_filt.sorted.bed -wo \
| awk '{sum += $NF} END {print sum}' \
> Cactus4_bp_overlaps_seixas_genes.txt


#bootsraps
# overlaps all genes
for i in {1..2000}; do
    bedtools intersect \
        -a circular_boostrap/bootstrap_cactus4_rep_${i}.bed \
        -b All_genes_no_loc.bed \
        -wo \
    | awk '{sum += $NF} END {print sum+0}' \
    >> Bootstrap_Cactus4_bp_overlaps_all_genes.txt
done

# overlaps Seixas genes
for i in {1..2000}; do
    bedtools intersect \
        -a circular_boostrap/bootstrap_cactus4_rep_${i}.bed \
        -b Seixas_genes_filt.sorted.bed \
        -wo \
    | awk '{sum += $NF} END {print sum+0}' \
    >> Bootstrap_Cactus4_bp_overlaps_seixas_genes.txt
done
