# Step 1: Intersect fragments with CDS regions
bedtools intersect -a cactus_iter12_cactus4_windows.bed \
-b CDS_autossomes.bed > cactus4_bp_in_cds.bed

# Step 2: Intersect fragments with non-CDS regions (everything else)
bedtools intersect -v \
-a cactus_iter12_cactus4_windows.bed \
-b CDS_autossomes.bed > cactus4_bp_in_non_cds.bed


# Step 3: Calculate total basepairs
echo "Basepairs in CDS regions:"
awk '{sum += $3 - $2} END {print sum}' cactus4_bp_in_cds.bed


echo "Basepairs in non-CDS regions:"
awk '{sum += $3 - $2} END {print sum}' cactus4_bp_in_non_cds.bed


# Step 4: Remove overlaps for non-redundant basepair counts
sort -k1,1 -k2,2n cactus4_bp_in_cds.bed | bedtools merge -i - \
> cactus4_bp_in_cds_no_overlap.bed
sort -k1,1 -k2,2n cactus4_bp_in_non_cds.bed | bedtools merge -i - \
> cactus4_bp_in_non_cds_no_overlap.bed

# Step 5: Recalculate basepairs without overlaps
echo "Basepairs in CDS regions (no overlaps):"
awk '{sum += $3 - $2} END {print sum}' cactus4_bp_in_cds_no_overlap.bed

echo "Basepairs in non-CDS regions (no overlaps):"
awk '{sum += $3 - $2} END {print sum}' cactus4_bp_in_non_cds_no_overlap.bed

