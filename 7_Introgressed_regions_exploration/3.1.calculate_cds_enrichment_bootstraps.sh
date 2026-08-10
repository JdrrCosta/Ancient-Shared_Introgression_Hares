#!/bin/bash

# Input files
CDS="CDS_autossomes.bed"

# Output header for counts file
echo -e "Bootstrap\tCDS\tNon_CDS" > CDS_cactus4_bootstrap_counts_bp_no_overlaps.tsv

# Loop over bootstrap BED files
for BEDFILE in ../circular_boostrap/bootstrap_cactus4_rep_{1..2000}.bed; do
    BASENAME=$(basename "$BEDFILE" .bed)
    
    # Intersect fragments with CDS regions (and merge overlapping intervals)
    bedtools intersect -a "$BEDFILE" -b $CDS | sort -k1,1 -k2,2n | bedtools merge -i - > "${BASENAME}_bp_in_cds.bed"
    
    # Intersect fragments with non-CDS regions (complement)
    bedtools intersect -v -a "$BEDFILE" -b $CDS | sort -k1,1 -k2,2n | bedtools merge -i - > "${BASENAME}_bp_in_non_cds.bed"
    
    # Count basepairs in each category
    CDS_COUNT=$(awk '{sum += $3 - $2} END {print sum}' < "${BASENAME}_bp_in_cds.bed")
    NON_CDS_COUNT=$(awk '{sum += $3 - $2} END {print sum}' < "${BASENAME}_bp_in_non_cds.bed")
    
    # Append counts to file
    echo -e "${BASENAME}\t${CDS_COUNT}\t${NON_CDS_COUNT}" >> CDS_bootstrap_counts_bp_no_overlaps.tsv
done
