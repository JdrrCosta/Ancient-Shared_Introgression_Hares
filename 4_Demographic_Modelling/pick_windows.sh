#!/bin/bash
# Sample windows per chromosome proportionally for BPP
# Ensures that each sampled window has length >= 1000 bp

INPUT_BED="non_genic_1000bp_windows_50K_spaced.bed"
OUTPUT_DIR="samples"
mkdir -p "$OUTPUT_DIR"

# Array of chromosomes
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23)

# Corresponding number of windows to sample
NUMS=(60 60 70 60 70 50 60 50 40 50 50 40 40 30 40 30 30 40 40 30 30 20 10)

# Loop over chromosomes
for i in "${!CHRS[@]}"; do
    chr=${CHRS[$i]}
    n=${NUMS[$i]}
    
    # Create a temporary file with only >=1000bp windows for this chromosome
    tmpfile=$(mktemp)
    awk -v c="$chr" '($1==c && ($3-$2)>=1000)' "$INPUT_BED" > "$tmpfile"
    
    # Sample windows and save to output
    bedtools sample -n "$n" -i "$tmpfile" > "$OUTPUT_DIR/${chr}_sample.bed"
    
    # Remove temporary file
    rm "$tmpfile"
done

# Optional: merge all sampled windows into one file
cat "$OUTPUT_DIR"/*_sample.bed > "$OUTPUT_DIR/all_sampled_windows.bed"

echo "Sampling complete. Check $OUTPUT_DIR/"

