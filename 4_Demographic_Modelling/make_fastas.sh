#make fastas for each locus
#!/bin/bash

BED="all_sampled_windows.bed"
REF="GCA_033115175.1_mLepEur2.pri_genomic_chr.fna"
VCF="High_cov_mono_filt.vcf.gz"
SAMPLES=("Lcas1" "Lcas2" "Lcas6" "Lcas7" "Lcas9" "Lcor10" "Lcor6" "Lcor7" "Lcor8" "Lcor9" "Ltim1" "Ltim2" "Ltim3" "Ltim4" "Ltim6")  # <-- replace with your individuals

while read chrom start end; do
    region="${chrom}_${start}_${end}"
    
    # Extract consensus for each sample
    for sample in "${SAMPLES[@]}"; do
        out="${region}_${sample}.fasta"
        samtools faidx $REF ${chrom}:${start}-${end} | \
        bcftools consensus -s $sample --haplotype 1pIu -I --missing N $VCF | \
        sed "s/^>.*$/>${sample}/" > $out
    done

    # Merge all samples into one FASTA per region
    cat ${region}_*.fasta > hares_${region}.fasta

    # Optional cleanup
    rm ${region}_*.fasta
done < $BED
