#intersect the cactus12 regions with genes
#relaxed
bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff -b cactus_iter12_cactus4_windows.bed | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_cactus4.txt


#write command for all 2000 files
for i in {1..2000}; do
  echo "Processing bootstrap ${i}..."
  bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
  -b ../circular_boostrap/bootstrap_cactus4_rep_${i}.bed | \
  awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
  uniq | wc -l >> cactus4_bootstrap_gene_number.txt
done
