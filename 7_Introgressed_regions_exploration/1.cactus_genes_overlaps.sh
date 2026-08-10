#genes in cactus4
bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff -b cactus_iter12_cactus4_windows.bed | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_cactus4.txt

#overlap between all genes with fragments in cactus4 and 
#Nuclear genes overlaps with introgression frequency outliers (Seixas et al.2018)
grep -w -f genes_in_cactus4.txt \
Nuclear_overlaps_outlier_intr_freq_seixas.txt

#Mitonuc genes with geographic introgression patterns similar to mtdna (Seixas et al.2018)
grep -w -f genes_in_cactus4.txt \
Mitonuc_genes_geographic_patterns_sim_mtdna.txt

#cds in cactus4
bedtools intersect \
-a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
-b cactus_iter12_cactus4_windows.bed \
| awk '$3 == "CDS"' | awk -F "gene=" '{print $2}' \
| awk -F ";" '{print $1}' | uniq > cds_in_cactus4_relaxed.txt

#mitogenes with cds in cactsu4
grep -i -w -F -f cds_in_cactus4_relaxed.txt \
mitogenes_mice.txt


#Seixas genes with cds in cactus4
grep -i -w -f \
Genes_CDS/cds_in_cactus4_relaxed.txt \
Seixas_overlaps/Nuclear_overlaps_outlier_intr_freq_seixas.txt

grep -i -w -f \
Genes_CDS/cds_in_cactus4_relaxed.txt \
Seixas_overlaps/Mitonuc_genes_geographic_patterns_sim_mtdna.txt


#genes in cactus12
bedtools intersect -a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff -b saguaro_iter12_cactus12_windows.bed | \
awk '$3 == "gene" {split($9, a, ";"); split(a[1], b, "-"); print b[2]}' | \
uniq > genes_in_cactus12.txt

#Overlap between all genes with fragments in cactus 12 and 
#Mitonuc genes with geographic introgression patterns similar to mtdna (Seixas et al.2025)
grep -w -f genes_in_cactus12.txt Mitonuc_genes_geographic_patterns_sim_mtdna.txt

#Nuclear genes overlaps woth introgression frequencyo utliers
grep -w -f genes_in_cactus12.txt \
Nuclear_overlaps_outlier_intr_freq_seixas.txt

#cds in cactus12
bedtools intersect \
-a GCF_033115175.1_mLepTim1.pri_genomic_chr.gff \
-b saguaro_iter12_cactus12_windows.bed -F 1 \
| awk '$3 == "CDS"' | awk -F "gene=" '{print $2}' \
| awk -F ";" '{print $1}' | uniq > cds_in_cactus12.txt

#mitogenes in cactus12 cds
grep -i -w -F -f cds_in_cactus12.txt mitogenes_mice.txt

#seixas genes with cds in cactus12
grep -i -w -f \
../genes/CDS_in_cactus12.txt \
Seixas_overlaps/Nuclear_overlaps_outlier_intr_freq_seixas.txt

grep -i -w -f \
../genes/CDS_in_cactus12.txt \
Seixas_overlaps/Mitonuc_genes_geographic_patterns_sim_mtdna.txt
