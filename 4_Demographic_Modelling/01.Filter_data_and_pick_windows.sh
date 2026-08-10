# Keep monomorphic SNPs and remove missing sites also filter by coverage to have failtul genotypes.
bcftools view -S samples.list -M2 Hares_Leur_ref_all.bcf.gz --threads 50 | \
bcftools filter -e 'F_MISSING>0 || QUAL<30' --threads 50 | \
bcftools filter -e 'FMT/DP<5 || FMT/DP>224' -Ob -o Lcas_Lcor_Ltim_bpp.bcf.gz --threads 50


#first get only high coverage Lcas Lcor and Ltim
bcftools view -S samples.list High_cov_mono_filt.vcf.gz --threads 50 | \
bcftools filter -e 'F_MISSING>0 || QUAL<30' --threads 50 | \
bcftools filter -e 'FMT/DP<5 || FMT/DP>224' -Ob -o Lcas_Lcor_Ltim_high_cov_bpp.bcf.gz --threads 50

#get bed file with genes
awk '$3=="gene"' GCF_033115175.1_mLepTim1.pri_genomic_chr.gff > genes.gff
awk -F " " '{print $1,$4,$5}' genes.gff > genes.bed

#get intergenic regions
sed -i 's/ \+/\t/g' genes.bed
grep "chr" genes.bed | grep -v "MT" | grep -v "Y" | grep -v "X" > genes_chr.bed
bedtools complement -i genes_chr.bed -g GCA_033115175.1_mLepEur2.pri_genomic_chr.fna.fai > non_genic.bed
grep "chr" non_genic.bed | grep -v "X" | grep -v "Y" > non_genic_chr.bed

#make windows
bedtools makewindows -b non_genic_chr.bed -w 1000 -s 50000 > non_genic_1000bp_windows_50K_spaced.bed

#check number of windows per chr
cut -f1 non_genic_1000bp_windows_50k_spaced.bed | sort | uniq -c

#pick 1000 windows with proportional representation of each chromossome
bash pick_windows.sh

#make fastas for each window
bash make_fastas.sh

#make bpp input
python make_bpp.py

#create tree
~/A_Lcas_Lcor_paper/software/bpp-4.8.0-linux-x86_64/bin/bpp --msci-create msci.txt
#run bpp
~/A_Lcas_Lcor_paper/software/bpp-4.8.0-linux-x86_64/bin/bpp --cfile
