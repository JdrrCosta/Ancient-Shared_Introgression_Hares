#Produce inputs for PCangsd
#Using samtools method for GL
#Major allele is the ref
angsd -GL 1 -out Lcas_GL -nThreads 15 -doGlf 2 -doMajorMinor 5 -anc ~/Z_refs/GCA_033115175.1_mLepEur2.pri_genomic_chr.fna \
-SNP_pval 1e-6 -doMaf 1 -minQ 30 -minMapQ 30 -minMAF 0.034 -uniqueOnly -remove_bads 1 -minInd 11 -doCounts 1 -setMaxDepth 132 \
-doSaf 1 -bam Lcas_bam.filelist

#only autossomes
zcat Lcas_GL.beagle.gz | awk '$1 ~ /^chr[0-9]+/' | gzip > Lcas_GL_chr_only.beagle.gz

#split per chr
for i in {1..23}; do
    zcat Lcas_GL_chr_only.beagle.gz | awk -F "_" -v chr="chr$i" '$1 == chr {print}' > Lcas_GL_chr$i.beagle
done

#Admixture
for run in {1..5}; do
  for k in {1..5}; do
    NGSadmix -likes Lcas_GL_chr_only.beagle.gz -K $k -o Lcas_NGSadmix_k${k}.${run} -P 20
  done
done


